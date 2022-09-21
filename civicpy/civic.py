import requests
from requests.packages.urllib3.util.retry import Retry
import importlib
import logging
import pandas as pd
import pickle
import os
from pathlib import Path
from collections import defaultdict, namedtuple
from civicpy import REMOTE_CACHE_URL, LOCAL_CACHE_PATH, CACHE_TIMEOUT_DAYS
import requests
from civicpy import exports
import deprecation
from civicpy.__version__ import __version__
from datetime import datetime, timedelta
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()
import re


CACHE = dict()

COORDINATE_TABLE = None
COORDINATE_TABLE_START = None
COORDINATE_TABLE_STOP = None
COORDINATE_TABLE_CHR = None

HPO_TERMS = dict()

FRESH_DELTA = timedelta(days=CACHE_TIMEOUT_DAYS)

MODULE = importlib.import_module('civicpy.civic')

API_URL = 'https://civicdb.org/api/graphql'

LINKS_URL = 'https://civicdb.org/links'

UNMARKED_PLURALS = {'evidence'}

CIVIC_TO_PYCLASS = {
    'evidence_items': 'evidence'
}


_CoordinateQuery = namedtuple('CoordinateQuery', ['chr', 'start', 'stop', 'alt', 'ref', 'build', 'key'])
_CoordinateQuery.__new__.__defaults__ = (None, None, "GRCh37", None)


class CoordinateQuery(_CoordinateQuery):  # Wrapping for documentation
    """
    A namedtuple with preset fields describing a genomic coordinate,
    for use with coordinate-based queries of CIViC Variants.

    :param str chr: A chromosome of value 1-23, X, Y
    :param int start: The chromosomal start position in base coordinates (1-based)
    :param int stop: The chromosomal stop position in base coordinates (1-based)
    :param str optional alt: The alternate nucleotide(s) at the designated coordinates
    :param str optional ref: The reference nucleotide(s) at the designated coordinates
    :param NCBI36,GRCh37,GRCh38 build: The reference build version of the coordinates
    :param Any optional key: A user-defined object linked to the coordinate
    """
    pass


def pluralize(string):
    if string in UNMARKED_PLURALS:
        return '{}_items'.format(string)
    if string.endswith('s'):
        return string
    return string + 's'


def singularize(string):
    string = string.rstrip('s')
    if string == 'evidence_item':
        string = 'evidence'
    return string


def search_url(element, use_search_meta):
    element = pluralize(element).lower()
    components = [API_URL, element]
    if use_search_meta:
        components.append('search')
    return '/'.join(components)


def snake_to_camel(snake_string):
    words = snake_string.split('_')
    cap_words = [x.capitalize() for x in words]
    return ''.join(cap_words)


def element_lookup_by_id(element_type, element_id):
    e = _request_by_ids(element_type, [int(element_id)])[0]
    e = _postprocess_response_element(e, element_type)
    return e


def get_class(element_type):
    e_string = singularize(element_type)
    class_string = snake_to_camel(e_string)
    cls = getattr(MODULE, class_string, CivicAttribute)
    return cls


def download_remote_cache(remote_cache_url=REMOTE_CACHE_URL, local_cache_path=LOCAL_CACHE_PATH):
    """
    Retrieve a remote cache file from URL and save to local filepath.

    :param remote_cache_url:    A URL string to a remote cache for retrieval.
                                This parameter defaults to REMOTE_CACHE_URL.

    :param local_cache_path:    A filepath destination string for the retrieved remote cache.
                                This parameter defaults to LOCAL_CACHE_PATH.

    :return:                    Returns True on success.
    """
    logging.warning(
        'Downloading remote cache from {}.'.format(remote_cache_url)
    )
    _make_local_cache_path_if_missing(local_cache_path)
    r = requests.get(remote_cache_url)
    r.raise_for_status()
    with open(local_cache_path, 'wb') as local_cache:
        local_cache.write(r.content)
    return True


def save_cache(local_cache_path=LOCAL_CACHE_PATH):
    """
    Save in-memory cache to local file.

    :param local_cache_path:    A filepath destination string for storing the cache.
                                This parameter defaults to LOCAL_CACHE_PATH.

    :return:                    Returns True on success.
    """
    with open(local_cache_path, 'wb') as pf:
        pickle.dump(CACHE, pf)
    return True


def cache_file_present(local_cache_path=LOCAL_CACHE_PATH):
    """
    Determines if a file exists at a given path.

    :param local_cache_path:    A filepath where cache is expected.
                                This parameter defaults to LOCAL_CACHE_PATH.

    :return:                    Returns True on success.
    """
    return os.path.isfile(local_cache_path)


def delete_local_cache(local_cache_path=LOCAL_CACHE_PATH):
    """
    Deletes local cache file.

    :param local_cache_path:    A filepath destination string for the cache to be deleted.
                                This parameter defaults to LOCAL_CACHE_PATH.

    :return:            Returns True on success.
    """
    return os.unlink(local_cache_path)


def load_cache(local_cache_path=LOCAL_CACHE_PATH, on_stale='auto'):
    """
    Load local file to in-memory cache.

    :param local_cache_path:    A filepath destination string for loading the cache.
                                This parameter defaults to LOCAL_CACHE_PATH.

    :param on_stale:    ['auto', 'reject', 'ignore', 'update']
                        auto:   If cache_path matches the filepath in
                                LOCAL_CACHE_PATH, follows 'update' behavior.
                                Otherwise, follows 'reject' behavior.
                        reject: Clear loaded cache from memory if stale.
                        ignore: Keep loaded cache in memory if stale.
                        update: Run update_cache and save fresh cache to
                                cache_path.
                        This parameter defaults to 'auto'.

    :return:            Returns True if content is loaded to in-memory cache.
    """
    downloaded_remote = False
    remote_url = REMOTE_CACHE_URL
    if local_cache_path == LOCAL_CACHE_PATH:
        if not cache_file_present():
            download_remote_cache(remote_cache_url=remote_url)
            downloaded_remote = True
    elif not cache_file_present(local_cache_path):
        raise FileNotFoundError("No cache found at {}".format(local_cache_path))
    with open(local_cache_path, 'rb') as pf:
        loaded_cache = pickle.load(pf)
    c = dict()
    variants = set()
    for k, v in loaded_cache.items():
        if isinstance(k, str):
            c[k] = v
        elif isinstance(k, int):
            c[hash(v)] = v
            if v.type == 'variant':
                variants.add(v)
        else:
            raise ValueError
    old_cache = MODULE.CACHE
    MODULE.CACHE = c
    for k, v in MODULE.CACHE.items():
        if isinstance(k, str):
            continue
        v.update()
    if _has_full_cached_fresh() or on_stale == 'ignore':
        _build_coordinate_table(variants)
        return True
    elif (on_stale == 'auto' and local_cache_path == LOCAL_CACHE_PATH) or on_stale == 'update':
        MODULE.CACHE = old_cache
        if downloaded_remote:
            logging.error(
                'Remote cache at {} is stale. Consider running `update_cache(from_remote_cache=False)` '
                "to create cache from API query (slow), or `load_cache(on_stale='ignore')` "
                "to load stale local cache (if present). "
                'Please create an issue at https://github.com/griffithlab/civicpy/issues '
                'if this is unexpected behavior.'.format(remote_url)
            )
            raise SystemError
        else:
            logging.warning(
                'Local cache at {} is stale, updating from remote.'.format(local_cache_path)
            )
            update_cache(local_cache_path=local_cache_path)
            return True
    elif on_stale == 'reject' or on_stale == 'auto':
        MODULE.CACHE = old_cache
        logging.warning(
            'Local cache at {} is stale and was not loaded. To load anyway, re-run '
            '`load_cache` with `on_stale` parameter set to desired behavior.'.format(local_cache_path)
        )
        return False
    raise NotImplementedError  # An unexpected condition occurred.


def update_cache(from_remote_cache=True, remote_cache_url=REMOTE_CACHE_URL,
                 local_cache_path=LOCAL_CACHE_PATH):
    """
    Update local cache file.

    :param from_remote_cache:   If set to True, update_cache will first download the
                                remote cache designated by REMOTE_CACHE_URL, store it
                                to LOCAL_CACHE_PATH, and then load the downloaded cache
                                into memory.
                                This parameter defaults to True.

    :param remote_cache_url:    A URL string to a remote cache for retrieval.
                                This parameter defaults to REMOTE_CACHE_URL.

    :param local_cache_path:    A filepath destination string for the retrieved remote cache.
                                This parameter defaults to LOCAL_CACHE_PATH.

    :return:                    Returns True on success.
    """
    _make_local_cache_path_if_missing(local_cache_path)
    if from_remote_cache:
        download_remote_cache(local_cache_path=local_cache_path, remote_cache_url=remote_cache_url)
        load_cache(local_cache_path=local_cache_path)
    else:
        genes = _get_elements_by_ids('gene', allow_cached=False, get_all=True)
        variants = _get_elements_by_ids('variant', allow_cached=False, get_all=True)
        evidence = _get_elements_by_ids('evidence', allow_cached=False, get_all=True)
        assertions = _get_elements_by_ids('assertion', allow_cached=False, get_all=True)
        variant_groups = _get_elements_by_ids('variant_group', allow_cached=False, get_all=True)
        for e in evidence:
            e.assertions = [a for a in assertions if a.id in e.assertion_ids]
            e._partial = False
            CACHE[hash(e)] = e
        for g in genes:
            g.variants = [v for v in variants if v.gene_id == g.id]
            g.assertions = [a for a in assertions if a.gene_id == g.id]
            g._partial = False
            CACHE[hash(g)] = g
        for v in variants:
            v.assertions = [a for a in assertions if a.variant_id == v.id]
            v.evidence_items = [e for e in evidence if e.variant_id == v.id]
            v.variant_groups = [vg for vg in variant_groups if v.id in vg.variant_ids]
            v._partial = False
            CACHE[hash(v)] = v
        for a in assertions:
            a.evidence_items = [e for e in evidence if e.id in a.evidence_ids]
            CACHE[hash(a)] = a
        for vg in variant_groups:
            vg.variants = [v for v in variants if v.id in vg.variant_ids]
            vg._partial = False
            CACHE[hash(vg)] = vg
        CACHE['full_cached'] = datetime.now()
        _build_coordinate_table(variants)
        save_cache(local_cache_path=local_cache_path)


def _make_local_cache_path_if_missing(local_cache_path):
    p = Path(local_cache_path)
    if not p.parent.is_dir():
        os.makedirs(p.parent)


class CivicRecord:
    """
    As a base class, :class:`CivicRecord` is used to define the characteristic of all records in CIViC. This class is not
    intended to be invoked directly by the end user, but provided for documentation of shared methods and variables in
    child classes.
    """

    _SIMPLE_FIELDS = {'id', 'type'}
    _COMPLEX_FIELDS = set()
    _OPTIONAL_FIELDS = set()

    def __init__(self, partial=False, **kwargs):
        """
        The record object may be initialized by the user, though the practice is discouraged. To do so, values for each
        of the object attributes (except ``type``) must be specified as keyword arguments, or the ``partial`` parameter must
        be set to **True**. If ``partial`` is set to **True**, the ``id`` keyword argument is still required.

        Users are encouraged to use the functions for `getting records`_ in lieu of directly initializing record
        objects.

        :param bool partial: Indicates whether the the set of object attributes passed is incomplete. If set to **True** the ``id`` keyword is required.
        """
        self._incomplete = set()
        self._partial = partial
        simple_fields = sorted(self._SIMPLE_FIELDS, reverse=True)
        simple_fields = sorted(simple_fields, key=lambda x: x in CivicRecord._SIMPLE_FIELDS, reverse=True)
        for field in simple_fields:
            try:
                self.__setattr__(field, kwargs[field])
            except KeyError:
                try:
                    object.__getattribute__(self, field)
                except AttributeError:
                    if (partial and field not in CivicRecord._SIMPLE_FIELDS) or field in self._OPTIONAL_FIELDS:
                        self._incomplete.add(field)     # Allow for incomplete data when partial flag set
                    else:
                        raise AttributeError('Expected {} attribute for {}, none found.'.format(field, self.type))

        for field in self._COMPLEX_FIELDS:
            try:
                v = kwargs[field]
                if v is None:
                    v = dict()
            except KeyError:
                if partial or field in self._OPTIONAL_FIELDS:
                    self._incomplete.add(field)
                    continue
                else:
                    raise AttributeError('Expected {} attribute for {}, none found.'.format(field, self.type))
            is_compound = isinstance(v, list)
            cls = get_class(field)
            if is_compound:
                result = list()
                for data in v:
                    if isinstance(data, dict):
                        data['type'] = data.get('type', singularize(field))
                        result.append(cls(partial=True, **data))
                    else:
                        result.append(data)
                self.__setattr__(field, result)
            else:
                t = v.get('type', field)
                v['type'] = CIVIC_TO_PYCLASS.get(t, t)
                if v.keys() == {'type'}:
                    self.__setattr__(field, {})
                else:
                    self.__setattr__(field, cls(partial=True, **v))

        self._partial = bool(self._incomplete)
        if not isinstance(self, CivicAttribute) and not self._partial and self.__class__.__name__ != 'CivicRecord':
            CACHE[hash(self)] = self

        self._include_status = ['accepted','submitted','rejected']

    def __dir__(self):
        return [attribute for attribute in super().__dir__() if not attribute.startswith('_')]

    def __repr__(self):
        return '<CIViC {} {}>'.format(self.type, self.id)

    def __getattr__(self, item):
        if self._partial and item in self._incomplete:
            self.update()
        return object.__getattribute__(self, item)

    def __hash__(self):
        return hash('{}:{}'.format(self.type, self.id))

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __setstate__(self, state):
        self.__dict__ = state

    def update(self, allow_partial=True, force=False, **kwargs):
        """
        Updates the record object from the cache or the server.
        Keyword arguments may be passed to ``kwargs``, which will update the corresponding attributes of the
        :class:`CivicRecord` instance.

        :param bool allow_partial: Flag to indicate whether the record will be updated according to the contents of CACHE, without requiring all attributes to be assigned.
        :param bool force: Flag to indicate whether to force an update from the server, even if a full record ecists in the cache.
        :return: True if record is complete after update, else False.
        """
        if kwargs:
            self.__init__(partial=allow_partial, force=force, **kwargs)
            return not self._partial

        if not force and CACHE.get(hash(self)):
            cached = CACHE[hash(self)]
            for field in self._SIMPLE_FIELDS | self._COMPLEX_FIELDS:
                v = getattr(cached, field)
                setattr(self, field, v)
            self._partial = False
            logging.info('Loading {} from cache'.format(str(self)))
            return True
        resp_dict = element_lookup_by_id(self.type, self.id)
        self.__init__(partial=False, **resp_dict)
        return True

    @property
    def site_link(self):
        """Returns a URL to the record on the CIViC web application."""
        return '/'.join([LINKS_URL, self.type, str(self.id)])


class Variant(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'allele_registry_id',
        'civic_actionability_score',
        'description',
        'gene_id',
        'name',
        'entrez_name',
        'entrez_id'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'assertions',
        'clinvar_entries',
        'coordinates',
        # 'errors',
        'evidence_items',
        'hgvs_expressions',
        #'lifecycle_actions',
        # 'provisional_values',
        'sources',
        'variant_aliases',
        'variant_groups',
        'variant_types'})

    def __init__(self, **kwargs):
        # Handle overloaded evidence_items from some advanced search views
        kwargs['type'] = 'variant'
        self._evidence_items = []
        self._assertions = []
        self._variant_groups = []
        coordinates = kwargs.get('coordinates')
        if coordinates:
            if coordinates.get('reference_bases') in ['', '-']:
                coordinates['reference_bases'] = None
            if coordinates.get('variant_bases') in ['', '-']:
                coordinates['variant_bases'] = None
        super().__init__(**kwargs)

    @property
    def evidence_sources(self):
        sources = set()
        for evidence in self.evidence_items:
            if evidence.source is not None:
                sources.add(evidence.source)
        return sources

    @property
    def aliases(self):
        return self.variant_aliases

    @property
    def groups(self):
        return self.variant_groups

    @property
    def types(self):
        return self.variant_types

    @property
    def summary(self):
        return self.description

    @summary.setter
    def summary(self, value):
        self.description = value

    @property
    def evidence(self):
        return self.evidence_items

    @property
    def evidence_items(self):
        return [e for e in self._evidence_items if e.status in self._include_status]

    @evidence_items.setter
    def evidence_items(self, value):
        self._evidence_items = value

    @property
    def assertions(self):
        return [a for a in self._assertions if a.status in self._include_status]

    @assertions.setter
    def assertions(self, value):
        self._assertions = value

    @property
    def variant_groups(self):
        return self._variant_groups

    @variant_groups.setter
    def variant_groups(self, value):
        self._variant_groups = value

    @property
    def gene(self):
        return _get_element_by_id('gene', self.gene_id)

    @property
    def is_insertion(self):
        ref = self.coordinates.reference_bases
        alt = self.coordinates.variant_bases
        return (ref is None and alt is not None) or (ref is not None and alt is not None and len(ref) < len(alt))

    @property
    def is_deletion(self):
        ref = self.coordinates.reference_bases
        alt = self.coordinates.variant_bases
        if alt is not None and (alt == '-' or alt == ''):
            alt = None
        return (ref is not None and alt is None) or (ref is not None and alt is not None and len(ref) > len(alt))

    def is_valid_for_vcf(self, emit_warnings=False):
        if self.coordinates.chromosome2 or self.coordinates.start2 or self.coordinates.stop2:
            warning = "Variant {} has a second set of coordinates. Skipping".format(self.id)
        elif self.coordinates.chromosome and self.coordinates.start and (self.coordinates.reference_bases or self.coordinates.variant_bases):
            if self._valid_ref_bases():
                if self._valid_alt_bases():
                    return True
                else:
                    warning = "Unsupported variant base(s) for variant {}. Skipping.".format(self.id)
            else:
                warning = "Unsupported reference base(s) for variant {}. Skipping.".format(self.id)
        else:
            warning = "Incomplete coordinates for variant {}. Skipping.".format(self.id)
        if emit_warnings:
            logging.warning(warning)
        return False

    def _valid_ref_bases(self):
        if self.coordinates.reference_bases is not None:
            return all([c.upper() in ['A', 'C', 'G', 'T', 'N'] for c in self.coordinates.reference_bases])
        else:
            return True

    def _valid_alt_bases(self):
        if self.coordinates.variant_bases is not None:
            return all([c.upper() in ['A', 'C', 'G', 'T', 'N'] for c in self.coordinates.variant_bases])
        else:
            return True

    def vcf_coordinates(self):
        ensembl_server = "https://grch37.rest.ensembl.org"
        if self.coordinates.reference_build != 'GRCh37':
            return
        if self.is_insertion:
            if not self.coordinates.representative_transcript:
                return
            else:
                start = self.coordinates.start
                ext = "/sequence/region/human/{}:{}-{}".format(self.coordinates.chromosome, start, start)
                r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "text/plain"})
                if self.coordinates.reference_bases == None or self.coordinates.reference_bases == '-' or self.coordinates.reference_bases == '':
                    ref = r.text
                else:
                    ref = "{}{}".format(r.text, self.coordinates.reference_bases)
                alt = "{}{}".format(r.text, self.coordinates.variant_bases)
        elif self.is_deletion:
            if not self.coordinates.representative_transcript:
                return
            else:
                start = self.coordinates.start - 1
                ext = "/sequence/region/human/{}:{}-{}".format(self.coordinates.chromosome, start, start)
                r = requests.get(ensembl_server+ext, headers={ "Content-Type" : "text/plain"})
                ref = "{}{}".format(r.text, self.coordinates.reference_bases)
                if self.coordinates.variant_bases == None or self.coordinates.variant_bases == '-' or self.coordinates.variant_bases == '':
                    alt = r.text
                else:
                    alt = "{}{}".format(r.text, self.coordinates.variant_bases)
        else:
            start = self.coordinates.start
            ref = self.coordinates.reference_bases
            alt = self.coordinates.variant_bases
        return (start, ref, alt)

    def csq_alt(self):
        if self.coordinates.reference_build != 'GRCh37':
            return
        if self.is_insertion:
            if not self.coordinates.representative_transcript:
                return
            else:
                return self.coordinates.variant_bases
        elif self.is_deletion:
            if not self.coordinates.representative_transcript:
                return
            else:
                return "-"
        else:
            return self.coordinates.variant_bases

    def hgvs_c(self):
        if self.coordinates.representative_transcript:
            hgvs_cs = [e for e in self.hgvs_expressions if (':c.' in e) and (self.coordinates.representative_transcript in e)]
            return hgvs_cs[0] if len(hgvs_cs) == 1 else ''
        else:
            return ''

    def hgvs_p(self):
        if self.coordinates.representative_transcript:
            hgvs_ps = [e for e in self.hgvs_expressions if (':p.' in e) and (self.coordinates.representative_transcript in e)]
            return hgvs_ps[0] if len(hgvs_ps) == 1 else ''
        else:
            return ''

    def sanitized_name(self):
        name = self.name
        regex = re.compile(r"^([A-Z]+)([0-9]+)(=)(.*)$")
        match = regex.match(name)
        if match is not None:
            name = "".join([match.group(1), match.group(2), match.group(1), match.group(4)])
        return name

    def csq(self, include_status=None):
        if self.csq_alt() is None:
            return []
        else:
            csq = []
            for evidence in self.evidence:
                if include_status is not None and evidence.status not in include_status:
                    continue
                special_character_table = str.maketrans(exports.VCFWriter.SPECIAL_CHARACTERS)
                csq.append('|'.join([
                    self.csq_alt(),
                    '&'.join(map(lambda t: t.name, self.variant_types)),
                    self.gene.name,
                    str(self.gene.entrez_id),
                    'transcript',
                    str(self.coordinates.representative_transcript),
                    self.hgvs_c(),
                    self.hgvs_p(),
                    self.sanitized_name(),
                    str(self.id),
                    '&'.join(map(lambda a: a.translate(special_character_table), self.variant_aliases)),
                    '&'.join(map(lambda e: e.translate(special_character_table), self.hgvs_expressions)),
                    str(self.allele_registry_id),
                    '&'.join(self.clinvar_entries),
                    str(self.civic_actionability_score),
                    "evidence",
                    str(evidence.id),
                    "https://civicdb.org/links/evidence/{}".format(evidence.id),
                    "{} ({})".format(evidence.source.citation_id, evidence.source.source_type),
                    str(evidence.variant_origin),
                    evidence.status,
                    str(evidence.clinical_significance or ''),
                    str(evidence.evidence_direction or ''),
                    str(evidence.disease),
                    '&'.join([str(drug) for drug in evidence.drugs]),
                    str(evidence.drug_interaction_type or ""),
                    '&'.join(["{} (HPO ID {})".format(phenotype.name, phenotype.hpo_id) for phenotype in evidence.phenotypes]),
                    evidence.evidence_level,
                    str(evidence.rating),
                    "",
                    "",
                    "",
                    "",
                    "",
                ]))
            for assertion in self.assertions:
                if include_status is not None and assertion.status not in include_status:
                    continue
                csq.append('|'.join([
                    self.csq_alt(),
                    '&'.join(map(lambda t: t.name, self.variant_types)),
                    self.gene.name,
                    str(self.gene.entrez_id),
                    'transcript',
                    str(self.coordinates.representative_transcript),
                    self.hgvs_c(),
                    self.hgvs_p(),
                    self.sanitized_name(),
                    str(self.id),
                    '&'.join(map(lambda a: a.translate(special_character_table), self.variant_aliases)),
                    '&'.join(map(lambda e: e.translate(special_character_table), self.hgvs_expressions)),
                    str(self.allele_registry_id),
                    '&'.join(self.clinvar_entries),
                    str(self.civic_actionability_score),
                    "assertion",
                    str(assertion.id),
                    "https://civicdb.org/links/assertion/{}".format(assertion.id),
                    "",
                    str(assertion.variant_origin),
                    assertion.status,
                    assertion.clinical_significance,
                    assertion.evidence_direction,
                    str(assertion.disease),
                    '&'.join([str(drug) for drug in assertion.drugs]),
                    str(assertion.drug_interaction_type or ''),
                    "",
                    "",
                    "",
                    "&".join([acmg_code.code for acmg_code in assertion.acmg_codes]),
                    "&".join([clingen_code.code for clingen_code in assertion.clingen_codes]),
                    str(assertion.amp_level or ''),
                    assertion.format_nccn_guideline(),
                    str(assertion.fda_regulatory_approval or ''),
                    str(assertion.fda_companion_test or ''),
                ]))
            return csq


class VariantGroup(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union(
        {'description', 'name', 'variant_ids'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        # 'errors',                 # TODO: Add support for these fields in advanced search endpoint
         # 'lifecycle_actions',
        # 'provisional_values',
        'sources',
        'variants'
    })

    def __init__(self, **kwargs):
        self._variants = []
        super().__init__(**kwargs)

    @property
    def variants(self):
        for variant in self._variants:
            variant._include_status = self._include_status
        return [v for v in self._variants if v.evidence]

    @variants.setter
    def variants(self, value):
        self._variants = value


class Gene(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union(
        {'description', 'entrez_id', 'name'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'aliases',
        # 'errors',                 # TODO: Add support for these fields in advanced search endpoint
        # /'lifecycle_actions',
        # 'provisional_values',
        'sources',
        'variants',
        'assertions',
    })

    def __init__(self, **kwargs):
        self._variants = []
        self._assertions = []
        super().__init__(**kwargs)

    @property
    def variants(self):
        for variant in self._variants:
            variant._include_status = self._include_status
        return [v for v in self._variants if v.evidence]

    @variants.setter
    def variants(self, value):
        self._variants = value

    @property
    def assertions(self):
        return [a for a in self._assertions if a.status in self._include_status]

    @assertions.setter
    def assertions(self, value):
        self._assertions = value


class Evidence(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'clinical_significance',
        'description',
        'drug_interaction_type',
        'evidence_direction',
        'evidence_level',
        'evidence_type',
        'gene_id',
        'name',
        # 'open_change_count',
        'rating',
        'status',
        'variant_id',
        'variant_origin',
        'assertion_ids',
    })
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'assertions',
        'disease',
        'drugs',
        # 'errors',
        # 'fields_with_pending_changes',
        #'lifecycle_actions',
        'phenotypes',
        'source'})

    def __init__(self, **kwargs):
        self._assertions = []
        super().__init__(**kwargs)

    @property
    def variant(self):
        return get_variant_by_id(self.variant_id)

    @property
    def gene(self):
        return get_gene_by_id(self.gene_id)

    @property
    def assertions(self):
        return [a for a in self._assertions if a.status in self._include_status]

    @assertions.setter
    def assertions(self, value):
        self._assertions = value

    @property
    def statement(self):
        return self.description

    @statement.setter
    def statement(self, value):
        self.description = value


class Assertion(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'amp_level',
        'clinical_significance',
        'description',
        'drug_interaction_type',
        'evidence_direction',
        # 'evidence_item_count',
        'evidence_type',
        'fda_companion_test',
        'fda_regulatory_approval',
        'name',
        'nccn_guideline',
        'nccn_guideline_version',
        # 'open_change_count',
        # 'pending_evidence_count',
        'status',
        'summary',
        'variant_origin',
        'gene_id',
        'variant_id',
        'evidence_ids',
    })

    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'acmg_codes',
        'clingen_codes',
        'disease',
        'drugs',
        'evidence_items',
        #'lifecycle_actions',
        'phenotypes',
    })

    def __init__(self, **kwargs):
        self._evidence_items = []
        super().__init__(**kwargs)

    @property
    def evidence(self):
        return self.evidence_items

    @property
    def evidence_items(self):
        return [e for e in self._evidence_items if e.status in self._include_status]

    @evidence_items.setter
    def evidence_items(self, value):
        self._evidence_items = value

    @property
    def hpo_ids(self):
        return [x.hpo_id for x in self.phenotypes if x.hpo_id]

    @property
    def variant(self):
        return get_variant_by_id(self.variant_id)

    @property
    def gene(self):
        return get_gene_by_id(self.gene_id)

    def format_nccn_guideline(self):
        if self.nccn_guideline is None:
            return ""
        else:
            return "{} (v{})".format(self.nccn_guideline, self.nccn_guideline_version)


class User(CivicRecord):

    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'name',
        'username',
        'role',
        'avatar_url',
        'area_of_expertise',
        'orcid',
        'display_name',
        'created_at',
        'url',
        'twitter_handle',
        'facebook_profile',
        'linkedin_profile',
        'bio',
        'featured_expert',
        # 'accepted_license',
        # 'signup_complete',
        # 'affiliation'
    })

    _OPTIONAL_FIELDS = CivicRecord._OPTIONAL_FIELDS.union({
        'country',
        'organization',
        'conflict_of_interest'
    })

    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union(_OPTIONAL_FIELDS)

    def __init__(self, **kwargs):
        self._created_at = None
        super().__init__(**kwargs)

    @property
    def created_at(self):
        assert self._created_at[-1] == 'Z'
        return datetime.fromisoformat(self._created_at[:-1])

    @created_at.setter
    def created_at(self, value):
        self._created_at = value


class Organization(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'name',
        'url',
        'description'
    })

    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'profile_image',
        'parent'
    })


class CivicAttribute(CivicRecord, dict):

    _SIMPLE_FIELDS = {'type'}
    _COMPLEX_FIELDS = set()

    def __repr__(self):
        try:
            _id = self.id
        except AttributeError:
            return '<CIViC Attribute {}>'.format(self.type)
        else:
            return '<CIViC Attribute {} {}>'.format(self.type, self.id)

    def __init__(self, **kwargs):
        kwargs['partial'] = False
        for k, v in kwargs.items():
            self.__setattr__(k, v)
        super().__init__(**kwargs)

    def __hash__(self):
        try:
            _id = self.id
        except AttributeError:
            raise NotImplementedError
        if _id is not None:
            return CivicRecord.__hash__(self)
        else:
            raise ValueError

    @property
    def site_link(self):
        return None

    def update(self):
        return NotImplementedError


class Drug(CivicAttribute):
    _SIMPLE_FIELDS = CivicAttribute._SIMPLE_FIELDS.union({'ncit_id', 'drug_url', 'name'})
    _COMPLEX_FIELDS = CivicAttribute._COMPLEX_FIELDS.union({'aliases'})

    def __str__(self):
        if self.ncit_id is None:
            return self.name
        else:
            return "{} (NCIt ID {})".format(self.name, self.ncit_id)


class Disease(CivicAttribute):
    _SIMPLE_FIELDS = CivicAttribute._SIMPLE_FIELDS.union({'display_name', 'doid', 'disease_url'})
    _COMPLEX_FIELDS = CivicAttribute._COMPLEX_FIELDS.union({'aliases'})

    def __str__(self):
        if self.doid is None:
            return self.name
        else:
            return "{} (DOID {})".format(self.name, self.doid)


class Phenotype(CivicAttribute):
    _SIMPLE_FIELDS = CivicAttribute._SIMPLE_FIELDS.union({'hpo_id', 'url', 'name'})

    def __str__(self):
        return "{} (HPO ID {})".format(self.name, self.hpo_id)


class Source(CivicAttribute):
    _SIMPLE_FIELDS = CivicAttribute._SIMPLE_FIELDS.union({
        'citation',
        'citation_id',
        'source_type',
        'abstract',
        'asco_abstract_id',
        'author_string',
        'full_journal_title',
        'journal',
        'pmc_id',
        'publication_date',
        'source_url',
        'title'
    })
    _COMPLEX_FIELDS = CivicAttribute._COMPLEX_FIELDS.union({'clinical_trials'})

    def __str__(self):
        return "{} ({} {})".format(self.citation, self.source_type, self.citation_id)


class Country(CivicAttribute):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({'iso', 'name'})


class LifecycleAction(CivicAttribute):
    _OPTIONAL_FIELDS = CivicAttribute._OPTIONAL_FIELDS.union({
        'submitted',
        'last_modified',
        'last_reviewed',
        'accepted'
    })
    _COMPLEX_FIELDS = CivicAttribute._COMPLEX_FIELDS.union(_OPTIONAL_FIELDS)


class BaseLifecycleAction(CivicAttribute):
    _SIMPLE_FIELDS = CivicAttribute._SIMPLE_FIELDS.union({
        'timestamp'
    })
    _COMPLEX_FIELDS = CivicAttribute._COMPLEX_FIELDS.union({
        'user'
    })

    def __init__(self, **kwargs):
        self._timestamp = None
        super().__init__(**kwargs)

    @property
    def timestamp(self):
        assert self._timestamp[-1] == 'Z'
        return datetime.fromisoformat(self._timestamp[:-1])

    @timestamp.setter
    def timestamp(self, value):
        self._timestamp = value


class Submitted(BaseLifecycleAction):
    pass


class LastModified(BaseLifecycleAction):
    pass


class LastReviewed(BaseLifecycleAction):
    pass


class Accepted(BaseLifecycleAction):
    pass


def get_cached(element_type, element_id):
    klass = get_class(element_type)
    r = klass(type=element_type, id=element_id, partial=True)
    return CACHE.get(hash(r), False)


def _has_full_cached_fresh(delta=FRESH_DELTA):
    s = 'full_cached'
    if CACHE.get(s, False):
        return CACHE[s] + delta > datetime.now()
    return False


def _get_elements_by_ids(element, id_list=[], allow_cached=True, get_all=False):
    if allow_cached:
        if not CACHE:
            load_cache()
        if not get_all:
            cached = [get_cached(element, element_id) for element_id in id_list]
            if all(cached):
                logging.info('Loading {} from cache'.format(pluralize(element)))
                return cached
        else:
            cached = [get_cached(element, element_id) for element_id in CACHE['{}_all_ids'.format(pluralize(element))]]
            logging.info('Loading {} from cache'.format(pluralize(element)))
            return cached
    if id_list and get_all:
        raise ValueError('Please pass list of ids or use the get_all flag, not both.')
    if get_all:
        logging.warning('Getting all {}. This may take a couple of minutes...'.format(pluralize(element)))
        response_elements = _request_all(element)
    else:
        response_elements = _request_by_ids(element, id_list)

    for e in response_elements:
        e = _postprocess_response_element(e, element)

    cls = get_class(element)
    elements = [cls(**x, partial=True) for x in response_elements]
    cache = [x['id'] for x in response_elements]
    CACHE['{}_all_ids'.format(pluralize(element))] = cache
    return elements


def _postprocess_response_element(e, element):
    e['type'] = element
    if element == 'assertion':
        e['gene_id'] = e['gene']['id']
        e['variant_id'] = e['variant']['id']
        e['evidence_ids'] = [evidence['id'] for evidence in e['evidenceItems']]
        e['status'] = e['status'].lower()
    elif element == 'evidence':
        e['gene_id'] = e['gene']['id']
        e['variant_id'] = e['variant']['id']
        e['assertion_ids'] = [a['id'] for a in e['assertions']]
        e['status'] = e['status'].lower()
    elif element == 'variant':
        e['gene_id'] = e['gene']['id']
        e['entrez_name'] = e['gene']['name']
        e['entrez_id'] = e['gene']['entrezId']
        build = e['referenceBuild']
        if build == 'GRCH37':
            build = 'GRCh37'
        elif build == 'GRCH38':
            build = 'GRCh38'
        e['coordinates'] = {
            'ensembl_version': e['ensemblVersion'],
            'reference_build': build,
            'reference_bases': e['referenceBases'],
            'variant_bases': e['variantBases'],
            'representative_transcript': None if e['primaryCoordinates'] is None else e['primaryCoordinates']['representativeTranscript'],
            'chromosome': None if e['primaryCoordinates'] is None else e['primaryCoordinates']['chromosome'],
            'start': None if e['primaryCoordinates'] is None else e['primaryCoordinates']['start'],
            'stop': None if e['primaryCoordinates'] is None else e['primaryCoordinates']['stop'],
            'representative_transcript2': None if e['secondaryCoordinates'] is None else e['secondaryCoordinates']['representativeTranscript'],
            'chromosome2': None if e['secondaryCoordinates'] is None else e['secondaryCoordinates']['chromosome'],
            'start2': None if e['secondaryCoordinates'] is None else e['secondaryCoordinates']['start'],
            'stop2': None if e['secondaryCoordinates'] is None else e['secondaryCoordinates']['stop'],
        }
    elif element == 'variant_group':
        e['variant_ids'] = [v['id'] for v in e['variants']['nodes']]
        del e['variants']
    return e


def _get_element_by_id(element, id, allow_cached=True):
    return _get_elements_by_ids(element, [id], allow_cached)[0]


def _request_by_ids(element, ids):
    payload_methods = {
        'evidence': _construct_get_evidence_payload,
        'gene': _construct_get_gene_payload,
        'variant': _construct_get_variant_payload,
        'assertion': _construct_get_assertion_payload,
        'variant_group': _construct_get_variant_group_payload,
    }
    payload_method = payload_methods[element]
    payload = payload_method()

    response_elements = []
    for i in ids:
        resp = requests.post(API_URL, json={'query': payload, 'variables': {'id': i}}, timeout=(10,200))
        resp.raise_for_status()
        response = resp.json()['data'][element]
        response_elements.append(response)
    return response_elements


def _request_all(element):
    payload_methods = {
        'evidence': _construct_get_all_evidence_payload,
        'gene': _construct_get_all_genes_payload,
        'variant': _construct_get_all_variants_payload,
        'assertion': _construct_get_all_assertions_payload,
        'variant_group': _construct_get_all_variant_groups_payload,
    }
    payload_method = payload_methods[element]
    payload = payload_method()

    after_cursor = None
    variables = { "after": after_cursor }
    resp = requests.post(API_URL, json={'query': payload, 'variables': variables}, timeout=(10,200))
    resp.raise_for_status()
    response = resp.json()['data'][pluralize(element)]
    response_elements = response['nodes']
    has_next_page = response['pageInfo']['hasNextPage']
    after_cursor = response['pageInfo']['endCursor']

    while has_next_page:
        variables = {
          "after": after_cursor
        }
        resp = requests.post(API_URL, json={'query': payload, 'variables': variables}, timeout=(10,200))
        resp.raise_for_status()
        response = resp.json()['data'][pluralize(element)]
        response_elements.extend(response['nodes'])
        has_next_page = response['pageInfo']['hasNextPage']
        after_cursor = response['pageInfo']['endCursor']

    return response_elements



def _construct_get_gene_payload():
    return """
        query gene($id: Int!) {
            gene(id: $id) {
                id
                name
                description
                entrez_id: entrezId
                aliases: geneAliases
                sources {
                    id
                    name
                    title
                    citation
                    citation_id: citationId
                    source_type: sourceType
                    abstract
                    asco_abstract_id: ascoAbstractId
                    author_string: authorString
                    full_journal_title: fullJournalTitle
                    journal
                    pmc_id: pmcId
                    publication_date: publicationDate
                    source_url: sourceUrl
                    clinical_trials: clinicalTrials {
                        id
                        name
                        description
                        nctId
                        url
                    }
                }
            }
        }"""


def _construct_get_all_genes_payload():
    return """
        query genes($after: String) {
            genes(after: $after) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                name
                description
                entrez_id: entrezId
                aliases: geneAliases
                sources {
                    id
                    name
                    title
                    citation
                    citation_id: citationId
                    source_type: sourceType
                    abstract
                    asco_abstract_id: ascoAbstractId
                    author_string: authorString
                    full_journal_title: fullJournalTitle
                    journal
                    pmc_id: pmcId
                    publication_date: publicationDate
                    source_url: sourceUrl
                    clinical_trials: clinicalTrials {
                        id
                        name
                        description
                        nctId
                        url
                    }
                }
              }
            }
        }"""

def _construct_get_variant_payload():
    return """
        query variant($id: Int!) {
            variant(id: $id) {
                id
                name
                allele_registry_id: alleleRegistryId
                civic_actionability_score: evidenceScore
                description
                gene {
                    id
                    name
                    entrezId
                }
                clinvar_entries: clinvarIds
                hgvs_expressions: hgvsDescriptions
                variant_aliases: variantAliases
                variant_types: variantTypes {
                    id
                    name
                    so_id: soid
                    description
                    url
                }
                sources {
                    id
                    name
                    title
                    citation
                    citation_id: citationId
                    source_type: sourceType
                    abstract
                    asco_abstract_id: ascoAbstractId
                    author_string: authorString
                    full_journal_title: fullJournalTitle
                    journal
                    pmc_id: pmcId
                    publication_date: publicationDate
                    source_url: sourceUrl
                    clinical_trials: clinicalTrials {
                        id
                        name
                        description
                        nctId
                        url
                    }
                }
                variantBases
                referenceBases
                referenceBuild
                primaryCoordinates {
                    chromosome
                    representativeTranscript
                    start
                    stop
                }
                secondaryCoordinates {
                    chromosome
                    representativeTranscript
                    start
                    stop
                }
                ensemblVersion
            }
        }"""


def _construct_get_all_variants_payload():
    return """
        query variants($after: String) {
            variants(after: $after, evidenceStatusFilter: ALL) {
                totalCount
                pageInfo {
                  hasNextPage
                  endCursor
                }
                nodes {
                    id
                    name
                    allele_registry_id: alleleRegistryId
                    civic_actionability_score: evidenceScore
                    description
                    gene {
                      id
                      name
                      entrezId
                    }
                    clinvar_entries: clinvarIds
                    hgvs_expressions: hgvsDescriptions
                    variant_aliases: variantAliases
                    variant_types: variantTypes {
                        id
                        name
                        so_id: soid
                        description
                        url
                    }
                    sources {
                        id
                        name
                        title
                        citation
                        citation_id: citationId
                        source_type: sourceType
                        abstract
                        asco_abstract_id: ascoAbstractId
                        author_string: authorString
                        full_journal_title: fullJournalTitle
                        journal
                        pmc_id: pmcId
                        publication_date: publicationDate
                        source_url: sourceUrl
                        clinical_trials: clinicalTrials {
                            id
                            name
                            description
                            nctId
                            url
                        }
                    }
                    variantBases
                    referenceBases
                    referenceBuild
                    primaryCoordinates {
                        chromosome
                        representativeTranscript
                        start
                        stop
                    }
                    secondaryCoordinates {
                        chromosome
                        representativeTranscript
                        start
                        stop
                    }
                    ensemblVersion
                }
            }
        }"""

def _construct_get_evidence_payload():
    return """
        query evidenceItem($id: Int!) {
            evidence: evidenceItem(id: $id) {
                id
                name
                clinical_significance: clinicalSignificance
                description
                drug_interaction_type: drugInteractionType
                evidence_direction: evidenceDirection
                evidence_level: evidenceLevel
                evidence_type: evidenceType
                status
                variant_origin: variantOrigin
                gene {
                  id
                }
                variant {
                  id
                }
                disease {
                  id
                  name
                  display_name: displayName
                  doid
                  disease_url: diseaseUrl
                  aliases: diseaseAliases
                }
                drugs {
                  id
                  name
                  ncit_id: ncitId
                  drug_url: drugUrl
                  aliases: drugAliases
                }
                phenotypes {
                  id
                  name
                  hpo_id: hpoId
                  url
                }
                assertions {
                  id
                }
                source {
                    id
                    name
                    title
                    citation
                    citation_id: citationId
                    source_type: sourceType
                    abstract
                    asco_abstract_id: ascoAbstractId
                    author_string: authorString
                    full_journal_title: fullJournalTitle
                    journal
                    pmc_id: pmcId
                    publication_date: publicationDate
                    source_url: sourceUrl
                    clinical_trials: clinicalTrials {
                        id
                        name
                        description
                        nctId
                        url
                    }
                }
                rating: evidenceRating
            }
        }"""



def _construct_get_all_evidence_payload():
    return """
        query evidenceItems($after: String) {
            evidence_items: evidenceItems(after: $after, status: ALL) {
                totalCount
                pageInfo {
                  hasNextPage
                  endCursor
                }
                nodes {
                    id
                    name
                    clinical_significance: clinicalSignificance
                    description
                    drug_interaction_type: drugInteractionType
                    evidence_direction: evidenceDirection
                    evidence_level: evidenceLevel
                    evidence_type: evidenceType
                    status
                    variant_origin: variantOrigin
                    gene {
                      id
                    }
                    variant {
                      id
                    }
                    disease {
                      id
                      name
                      display_name: displayName
                      doid
                      disease_url: diseaseUrl
                      aliases: diseaseAliases
                    }
                    drugs {
                      id
                      name
                      ncit_id: ncitId
                      drug_url: drugUrl
                      aliases: drugAliases
                    }
                    phenotypes {
                      id
                      name
                      hpo_id: hpoId
                      url
                    }
                    assertions {
                      id
                    }
                    source {
                        id
                        name
                        title
                        citation
                        citation_id: citationId
                        source_type: sourceType
                        abstract
                        asco_abstract_id: ascoAbstractId
                        author_string: authorString
                        full_journal_title: fullJournalTitle
                        journal
                        pmc_id: pmcId
                        publication_date: publicationDate
                        source_url: sourceUrl
                        clinical_trials: clinicalTrials {
                            id
                            name
                            description
                            nctId
                            url
                        }
                    }
                    rating: evidenceRating
                }
            }
        }"""



def _construct_get_assertion_payload():
    return """
        query assertion($id: Int!) {
            assertion(id: $id) {
                id
                name
                amp_level: ampLevel
                clinical_significance: clinicalSignificance
                description
                drug_interaction_type: drugInteractionType
                evidence_direction: assertionDirection
                evidence_type: assertionType
                fda_companion_test: fdaCompanionTest
                fda_regulatory_approval: regulatoryApproval
                name
                nccn_guideline: nccnGuideline {
                  name
                }
                nccn_guideline_version: nccnGuidelineVersion
                status
                summary
                variant_origin: variantOrigin
                gene {
                  id
                }
                variant {
                  id
                }
                acmg_codes: acmgCodes {
                  id
                  code
                  description
                }
                clingen_codes: clingenCodes {
                    id
                    code
                    description
                }
                disease {
                  id
                  name
                  display_name: displayName
                  doid
                  disease_url: diseaseUrl
                  aliases: diseaseAliases
                }
                drugs {
                  id
                  name
                  ncit_id: ncitId
                  drug_url: drugUrl
                  aliases: drugAliases
                }
                evidenceItems {
                  id
                }
                phenotypes {
                  id
                  name
                  hpo_id: hpoId
                  url
                }
            }
        }"""


def _construct_get_all_assertions_payload():
    return """
        query assertions($after: String) {
            assertions(after: $after, status: ALL) {
                totalCount
                pageInfo {
                  hasNextPage
                  endCursor
                }
                nodes {
                    id
                    name
                    amp_level: ampLevel
                    clinical_significance: clinicalSignificance
                    description
                    drug_interaction_type: drugInteractionType
                    evidence_direction: assertionDirection
                    evidence_type: assertionType
                    fda_companion_test: fdaCompanionTest
                    fda_regulatory_approval: regulatoryApproval
                    name
                    nccn_guideline: nccnGuideline {
                      name
                    }
                    nccn_guideline_version: nccnGuidelineVersion
                    status
                    summary
                    variant_origin: variantOrigin
                    gene {
                      id
                    }
                    variant {
                      id
                    }
                    acmg_codes: acmgCodes {
                      id
                      code
                      description
                    }
                    clingen_codes: clingenCodes {
                        id
                        code
                        description
                    }
                    disease {
                      id
                      name
                      display_name: displayName
                      doid
                      disease_url: diseaseUrl
                      aliases: diseaseAliases
                    }
                    drugs {
                      id
                      name
                      ncit_id: ncitId
                      drug_url: drugUrl
                      aliases: drugAliases
                    }
                    evidenceItems {
                      id
                    }
                    phenotypes {
                      id
                      name
                      hpo_id: hpoId
                      url
                    }
                }
            }
        }"""


def _construct_get_variant_group_payload():
    return """
        query variantGroup($id: Int!) {
            variant_group: variantGroup(id: $id) {
                id
                name
                description
                variants(first: 100) {
                  nodes {
                    id
                  }
                }
                sources {
                    id
                    name
                    title
                    citation
                    citation_id: citationId
                    source_type: sourceType
                    abstract
                    asco_abstract_id: ascoAbstractId
                    author_string: authorString
                    full_journal_title: fullJournalTitle
                    journal
                    pmc_id: pmcId
                    publication_date: publicationDate
                    source_url: sourceUrl
                    clinical_trials: clinicalTrials {
                        id
                        name
                        description
                        nctId
                        url
                    }
                }
            }
        }"""

def _construct_get_all_variant_groups_payload():
    return """
        query variantGroups($after: String) {
            variant_groups: variantGroups(after: $after) {
              totalCount
              pageInfo {
                hasNextPage
                endCursor
              }
              nodes {
                id
                name
                description
                variants(first: 100) {
                  nodes {
                    id
                  }
                }
                sources {
                    id
                    name
                    title
                    citation
                    citation_id: citationId
                    source_type: sourceType
                    abstract
                    asco_abstract_id: ascoAbstractId
                    author_string: authorString
                    full_journal_title: fullJournalTitle
                    journal
                    pmc_id: pmcId
                    publication_date: publicationDate
                    source_url: sourceUrl
                    clinical_trials: clinicalTrials {
                        id
                        name
                        description
                        nctId
                        url
                    }
                }
              }
            }
        }"""


def _construct_get_all_payload(page):
    queries = [
        {
            'field': 'id',
            'condition': {
                'name': 'is_greater_than',
                'parameters': [
                    -1
                ]
            }
        }
    ]
    payload = {
        'operator': 'OR',
        'queries': queries,
        'page': page,
        'count': 500,
    }
    return payload


def _construct_query_payload(id_list):
    queries = list()
    for element_id in id_list:
        query = {
            'field': 'id',
            'condition': {
                'name': 'is_equal_to',
                'parameters': [
                    element_id
                ]
            }
        }
        queries.append(query)
    payload = {
        'operator': 'OR',
        'queries': queries
    }
    return payload


def get_evidence_by_ids(evidence_id_list):
    logging.info('Getting evidence...')
    evidence = _get_elements_by_ids('evidence', evidence_id_list)
    logging.info('Caching variant details...')
    variant_ids = [x.variant.id for x in evidence]    # Add variants to cache
    _get_elements_by_ids('variant', variant_ids)
    logging.info('Caching gene details...')
    gene_ids = [x.gene.id for x in evidence]          # Add genes to cache
    _get_elements_by_ids('gene', gene_ids)
    for e in evidence:                        # Load from cache
        e.variant.update()
        e.gene.update()
    return evidence


def get_evidence_by_id(evidence_id):
    return get_evidence_by_ids([evidence_id])[0]


def get_assertions_by_ids(assertion_id_list=[], get_all=False):
    logging.info('Getting assertions...')
    assertions = _get_elements_by_ids('assertion', assertion_id_list, get_all=get_all)
    logging.info('Caching variant details...')
    variant_ids = [x.variant.id for x in assertions]    # Add variants to cache
    _get_elements_by_ids('variant', variant_ids)
    logging.info('Caching gene details...')
    gene_ids = [x.gene.id for x in assertions]          # Add genes to cache
    _get_elements_by_ids('gene', gene_ids)
    for assertion in assertions:                        # Load from cache
        assertion.variant.update()
        assertion.gene.update()
    return assertions


def get_assertion_by_id(assertion_id):
    return get_assertions_by_ids([assertion_id])[0]


def get_all_assertions(include_status=['accepted','submitted','rejected'], allow_cached=True):
    assertions = _get_elements_by_ids('assertion', allow_cached=allow_cached, get_all=True)
    return [a for a in assertions if a.status in include_status]


def search_assertions_by_coordinates(coordinates, search_mode='any'):
    variants = search_variants_by_coordinates(coordinates, search_mode=search_mode)
    assertions = set()
    for v in variants:
        if v.assertions:
            assertions.update(v.assertions)
    return list(assertions)


def get_variants_by_ids(variant_id_list):
    logging.info('Getting variants...')
    variants = _get_elements_by_ids('variant', variant_id_list)
    gene_ids = set()
    for variant in variants:
        gene_ids.add(variant.gene_id)
    if gene_ids:
        logging.info('Caching gene details...')
        _get_elements_by_ids('gene', gene_ids)
    return variants


def get_variant_by_id(variant_id):
    return get_variants_by_ids([variant_id])[0]


def get_variant_groups_by_ids(variant_group_id_list):
    logging.info('Getting variant groups...')
    return _get_elements_by_ids('variant_group', variant_group_id_list)


def get_variant_group_by_id(variant_group_id):
    return get_variant_groups_by_ids([variant_group_id])[0]


def _build_coordinate_table(variants):
    variant_records = list()
    for v in variants:
        c = v.coordinates
        start = getattr(c, 'start', None)
        stop = getattr(c, 'stop', None)
        chr = getattr(c, 'chromosome', None)
        alt = getattr(c, 'variant_bases', None)
        ref = getattr(c, 'reference_bases', None)
        if all([start, stop, chr]):
            variant_records.append([chr, start, stop, alt, ref, hash(v)])
        else:
            continue
        start = getattr(c, 'start2', None)
        stop = getattr(c, 'stop2', None)
        chr = getattr(c, 'chromosome2', None)
        if all([start, stop, chr]):
            variant_records.append([chr, start, stop, None, None, hash(v)])
    df = pd.DataFrame.from_records(
        variant_records,
        columns=['chr', 'start', 'stop', 'alt', 'ref', 'v_hash']
    ).sort_values(by=['chr', 'start', 'stop', 'alt', 'ref'])
    MODULE.COORDINATE_TABLE = df
    MODULE.COORDINATE_TABLE_START = df.start.sort_values()
    MODULE.COORDINATE_TABLE_STOP = df.stop.sort_values()
    MODULE.COORDINATE_TABLE_CHR = df.chr.sort_values()


def get_all_variants(include_status=['accepted','submitted','rejected'], allow_cached=True):
    variants = _get_elements_by_ids('variant', allow_cached=allow_cached, get_all=True)
    if include_status:
        assert CACHE.get('evidence_items_all_ids', False)
        resp = list()
        for v in variants:
            v._include_status = include_status
            if v.evidence:
                resp.append(v)
        return resp
    else:
        return variants


def get_all_variant_groups(allow_cached=True):
    variant_groups = _get_elements_by_ids('variant_group', allow_cached=allow_cached, get_all=True)
    return variant_groups

def search_variants_by_allele_registry_id(caid):
    """
    Search the cache for variants matching the queried Allele Registry ID (CAID)

    :param String caid: Allele Registry ID to query

    :return: Returns a list of variant hashes matching the Allele Registry ID
    """
    return search_variants_by_attribute('allele_registry_id', caid)

def search_variants_by_name(name):
    """
    Search the cache for variants matching the queried name

    :param String name: Variant name to query

    :return: Returns a list of variant hashes matching the name
    """
    return search_variants_by_attribute('name', name)

def search_variants_by_hgvs(hgvs):
    """
    Search the cache for variants matching the queried HGVS expression

    :param String name: HGVS expression to query

    :return: Returns a list of variant hashes matching the HGVS expression
    """
    return search_variants_by_list_field('hgvs_expressions', hgvs)

def search_variants_by_attribute(attribute, value):
    variants = get_all_variants()
    return [v for v in variants if getattr(v, attribute) == value]

def search_variants_by_list_field(field, value):
    variants = get_all_variants()
    matched_variants = []
    return [v for v in variants if value in getattr(v, field)]

def search_variants_by_coordinates(coordinate_query, search_mode='any'):
    """
    Search the cache for variants matching provided coordinates using the corresponding search mode.

    :param CoordinateQuery coordinate_query: Coordinates to query

    :param any,query_encompassing,variant_encompassing,exact search_mode:
                *any* : any overlap between a query and a variant is a match\n
                *query_encompassing* : CIViC variant records must fit within the coordinates of the query\n
                *record_encompassing* : CIViC variant records must encompass the coordinates of the query\n
                *exact* : variants must match coordinates precisely, as well as reference allele(s) and alternate allele(s).
                Use ``'*'`` in the coordinate_query as a wildcard for reference and/or alternate alleles.
                Using ``None`` in the coordinate_query for reference or alternate alleles will only match
                variants that have no reference or alternate alleles, respectively (e.g. indels) \n
                search_mode is *any* by default

    :return:    Returns a list of variant hashes matching the coordinates and search_mode
    """
    get_all_variants()
    if coordinate_query.build == 'GRCh37':
        ct = COORDINATE_TABLE
        start_idx = COORDINATE_TABLE_START
        stop_idx = COORDINATE_TABLE_STOP
        chr_idx = COORDINATE_TABLE_CHR
        start = int(coordinate_query.start)
        stop = int(coordinate_query.stop)
        chromosome = str(coordinate_query.chr)
        # overlapping = (start <= ct.stop) & (stop >= ct.start)
        left_idx = chr_idx.searchsorted(chromosome)
        right_idx = chr_idx.searchsorted(chromosome, side='right')
        chr_ct_idx = chr_idx[left_idx:right_idx].index
        right_idx = start_idx.searchsorted(stop, side='right')
        start_ct_idx = start_idx[:right_idx].index
        left_idx = stop_idx.searchsorted(start)
        stop_ct_idx = stop_idx[left_idx:].index
        match_idx = chr_ct_idx & start_ct_idx & stop_ct_idx
        m_df = ct.loc[match_idx, ]
        if search_mode == 'any':
            var_digests = m_df.v_hash.to_list()
            return [CACHE[v] for v in var_digests]
        elif search_mode == 'query_encompassing':
            match_idx = (start <= m_df.start) & (stop >= m_df.stop)
        elif search_mode == 'variant_encompassing':
            match_idx = (start >= m_df.start) & (stop <= m_df.stop)
        elif search_mode == 'exact':
            match_idx = (start == m_df.start) & (stop == m_df.stop)
            if coordinate_query.alt is not None and coordinate_query.alt != '*':
                if coordinate_query.alt == '-':
                    raise ValueError("Unexpected alt `-` in coordinate query. Did you mean `None`?")
                match_idx = match_idx & (coordinate_query.alt == m_df.alt)
            elif coordinate_query.alt is None:
                match_idx = match_idx & pd.isnull(m_df.alt)
            if (coordinate_query.ref is not None and coordinate_query.ref != '*'):
                if coordinate_query.ref == '-':
                    raise ValueError("Unexpected ref `-` in coordinate query. Did you mean `None`?")
                match_idx = match_idx & (coordinate_query.ref == m_df.ref)
            elif coordinate_query.ref is None:
                match_idx = match_idx & pd.isnull(m_df.ref)
        else:
            raise ValueError("unexpected search mode")
        var_digests = m_df.loc[match_idx,].v_hash.to_list()
        return [CACHE[v] for v in var_digests]
    else:
        if search_mode == 'exact':
            if coordinate_query.alt or coordinate_query.ref:
                if coordinate_query.alt == '*' or coordinate_query.ref == '*':
                    raise ValueError("Can't use wildcard when searching for non-GRCh37 coordinates")
                if coordinate_query.alt == '-':
                    raise ValueError("Unexpected alt `-` in coordinate query. Did you mean `None`?")
                if coordinate_query.ref == '-':
                    raise ValueError("Unexpected ref `-` in coordinate query. Did you mean `None`?")
                hgvs = _construct_hgvs_for_coordinate_query(coordinate_query)
                if hgvs is not None:
                    s = requests.Session()
                    retry = Retry(
                        total=5,
                        read=5,
                        connect=5,
                        backoff_factor=0.3,
                        status_forcelist=(500, 502, 504),
                    )
                    adapter = requests.adapters.HTTPAdapter(max_retries=retry)
                    s.mount('http://', adapter)
                    r = s.get(url=_allele_registry_url(), params={'hgvs': hgvs})
                    data = r.json()
                    if '@id' in data:
                        allele_registry_id = data['@id'].split('/')[-1]
                        if not allele_registry_id == '_:CA':
                            return search_variants_by_allele_registry_id(allele_registry_id)
            else:
                raise ValueError("alt or ref required for non-GRCh37 coordinate queries")
        else:
            raise ValueError("Only exact search mode is supported for non-GRCh37 coordinate queries")

def _allele_registry_url():
    return "http://reg.genome.network/allele"

def _construct_hgvs_for_coordinate_query(coordinate_query):
    if coordinate_query.build == 'GRCh38':
        chromosome = _refseq_sequence_b38(coordinate_query.chr)
    elif coordinate_query.build == 'NCBI36':
        chromosome = _refseq_sequence_b36(coordinate_query.chr)
    else:
        raise ValueError("unexpected reference build")
    if chromosome is None:
        return None
    base_hgvs = "{}:g.{}".format(chromosome, coordinate_query.start)
    variant_type = _variant_type(coordinate_query)
    if variant_type == "deletion":
        if len(coordinate_query.ref) > 1:
            return"{}_{}del".format(base_hgvs, coordinate_query.stop)
        else:
            return "{}del".format(base_hgvs)
    elif variant_type == "substitution":
        return "{}{}>{}".format(base_hgvs, coordinate_query.ref, coordinate_query.alt)
    elif variant_type == "insertion":
        return "{}_{}ins{}".format(base_hgvs, coordinate_query.stop, coordinate_query.alt)
    elif variant_type == "indel":
        if len(coordinate_query.ref) > 1:
          return "{}_{}delins{}".format(base_hgvs, coordinate_query.stop, coordinate_query.alt)
        else:
          return "{}delins{}".format(base_hgvs, coordinate_query.alt)
    else:
        return None

def _variant_type(coordinate_query):
    if not coordinate_query.ref and not coordinate_query.alt:
        return None
    elif coordinate_query.ref and not coordinate_query.alt:
        return "deletion"
    elif not coordinate_query.ref and coordinate_query.alt:
        return "insertion"
    elif len(coordinate_query.ref) == 1 and len(coordinate_query.alt) == 1:
        return "substitution"
    elif len(coordinate_query.ref) > 1 and len(coordinate_query.alt) > 1:
        return "indel"
    else:
        return None

def _refseq_sequence_b36(chromosome):
    chromosome = chromosome.replace('chr', '')
    sequences = {
      '1' : 'NC_000001.9',
      '2' : 'NC_000002.10',
      '3' : 'NC_000003.10',
      '4' : 'NC_000004.10',
      '5' : 'NC_000005.8',
      '6' : 'NC_000006.10',
      '7' : 'NC_000007.12',
      '8' : 'NC_000008.9',
      '9' : 'NC_000009.10',
      '10' : 'NC_000010.9',
      '11' : 'NC_000011.8',
      '12' : 'NC_000012.10',
      '13' : 'NC_000013.9',
      '14' : 'NC_000014.7',
      '15' : 'NC_000015.8',
      '16' : 'NC_000016.8',
      '17' : 'NC_000017.9',
      '18' : 'NC_000018.8',
      '19' : 'NC_000019.8',
      '20' : 'NC_000020.9',
      '21' : 'NC_000021.7',
      '22' : 'NC_000022.9',
      'X' : 'NC_000023.9',
      'Y' : 'NC_000024.8',
    }
    if chromosome not in sequences:
        return None
    return sequences[chromosome]

def _refseq_sequence_b38(chromosome):
    chromosome = chromosome.replace('chr', '')
    sequences = {
      '1' : 'NC_000001.11',
      '2' : 'NC_000002.12',
      '3' : 'NC_000003.12',
      '4' : 'NC_000004.12',
      '5' : 'NC_000005.10',
      '6' : 'NC_000006.12',
      '7' : 'NC_000007.14',
      '8' : 'NC_000008.11',
      '9' : 'NC_000009.12',
      '10' : 'NC_000010.11',
      '11' : 'NC_000011.10',
      '12' : 'NC_000012.12',
      '13' : 'NC_000013.11',
      '14' : 'NC_000014.9',
      '15' : 'NC_000015.10',
      '16' : 'NC_000016.10',
      '17' : 'NC_000017.11',
      '18' : 'NC_000018.10',
      '19' : 'NC_000019.10',
      '20' : 'NC_000020.11',
      '21' : 'NC_000021.9',
      '22' : 'NC_000022.11',
      'X' : 'NC_000023.11',
      'Y' : 'NC_000024.10',
    }
    if chromosome not in sequences:
        return None
    return sequences[chromosome]

# TODO: Refactor this method
def bulk_search_variants_by_coordinates(sorted_queries, search_mode='any'):
    """
    An interator to search the cache for variants matching the set of sorted coordinates and yield
    matches corresponding to the search mode.

    :param list[CoordinateQuery] sorted_queries: Sorted list of coordinates to query

    :param any,query_encompassing,variant_encompassing,exact search_mode:
                *any* : any overlap between a query and a variant is a match\n
                *query_encompassing* : CIViC variant records must fit within the coordinates of the query\n
                *record_encompassing* : CIViC variant records must encompass the coordinates of the query\n
                *exact* : variants must match coordinates precisely, as well as reference allele(s) and alternate allele(s).
                Use ``'*'`` in the coordinate_query as a wildcard for reference and/or alternate alleles.
                Using ``None`` in the coordinate_query for reference or alternate alleles will only match
                variants that have no reference or alternate alleles, respectively (e.g. indels) \n
                search_mode is *any* by default

    :return:    returns a dictionary of Match lists, keyed by query
    """

    def is_sorted(prev_q, current_q):
        if prev_q['chr'] < current_q['chr']:
            return True
        if prev_q['chr'] > current_q['chr']:
            return False
        if prev_q['start'] < current_q['start']:
            return True
        if prev_q['start'] > current_q['start']:
            return False
        if prev_q['stop'] < current_q['stop']:
            return True
        if prev_q['stop'] > current_q['stop']:
            return False
        return True

    ct_pointer = 0
    query_pointer = 0
    last_query_pointer = -1
    match_start = None
    ct = MODULE.COORDINATE_TABLE
    matches = defaultdict(list)
    Match = namedtuple('Match', ct.columns)

    def append_match(matches_list, query, ct_row):
        matches_list[query].append(Match(**ct_row.to_dict()))

    while query_pointer < len(sorted_queries) and ct_pointer < len(ct):
        if last_query_pointer != query_pointer:
            q = sorted_queries[query_pointer]
            if q.build != 'GRCh37':
                raise ValueError("Bulk coordinate search only supports build GRCh37")
            if match_start is not None:
                ct_pointer = match_start
                match_start = None
            last_query_pointer = query_pointer
        c = ct.iloc[ct_pointer]
        q_chr = str(q.chr)
        c_chr = c.chr
        if q_chr < c_chr:
            query_pointer += 1
            continue
        if q_chr > c_chr:
            ct_pointer += 1
            continue
        q_start = int(q.start)
        c_start = c.start
        q_stop = int(q.stop)
        c_stop = c.stop
        if q_start > c_stop:
            ct_pointer += 1
            continue
        if q_stop < c_start:
            query_pointer += 1
            continue
        if search_mode == 'any':
            append_match(matches, q, c)
        elif search_mode == 'exact' and q_start == c_start and q_stop == c_stop:
            q_alt = q.alt
            c_alt = c.alt
            q_ref = q.ref
            c_ref = c.ref
            if q_alt == '-':
                raise ValueError("Unexpected alt `-` in coordinate query. Did you mean `None`?")
            if q_ref == '-':
                raise ValueError("Unexpected ref `-` in coordinate query. Did you mean `None`?")
            if (not (q_alt != '*' and q_alt != c_alt)) and (not (q_ref != '*' and q_ref != c_ref)):
                append_match(matches, q, c)
        elif search_mode == 'query_encompassing' and q_start <= c_start and q_stop >= c_stop:
            append_match(matches, q, c)
        elif search_mode == 'record_encompassing' and c_start <= q_start and c_stop >= q_stop:
            append_match(matches, q, c)
        if match_start is None:
            match_start = ct_pointer
        ct_pointer += 1
    return dict(matches)


def get_genes_by_ids(gene_id_list):
    """
    :param list gene_id_list: A list of CIViC gene IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Gene` objects.
    """
    logging.info('Getting genes...')
    genes = _get_elements_by_ids('gene', gene_id_list)  # Advanced search results are incomplete
    variant_ids = set()
    for gene in genes:
        for variant in gene.variants:
            variant_ids.add(variant.id)
    if variant_ids:
        logging.info('Caching variant details...')
        _get_elements_by_ids('variant', variant_ids)
    for gene in genes:
        for variant in gene.variants:
            variant.update()
    return genes


def get_gene_by_id(gene_id):
    """
    :param int gene_id: A single CIViC gene ID.
    :returns: A :class:`Gene` object.
    """
    return get_genes_by_ids([gene_id])[0]


def get_all_genes(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all genes.

    :param list include_status: A list of statuses. Only genes and their associated entities matching the given statuses will be returned.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Gene` objects.
    """
    genes = _get_elements_by_ids('gene', get_all=True, allow_cached=allow_cached)
    if include_status:
        assert CACHE.get('variants_all_ids', False)
        assert CACHE.get('evidence_items_all_ids', False)
        resp = list()
        for g in genes:
            g._include_status = include_status
            if g.variants:
                resp.append(g)
        return resp
    else:
        return genes


def get_all_evidence(include_status=['accepted','submitted','rejected'], allow_cached=True):
    evidence = _get_elements_by_ids('evidence', get_all=True, allow_cached=allow_cached)
    return [e for e in evidence if e.status in include_status]


def get_HPO_terms_by_ids(hpo_id_list):
    if not HPO_TERMS:
        _load_HPO()
    return [HPO_TERMS[x] for x in hpo_id_list]


def _load_HPO():
    url = 'https://civicdb.org/api/phenotypes?count=100000'
    resp = requests.get(url)
    resp.raise_for_status()
    for h in resp.json():
        HPO_TERMS[h['id']] = h

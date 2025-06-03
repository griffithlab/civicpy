import requests
from requests.packages.urllib3.util.retry import Retry
import importlib
import logging
import pandas as pd
import pickle
import os
from pathlib import Path
from collections import defaultdict, namedtuple
import requests
import deprecation
from datetime import datetime, timedelta
from backports.datetime_fromisoformat import MonkeyPatch
MonkeyPatch.patch_fromisoformat()
import re

from civicpy import REMOTE_CACHE_URL, LOCAL_CACHE_PATH, CACHE_TIMEOUT_DAYS
from civicpy.__version__ import __version__
from civicpy import exports
from civicpy import graphql_payloads
from civicpy import utils


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


CIVIC_TO_PYCLASS = {
    'evidence_items': 'evidence',
    'five_prime_coordinates': 'coordinate',
    'three_prime_coordinates': 'coordinate',
    'five_prime_start_exon_coordinates': 'exon_coordinate',
    'five_prime_end_exon_coordinates': 'exon_coordinate',
    'three_prime_start_exon_coordinates': 'exon_coordinate',
    'three_prime_end_exon_coordinates': 'exon_coordinate',
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


def element_lookup_by_id(element_type, element_id):
    e = _request_by_ids(element_type, [int(element_id)])[0]
    e = _postprocess_response_element(e, element_type)
    return e


def get_class(element_type):
    e_string = utils.singularize(element_type)
    class_string = utils.snake_to_camel(e_string)
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
    variants_with_coords = set()
    for k, v in loaded_cache.items():
        if isinstance(k, str):
            c[k] = v
        elif isinstance(k, int):
            c[hash(v)] = v
            if isinstance(v, GeneVariant) or isinstance(v, FusionVariant):
                variants_with_coords.add(v)
        else:
            raise ValueError
    old_cache = MODULE.CACHE
    MODULE.CACHE = c
    for k, v in MODULE.CACHE.items():
        if isinstance(k, str):
            continue
        v.update()
    if _has_full_cached_fresh() or on_stale == 'ignore':
        _build_coordinate_table(variants_with_coords)
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
        molecular_profiles = _get_elements_by_ids('molecular_profile', allow_cached=False, get_all=True)
        genes = _get_elements_by_ids('gene', allow_cached=False, get_all=True)
        factors = _get_elements_by_ids('factor', allow_cached=False, get_all=True)
        fusions = _get_elements_by_ids('fusion', allow_cached=False, get_all=True)
        variants = _get_elements_by_ids('variant', allow_cached=False, get_all=True)
        evidence = _get_elements_by_ids('evidence', allow_cached=False, get_all=True)
        assertions = _get_elements_by_ids('assertion', allow_cached=False, get_all=True)
        variant_groups = _get_elements_by_ids('variant_group', allow_cached=False, get_all=True)
        sources = _get_elements_by_ids('source', allow_cached=False, get_all=True)
        diseases = _get_elements_by_ids('disease', allow_cached=False, get_all=True)
        therapies = _get_elements_by_ids('therapy', allow_cached=False, get_all=True)
        phenotypes = _get_elements_by_ids('phenotype', allow_cached=False, get_all=True)
        organizations = _get_elements_by_ids('organization', allow_cached=False, get_all=True)
        endorsements = _get_elements_by_ids('endorsement', allow_cached=False, get_all=True)
        for e in evidence:
            e.assertions = [a for a in assertions if a.id in e.assertion_ids]
            e.therapies = [t for t in therapies if t.id in e.therapy_ids]
            e._partial = False
            CACHE[hash(e)] = e
        for g in genes:
            g.sources = [s for s in sources if s.id in g.source_ids]
            g.variants = [v for v in variants if v.feature_id == g.id]
            g._partial = False
            CACHE[hash(g)] = g
        for f in factors:
            f.sources = [s for s in sources if s.id in f.source_ids]
            f.variants = [v for v in variants if v.feature_id == f.id]
            f._partial = False
            CACHE[hash(f)] = f
        for f in fusions:
            f.sources = [s for s in sources if s.id in f.source_ids]
            f.variants = [v for v in variants if v.feature_id == f.id]
            f._partial = False
            CACHE[hash(f)] = f
        for v in variants:
            v.variant_groups = [vg for vg in variant_groups if v.id in vg.variant_ids]
            v.molecular_profiles = [mp for mp in molecular_profiles if v.id in mp.variant_ids]
            v._partial = False
            CACHE[hash(v)] = v
        for a in assertions:
            a.evidence_items = [e for e in evidence if e.id in a.evidence_ids]
            a.therapies = [t for t in therapies if t.id in a.therapy_ids]
            a.endorsements = [e for e in endorsements if e.id in a.endorsement_ids]
            a._partial = False
            CACHE[hash(a)] = a
        for vg in variant_groups:
            vg.sources = [s for s in sources if s.id in vg.source_ids]
            vg.variants = [v for v in variants if v.id in vg.variant_ids]
            vg._partial = False
            CACHE[hash(vg)] = vg
        for mp in molecular_profiles:
            mp.sources = [s for s in sources if s.id in mp.source_ids]
            mp.evidence_items = [e for e in evidence if e.molecular_profile_id == mp.id]
            mp.variants = [v for v in variants if v.id in mp.variant_ids]
            mp.assertions = [a for a in assertions if a.molecular_profile_id == mp.id]
            updated_parsed_name = []
            for pn in mp.parsed_name:
                if pn.type == 'Feature':
                    if pn.featureType == 'GENE':
                        pn = [g for g in genes if g.id == pn.id][0]
                    elif pn.featureType == 'FACTOR':
                        pn = [f for f in factors if f.id == pn.id][0]
                    elif pn.featureType == 'FUSION':
                        pn = [f for f in fusions if f.id == pn.id][0]
                elif pn.type == 'Variant':
                    pn = [v for v in variants if v.id == pn.id][0]
                else:
                    pn = pn.text
                updated_parsed_name.append(pn)
            mp.parsed_name = updated_parsed_name
            mp._partial = False
            CACHE[hash(mp)] = mp
        for s in sources:
            s.evidence_items = [e for e in evidence if s.id == e.source_id]
            s.genes = [g for g in genes if s.id in g.source_ids]
            s.factors = [f for f in factors if s.id in f.source_ids]
            s.fusions = [f for f in fusions if s.id in f.source_ids]
            s.molecular_profiles = [m for m in molecular_profiles if s.id in m.source_ids]
            s._partial = False
            CACHE[hash(s)] = s
        for d in diseases:
            d.evidence_items = [e for e in evidence if d.id == e.disease_id]
            d.assertions = [a for a in assertions if d.id == a.disease_id]
            d._partial = False
            CACHE[hash(d)] = d
        for t in therapies:
            t.evidence_items = [e for e in evidence if t.id in e.therapy_ids]
            t.assertions = [a for a in assertions if t.id in a.therapy_ids]
            t._partial = False
            CACHE[hash(t)] = t
        for p in phenotypes:
            p.evidence_items = [e for e in evidence if p.id in e.phenotype_ids]
            p.assertions = [a for a in assertions if p.id in a.phenotype_ids]
            p._partial = False
            CACHE[hash(p)] = p
        for o in organizations:
            o.endorsements = [e for e in endorsements if e.organization_id == o.id]
            o._partial = False
            CACHE[hash(o)] = o
        for e in endorsements:
            e._partial = False
            CACHE[hash(e)] = e
        CACHE['full_cached'] = datetime.now()
        _build_coordinate_table(variants)
        save_cache(local_cache_path=local_cache_path)


def _make_local_cache_path_if_missing(local_cache_path):
    p = Path(local_cache_path)
    if not p.parent.is_dir():
        os.makedirs(p.parent)


def _build_coordinate_table(variants):
    variant_records = list()
    for v in variants:
        if isinstance(v, GeneVariant):
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
        if isinstance(v, FusionVariant):
            c = v.five_prime_coordinates
            start = getattr(c, 'start', None)
            stop = getattr(c, 'stop', None)
            chr = getattr(c, 'chromosome', None)
            alt = getattr(c, 'variant_bases', None)
            ref = getattr(c, 'reference_bases', None)
            if all([start, stop, chr]):
                variant_records.append([chr, start, stop, alt, ref, hash(v)])
            else:
                continue
            c = v.three_prime_coordinates
            start = getattr(c, 'start', None)
            stop = getattr(c, 'stop', None)
            chr = getattr(c, 'chromosome', None)
            alt = getattr(c, 'variant_bases', None)
            ref = getattr(c, 'reference_bases', None)
            if all([start, stop, chr]):
                variant_records.append([chr, start, stop, alt, ref, hash(v)])
            else:
                continue
    df = pd.DataFrame.from_records(
        variant_records,
        columns=['chr', 'start', 'stop', 'alt', 'ref', 'v_hash']
    ).sort_values(by=['chr', 'start', 'stop', 'alt', 'ref'])
    MODULE.COORDINATE_TABLE = df
    MODULE.COORDINATE_TABLE_START = df.start.sort_values()
    MODULE.COORDINATE_TABLE_STOP = df.stop.sort_values()
    MODULE.COORDINATE_TABLE_CHR = df.chr.sort_values()


def _is_valid(warnings: list[str], emit_warnings: bool = False) -> bool:
    """Determine whether object is valid

    :param warnings: List of warnings. If warnings exist, then object is invalid
    :param emit_warnings: Whether to log warnings, defaults to False
    :return: ``True`` if object is valid. ``False`` otherwise
    """
    if not warnings:
        return True

    if emit_warnings:
        for warning in warnings:
            logging.warning(warning)
    return False


def _is_valid_for_gks_json(cls, emit_warnings: bool = False) -> bool:
    """Determine whether ``cls`` is able to be represented as GKS model

    :param cls: The instance to validate
    :param emit_warnings: Whether to log warnings, defaults to False
    :raises TypeError: If ``cls`` is not of type Assertion or Evidence
    :return: ``True`` if ``cls`` is able to be represented as GKS model. ``False``
        otherwise
    """
    if not isinstance(cls, (Assertion, Evidence)):
        raise TypeError(f"{type(cls)} is not an instance of `Assertion` or `Evidence`.")

    prefix = f"{cls.__class__.__name__} {cls.id}"
    warnings = []

    if cls.status != "accepted":
        warnings.append(f"{prefix} does not have 'accepted' status. Skipping")

    record_type = cls.evidence_type if isinstance(cls, Evidence) else cls.assertion_type
    if record_type not in ("DIAGNOSTIC", "PREDICTIVE", "PROGNOSTIC"):
        warnings.append(
            f"{prefix} type is not one of: 'DIAGNOSTIC', 'PREDICTIVE', or 'PROGNOSTIC'. Skipping"
        )

    len_mp_variants = len(cls.molecular_profile.variants)
    if len_mp_variants > 1:
        warnings.append(f"{prefix} has a complex molecular profile. Skipping")
    elif len_mp_variants == 1:
        if not isinstance(cls.molecular_profile.variants[0], GeneVariant):
            warnings.append(f"{prefix} variant is not a ``GeneVariant``. Skipping")
    else:
        warnings.append(f"{prefix} has no variants. Skipping")

    return _is_valid(warnings, emit_warnings)


class CivicRecord:
    """
    As a base class, :class:`CivicRecord` is used to define the characteristic of all records in CIViC. This class is not
    intended to be invoked directly by the end user, but provided for documentation of shared methods and variables in
    child classes.
    """

    _SIMPLE_FIELDS = {'id', 'type'}
    _COMPLEX_FIELDS = set()
    _NULLABLE_COMPLEX_FIELDS = set()
    _OPTIONAL_FIELDS = set()

    def __init__(self, partial=False, **kwargs):
        """
        The record object may be initialized by the user, though the practice is discouraged. To do so, values for each
        of the object attributes (except ``type``) must be specified as keyword arguments, or the ``partial`` parameter must
        be set to **True**. If ``partial`` is set to **True**, the ``id`` keyword argument is still required.

        Users are encouraged to use the functions for :ref:`getting_records` in lieu of directly initializing record
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
            cls = get_class(CIVIC_TO_PYCLASS.get(field, field))
            if is_compound:
                result = list()
                for data in v:
                    if isinstance(data, dict):
                        data['type'] = data.get('type', utils.singularize(field))
                        result.append(cls(partial=True, **data))
                    else:
                        result.append(data)
                self.__setattr__(field, result)
            else:
                t = v.get('type', field)
                v['type'] = CIVIC_TO_PYCLASS.get(t, t)
                if v.keys() == {'type'}:
                    if field in self._NULLABLE_COMPLEX_FIELDS:
                        self.__setattr__(field, None)
                    else:
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
        :param bool force: Flag to indicate whether to force an update from the server, even if a full record exists in the cache.
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


class MolecularProfile(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'description',
        'molecular_profile_score',
        'name',
        'variant_ids',
        'source_ids',
    })
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'aliases',
        'assertions',
        'evidence_items',
        'sources',
        'variants',
        'parsed_name',
    })

    def __init__(self, **kwargs):
        self._evidence_items = []
        self._assertions = []
        self._variants = []
        self._sources = []
        super().__init__(**kwargs)

    @property
    def evidence_sources(self):
        """
        A list of :class:`Source` records associated with all the :class:`Evidence` records under this molecular profile.
        """
        sources = set()
        for evidence in self.evidence_items:
            if evidence.source is not None:
                sources.add(evidence.source)
        return sources

    @property
    def summary(self):
        """
        A shorthand for the description.
        """
        return self.description

    #@summary.setter
    #def summary(self, value):
    #    self.description = value

    @property
    def evidence(self):
        """
        A shorthand for evidence_items.
        """
        return self.evidence_items

    @property
    def evidence_items(self):
        """
        A list of :class:`Evidence` records associated with this molecular profile.
        """
        return [e for e in self._evidence_items if e.status in self._include_status]

    @evidence_items.setter
    def evidence_items(self, value):
        self._evidence_items = value

    @property
    def assertions(self):
        """
        A list of :class:`Assertion` records associated with this molecular profile.
        """
        return [a for a in self._assertions if a.status in self._include_status]

    @assertions.setter
    def assertions(self, value):
        self._assertions = value

    @property
    def variants(self):
        """
        A list :class:`Variant` objects involved in this molecular profile.
        """
        return self._variants

    @variants.setter
    def variants(self, value):
        self._variants = value

    @property
    def sources(self):
        """
        A list :class:`Source` objects involved in this molecular profile.
        """
        return self._sources

    @sources.setter
    def sources(self, value):
        self._sources = value

    def sanitized_name(self):
        name = self.name
        words = []
        for word in name.split(' '):
            regex = re.compile(r"^([A-Z]+)([0-9]+)(=)(.*)$")
            match = regex.match(word)
            if match is not None:
                word = "".join([match.group(1), match.group(2), match.group(1), match.group(4)])
            words.append(word)
        return ' '.join(words)


class Variant(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'subtype',
        'feature_id',
        'name',
        'single_variant_molecular_profile_id',
    })
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        # 'errors',
        #'lifecycle_actions',
        # 'provisional_values',
        'variant_aliases',
        'variant_groups',
        'variant_types'
    })

    def __init__(self, **kwargs):
        kwargs['type'] = 'variant'
        self._variant_groups = []
        self._molecular_profiles = []
        super().__init__(**kwargs)

    def __repr__(self):
        return '<CIViC {} ({}) {}>'.format(self.type, self.subtype, self.id)

    @property
    def aliases(self):
        """
        A curated list of aliases by which this variant is references. Shorthand for the ``variant_aliases`` attribute.
        """
        return self.variant_aliases

    @property
    def groups(self):
        """
        A list of :class:`VariantGroup` records to which this variant belongs. Shorthand for the ``variant_groups`` attribute.
        """
        return self.variant_groups

    @property
    def types(self):
        """
        A list of :class:`CivicAttribute` objects describing `variant types`_ from the `Sequence Ontology`_.

        .. _variant types: https://docs.civicdb.org/en/latest/model/variants/types.html
        .. _Sequence Ontology: http://www.sequenceontology.org/
        """
        return self.variant_types

    @property
    def variant_groups(self):
        """
        A list of :class:`VariantGroup` records to which this variant belongs.
        """
        return self._variant_groups

    @variant_groups.setter
    def variant_groups(self, value):
        self._variant_groups = value

    @property
    def molecular_profiles(self):
        """
        A list of :class:`MolecularProfile` records involving this variant.
        """
        for mp in self._molecular_profiles:
            mp._include_status = self._include_status
        return [m for m in self._molecular_profiles if m.evidence_items or m.assertions]

    @molecular_profiles.setter
    def molecular_profiles(self, value):
        self._molecular_profiles = value

    @property
    def single_variant_molecular_profile(self):
        """
        The :class:`MolecularProfile` record representing the single variant on its own.
        """
        mp = _get_element_by_id('molecular_profile', self.single_variant_molecular_profile_id)
        mp._include_status = self._include_status
        return mp

    def is_valid_for_vcf(self, emit_warnings=False):
        return False


class GeneVariant(Variant):
    _SIMPLE_FIELDS = Variant._SIMPLE_FIELDS.union({
        'allele_registry_id',
        'entrez_name',
        'entrez_id'
    })
    _COMPLEX_FIELDS = Variant._COMPLEX_FIELDS.union({
        'clinvar_entries',
        'coordinates',
        'hgvs_expressions',
        #'lifecycle_actions',
        # 'provisional_values',
    })

    @property
    def gene(self):
        """
        The :class:`Gene` record this variant belongs to.
        """
        return _get_element_by_id('gene', self.feature_id)

    @property
    def feature(self):
        """
        The :class:`Gene` feature this variant belongs to.
        """
        return self.gene

    @property
    def is_insertion(self):
        """
        Based on the coordiantes, True if the variant is an insertion, else False.
        """
        ref = self.coordinates.reference_bases
        alt = self.coordinates.variant_bases
        return (ref is None and alt is not None) or (ref is not None and alt is not None and len(ref) < len(alt))

    @property
    def is_deletion(self):
        """
        Based on the coordiantes, True if the variant is a deletion, else False.
        """
        ref = self.coordinates.reference_bases
        alt = self.coordinates.variant_bases
        if alt is not None and (alt == '-' or alt == ''):
            alt = None
        return (ref is not None and alt is None) or (ref is not None and alt is not None and len(ref) > len(alt))

    def is_valid_for_vcf(self, emit_warnings=False):
        warnings = []
        if self.coordinates is None:
            warnings.append("Variant {} has no coordinates. Skipping".format(self.id))
        if self.coordinates.reference_build != 'GRCh37':
            warnings.append("Variant coordinate reference build is not GRCh37 for variant {}. Skipping.".format(self.id))
        if (self.is_insertion or self.is_deletion) and self.coordinates.representative_transcript is None:
            warnings.append("Variant {} is an indel but coordinates are missing a representative transcript. Skipping.".format(self.id))
        if not self._valid_alt_bases():
            warnings.append("Unsupported variant base(s) for variant {}. Skipping.".format(self.id))
        if not self._valid_ref_bases():
            warnings.append("Unsupported reference base(s) for variant {}. Skipping.".format(self.id))
        if not(self.coordinates.chromosome and self.coordinates.start and (self.coordinates.reference_bases or self.coordinates.variant_bases)):
            warnings.append("Incomplete coordinates for variant {}. Skipping.".format(self.id))

        return _is_valid(warnings, emit_warnings)

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



class FactorVariant(Variant):
    _SIMPLE_FIELDS = Variant._SIMPLE_FIELDS.union({
        'ncit_id',
    })

    @property
    def factor(self):
        """
        The :class:`Factor` record this variant belongs to.
        """
        return _get_element_by_id('factor', self.feature_id)

    @property
    def feature(self):
        """
        The :class:`Factor` feature this variant belongs to.
        """
        return self.factor


class FusionVariant(Variant):
    _SIMPLE_FIELDS = Variant._SIMPLE_FIELDS.union({
        'vicc_compliant_name',
    })
    _COMPLEX_FIELDS = Variant._COMPLEX_FIELDS.union({
        'five_prime_coordinates',
        'three_prime_coordinates',
        'five_prime_start_exon_coordinates',
        'five_prime_end_exon_coordinates',
        'three_prime_start_exon_coordinates',
        'three_prime_end_exon_coordinates',
    })
    _NULLABLE_COMPLEX_FIELDS = Variant._NULLABLE_COMPLEX_FIELDS.union({
        'five_prime_coordinates',
        'three_prime_coordinates',
        'five_prime_start_exon_coordinates',
        'five_prime_end_exon_coordinates',
        'three_prime_start_exon_coordinates',
        'three_prime_end_exon_coordinates',
    })

    @property
    def fusion(self):
        """
        The :class:`Fusion` record this variant belongs to.
        """
        return _get_element_by_id('fusion', self.feature_id)

    @property
    def feature(self):
        """
        The :class:`Fusion` feature this variant belongs to.
        """
        return self.fusion


class VariantGroup(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union(
        {'description', 'name', 'variant_ids', 'source_ids'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        # 'errors',                 # TODO: Add support for these fields in advanced search endpoint
         # 'lifecycle_actions',
        # 'provisional_values',
        'sources',
        'variants'
    })

    def __init__(self, **kwargs):
        self._variants = []
        self._sources = []
        super().__init__(**kwargs)

    @property
    def variants(self):
        for variant in self._variants:
            variant._include_status = self._include_status
        return self._variants

    @variants.setter
    def variants(self, value):
        self._variants = value

    @property
    def sources(self):
        return self._sources

    @sources.setter
    def sources(self, value):
        self._sources = value


class Gene(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union(
        {'description', 'entrez_id', 'name', 'source_ids'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'aliases',
        # 'errors',                 # TODO: Add support for these fields in advanced search endpoint
        # /'lifecycle_actions',
        # 'provisional_values',
        'sources',
        'variants',
    })

    def __init__(self, **kwargs):
        self._variants = []
        self._sources = []
        super().__init__(**kwargs)

    @property
    def variants(self):
        """
        A list of :class:`Variant` records associated with this gene.
        """
        for variant in self._variants:
            variant._include_status = self._include_status
        return [v for v in self._variants if v.molecular_profiles]

    @variants.setter
    def variants(self, value):
        self._variants = value

    @property
    def sources(self):
        """
        A list of :class:`Source` records associated with the gene description.
        """
        return self._sources

    @sources.setter
    def sources(self, value):
        self._sources = value


class Factor(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union(
        {'description', 'ncit_id', 'name', 'full_name', 'source_ids'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'aliases',
        # 'errors',                 # TODO: Add support for these fields in advanced search endpoint
        # /'lifecycle_actions',
        # 'provisional_values',
        'sources',
        'variants',
    })

    def __init__(self, **kwargs):
        self._variants = []
        self._sources = []
        super().__init__(**kwargs)

    @property
    def variants(self):
        """
        A list of :class:`Variant` records associated with this factor.
        """
        for variant in self._variants:
            variant._include_status = self._include_status
        return [v for v in self._variants if v.molecular_profiles]

    @variants.setter
    def variants(self, value):
        self._variants = value

    @property
    def sources(self):
        """
        A list of :class:`Source` records associated with the factor description.
        """
        return self._sources

    @sources.setter
    def sources(self, value):
        self._sources = value


class Fusion(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union(
        {'description', 'name', 'five_prime_partner_status', 'three_prime_partner_status', 'five_prime_gene_id', 'three_prime_gene_id', 'source_ids'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'aliases',
        # 'errors',                 # TODO: Add support for these fields in advanced search endpoint
        # /'lifecycle_actions',
        # 'provisional_values',
        'sources',
        'variants',
    })

    def __init__(self, **kwargs):
        self._variants = []
        self._sources = []
        super().__init__(**kwargs)

    @property
    def variants(self):
        """
        A list of :class:`Variant` records associated with this fusion.
        """
        for variant in self._variants:
            variant._include_status = self._include_status
        return [v for v in self._variants if v.molecular_profiles]

    @variants.setter
    def variants(self, value):
        self._variants = value

    @property
    def sources(self):
        """
        A list of :class:`Source` records associated with the fusion description.
        """
        return self._sources

    @sources.setter
    def sources(self, value):
        self._sources = value

    @property
    def five_prime_gene(self):
        """
        The :class:`Gene` record of the 5' fusion partner if that partner is ``KNOWN``.
        """
        if self.five_prime_gene_id:
            return get_gene_by_id(self.five_prime_gene_id)
        else:
            return None

    @property
    def three_prime_gene(self):
        """
        The :class:`Gene` record of the 3' fusion partner if that partner is ``KNOWN``.
        """
        if self.three_prime_gene_id:
            return get_gene_by_id(self.three_prime_gene_id)
        else:
            return None


class Evidence(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'assertion_ids',
        'description',
        'disease_id',
        'evidence_direction',
        'evidence_level',
        'evidence_type',
        'molecular_profile_id',
        'name',
        'phenotype_ids',
        'rating',
        'significance',
        'source_id',
        'status',
        'therapy_ids',
        'therapy_interaction_type',
        'variant_origin',
    })
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'assertions',
        'phenotypes',
        'therapies',
    })

    def __init__(self, **kwargs):
        self._assertions = []
        self._therapies = []
        self._phenotypes = []
        super().__init__(**kwargs)

    @property
    def molecular_profile(self):
        """
        The :class:`MolecularProfile` object this evidence item belongs to.
        """
        return get_molecular_profile_by_id(self.molecular_profile_id)

    @property
    def source(self):
        """
        A :class:`Source` object from which this evidence was derived.
        """
        return get_source_by_id(self.source_id)

    @property
    def assertions(self):
        """
        CIViC :class:`Assertion` records containing this evidence.
        """
        return [a for a in self._assertions if a.status in self._include_status]

    @assertions.setter
    def assertions(self, value):
        self._assertions = value

    @property
    def disease(self):
        """
        The :class:`Disease` record of the cancer or cancer subtype context for the evidence record. **None** for functional evidence_type.
        """
        if self.disease_id is not None:
            return get_disease_by_id(self.disease_id)
        else:
            return None


    @property
    def therapies(self):
        """
        Zero or more :class:`Therapy` records, linked to corresponding NCIt terms when applicable. Only used with therapeutic response predictive evidence_type.
        """
        return self._therapies

    @therapies.setter
    def therapies(self, value):
        self._therapies = value

    @property
    def phenotypes(self):
        """
        Zero or more :class:`Phenotype` records, linked to corresponding `Human Phenotype Ontology (HPO)`_ terms when applicable.

        .. _Human Phenotype Ontology (HPO): https://hpo.jax.org/
        """
        return self._phenotypes

    @phenotypes.setter
    def phenotypes(self, value):
        self._phenotypes = value

    @property
    def statement(self):
        """
        A shorthand for the evidence ``description``.
        """
        return self.description

    @property
    def variants(self):
        """
        One or more :class:`Variant` objects associated with this evidence's molecular profile.
        """
        return self.molecular_profile.variants


    def is_valid_for_gks_json(self, emit_warnings: bool = False) -> bool:
        """Determine whether Evidence is able to be represented as GKS model

        :param emit_warnings: Whether to log warnings, defaults to False
        :return: ``True`` if Evidence is able to be represented as GKS model. ``False``
            otherwise
        """
        return _is_valid_for_gks_json(self, emit_warnings)


class Assertion(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'amp_level',
        'assertion_direction',
        'assertion_type',
        'description',
        'disease_id',
        'endorsement_ids',
        'evidence_ids',
        'fda_companion_test',
        'fda_regulatory_approval',
        'molecular_profile_id',
        'name',
        'nccn_guideline',
        'nccn_guideline_version',
        'phenotype_ids',
        'significance',
        'status',
        'summary',
        'therapy_ids',
        'therapy_interaction_type',
        'variant_origin',
    })

    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'acmg_codes',
        'clingen_codes',
        'endorsements',
        'evidence_items',
        'phenotypes',
        'therapies',
    })

    def __init__(self, **kwargs):
        self._evidence_items = []
        self._therapies = []
        self._phenotypes = []
        self._endorsements = []
        super().__init__(**kwargs)

    @property
    def evidence(self):
        """
        A shortcut for the :attr:`evidence_items` property.
        """
        return self.evidence_items

    @property
    def evidence_items(self):
        """
        A list of :class:`Evidence` records supporting this assertion.
        """
        return [e for e in self._evidence_items if e.status in self._include_status]

    @evidence_items.setter
    def evidence_items(self, value):
        self._evidence_items = value

    @property
    def disease(self):
        """
        The :class:`Disease` record of the cancer or cancer subtype context for the assertion, linked to a corresponding `Disease Ontology`_ term when applicable.

        .. _Disease Ontology: http://disease-ontology.org/
        """
        return get_disease_by_id(self.disease_id)

    @property
    def endorsements(self):
        """
        Zero or more :class:`Endorsement` records representing organizations endorsing this assertion.
        """
        return self._endorsements

    @endorsements.setter
    def endorsements(self, value):
        self._endorsements = value

    @property
    def therapies(self):
        """
        Zero or more :class:`Therapy` records, linked to corresponding `NCIt`_ terms when applicable. Only used with therapeutic response predictive evidence_type.

        .. _NCIt: https://ncit.nci.nih.gov/ncitbrowser/
        """
        return self._therapies

    @therapies.setter
    def therapies(self, value):
        self._therapies = value

    @property
    def phenotypes(self):
        """
        Zero or more :class:`Phenotype` records associated with the assertion, linked to corresponding `Human Phenotype Ontology (HPO)`_ terms.

        .. _Human Phenotype Ontology (HPO): https://hpo.jax.org/
        """
        return self._phenotypes

    @phenotypes.setter
    def phenotypes(self, value):
        self._phenotypes = value

    @property
    def molecular_profile(self):
        """
        The :class:`MolecularProfile` object this assertion belongs to.
        """
        return get_molecular_profile_by_id(self.molecular_profile_id)

    @property
    def hpo_ids(self):
        """
        A list of `HPO`_ IDs of the :attr:`phenotypes` associated with this assertion

        .. _HPO: https://hpo.jax.org/
        """
        return [x.hpo_id for x in self.phenotypes if x.hpo_id]

    def format_nccn_guideline(self):
        if self.nccn_guideline is None:
            return ""
        else:
            return "{} (v{})".format(self.nccn_guideline, self.nccn_guideline_version)

    @property
    def variants(self):
        """
        One or more :class:`Variant` objects associated with this assertion's molecular profile.
        """
        return self.molecular_profile.variants

    def is_valid_for_gks_json(self, emit_warnings: bool = False) -> bool:
        """Determine whether Assertion is able to be represented as GKS model

        :param emit_warnings: Whether to log warnings, defaults to False
        :return: ``True`` if Assertion is able to be represented as GKS model. ``False``
            otherwise
        """
        return _is_valid_for_gks_json(self, emit_warnings)

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
        'description',
    })

    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        #'profile_image',
        'endorsements',
    })

    def __init__(self, **kwargs):
        self._endorsements = []
        super().__init__(**kwargs)

    @property
    def endorsements(self):
        """
        Zero or more :class:`Endorsement` records representing assertions endorsed by this organization.
        """
        return self._endorsements

    @endorsements.setter
    def endorsements(self, value):
        self._endorsements = value


class Therapy(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'ncit_id',
        'therapy_url',
        'name'
    })
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'aliases'
    })

    def __init__(self, **kwargs):
        self._evidence_items = []
        self._assertions = []
        super().__init__(**kwargs)

    def __str__(self):
        if self.ncit_id is None:
            return self.name
        else:
            return "{} (NCIt ID {})".format(self.name, self.ncit_id)

    @property
    def evidence(self):
        """
        A shortcut for the :attr:`evidence_items` property.
        """
        return self.evidence_items

    @property
    def evidence_items(self):
        """
        A list of :class:`Evidence` records linked to this therapy.
        """
        return [e for e in self._evidence_items if e.status in self._include_status]

    @evidence_items.setter
    def evidence_items(self, value):
        self._evidence_items = value

    @property
    def assertions(self):
        """
        A list of :class:`Assertion` records linked to this therapy.
        """
        return [a for a in self._assertions if a.status in self._include_status]

    @assertions.setter
    def assertions(self, value):
        self._assertions = value


class Disease(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'name',
        'doid',
        'disease_url'
    })
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'aliases'
    })

    def __init__(self, **kwargs):
        self._evidence_items = []
        self._assertions = []
        super().__init__(**kwargs)

    def __str__(self):
        if self.doid is None:
            return self.name
        else:
            return "{} (DOID {})".format(self.name, self.doid)

    @property
    def evidence(self):
        """
        A shortcut for the :attr:`evidence_items` property.
        """
        return self.evidence_items

    @property
    def evidence_items(self):
        """
        A list of :class:`Evidence` records linked to this disease.
        """
        return [e for e in self._evidence_items if e.status in self._include_status]

    @evidence_items.setter
    def evidence_items(self, value):
        self._evidence_items = value

    @property
    def assertions(self):
        """
        A list of :class:`Assertion` records linked to this disease.
        """
        return [a for a in self._assertions if a.status in self._include_status]

    @assertions.setter
    def assertions(self, value):
        self._assertions = value


class Phenotype(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'hpo_id',
        'phenotype_url',
        'name'
    })

    def __init__(self, **kwargs):
        self._evidence_items = []
        self._assertions = []
        super().__init__(**kwargs)

    def __str__(self):
        return "{} (HPO ID {})".format(self.name, self.hpo_id)

    @property
    def evidence(self):
        """
        A shortcut for the :attr:`evidence_items` property.
        """
        return self.evidence_items

    @property
    def evidence_items(self):
        """
        A list of :class:`Evidence` records linked to this phenotype.
        """
        return [e for e in self._evidence_items if e.status in self._include_status]

    @evidence_items.setter
    def evidence_items(self, value):
        self._evidence_items = value

    @property
    def assertions(self):
        """
        A list of :class:`Assertion` records linked to this phenotype.
        """
        return [a for a in self._assertions if a.status in self._include_status]

    @assertions.setter
    def assertions(self, value):
        self._assertions = value


class Source(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
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
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'clinical_trials'
    })

    def __init__(self, **kwargs):
        self._evidence_items = []
        self._genes = []
        self._factors = []
        self._fusions = []
        self._molecular_profiles = []
        super().__init__(**kwargs)

    def __str__(self):
        return "{} ({} {})".format(self.citation, self.source_type, self.citation_id)

    @property
    def evidence(self):
        """
        A shortcut for the :attr:`evidence_items` property.
        """
        return self.evidence_items

    @property
    def evidence_items(self):
        """
        A list of :class:`Evidence` records linked to this source.
        """
        return [e for e in self._evidence_items if e.status in self._include_status]

    @evidence_items.setter
    def evidence_items(self, value):
        self._evidence_items = value

    @property
    def genes(self):
        """
        A list of :class:`Gene` records supported by this source.
        """
        return self._genes

    @genes.setter
    def genes(self, value):
        self._genes = value

    @property
    def fusions(self):
        """
        A list of :class:`Fusion` records supported by this source.
        """
        return self._fusions

    @fusions.setter
    def fusions(self, value):
        self._fusions = value

    @property
    def factors(self):
        """
        A list of :class:`Factor` records supported by this source.
        """
        return self._factors

    @factors.setter
    def factors(self, value):
        self._factors = value

    @property
    def molecular_profiles(self):
        """
        A list of :class:`MolecularProfile` records supported by this source.
        """
        return self._molecular_profiles

    @molecular_profiles.setter
    def molecular_profiles(self, value):
        self._molecular_profiles = value


class Endorsement(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'assertion_id',
        'organization_id',
        'status',
        'last_reviewed',
        'ready_for_clinvar_submission',
    })

    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        #'profile_image',
    })

    @property
    def assertion(self):
        """
        The :class:`Assertion` object this endorsement endorses.
        """
        return get_assertion_by_id(self.assertion_id)

    @property
    def organization(self):
        """
        The :class:`Organization` object this endorsement was made on behalf of.
        """
        return get_organization_by_id(self.organization_id)


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


class Coordinate(CivicAttribute):
    _SIMPLE_FIELDS = CivicAttribute._SIMPLE_FIELDS.union({
         'chromosome',
         'start',
         'stop',
         'reference_bases',
         'variant_bases',
         'ensembl_version',
         'representative_transcript',
         'reference_build',
    })

    def __repr__(self):
        return '<CIViC Coordinate>'.format(self.type)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class ExonCoordinate(CivicAttribute):
    _SIMPLE_FIELDS = CivicAttribute._SIMPLE_FIELDS.union({
         'chromosome',
         'ensembl_id',
         'ensembl_version',
         'exon',
         'exon_offset',
         'exon_offset_direction',
         'reference_build',
         'representative_transcript',
         'start',
         'stop',
         'strand',
    })

    def __repr__(self):
        return '<CIViC ExonCoordinate>'.format(self.type)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


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
                logging.info('Loading {} from cache'.format(utils.pluralize(element)))
                return cached
        else:
            cached = [get_cached(element, element_id) for element_id in CACHE['{}_all_ids'.format(utils.pluralize(element))]]
            logging.info('Loading {} from cache'.format(utils.pluralize(element)))
            return cached
    if id_list and get_all:
        raise ValueError('Please pass list of ids or use the get_all flag, not both.')
    if get_all:
        logging.warning('Getting all {}. This may take a couple of minutes...'.format(utils.pluralize(element)))
        response_elements = _request_all(element)
    else:
        response_elements = _request_by_ids(element, id_list)

    elements = []
    ids = []
    for e in response_elements:
        e = _postprocess_response_element(e, element)
        if element == 'variant':
            cls = get_class(e['subtype'])
        else:
            cls = get_class(e['type'])
        partial_element = cls(**e, partial=True)
        ids.append(e['id'])
        elements.append(partial_element)

    CACHE['{}_all_ids'.format(utils.pluralize(element))] = ids
    return elements


def _postprocess_response_element(e, element):
    if e is None:
        raise Exception("{} not found".format(element.title()))
    e['type'] = element
    if element == 'assertion':
        e['molecular_profile_id'] = e['molecular_profile']['id']
        e['evidence_ids'] = [evidence['id'] for evidence in e['evidenceItems']]
        e['disease_id'] = e['disease']['id'] if e['disease'] is not None else None
        e['endorsement_ids'] = [endorsement['id'] for endorsement in e['endorsements']['nodes']]
        e['therapy_ids'] = [t['id'] for t in e['therapies']]
        e['phenotype_ids'] = [p['id'] for p in e['phenotypes']]
        e['status'] = e['status'].lower()
        del e['therapies']
        del e['endorsements']
    elif element == 'endorsement':
        e['assertion_id'] = e['assertion']['id']
        e['organization_id'] = e['organization']['id']
    elif element == 'evidence':
        e['source_id'] = e['source']['id']
        e['molecular_profile_id'] = e['molecular_profile']['id']
        e['assertion_ids'] = [a['id'] for a in e['assertions']]
        e['disease_id'] = e['disease']['id'] if e['disease'] is not None else None
        e['therapy_ids'] = [t['id'] for t in e['therapies']]
        e['phenotype_ids'] = [p['id'] for p in e['phenotypes']]
        e['status'] = e['status'].lower()
        del e['therapies']
        del e['phenotypes']
    elif element == 'factor':
        e['source_ids'] = [v['id'] for v in e['sources']]
        del e['sources']
    elif element == 'fusion':
        e['source_ids'] = [v['id'] for v in e['sources']]
        del e['sources']
        if e['threePrimeGene']:
            e['three_prime_gene_id'] = e['threePrimeGene']['id']
        else:
            e['three_prime_gene_id'] = None
        if e['fivePrimeGene']:
            e['five_prime_gene_id'] = e['fivePrimeGene']['id']
        else:
            e['five_prime_gene_id'] = None
    elif element == 'gene':
        e['source_ids'] = [v['id'] for v in e['sources']]
        del e['sources']
    elif element == 'molecular_profile':
        e['source_ids'] = [s['id'] for s in e['sources']]
        del e['sources']
        e['variant_ids'] = [v['id'] for v in e['variants']]
        del e['variants']
    elif element == 'variant':
        e['feature_id'] = e['feature']['id']
        if e['__typename'] == 'GeneVariant':
            e['subtype'] = 'gene_variant'
            e['entrez_id'] = e['feature']['featureInstance']['entrezId']
            e['entrez_name'] = e['feature']['name']
            build = e['coordinates']['reference_build']
            if build == 'GRCH37':
                build = 'GRCh37'
            elif build == 'GRCH38':
                build = 'GRCh38'
            e['coordinates']['reference_build'] = build
            if e['coordinates']['reference_bases'] in ['', '-']:
                e['coordinates']['reference_bases'] = None
            if e['coordinates']['variant_bases'] in ['', '-']:
                e['coordinates']['variant_bases'] = None
        elif e['__typename'] == 'FactorVariant':
            e['subtype'] = 'factor_variant'
        elif e['__typename'] == 'FusionVariant':
            if e['five_prime_start_exon_coordinates'] and e['five_prime_start_exon_coordinates']['exon_offset'] is None:
                e['five_prime_start_exon_coordinates']['exon_offset'] = 0
            if e['five_prime_end_exon_coordinates'] and e['five_prime_end_exon_coordinates']['exon_offset'] is None:
                e['five_prime_end_exon_coordinates']['exon_offset'] = 0
            if e['three_prime_start_exon_coordinates'] and e['three_prime_start_exon_coordinates']['exon_offset'] is None:
                e['three_prime_start_exon_coordinates']['exon_offset'] = 0
            if e['three_prime_end_exon_coordinates'] and e['three_prime_end_exon_coordinates']['exon_offset'] is None:
                e['three_prime_end_exon_coordinates']['exon_offset'] = 0
            e['subtype'] = 'fusion_variant'
        else:
            raise Exception("Variant type {} not supported yet".format(e['__typename']))
    elif element == 'variant_group':
        e['source_ids'] = [v['id'] for v in e['sources']]
        del e['sources']
        e['variant_ids'] = [v['id'] for v in e['variants']['nodes']]
        del e['variants']
    return e


def _get_element_by_id(element, id, allow_cached=True):
    return _get_elements_by_ids(element, [id], allow_cached)[0]


def _request_by_ids(element, ids):
    payload_methods = {
        'evidence': graphql_payloads._construct_get_evidence_payload,
        'gene': graphql_payloads._construct_get_gene_payload,
        'factor': graphql_payloads._construct_get_factor_payload,
        'fusion': graphql_payloads._construct_get_fusion_payload,
        'variant': graphql_payloads._construct_get_variant_payload,
        'assertion': graphql_payloads._construct_get_assertion_payload,
        'variant_group': graphql_payloads._construct_get_variant_group_payload,
        'molecular_profile': graphql_payloads._construct_get_molecular_profile_payload,
        'source': graphql_payloads._construct_get_source_payload,
        'disease': graphql_payloads._construct_get_disease_payload,
        'therapy': graphql_payloads._construct_get_therapy_payload,
        'phenotype': graphql_payloads._construct_get_phenotype_payload,
        'organization': graphql_payloads._construct_get_organization_payload,
        'endorsement': graphql_payloads._construct_get_endorsement_payload,
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
        'evidence': graphql_payloads._construct_get_all_evidence_payload,
        'gene': graphql_payloads._construct_get_all_genes_payload,
        'factor': graphql_payloads._construct_get_all_factors_payload,
        'fusion': graphql_payloads._construct_get_all_fusions_payload,
        'variant': graphql_payloads._construct_get_all_variants_payload,
        'assertion': graphql_payloads._construct_get_all_assertions_payload,
        'variant_group': graphql_payloads._construct_get_all_variant_groups_payload,
        'molecular_profile': graphql_payloads._construct_get_all_molecular_profiles_payload,
        'source': graphql_payloads._construct_get_all_sources_payload,
        'disease': graphql_payloads._construct_get_all_diseases_payload,
        'therapy': graphql_payloads._construct_get_all_therapies_payload,
        'phenotype': graphql_payloads._construct_get_all_phenotypes_payload,
        'organization': graphql_payloads._construct_get_all_organizations_payload,
        'endorsement': graphql_payloads._construct_get_all_endorsements_payload,
    }
    payload_method = payload_methods[element]
    payload = payload_method()

    after_cursor = None
    variables = { "after": after_cursor }
    resp = requests.post(API_URL, json={'query': payload, 'variables': variables}, timeout=(10,200))
    resp.raise_for_status()
    response = resp.json()['data'][utils.pluralize(element)]
    response_elements = response['nodes']
    has_next_page = response['pageInfo']['hasNextPage']
    after_cursor = response['pageInfo']['endCursor']

    while has_next_page:
        variables = {
          "after": after_cursor
        }
        resp = requests.post(API_URL, json={'query': payload, 'variables': variables}, timeout=(10,200))
        resp.raise_for_status()
        response = resp.json()['data'][utils.pluralize(element)]
        response_elements.extend(response['nodes'])
        has_next_page = response['pageInfo']['hasNextPage']
        after_cursor = response['pageInfo']['endCursor']

    return response_elements

#########################
# Get Entities By ID(s) #
#########################

# Evidence

def get_evidence_by_ids(evidence_id_list):
    """
    :param list evidence_id_list: A list of CIViC evidence item IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`EvidenceItem` objects.
    """
    logging.info('Getting evidence...')
    evidence = _get_elements_by_ids('evidence', evidence_id_list)
    logging.info('Caching evidence details...')
    for e in evidence:
        e._include_status = ['accepted', 'submitted', 'rejected']
    mp_ids = [x.molecular_profile.id for x in evidence]    # Add molecular profiles to cache
    _get_elements_by_ids('molecular_profile', mp_ids)
    for e in evidence:                        # Load from cache
        e.molecular_profile.update()
    return evidence


def get_evidence_by_id(evidence_id):
    """
    :param int phenotype_id: A single CIViC evidence item ID.
    :returns: A :class:`EvidenceItem` object.
    """
    return get_evidence_by_ids([evidence_id])[0]


# Molecular Profile

def get_molecular_profiles_by_ids(mp_id_list):
    """
    :param list mp_id_list: A list of CIViC molecular profile IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`MolecularProfile` objects.
    """
    logging.info('Getting molecular profiles...')
    mps = _get_elements_by_ids('molecular_profile', mp_id_list)
    for mp in mps:
        mp._include_status = ['accepted', 'submitted', 'rejected']
    #logging.info('Caching molecular profile details...')
    return mps


def get_molecular_profile_by_id(mp_id):
    """
    :param int mp_id: A single CIViC molecular profile ID.
    :returns: A :class:`MolecularProfile` object.
    """
    return get_molecular_profiles_by_ids([mp_id])[0]


# Assertion

def get_assertions_by_ids(assertion_id_list):
    """
    :param list assertion_id_list: A list of CIViC assertion IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Assertion` objects.
    """
    logging.info('Getting assertions...')
    assertions = _get_elements_by_ids('assertion', assertion_id_list)
    for a in assertions:
        a._include_status = ['accepted', 'submitted', 'rejected']
    logging.info('Caching variant details...')
    mp_ids = [x.molecular_profile.id for x in assertions]    # Add molecular profile to cache
    _get_elements_by_ids('molecular_profile', mp_ids)
    for assertion in assertions:                        # Load from cache
        assertion.molecular_profile.update()
    return assertions


def get_assertion_by_id(assertion_id):
    """
    :param int assertion_id: A single CIViC assertion ID.
    :returns: A :class:`Assertion` object.
    """
    return get_assertions_by_ids([assertion_id])[0]


# Variant

def get_variants_by_ids(variant_id_list):
    """
    :param list variant_id_list: A list of CIViC variant IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Variant` objects.
    """
    logging.info('Getting variants...')
    variants = _get_elements_by_ids('variant', variant_id_list)
    gene_ids = set()
    factor_ids = set()
    fusion_ids = set()
    for variant in variants:
        if isinstance(variant, GeneVariant):
            gene_ids.add(variant.feature_id)
        elif isinstance(variant, FactorVariant):
            factor_ids.add(variant.feature_id)
        elif isinstance(variant, FusionVariant):
            fusion_ids.add(variant.feature_id)
        variant._include_status = ['accepted', 'submitted', 'rejected']
    if gene_ids:
        logging.info('Caching gene details...')
        _get_elements_by_ids('gene', gene_ids)
    if factor_ids:
        logging.info('Caching factor details...')
        _get_elements_by_ids('factor', factor_ids)
    if fusion_ids:
        logging.info('Caching fusion details...')
        _get_elements_by_ids('fusion', fusion_ids)
    return variants


def get_variant_by_id(variant_id):
    """
    :param int variant_id: A single CIViC variant ID.
    :returns: A :class:`Variant` object.
    """
    return get_variants_by_ids([variant_id])[0]


# Variant Group

def get_variant_groups_by_ids(variant_group_id_list):
    """
    :param list variant_group_id_list: A list of CIViC variant group IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`VariantGroup` objects.
    """
    logging.info('Getting variant groups...')
    vgs = _get_elements_by_ids('variant_group', variant_group_id_list)
    for vg in vgs:
        vg._include_status = ['accepted', 'submitted', 'rejected']
    return vgs


def get_variant_group_by_id(variant_group_id):
    """
    :param int variant_group_id: A single CIViC variant group ID.
    :returns: A :class:`VariantGroup` object.
    """
    return get_variant_groups_by_ids([variant_group_id])[0]


# Feature

def get_features_by_ids(feature_id_list):
    """
    :param list feature_id_list: A list of CIViC feature IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Gene`, `Fusion`, and/or `Factor` objects.
    """
    logging.info('Getting features...')
    features = []
    for feature_id in feature_id_list:
        feature = None
        try:
            feature = _get_element_by_id('gene', feature_id)
        except:
            pass
        try:
            feature = _get_element_by_id('fusion', feature_id)
        except:
            pass
        try:
            feature = _get_element_by_id('factor', feature_id)
        except:
            pass
        if feature is None:
            raise Exception("Feature {} not found".format(feature_id))
        else:
            features.append(feature)
    variant_ids = set()
    for feature in features:
        feature._include_status = ['accepted', 'submitted', 'rejected']
        for variant in feature.variants:
            variant_ids.add(variant.id)
    if variant_ids:
        logging.info('Caching variant details...')
        _get_elements_by_ids('variant', variant_ids)
    for feature in features:
        for variant in feature.variants:
            variant.update()
    return features


def get_feature_by_id(feature_id):
    """
    :param int gene_id: A single CIViC feature ID.
    :returns: A :class:`Gene`, `Fusion`, or `Factor` object.
    """
    return get_features_by_ids([feature_id])[0]


def get_genes_by_ids(gene_id_list):
    """
    :param list gene_id_list: A list of CIViC gene feature IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Gene` objects.
    """
    logging.info('Getting genes...')
    genes = _get_elements_by_ids('gene', gene_id_list)
    variant_ids = set()
    for gene in genes:
        gene._include_status = ['accepted', 'submitted', 'rejected']
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
    :param int gene_id: A single CIViC gene feature ID.
    :returns: A :class:`Gene` object.
    """
    return get_genes_by_ids([gene_id])[0]


def get_fusions_by_ids(fusion_id_list):
    """
    :param list fusion_id_list: A list of CIViC fusion feature IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Fusion` objects.
    """
    logging.info('Getting fusions...')
    fusions = _get_elements_by_ids('fusion', fusion_id_list)
    variant_ids = set()
    for fusion in fusions:
        fusion._include_status = ['accepted', 'submitted', 'rejected']
        for variant in fusion.variants:
            variant_ids.add(variant.id)
    if variant_ids:
        logging.info('Caching variant details...')
        _get_elements_by_ids('variant', variant_ids)
    for fusion in fusions:
        for variant in fusion.variants:
            variant.update()
    return fusions


def get_fusion_by_id(fusion_id):
    """
    :param int fusion_id: A single CIViC fusion feature ID.
    :returns: A :class:`Fusion` object.
    """
    return get_fusions_by_ids([fusion_id])[0]


def get_factors_by_ids(factor_id_list):
    """
    :param list factor_id_list: A list of CIViC factor feature IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Factor` objects.
    """
    logging.info('Getting factors...')
    factors = _get_elements_by_ids('factor', factor_id_list)
    variant_ids = set()
    for factor in factors:
        factor._include_status = ['accepted', 'submitted', 'rejected']
        for variant in factor.variants:
            variant_ids.add(variant.id)
    if variant_ids:
        logging.info('Caching variant details...')
        _get_elements_by_ids('variant', variant_ids)
    for factor in factors:
        for variant in factor.variants:
            variant.update()
    return factors


def get_factor_by_id(factor_id):
    """
    :param int factor_id: A single CIViC factor feature ID.
    :returns: A :class:`Factor` object.
    """
    return get_factors_by_ids([factor_id])[0]


# Source

def get_sources_by_ids(source_id_list):
    """
    :param list source_id_list: A list of CIViC source IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Source` objects.
    """
    logging.info('Getting sources...')
    sources = _get_elements_by_ids('source', source_id_list)
    return sources


def get_source_by_id(source_id):
    """
    :param int source_id: A single CIViC source ID.
    :returns: A :class:`Source` object.
    """
    return get_sources_by_ids([source_id])[0]


# Disease

def get_diseases_by_ids(disease_id_list):
    """
    :param list disease_id_list: A list of CIViC disease IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Disease` objects.
    """
    logging.info('Getting diseases...')
    diseases = _get_elements_by_ids('disease', disease_id_list)
    return diseases


def get_disease_by_id(disease_id):
    """
    :param int disease_id: A single CIViC disease ID.
    :returns: A :class:`Disease` object.
    """
    return get_diseases_by_ids([disease_id])[0]


# Therapy

def get_therapies_by_ids(therapy_id_list):
    """
    :param list therapy_id_list: A list of CIViC therapy IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Therapy` objects.
    """
    logging.info('Getting therapies...')
    therapies = _get_elements_by_ids('therapy', therapy_id_list)
    return therapies


def get_therapy_by_id(therapy_id):
    """
    :param int therapy_id: A single CIViC therapy ID.
    :returns: A :class:`Therapy` object.
    """
    return get_therapies_by_ids([therapy_id])[0]


# Phenotype

def get_phenotypes_by_ids(phenotype_id_list):
    """
    :param list phenotype_id_list: A list of CIViC phenotype IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Phenotype` objects.
    """
    logging.info('Getting phenotypes...')
    phenotypes = _get_elements_by_ids('phenotype', phenotype_id_list)
    return phenotypes


def get_phenotype_by_id(phenotype_id):
    """
    :param int phenotype_id: A single CIViC phenotype ID.
    :returns: A :class:`Phenotype` object.
    """
    return get_phenotypes_by_ids([phenotype_id])[0]


# Organization

def get_organizations_by_ids(organization_id_list):
    """
    :param list organization_id_list: A list of CIViC organization IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Organization` objects.
    """
    logging.info('Getting organizations...')
    organizations = _get_elements_by_ids('organization', organization_id_list)
    return organizations

def get_organization_by_id(organization_id):
    """
    :param int organization_id: A single CIViC organization ID.
    :returns: A :class:`Organization` object.
    """
    return get_organizations_by_ids([organization_id])[0]


# Endorsement

def get_endorsements_by_ids(endorsement_id_list):
    """
    :param list endorsement_id_list: A list of CIViC endorsement IDs to query against to cache and (as needed) CIViC.
    :returns: A list of :class:`Endorsement` objects.
    """
    logging.info('Getting endorsements...')
    endorsements = _get_elements_by_ids('endorsement', endorsement_id_list)
    return endorsements

def get_endorsement_by_id(endorsement_id):
    """
    :param int endorsement_id: A single CIViC endorsement ID.
    :returns: A :class:`Endorsement` object.
    """
    return get_endorsements_by_ids([endorsement_id])[0]


###########
# Get All #
###########

# Assertion

def get_all_assertions(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all assertions.

    :param list include_status: A list of statuses. Only assertions matching the given statuses will be returned.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Assertion` objects.
    """
    assertions = _get_elements_by_ids('assertion', allow_cached=allow_cached, get_all=True)
    return [a for a in assertions if a.status in include_status]

def get_all_assertions_ready_for_clinvar_submission_for_org(organization_id, allow_cached=True):
    """
    Queries CIViC for all assertions endorsed by a specific organization that are ready for submission to ClinVar.

    :param int organization_id: The CIViC organization ID that endorsed the assertion(s) for submission to ClinVar.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Assertion` objects endorsed by a specific organization that are ready for submission to ClinVar.
    """
    endorsements = get_all_endorsements(include_status=["accepted"], allow_cached=allow_cached)
    assertions = []
    for e in endorsements:
        if e.organization_id == organization_id and e.ready_for_clinvar_submission:
            assertions.append(e.assertion)
    return assertions


# Molecular Profile

def get_all_molecular_profiles(include_status=['accepted', 'submitted', 'rejected'], allow_cached=True):
    """
    Queries CIViC for all molecular profiles.

    :param list include_status: A list of statuses. Only molecular profiles and their associated entities matching the given statuses will be returned. Use **None** to include molecular profiles without any associated entities.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`MolecularProfile` objects.
    """
    mps = _get_elements_by_ids('molecular_profile', allow_cached=allow_cached, get_all=True)
    if include_status:
        assert CACHE.get('evidence_items_all_ids', False)
        resp = list()
        for mp in mps:
            mp._include_status = include_status
            if mp.evidence:
                resp.append(mp)
        return resp
    else:
        return mps


# Variant

def get_all_variants(include_status=['accepted', 'submitted', 'rejected'], allow_cached=True):
    """
    Queries CIViC for all variants.

    :param list include_status: A list of statuses. Only variants and their associated entities matching the given statuses will be returned. Use **None** to include variants without any associated entities.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Variant` objects.
    """
    variants = _get_elements_by_ids('variant', allow_cached=allow_cached, get_all=allow_cached)
    if include_status:
        assert CACHE.get('evidence_items_all_ids', False)
        assert CACHE.get('assertions_all_ids', False)
        resp = list()
        for v in variants:
            v._include_status = include_status
            if v.molecular_profiles:
                resp.append(v)
        return resp
    else:
        return variants


def get_all_gene_variants(include_status=['accepted', 'submitted', 'rejected'], allow_cached=True):
    """
    Queries CIViC for all gene variants.

    :param list include_status: A list of statuses. Only variants and their associated entities matching the given statuses will be returned. Use **None** to include variants without any associated entities.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Variant` objects of **subtype** **gene_variant**.
    """
    variants = get_all_variants(include_status=include_status, allow_cached=allow_cached)
    return [v for v in variants if v.subtype == 'gene_variant']


def get_all_fusion_variants(include_status=['accepted', 'submitted', 'rejected'], allow_cached=True):
    """
    Queries CIViC for all fusion variants.

    :param list include_status: A list of statuses. Only variants and their associated entities matching the given statuses will be returned. Use **None** to include variants without any associated entities.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Variant` objects of **subtype** **fusion_variant**.
    """
    variants = get_all_variants(include_status=include_status, allow_cached=allow_cached)
    return [v for v in variants if v.subtype == 'fusion_variant']


def get_all_factor_variants(include_status=['accepted', 'submitted', 'rejected'], allow_cached=True):
    """
    Queries CIViC for all factor variants.

    :param list include_status: A list of statuses. Only variants and their associated entities matching the given statuses will be returned. Use **None** to include variants without any associated entities.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Variant` objects of **subtype** **factor_variant**.
    """
    variants = get_all_variants(include_status=include_status, allow_cached=True)
    return [v for v in variants if v.subtype == 'factor_variant']


# Variant Group

def get_all_variant_groups(allow_cached=True):
    """
    Queries CIViC for all variant groups.

    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`VariantGroup` objects.
    """
    variant_groups = _get_elements_by_ids('variant_group', allow_cached=allow_cached, get_all=True)
    return variant_groups


# Feature

def get_all_features(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all features.

    :param list include_status: A list of statuses. Only features and their associated entities matching the given statuses will be returned. Use **None** to include features without any associated entities.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Gene`, :class:`Fusion`, and/or :class:`Factor` objects.
    """
    genes = _get_elements_by_ids('gene', get_all=True, allow_cached=allow_cached)
    fusions = _get_elements_by_ids('fusion', get_all=True, allow_cached=allow_cached)
    factors = _get_elements_by_ids('factor', get_all=True, allow_cached=allow_cached)
    features = []
    features.extend(genes)
    features.extend(fusions)
    features.extend(factors)
    if include_status:
        assert CACHE.get('variants_all_ids', False)
        assert CACHE.get('evidence_items_all_ids', False)
        resp = list()
        for f in features:
            f._include_status = include_status
            if f.variants:
                resp.append(f)
        return resp
    else:
        return features



def get_all_genes(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all gene features.

    :param list include_status: A list of statuses. Only genes and their associated entities matching the given statuses will be returned. Use **None** to include genes without any associated entities.
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


def get_all_fusions(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all fusion features.

    :param list include_status: A list of statuses. Only fusions and their associated entities matching the given statuses will be returned. Use **None** to include fusions without any associated entities.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Fusion` objects.
    """
    fusions = _get_elements_by_ids('fusion', get_all=True, allow_cached=allow_cached)
    if include_status:
        assert CACHE.get('variants_all_ids', False)
        assert CACHE.get('evidence_items_all_ids', False)
        resp = list()
        for f in fusions:
            f._include_status = include_status
            if f.variants:
                resp.append(f)
        return resp
    else:
        return fusions


def get_all_factors(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all factor features.

    :param list include_status: A list of statuses. Only factors and their associated entities matching the given statuses will be returned. Use **None** to include factors without any associated entities.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Factor` objects.
    """
    factors = _get_elements_by_ids('factor', get_all=True, allow_cached=allow_cached)
    if include_status:
        assert CACHE.get('variants_all_ids', False)
        assert CACHE.get('evidence_items_all_ids', False)
        resp = list()
        for f in factors:
            f._include_status = include_status
            if f.variants:
                resp.append(f)
        return resp
    else:
        return factors


# Evidence

def get_all_evidence(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all evidence items.

    :param list include_status: A list of statuses. Only evidence items matching the given statuses will be returned.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`EvidenceItem` objects.
    """
    evidence = _get_elements_by_ids('evidence', get_all=True, allow_cached=allow_cached)
    return [e for e in evidence if e.status in include_status]


# Source

def get_all_sources(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all sources.

    :param list include_status: A list of statuses. Only sources and their associated entities matching the given statuses will be returned.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Source` objects.
    """
    sources = _get_elements_by_ids('source', get_all=True, allow_cached=allow_cached)
    if include_status:
        assert CACHE.get('evidence_items_all_ids', False)
        resp = list()
        for s in sources:
            s._include_status = include_status
            if s.evidence_items:
                resp.append(s)
        return resp
    else:
        return sources


# Disease

def get_all_diseases(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all diseases.

    :param list include_status: A list of statuses. Only diseases and their associated entities matching the given statuses will be returned.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Disease` objects.
    """
    diseases = _get_elements_by_ids('disease', get_all=True, allow_cached=allow_cached)
    if include_status:
        assert CACHE.get('evidence_items_all_ids', False)
        assert CACHE.get('assertions_all_ids', False)
        resp = list()
        for d in diseases:
            d._include_status = include_status
            if d.evidence_items or d.assertions:
                resp.append(d)
        return resp
    else:
        return diseases


# Therapy

def get_all_therapies(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all therapies.

    :param list include_status: A list of statuses. Only therapies and their associated entities matching the given statuses will be returned.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Therapy` objects.
    """
    therapies = _get_elements_by_ids('therapy', get_all=True, allow_cached=allow_cached)
    if include_status:
        assert CACHE.get('evidence_items_all_ids', False)
        assert CACHE.get('assertions_all_ids', False)
        resp = list()
        for t in therapies:
            t._include_status = include_status
            if t.evidence_items or t.assertions:
                resp.append(t)
        return resp
    else:
        return therapies


# Phenotype

def get_all_phenotypes(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all phenotypes.

    :param list include_status: A list of statuses. Only phenotypes and their associated entities matching the given statuses will be returned.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Phenotype` objects.
    """
    phenotypes = _get_elements_by_ids('phenotype', get_all=True, allow_cached=allow_cached)
    if include_status:
        assert CACHE.get('evidence_items_all_ids', False)
        assert CACHE.get('assertions_all_ids', False)
        resp = list()
        for p in phenotypes:
            p._include_status = include_status
            if p.evidence_items or p.assertions:
                resp.append(p)
        return resp
    else:
        return phenotypes


# Organization

def get_all_organizations(allow_cached=True):
    """
    Queries CIViC for all endorsements.

    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Organization` objects.
    """
    organizations = _get_elements_by_ids('organization', get_all=True, allow_cached=allow_cached)
    return organizations


# Endorsement

def get_all_endorsements(include_status=['accepted','submitted','rejected'], allow_cached=True):
    """
    Queries CIViC for all endorsements.

    :param list include_status: A list of statuses. Only endorsements for assertions matching the given statuses will be returned.
    :param bool allow_cached: Indicates whether or not object retrieval from CACHE is allowed. If **False** it will query the CIViC database directly.
    :returns: A list of :class:`Endorsement` objects.
    """
    endorsements = _get_elements_by_ids('endorsement', get_all=True, allow_cached=allow_cached)
    if include_status:
        assert CACHE.get('assertions_all_ids', False)
        resp = list()
        for e in endorsements:
            if e.assertion.status in include_status:
                resp.append(e)
        return resp
    else:
        return endorsements


#########################
# Search by Coordinates #
#########################

def search_evidence_by_coordinates(coordinates, search_mode='any'):
    """
    Search the cache for variants matching provided coordinates using the corresponding search mode and return all evidence items linked to any molecular profile involving those variants.

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

    :return:    A list of :class:`EvidenceItem` objects linked to molecular profiles involving variants matching the coordinates and search_mode
    """
    variants = search_variants_by_coordinates(coordinates, search_mode=search_mode)
    evidence = set()
    for v in variants:
        for mp in v.molecular_profiles:
            if mp.evidence:
                evidence.update(mp.evidence)
    return list(evidence)


def search_assertions_by_coordinates(coordinates, search_mode='any'):
    """
    Search the cache for variants matching provided coordinates using the corresponding search mode and return all assertions linked to any molecular profile involving those variants.

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

    :return:    A list of :class:`Assertion` objects linked to molecular profiles involving variants matching the coordinates and search_mode
    """
    variants = search_variants_by_coordinates(coordinates, search_mode=search_mode)
    assertions = set()
    for v in variants:
        for mp in v.molecular_profiles:
            if mp.assertions:
                assertions.update(mp.assertions)
    return list(assertions)


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
        match_idx = list(set(chr_ct_idx) & set(start_ct_idx) & set(stop_ct_idx))
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


##################################
# Search/get by other attributes #
##################################

# Genes

def get_gene_by_entrez_id(entrez_id):
    """
    :param str entrez_id: A gene `Entrez ID`_.
    :returns: A :class:`Gene` object.

    .. _Entrez ID: https://www.ncbi.nlm.nih.gov/gene/
    """
    genes = _get_elements_by_ids('gene', get_all=True)
    matching_genes = [g for g in genes if g.entrez_id == entrez_id]
    if len(matching_genes) == 0:
        raise Exception("No Gene with Entrez ID: {}".format(entrez_id))
    return matching_genes[0]


def get_gene_by_name(name):
    """
    :param str name: A `HGNC Gene Symbol`_.
    :returns: A :class:`Gene` object.

    .. _HGNC Gene Symbol: https://www.genenames.org/
    """
    genes = _get_elements_by_ids('gene', get_all=True)
    matching_genes = [g for g in genes if g.name == name]
    if len(matching_genes) == 0:
        raise Exception("No Gene with HGNC Gene Symbol: {}".format(name))
    return matching_genes[0]


# Factors

def get_factor_by_ncit_id(ncit_id):
    """
    :param str ncit_id: A factor `NCIthesaurus ID`_.
    :returns: A :class:`Factor` object.

    .. _NCIthesaurus ID: https://ncithesaurus.nci.nih.gov/ncitbrowser/
    """
    factors = _get_elements_by_ids('factor', get_all=True)
    matching_factors = [f for f in factors if f.ncit_id == ncit_id]
    if len(matching_factors) == 0:
        raise Exception("No Factor with NCIt ID: {}".format(ncit_id))
    return matching_factors[0]


def get_factor_by_name(name):
    """
    :param str name: A factor name or full name.
    :returns: A :class:`Factor` object.
    """
    factors = _get_elements_by_ids('factor', get_all=True)
    matching_factors = [f for f in factors if f.name == name or f.full_name == name]
    if len(matching_factors) == 0:
        raise Exception("No Factor with name or full name: {}".format(name))
    return matching_factors[0]


# Fusion

def get_fusion_by_name(name):
    """
    :param str name: A fusion name.
    :returns: A :class:`Fusion` object.
    """
    fusions = _get_elements_by_ids('fusion', get_all=True)
    matching_fusions = [f for f in fusions if f.name == name]
    if len(matching_fusions) == 0:
        raise Exception("No Fusion with name: {}".format(name))
    return matching_fusions[0]


def search_fusions_by_partner_gene_id(partner_gene_id):
    """
    :param int partner_gene_id: A CIViC ID of one of the gene partners.
    :returns: A list of :class:`Fusion` object.
    """
    fusions = _get_elements_by_ids('fusion', get_all=True)
    matching_fusions = [f for f in fusions if f.five_prime_gene_id == partner_gene_id or f.three_prime_gene_id == partner_gene_id]
    return matching_fusions


# Variants

def search_variants_by_allele_registry_id(caid):
    """
    Search the cache for variants matching the queried Allele Registry ID (CAID)

    :param str caid: `Allele Registry ID`_ to query
    :return: Returns a list of variant hashes matching the Allele Registry ID

    .. _Allele Registry ID: https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/landing
    """
    return search_variants_by_attribute('allele_registry_id', caid)


def search_variants_by_name(name):
    """
    Search the cache for variants matching the queried name

    :param str name: Variant name to query
    :return: Returns a list of variant hashes matching the name
    """
    return search_variants_by_attribute('name', name)


def search_variants_by_hgvs(hgvs):
    """
    Search the cache for variants matching the queried HGVS expression

    :param str name: HGVS expression to query
    :return: Returns a list of variant hashes matching the HGVS expression
    """
    return search_variants_by_list_field('hgvs_expressions', hgvs)


def search_variants_by_attribute(attribute, value):
    variants = get_all_variants()
    return [v for v in variants if hasattr(v, attribute) and getattr(v, attribute) == value]


def search_variants_by_list_field(field, value):
    variants = get_all_variants()
    matched_variants = []
    return [v for v in variants if hasattr(v, field) and value in getattr(v, field)]


# Source

def get_pubmed_source_by_id(pmid):
    """
    :param str pmid: A PubMed ID.
    :returns: A :class:`Source` object.
    """
    sources = _get_elements_by_ids('source', get_all=True)
    matching_sources = [s for s in sources if s.citation_id == pmid and s.source_type == 'PUBMED']
    if len(matching_sources) == 0:
        raise Exception("No PubMed sources with PMID: {}".format(pmid))
    return matching_sources[0]


def get_ash_source_by_doi(doi):
    """
    :param str doi: A ASH abstract DOI.
    :returns: A :class:`Source` object.
    """
    sources = _get_elements_by_ids('source', get_all=True)
    matching_sources = [s for s in sources if s.citation_id == doi and s.source_type == 'ASH']
    if len(matching_sources) == 0:
        raise Exception("No ASH sources with DOI: {}".format(doi))
    return matching_sources[0]


def get_asco_source_by_id(asco_id):
    """
    :param str asco_id: A ASCO Web ID. This is the identification number found in the URL of the abstract.
    :returns: A :class:`Source` object.
    """
    sources = _get_elements_by_ids('source', get_all=True)
    matching_sources = [s for s in sources if s.citation_id == asco_id and s.source_type == 'ASCO']
    if len(matching_sources) == 0:
        raise Exception("No ASCO sources with ID: {}".format(asco_id))
    return matching_sources[0]


# Disease

def get_disease_by_doid(doid):
    """
    :param str doid: A single `Disease Ontology ID`_.
    :returns: A :class:`Disease` object.

    .. _Disease Ontology ID: https://disease-ontology.org/
    """
    diseases = _get_elements_by_ids('disease', get_all=True)
    matching_diseases = [d for d in diseases if d.doid == doid]
    if len(matching_diseases) == 0:
        raise Exception("No diseases with DO ID: {}".format(doid))
    return matching_diseases[0]


def get_disease_by_name(name):
    """
    :param str name: A single `Disease Ontology`_ name.
    :returns: A :class:`Disease` object.

    .. _Disease Ontology: https://disease-ontology.org/
    """
    diseases = _get_elements_by_ids('disease', get_all=True)
    matching_diseases = [d for d in diseases if d.name == name]
    if len(matching_diseases) == 0:
        raise Exception("No diseases with DO name: {}".format(name))
    return matching_diseases[0]


# Therapy

def get_therapy_by_ncit_id(ncit_id):
    """
    :param str ncit_id: A single `NCIthesaurus ID`_.
    :returns: A :class:`Therapy` object.

    .. _NCIthesaurus ID: https://ncithesaurus.nci.nih.gov/ncitbrowser/
    """
    therapies = _get_elements_by_ids('therapy', get_all=True)
    matching_therapies = [t for t in therapies if t.ncit_id == ncit_id]
    if len(matching_therapies) == 0:
        raise Exception("No therapies with NCIt ID: {}".format(ncit_id))
    return matching_therapies[0]


def get_therapy_by_name(name):
    """
    :param str name: A single `NCIthesaurus`_ name.
    :returns: A :class:`Therapy` object.

    .. _NCIthesaurus: https://ncithesaurus.nci.nih.gov/ncitbrowser/
    """
    therapies = _get_elements_by_ids('therapy', get_all=True)
    matching_therapies = [t for t in therapies if t.name == name]
    if len(matching_therapies) == 0:
        raise Exception("No therapies with NCIt name: {}".format(name))
    return matching_therapies[0]


# Phenotype

def get_phenotype_by_hpo_id(hpo_id):
    """
    :param str hpo_id: A single `Human Phenotype Ontology ID`_.
    :returns: A :class:`Phenotype` object.

    .. _Human Phenotype Ontology ID: https://hpo.jax.org/
    """
    phenotypes = _get_elements_by_ids('phenotype', get_all=True)
    matching_phenotypes = [p for p in phenotypes if p.hpo_id == hpo_id]
    if len(matching_phenotypes) == 0:
        raise Exception("No phenotypes with HPO ID: {}".format(hpo_id))
    return matching_phenotypes[0]


def get_phenotype_by_name(name):
    """
    :param str name: A single `Human Phenotype Ontology`_ name (sometimes also referred to as HPO class).
    :returns: A :class:`Phenotype` object.

    .. _Human Phenotype Ontology: https://hpo.jax.org/
    """
    phenotypes = _get_elements_by_ids('phenotype', get_all=True)
    matching_phenotypes = [p for p in phenotypes if p.name == name]
    if len(matching_phenotypes) == 0:
        raise Exception("No phenotypes with name: {}".format(name))
    return matching_phenotypes[0]

# Endorsement

def search_endorsements_by_organization_id(organization_id):
    """
    :param int organization_id: A CIViC :class:`Organization` ID.
    :returns: A list of :class:`Endorsement` objects.
    """
    endorsements = _get_elements_by_ids('endorsement', get_all=True)
    matching_endorsements = [e for e in endorsements if e.organization_id == organization_id]
    return matching_endorsements

def search_endorsements_by_assertion_id(assertion_id):
    """
    :param int assertion_id: A CIViC :class:`Assertion` ID.
    :returns: A list of :class:`Endorsement` objects.
    """
    endorsements = _get_elements_by_ids('endorsement', get_all=True)
    matching_endorsements = [e for e in endorsements if e.assertion_id == assertion_id]
    return matching_endorsements


def get_all_endorsements_ready_for_clinvar_submission_for_org(
    organization_id: int,
) -> list[Endorsement]:
    """
    Queries CIViC for all endorsements by a specific organization that are ready for submission to ClinVar.

    :param int organization_id: A CIViC :class:`Organization` ID.
    :returns: A list of :class:`Endorsement` objects endorsed by a specific organization
        that are ready for submission to ClinVar.
    """
    endorsements = search_endorsements_by_organization_id(organization_id)
    return [e for e in endorsements if e.ready_for_clinvar_submission]

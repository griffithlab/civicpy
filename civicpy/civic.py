import requests
import importlib
import logging
import datetime
import pandas as pd
import pickle
from civicpy import DATA_ROOT


CACHE = dict()

CACHE_FILE = DATA_ROOT / 'CACHE.pkl'

COORDINATE_TABLE = None

HPO_TERMS = dict()

FRESH_DELTA = datetime.timedelta(days=7)

MODULE = importlib.import_module('civicpy.civic')

API_URL = 'https://civicdb.org/api'

LINKS_URL = 'https://civicdb.org/links'

UNMARKED_PLURALS = {'evidence'}

CIVIC_TO_PYCLASS = {
    'evidence_items': 'evidence'
}


def pluralize(string):
    if string in UNMARKED_PLURALS:
        return f'{string}_items'
    if string.endswith('s'):
        return string
    return string + 's'


def singularize(string):
    string = string.rstrip('s')
    if string == 'evidence_item':
        string = 'evidence'
    return string


def search_url(element):
    element = pluralize(element).lower()
    return '/'.join([API_URL, element, 'search'])


def snake_to_camel(snake_string):
    words = snake_string.split('_')
    cap_words = [x.capitalize() for x in words]
    return ''.join(cap_words)


def element_lookup_by_id(element_type, element_id):
    e_string = pluralize(element_type.lower())
    url = '/'.join([API_URL, e_string, str(element_id)])
    resp = requests.get(url)
    resp.raise_for_status()
    resp_dict = resp.json()
    return resp_dict


def get_class(element_type):
    e_string = singularize(element_type)
    class_string = snake_to_camel(e_string)
    cls = getattr(MODULE, class_string, CivicAttribute)
    return cls


def save_cache():
    with open(CACHE_FILE, 'wb') as pf:
        pickle.dump(CACHE, pf)


def load_cache():
    with open(CACHE_FILE, 'rb') as pf:
        old_cache = pickle.load(pf)
    c = dict()
    variants = set()
    for k, v in old_cache.items():
        if isinstance(k, str):
            c[k] = v
        elif isinstance(k, int):
            c[hash(v)] = v
            if v.type == 'variant':
                variants.add(v)
        else:
            raise ValueError
    MODULE.CACHE = c
    _build_coordinate_table(variants)

class CivicRecord:

    _SIMPLE_FIELDS = {'id', 'type'}
    _COMPLEX_FIELDS = set()

    def __init__(self, partial=False, **kwargs):
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
                    if partial and field not in CivicRecord._SIMPLE_FIELDS:
                        self._incomplete.add(field)     # Allow for incomplete data when partial flag set
                    else:
                        raise AttributeError(f'Expected {field} attribute for {self.type}, none found.')

        for field in self._COMPLEX_FIELDS:
            try:
                v = kwargs[field]
            except KeyError:
                if partial:
                    self._incomplete.add(field)
                    continue
                else:
                    raise AttributeError(f'Expected {field} attribute for {self.type}, none found.')
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
                self.__setattr__(field, cls(partial=True, **v))

        self._partial = bool(self._incomplete)
        if not isinstance(self, CivicAttribute) and not self._partial and self.__class__.__name__ != 'CivicRecord':
            CACHE[hash(self)] = self

    def __repr__(self):
        return f'<CIViC {self.type} {self.id}>'

    def __getattr__(self, item):
        if self._partial and item in self._incomplete:
            self.update()
        return object.__getattribute__(self, item)

    def __hash__(self):
        return hash(f'{self.type}:{self.id}')

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __setstate__(self, state):
        self.__dict__ = state

    def update(self, allow_partial=True, force=False, **kwargs):
        """Updates record and returns True if record is complete after update, else False."""
        if kwargs:
            self.__init__(partial=allow_partial, force=force, **kwargs)
            return not self._partial

        if not force and CACHE.get(hash(self)):
            cached = CACHE[hash(self)]
            for field in self._SIMPLE_FIELDS | self._COMPLEX_FIELDS:
                v = getattr(cached, field)
                setattr(self, field, v)
            self._partial = False
            logging.info(f'Loading {str(self)} from cache')
            return True
        resp_dict = element_lookup_by_id(self.type, self.id)
        self.__init__(partial=False, **resp_dict)
        return True

    @property
    def site_link(self):
        return '/'.join([LINKS_URL, self.type, str(self.id)])


class Variant(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'allele_registry_id',
        'civic_actionability_score',
        'description',
        'entrez_id',
        'entrez_name',
        'gene_id',
        'name'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'assertions',
        'clinvar_entries',
        'coordinates',
        # 'errors',
        'evidence_items',
        'hgvs_expressions',
        # 'lifecycle_actions',
        # 'provisional_values',
        'sources',
        'variant_aliases',
        'variant_groups',
        'variant_types'})

    def __init__(self, **kwargs):
        # Handle overloaded evidence_items from some advanced search views
        evidence_items = kwargs.get('evidence_items')
        kwargs['type'] = 'variant'
        if evidence_items and not isinstance(evidence_items, list):
                del(kwargs['evidence_items'])
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
    def evidence(self):
        return self.evidence_items

    @property
    def gene(self):
        return _get_element_by_id('gene', self.gene_id)


class Gene(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union(
        {'description', 'entrez_id', 'name'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'aliases',
        # 'errors',                 # TODO: Add support for these fields in advanced search endpoint
        # 'lifecycle_actions',
        # 'provisional_values',
        # 'sources',
        'variants'
    })


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
        'open_change_count',
        'rating',
        'status',
        'variant_id',
        'variant_origin'})
    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'assertions',
        'disease',
        'drugs',
        'errors',
        'fields_with_pending_changes',
        'lifecycle_actions',
        'phenotypes',
        'source'})

    @property
    def variant(self):
        return get_variants_by_ids([self.variant_id])[0]


class Assertion(CivicRecord):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({
        'allele_registry_id',
        'amp_level',
        'clinical_significance',
        'description',
        'drug_interaction_type',
        'evidence_direction',
        'evidence_item_count',
        'evidence_type',
        'fda_companion_test',
        'fda_regulatory_approval',
        'name',
        'nccn_guideline',
        'nccn_guideline_version',
        'open_change_count',
        'pending_evidence_count',
        'status',
        'summary',
        'variant_origin'
    })

    _COMPLEX_FIELDS = CivicRecord._COMPLEX_FIELDS.union({
        'acmg_codes',
        'disease',
        'drugs',
        'evidence_items',
        'gene',
        'lifecycle_actions',
        'phenotypes',
        'variant'
    })

    @property
    def evidence(self):
        return self.evidence_items

    @property
    def hpo_ids(self):
        return [x.hpo_id for x in self.phenotypes if x.hpo_id]


class CivicAttribute(CivicRecord, dict):

    _SIMPLE_FIELDS = {'type'}
    _COMPLEX_FIELDS = set()

    def __repr__(self):
        try:
            _id = self.id
        except AttributeError:
            return f'<CIViC Attribute {self.type}>'
        else:
            return f'<CIViC Attribute {self.type} {self.id}>'

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
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({'pubchem_id'})


class Disease(CivicAttribute):
    _SIMPLE_FIELDS = CivicRecord._SIMPLE_FIELDS.union({'display_name', 'doid', 'url'})


def get_cached(element_type, element_id):
    klass = get_class(element_type)
    r = klass(type=element_type, id=element_id, partial=True)
    return CACHE.get(hash(r), False)


def _has_all_cached_fresh(element):
    s = '{}_all_cached'.format(pluralize(element))
    if CACHE.get(s, False):
        return CACHE[s] + FRESH_DELTA > datetime.datetime.now()
    return False


def _get_elements_by_ids(element, id_list=[], allow_cached=True, get_all=False):
    if allow_cached and not get_all:
        cached = [get_cached(element, element_id) for element_id in id_list]
        if all(cached):
            logging.info(f'Loading {pluralize(element)} from cache')
            return cached
    elif allow_cached and _has_all_cached_fresh(element):
        cached = [get_cached(element, element_id) for element_id in CACHE['{}_all_ids'.format(pluralize(element))]]
        logging.info(f'Loading {pluralize(element)} from cache')
        return cached
    if id_list and get_all:
        raise ValueError('Please pass list of ids or use the get_all flag, not both.')
    if get_all:
        payload = _construct_get_all_payload()
        logging.warning('Getting all {}. This may take a couple of minutes...'.format(pluralize(element)))
    else:
        payload = _construct_query_payload(id_list)
    url = search_url(element)
    response = requests.post(url, json=payload)
    response.raise_for_status()
    cls = get_class(element)
    elements = [cls(**x) for x in response.json()['results']]
    CACHE['{}_all_cached'.format(pluralize(element))] = datetime.datetime.now()
    CACHE['{}_all_ids'.format(pluralize(element))] = [x['id'] for x in response.json()['results']]
    return elements


def _get_element_by_id(element, id, allow_cached=True):
    return _get_elements_by_ids(element, [id], allow_cached)[0]


def _construct_get_all_payload():
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
        'queries': queries
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


def get_all_assertion_ids():
    return _get_all_element_ids('assertions')


def get_all_assertions():
    return get_assertions_by_ids(get_all=True)


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


def get_all_variants(allow_cached=True):
    precached = _has_all_cached_fresh('variants')
    variants = _get_all_genes_and_variants(allow_cached)['variants']
    if not (precached and allow_cached):
        _build_coordinate_table(variants)
    return variants



def get_all_variant_ids():
    return _get_all_element_ids('variants')


def _get_all_genes_and_variants(allow_cached=True):
    variants = _get_elements_by_ids('variants', get_all=True, allow_cached=allow_cached)
    genes = _get_elements_by_ids('gene', get_all=True, allow_cached=allow_cached)
    for variant in variants:
        variant.gene.update()
    return {'genes': genes, 'variants': variants}


def _get_all_element_ids(element):
    url = f'https://civicdb.org/api/{element}?count=100000'
    resp = requests.get(url)
    resp.raise_for_status()
    return [x['id'] for x in resp.json()['records']]


def get_genes_by_ids(gene_id_list):
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
    return get_genes_by_ids([gene_id])[0]


def get_all_gene_ids():
    return _get_all_element_ids('genes')


def get_all_genes():
    return _get_all_genes_and_variants()['genes']


def get_all_evidence_ids():
    return _get_all_element_ids('evidence_items')


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

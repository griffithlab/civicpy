import requests
import importlib
import logging

CACHE = dict()

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
    if e_string == 'evidence':
        e_string = 'evidence_items'
    url = '/'.join([API_URL, e_string, str(element_id)])
    resp = requests.get(url)
    resp.raise_for_status()
    resp_dict = resp.json()
    return resp_dict


def get_class(element_type):
    e_string = singularize(element_type)
    class_string = snake_to_camel(e_string)
    cls = getattr(MODULE, class_string, Attribute)
    return cls


class CivicRecord:

    SIMPLE_FIELDS = {'id', 'type'}
    COMPLEX_FIELDS = set()

    def __init__(self, partial=False, **kwargs):
        self._incomplete = set()
        self.partial = partial
        simple_fields = sorted(self.SIMPLE_FIELDS, reverse=True)
        simple_fields = sorted(simple_fields, key=lambda x: x in CivicRecord.SIMPLE_FIELDS, reverse=True)
        for field in simple_fields:
            try:
                self.__setattr__(field, kwargs[field])
            except KeyError:
                try:
                    object.__getattribute__(self, field)
                except AttributeError:
                    if partial and field not in CivicRecord.SIMPLE_FIELDS:
                        self._incomplete.add(field)     # Allow for incomplete data when partial flag set
                    else:
                        raise AttributeError(f'Expected {field} attribute for {self.type}, none found.')

        for field in self.COMPLEX_FIELDS:
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

        self.partial = bool(self._incomplete)
        if not isinstance(self, Attribute) and not self.partial and self.__class__.__name__ != 'CivicRecord':
            CACHE[hash(self)] = self

    def __repr__(self):
        return f'<CIViC {self.type} {self.id}>'

    def __getattr__(self, item):
        if self.partial and item in self._incomplete:
            self.update()
        return object.__getattribute__(self, item)

    def __hash__(self):
        return hash(f'{self.type}:{self.id}')

    def update(self, allow_partial=True, **kwargs):
        """Updates record and returns True if record is complete after update, else False."""
        if kwargs:
            self.__init__(partial=allow_partial, **kwargs)
            return not self.partial

        if CACHE.get(hash(self)):
            cached = CACHE[hash(self)]
            for field in self.SIMPLE_FIELDS | self.COMPLEX_FIELDS:
                v = getattr(cached, field)
                setattr(self, field, v)
            self.partial=False
            logging.info(f'Loading {str(self)} from cache')
            return True
        resp_dict = element_lookup_by_id(self.type, self.id)
        self.__init__(partial=False, **resp_dict)
        return True

    @property
    def site_link(self):
        return '/'.join([LINKS_URL, self.type, str(self.id)])


class Variant(CivicRecord):
    SIMPLE_FIELDS = {
        'allele_registry_id',
        'civic_actionability_score',
        'description',
        'entrez_id',
        'entrez_name',
        'gene_id',
        'id',
        'name',
        'type'}
    COMPLEX_FIELDS = {
        'assertions',
        'clinvar_entries',
        'coordinates',
        'errors',
        'evidence_items',
        'hgvs_expressions',
        'lifecycle_actions',
        'provisional_values',
        'sources',
        'variant_aliases',
        'variant_groups',
        'variant_types'}

    def __init__(self, **kwargs):
        # Handle overloaded evidence_items from some advanced search views
        evidence_items = kwargs.get('evidence_items')
        if evidence_items and not isinstance(evidence_items, list):
                del(kwargs['evidence_items'])
        super().__init__(**kwargs)

    @property
    def source_ids(self):
        ids = set()
        for evidence in self.evidence_items:
            if evidence.source.pubmed_id:
                ids.add(f'PMID:{evidence.source.pubmed_id}')
            else:
                raise NotImplementedError('No method (yet) to support evidence other than PMID')
        return ids

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


class Gene(CivicRecord):
    SIMPLE_FIELDS = {'description', 'entrez_id', 'id', 'name', 'type'}
    COMPLEX_FIELDS = {
        'aliases',
        # 'errors',                 # TODO: Add support for these fields in advanced search endpoint
        # 'lifecycle_actions',
        # 'provisional_values',
        # 'sources',
        'variants'
    }


class Evidence(CivicRecord):
    SIMPLE_FIELDS = {
        'clinical_significance',
        'description',
        'drug_interaction_type',
        'evidence_direction',
        'evidence_level',
        'evidence_type',
        'gene_id',
        'id',
        'name',
        'open_change_count',
        'rating',
        'status',
        'type',
        'variant_id',
        'variant_origin'}
    COMPLEX_FIELDS = {
        'assertions',
        'disease',
        'drugs',
        'errors',
        'fields_with_pending_changes',
        'lifecycle_actions',
        'phenotypes',
        'source'}


class Assertion(CivicRecord):
    SIMPLE_FIELDS = {
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
        'id',
        'name',
        'nccn_guideline',
        'nccn_guideline_version',
        'open_change_count',
        'pending_evidence_count',
        'status',
        'summary',
        'type',
        'variant_origin'
    }

    COMPLEX_FIELDS = CivicRecord.COMPLEX_FIELDS.union({
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


class Attribute(CivicRecord):

    SIMPLE_FIELDS = {'type'}
    COMPLEX_FIELDS = set()

    def __repr__(self):
        return f'<CIViC Attribute {self.type}>'

    def __init__(self, **kwargs):
        kwargs['partial'] = False
        for k, v in kwargs.items():
            self.__setattr__(k, v)
        super().__init__(**kwargs)

    def __hash__(self):
        return object.__hash__(self)

    @property
    def site_link(self):
        return None


class Drug(Attribute):
    SIMPLE_FIELDS = CivicRecord.SIMPLE_FIELDS.union({'pubchem_id'})


class Disease(Attribute):
    SIMPLE_FIELDS = CivicRecord.SIMPLE_FIELDS.union({'display_name', 'doid', 'url'})


def get_cached(element_type, element_id):
    r = Variant(type=element_type, id=element_id, partial=True)
    return CACHE.get(hash(r), False)


def get_elements_by_ids(element, id_list, allow_cached=True):
    if allow_cached:
        cached = [get_cached(element, element_id) for element_id in id_list]
        if all(cached):
            logging.info(f'Loading {pluralize(element)} from cache')
            return cached
    payload = _construct_query_payload(id_list)
    url = search_url(element)
    response = requests.post(url, json=payload)
    response.raise_for_status()
    cls = get_class(element)
    elements = [cls(**x) for x in response.json()['results']]
    return elements


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


def get_assertions_by_ids(id_list):
    logging.info('Getting assertions...')
    assertions = get_elements_by_ids('assertion', id_list)
    logging.info('Caching variant details...')
    variant_ids = [x.variant.id for x in assertions]    # Add variants to cache
    get_elements_by_ids('variant', variant_ids)
    logging.info('Caching gene details...')
    gene_ids = [x.gene.id for x in assertions]          # Add genes to cache
    get_elements_by_ids('gene', gene_ids)
    for assertion in assertions:                        # Load from cache
        assertion.variant.update()
        assertion.gene.update()
    return assertions


def get_variants_by_ids(id_list):
    logging.info('Getting variants...')
    return get_elements_by_ids('variant', id_list)


def get_genes_by_ids(id_list):
    logging.info('Getting genes...')
    genes = get_elements_by_ids('gene', id_list)  # Advanced search results are incomplete
    variant_ids = set()
    for gene in genes:
        for variant in gene.variants:
            variant_ids.add(variant.id)
    if variant_ids:
        logging.info('Caching variant details...')
        get_elements_by_ids('variant', variant_ids)
    for gene in genes:
        for variant in gene.variants:
            variant.update()
    return genes

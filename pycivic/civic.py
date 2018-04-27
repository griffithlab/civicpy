import requests
import importlib

MODULE = importlib.import_module('pycivic.civic')

API_URL = 'https://civicdb.org/api'


def search_url(element):
    return '/'.join([API_URL, element, 'search'])


class CivicRecord:

    SIMPLE_FIELDS = ['id', 'type', 'name']
    COMPLEX_FIELDS = []

    def __init__(self, partial=False, **kwargs):
        self._incomplete = set()
        self.partial = partial
        for field in self.SIMPLE_FIELDS:
            if partial:
                try:
                    self.__setattr__(field, kwargs[field])
                except KeyError:
                    self._incomplete.add(field)      # Allow for incomplete data when partial flag set
            else:
                self.__setattr__(field, kwargs[field])  # Intentionally raise error if missing expected fields
                                                        #  w/o partial flag
        for field in self.COMPLEX_FIELDS:
            v = kwargs[field]                           # Will error out if missing expected complex field
            is_compound = isinstance(v, list)
            if is_compound:
                klass_string = field.rstrip('s').capitalize()
                klass = getattr(MODULE, klass_string, CivicRecord)
                result = list()
                for data in v:
                    result.append(klass(partial=True, **data))
                self.__setattr__(field, result)
            else:
                klass = getattr(MODULE, field.capitalize(), CivicRecord)
                self.__setattr__(field, klass(partial=True, **v))

        self.partial = bool(self._incomplete)

    def __repr__(self):
        return f'[<CIViC {self.type}>]: {self.id}'

    def __getattr__(self, item):
        if self.partial and item in self._incomplete:
            self.update()
        return object.__getattr__(self, item)

    def update(self, data=None):
        raise NotImplementedError


class Variant(CivicRecord):
    SIMPLE_FIELDS = CivicRecord.SIMPLE_FIELDS + ['description']


class Assertion(CivicRecord):

    SIMPLE_FIELDS = CivicRecord.SIMPLE_FIELDS + \
      [
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
      ]

    COMPLEX_FIELDS = [
        'acmg_codes',
        'disease',
        'drugs',
        # 'errors',  TODO: evaluate the expected contents of errors and consider adding support
        'evidence_items',
        'gene',
        'lifecycle_actions',
        'phenotypes',
        'provisional_values',
        # 'state_params',
        'variant'
    ]


def get_assertions(assertion_id_list):
    queries = list()
    for assertion_id in assertion_id_list:
        query = {
            'field': 'id',
            'condition': {
                'name': 'is_equal_to',
                'parameters': [
                    assertion_id
                ]
            }
        }
        queries.append(query)
    payload = {
        'operator': 'OR',
        'queries': queries
    }
    url = search_url('assertions')
    response = requests.post(url, json=payload)
    response.raise_for_status()
    # return response.json()['results']
    assertions = [Assertion(**x) for x in response.json()['results']]
    return assertions

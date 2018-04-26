import requests


class RecordSet(set):

    pass


class CivicRecord:

    def __init__(self, **kwargs):
        self.id = kwargs['id']
        self.type = kwargs['type']
        self.name = kwargs['name']

    def __repr__(self):
        return f'[<CIViC {self.type}>]: {self.id}'


class Assertion(CivicRecord):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # TODO: other Assertion fields


API_URL = 'https://civicdb.org/api'
ASSERTION_SEARCH_URL = '/'.join([API_URL, 'assertions', 'search'])
VARIANT_SEARCH_URL   = '/'.join([API_URL, 'variants',   'search'])
EVIDENCE_SEARCH_URL  = '/'.join([API_URL, 'evidence',   'search'])


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
    response = requests.post(ASSERTION_SEARCH_URL, json=payload)
    response.raise_for_status()
    assertions = [Assertion(**x) for x in response.json()['results']]
    return assertions


def search_assertions(queries):
    raise NotImplementedError
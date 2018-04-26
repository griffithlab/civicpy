import requests


class CivicRecord:

    def __init__(self, **kwargs):
        self.id = kwargs['id']
        self.type = kwargs['type']


class Assertion(CivicRecord):

    def __init__(self, **kwargs):
        super().__init__(kwargs)
        # TODO: other Assertion fields


class CivicDb:

    API_URL = 'https://civicdb.org/api'

    def __init__(self):
        self._assertions = dict()
        self._variants = dict()
        self._evidence = dict()

    def get_assertions(self, assertion_id_list, update=False):
        out = set()
        for assertion_id in assertion_id_list:
            out.add(self.get_assertion(assertion_id, update))
        return out

    def get_assertion(self, assertion_id, update=False):
        assertion_id_int = int(assertion_id)
        if not self._assertions.get(assertion_id_int, False) or update:
            url = '/'.join(CivicDb.API_URL, assertion_id_int)
            resp = requests.get(url)
            resp.raise_for_status()
            self._assertions[assertion_id_int] = Assertion(resp.json())
        return self._assertions[assertion_id_int]

    def get_all_assertions(self):
        raise NotImplementedError

    def list_loaded_assertions(self):
        return list(self._assertions.values())
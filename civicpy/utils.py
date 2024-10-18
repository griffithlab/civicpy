UNMARKED_PLURALS = {'evidence'}

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
    if '_' not in snake_string and any(ele.isupper() for ele in snake_string):
        #this is already camel case
        return snake_string
    else:
        words = snake_string.split('_')
        cap_words = [x.capitalize() for x in words]
        return ''.join(cap_words)

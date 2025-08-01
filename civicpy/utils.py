UNMARKED_PLURALS = {'evidence'}

def pluralize(string):
    if string == 'therapy':
        return 'therapies'
    if string in UNMARKED_PLURALS:
        return '{}_items'.format(string)
    if string.endswith('s'):
        return string
    return string + 's'


def singularize(string):
    string = string.rstrip('s')
    if string == 'evidence_item':
        string = 'evidence'
    elif string == 'therapie':
        string = 'therapy'
    return string


def snake_to_camel(snake_string):
    words = snake_string.split('_')
    cap_words = [x.capitalize() for x in words]
    return ''.join(cap_words)

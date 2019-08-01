"""

F90wrap Plugin for QUIP documentation generation

Takes the lines of a subroutine, extracts the argument information from the param_register calls and makes a table to
be added to the docstring.
Adds nothing if there are no param_register instances found

Takes:
subroutine_lines: list of strings, as the lines of the subroutine
doc: the docstring to append to (out.doc in parser.py)

"""


import re

args_str_re = re.compile(r"""^\s*call param_register\([a-zA-Z][a-zA-Z0-9_%]*\s*,\s*(['"])([a-zA-Z][a-zA-Z0-9_]*)\1\s*,\s*(['"]?)(.*?)\3\s*,\s*([a-zA-Z_][a-zA-Z0-9_%]*)\s*,.+?help_string=(['"])(.+?)\6\)""")


def find_params(lines):
    """
    Finds the parameters line by line
    """

    spec = []

    for l in lines:
        m = args_str_re.search(l)
        if m:
            arg_data = dict(name=m.group(2),
                            value=m.group(4),
                            var=m.group(5).lower(),
                            doc=m.group(7))

            arg_data['type'] = infer_type(arg_data['value'], arg_data['var'])
            spec.append(arg_data)

    if spec:
        return spec
    else:
        # so there was nothing added ot it
        return None


def magic_table(spec):

    if len(spec) == 0:
        print('oh shit, zero len')
        return None

    name_list = ['Name']
    type_list = ['Type']
    value_list = ['Value']
    # var_list = ['Var']
    doc_list = ['Doc']

    for arg in spec:
        name_list.append(arg['name'])
        type_list.append(arg['type'])
        value_list.append(arg['value'])
        # var_list.append(arg['var'])
        doc_list.append(arg['doc'])

    max_name_len = max(len(name) for name in name_list)
    max_type_len = max(len(typ) for typ in type_list)
    max_value_len = max(len(value) for value in value_list)
    # max_var_len = max(len(var) for var in var_list)

    # cols = (max_name_len, max_value_len, max_var_len, 40)
    cols = (max_name_len, max_type_len, max_value_len, 40)

    args_str_lines = ['.. rubric:: args_str options', '']
    fmt = "%-{:d}s %-{:d}s %-{:d}s %-{:d}s".format(*cols)

    for i, (name, type_, default, doc) in enumerate(zip(name_list, type_list, value_list, doc_list)):
        if i == 0:
            args_str_lines.append(fmt % ('=' * cols[0], '=' * cols[1], '=' * cols[2], '=' * cols[3]))
        doc_words = doc.split()
        while doc_words:
            doc_line = ''
            while doc_words:
                word = doc_words.pop(0)
                if len(doc_line) + 1 + len(word) > cols[3]:
                    doc_words.insert(0, word)  # put it back
                    break
                else:
                    doc_line = doc_line + ' ' + word

            args_str_lines.append(fmt % (name, type_, default, doc_line.strip()))
            name = type_ = default = ''
        if i == 0 or i == len(name_list) - 1:
            args_str_lines.append(fmt % ('=' * cols[0], '=' * cols[1], '=' * cols[2], '=' * cols[3]))

    args_str_lines.extend(['', ''])

    return [line + '\n' for line in args_str_lines]


def infer_type(value, variable=None):
    """
    Tries to infer the type of a variable, based on the default value if given.

    """

    if value in ['T', 'F']:
        return 'bool'

    try:
        throwaway = int(value)
        return 'int'
    except ValueError:
        pass

    try:
        throwaway = float(value)
        return 'float'
    except ValueError:
        pass

    return 'None'


def doc_plugin(subroutine_lines, name, run_type=None):
    """
    F90wrap Plugin for QUIP documentation generation

    Takes the lines of a subroutine, extracts the argument information from the param_register calls and makes a table to
    be added to the docstring.
    Adds nothing if there are no param_register instances found

    Takes:
    subroutine_lines: list of strings, as the lines of the subroutine
    doc: the docstring to append to (out.doc in parser.py)

    """

    spec = find_params(subroutine_lines)
    if spec is None:
        print('QUIP_doc_plugin: No args found in ', name)
        table_string = []
    else:
        table_string = magic_table(spec)
        # doc.append(table_string)
        print('QUIP_doc_plugin: args found in {}s. Table added to doc as follows'.format(name))
        print(table_string)

    return table_string

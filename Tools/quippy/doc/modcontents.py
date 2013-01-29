import sphinx

if sphinx.__version__ < '1.0.1':
    raise RuntimeError("Sphinx 1.0.1 or newer is required")

import inspect
import pydoc

def process_docstring(app, what, name, obj, options, lines):
    if what == 'module':
        print 'Adding contents listing for module %s' % name

        title = 'Module contents for :mod:`%s`:' % name
        lines.append(title)

        classes = module_classes(obj)
        if classes:
            lines.append("""
.. rubric:: Classes

.. autosummary::
    
%s

""" % '\n'.join(['    %s' % cls for cls in classes]))

        funcs = module_functions(obj)
        if funcs:
            lines.append("""
.. rubric:: Functions

.. autosummary::

%s
            
""" % '\n'.join(['    %s' % func for func in funcs]))

        # FIXME add documentation for module attributes
        attributes = module_attributes(obj)
        if attributes:
            lines.append("""
.. rubric:: Attributes

%s

""" % attributes_table(obj, attributes))



def attributes_table(mod, attributes):
    MAX_LEN=100

    lines = ['']

    values = [ str(getattr(mod, attr)) for attr in attributes ]
    for i, value in enumerate(values):
        if len(value) > MAX_LEN:
            values[i] = '---'

    attributes = [':attr:`%s`' % attr for attr in attributes ]

    col1 = max(len(attr) for attr in attributes)
    col2 = max(len(value) for value in values)
    
    fmt = '%%-%ds %%-%ds' % (col1, col2)
    head = fmt % ('='*col1, '='*col2)

    lines.append(head)
    lines.append(fmt % ('Name', 'Value'))
    lines.append(head)
    for attr, value in zip(attributes, values):
        lines.append(fmt % (attr, value))
        
    lines.append(head)
    lines.append('')
    return '\n'.join(lines)


def module_functions(mod):
    if hasattr(mod, '__alldoc__'):
        allsymbols = mod.__alldoc__
    elif hasattr(mod, '__all__'):
        allsymbols = mod.__all__
    else:
        allsymbols = [name for (name, obj) in inspect.getmembers(mod)]
    return [name for name in allsymbols if
            not name.startswith('_') and
            inspect.isfunction(getattr(mod, name)) and
            pydoc.getdoc(getattr(mod, name))]

def module_classes(mod):
    if hasattr(mod, '__alldoc__'):
        allsymbols = mod.__alldoc__
    elif hasattr(mod, '__all__'):
        allsymbols = mod.__all__
    else:
        allsymbols = [name for (name, obj) in inspect.getmembers(mod)]
    return [name for name in allsymbols if
            not name.startswith('_') and
            inspect.isclass(getattr(mod, name)) and
            pydoc.getdoc(getattr(mod, name))]

def module_attributes(mod):
    if hasattr(mod, '__alldoc__'):
        allsymbols = mod.__alldoc__
    elif hasattr(mod, '__all__'):
        allsymbols = mod.__all__
    else:
        allsymbols = [name for (name, obj) in inspect.getmembers(mod)]
    return [name for name in allsymbols if
            not name.startswith('_') and
            not callable(getattr(mod, name)) and
            pydoc.getdoc(getattr(mod, name))]


def setup(app):
    print('modcontents extension loading')
    app.connect('autodoc-process-docstring', process_docstring)

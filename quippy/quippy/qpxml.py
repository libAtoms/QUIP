"""
quippy.xml

Functions to manipulate GAP xml files.

Works best if lxml is installed, but will fall back to the
standard library module.

Available in script form as combine_gap_xml.py

"""

import re
from collections import defaultdict
from copy import deepcopy

# Needs extra options for lxml to allow it to be
# global and overridden
parser = None

try:
    import lxml.etree as ElementTree
    parser = ElementTree.XMLParser(huge_tree=True, strip_cdata=False)
except ImportError:
    from xml.etree import ElementTree


def get_gap_version(tree):
    """Extract the ``gap_version`` from an xml tree
    and return as an integer. Will return the maximum
    if ``gap_version`` appears more than once.

    Parameters
    ----------
    tree : Element
        An ElementTree of a GAP xml
    """
    # Use xpath specification to get all occurances
    # Pull only values in a 'GAP_params'
    return max([int(elem.attrib['gap_version']) for elem in
                tree.findall('./GAP_params/[@gap_version]')])


def set_max_gap_version(base, extras):
    """Set the ``gap_version`` attribute in ``base`` to be the
    maximum of all ``gap_version`` found. Modify ``base`` in place
    and return nothing.

    Parameters
    ----------
    base : Element
        The ElementTree to be updated with maximum gap_version
    extras : list of Element
        All other ElementTrees to be scanned for gap_version

    """
    max_gap_version = max([get_gap_version(tree)
                           for tree in [base] + list(extras)])

    for elem in base.findall('./GAP_params/[@gap_version]'):
        elem.set('gap_version', str(max_gap_version))


def sum_e0(base, extras):
    """Update the values of ``e0`` in the base tree with the sum
    of all values found in the given potentials. Modify ``base`` in
    place and return nothing.

    Parameters
    ----------
    base : Element
        The ElementTree to be modified with the sum of ``e0`` values
    extras : list of Element
        All other ElementTrees to added to the summed ``e0`` values

    """
    # Start at 0.0 for uninitialised values
    total = defaultdict(float)
    # add values for the potentials being combined
    for tree in list(extras) + [base]:
        for elem in tree.findall('./GAP_params/GAP_data/e0'):
            total[elem.attrib['Z']] += float(elem.attrib['value'])

    # pop the values so we can tell if any are unused
    for e0 in base.findall('./GAP_params/GAP_data/e0'):
        e0.attrib['value'] = '{:.21g}'.format(total.pop(e0.attrib['Z']))

    # create new elements for anything not already present in base
    for z, e0 in total.items():
        ElementTree.SubElement(
            base.find('./GAP_params/GAP_data'),  # End of GAP data
            'e0', attrib={'Z': str(z), 'value': '{:.21g}'.format(e0)})


def merge_descriptors(base, extras, remove_xyz=True, label=None):
    """Extract all GAP descriptors from the potentials in extras and
    merge them into base, updating the lables and number of
    coordinates. Modify ``base`` in place and return nothing.

    Parameters
    ----------
    base : Element
        The ElementTree in which all the descriptors will be combined
    extras : list of Element
        All other ElementTrees to extract descriptors from
    remove_xyz : bool
        If True, remove any ``XYZ_data`` elements; If False then leave
        anything found in base (XYZ_data is not combined)
    label : str or None
        Use the string as the label for the new potential, otherwise
        use the label from base.

    """
    # Make sure the labels for everything are aligned
    if label is None:
        # Use the base as the basis for everything
        label = base.find('GAP_params').attrib['label']
    else:
        # Set custom labels in the Potential so it will properly init
        base.find('Potential').attrib['label'] = label
        base.find('Potential').attrib['init_args'] = re.sub(
            r'label\s*=\s*\w*', 'label={}'.format(label),
            base.find('Potential').attrib['init_args'])
        base.find('GAP_params').attrib['label'] = label

    if remove_xyz:
        # getparent only avaiable in lxml so find all parents first
        for xyz_data_parent in base.findall('GAP_params/XYZ_data/..'):
            for xyz_data in xyz_data_parent.findall('XYZ_data'):
                xyz_data_parent.remove(xyz_data)

    base_sparse = base.find('GAP_params/gpSparse')

    for extra in extras:
        for new_coordinate in extra.findall(
                'GAP_params/gpSparse/gpCoordinates'):
            add_coordinate = deepcopy(new_coordinate)
            base_sparse.append(add_coordinate)

    # update labels and total coordinates;
    # ensures labels are the same even with custom label
    # 1-based indexing
    n_coordinate = 0  # in case there are none!
    for n_coordinate, coordinate in enumerate(
            base_sparse.findall('gpCoordinates'), 1):
        coordinate.attrib['label'] = '{}{}'.format(label, n_coordinate)

    base_sparse.attrib['label'] = label
    base_sparse.attrib['n_coordinate'] = str(n_coordinate)


def combine_xml(base_filename, extra_filenames, remove_xyz=True, label=None):
    """Merge the GAP potentials in the xml files given in the arguments.
    Return a string containing the combined potential.

    Use this method to either write a new xml file or pass directly into
    ``quippy.Potential``:

        pot = Potential("IP GAP",
                        param_str=combine_xml('pot1.xml', ['pot2.xml']))

    Parameters
    ----------
    base_filename : str
        Path to xml file to use as the basis for combined xml file
    extra_filenames : list of str
        Paths to xml files to combine into base
    remove_xyz : bool
        If True, any XYZ_data will be removed from the output. If lxml is not
        avaiable XYZ_data will lose CDATA status and slow parsing. There is
        also no effort made to combine XYZ_data from the extra xml files
    label : str or None
        If not None, the string is used to label the combined potential in
        the output, otherwise the label from base is used.

    Returns
    -------
    combined_xml : str
        String containing the combined xml contents

    """

    base_xml = ElementTree.parse(base_filename, parser).getroot()
    extra_xml = [ElementTree.parse(extra, parser).getroot() for extra in
                 extra_filenames]

    # set max gap_version
    set_max_gap_version(base_xml, extra_xml)
    # sum all the zeros
    sum_e0(base_xml, extra_xml)
    # combine everything into base
    merge_descriptors(base_xml, extra_xml, remove_xyz=remove_xyz, label=label)
    # Done, send it back
    return ElementTree.tostring(base_xml)

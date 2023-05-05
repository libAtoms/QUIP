from ase.calculators.calculator import all_changes
import numpy as np
from xml.etree.ElementTree import parse, fromstring, tostring, ElementTree
import re
import os
from quippy.descriptors import Descriptor
from quippy.potential import Potential
from ase.data import chemical_symbols
from scipy.stats import norm
from scipy.linalg import solve_triangular
from ase.calculators.mixing import AverageCalculator


class DescXMLWrapper():
    '''
    Small shell class for storing key descriptor info from a GAP xml read

    '''
    _Z_regex = "(Z|z)[1-9]*\s?=\s?([1-9]+)"  # RegEx to search command line for "Z=, Z1 = or z2= style args"
    _Z_regex = re.compile(_Z_regex)

    # RegEx to find the dot_product exponents from command line args
    _exponent_regex = "exponents\s?=\s?.((\s?(-\d+)\s?)+)."
    _exponent_regex = re.compile(_exponent_regex)

    def __init__(self, desc_xml):
        self.weights = np.array([float(child.attrib["alpha"])
                                for child in desc_xml if "sparseX" in child.tag])
        self.sparse_cuts = np.array([float(
            child.attrib["sparseCutoff"]) for child in desc_xml if "sparseX" in child.tag])

        self.nsparse = self.weights.shape[0]

        self.delta = float(desc_xml.attrib["signal_variance"])
        self._int_cov_type = int(desc_xml.attrib["covariance_type"])

        if self._int_cov_type == 1:  # ard_se
            self.cov_prop = float([child.text.split()[0]
                                  for child in desc_xml if "theta" in child.tag][0])
            self.cov_type = "ard_se"
        elif self._int_cov_type == 2:  # dot_product
            self.cov_prop = float(desc_xml.attrib["zeta"])
            self.cov_type = "dot_product"
        else:
            self.cov_prop = []
            self.cov_type = "UNKNOWN"

        _desc_children = [
            child for child in desc_xml if "descriptor" in child.tag]

        self._cmd = _desc_children[0].text

        self.desc_type = self._cmd.split(" ")[0]

        Zs = self._Z_regex.findall(self._cmd)

        self.Zs = [Z[1] for Z in Zs]

        self.name = "".join([chemical_symbols[int(Z)]
                            for Z in self.Zs]) + " " + self.desc_type

        self.quip_desc = Descriptor(self._cmd)


class GAPXMLWrapper():
    '''
    Small shell class to store key GAP information
    '''

    def __init__(self, xml, mean_weights=None):
        '''
        Read XML ETree for a GAP potential

        Useful Attributes:
        self.gap_label : Label of the GAP potential
        self.isolated_atom_energies : dict of E0s; E0 for species Z is self.isolated_atom_energies[Z]
        self.num_desc : Number of descriptors
        '''
        self._xml_tree = xml
        root = self._xml_tree.getroot()

        self.gap_label = root[0].attrib["label"]

        iso_atom_energies = root[1][0][:]

        self.isolated_atom_energies = {}

        for item in iso_atom_energies:
            self.isolated_atom_energies[int(item.attrib["Z"])] = float(
                item.attrib["value"])

        descriptors = root[1][1][:]

        self.num_desc = len(descriptors)

        self.descriptors = [DescXMLWrapper(desc_xml)
                            for desc_xml in descriptors]

        self.weights = np.block([desc.weights for desc in self.descriptors])

        self.total_nsparse = self.weights.shape[0]

        if mean_weights is not None:
            self.mean_weights = mean_weights
        else:
            self.mean_weights = self.weights.copy()

    def save(self, fname):
        '''
        Save internal XML tree back to file
        '''

        with open(fname, "wb") as f:
            self._xml_tree.write(f)

    def as_potential(self):
        pot = Potential(param_str=tostring(self._xml_tree.getroot()))
        pot.xml = self
        return pot

    def _posterior_sample(self):

        if self.R is not None:
            z = norm.rvs(size=self.total_nsparse)

            return self.mean_weights + solve_triangular(self.R, z)
        else:
            raise FileNotFoundError(
                f"R matrix not found. {self.R_fname} does not exist")

    def _xml_sample(self):
        new_weights = self._posterior_sample()

        new_root = fromstring(tostring(self._xml_tree.getroot()))

        new_tree = ElementTree(new_root)

        new_descs = new_root[1][1][:]

        isparse = 0

        for i in range(self.num_desc):

            desc_nsparse = self.descriptors[i].nsparse

            desc_wts = [child for child in new_descs[i]
                        if "sparseX" in child.tag]

            for j in range(desc_nsparse):
                desc_wts[j].attrib["alpha"] = str(new_weights[isparse])
                isparse += 1

        XYZ_data = [child for child in new_root[1] if child.tag == "XYZ_data"]

        if len(XYZ_data):  # Parent xml has XYZ data
            for dat in XYZ_data:
                new_root[1].remove(dat)

        return new_tree

    def draw_posterior_samples(self, num_samples=1):
        if num_samples == 1:
            return GAPXMLWrapper(self._xml_sample(), mean_weights=self.mean_weights)
        else:
            return [self.draw_posterior_samples() for i in range(num_samples)]


def read_xml(path_to_xml):

    gap = GAPXMLWrapper(parse(path_to_xml))

    R_fname = path_to_xml + ".R." + gap.gap_label

    if os.path.exists(R_fname):
        gap.R = np.loadtxt(R_fname).reshape(
            (gap.total_nsparse, gap.total_nsparse)).T
    else:
        gap.R = None

    return gap


def get_xml_committee(path_to_xml, num_committors, return_core_wrapper=False):
    gap_wrapper = read_xml(path_to_xml)

    committee = gap_wrapper.draw_posterior_samples(num_committors)

    if return_core_wrapper:
        return committee, gap_wrapper
    else:
        return committee


def get_calc_committee(path_to_xml, num_committors, return_core_wrapper=False):
    xml_committee, gap_wrapper = get_xml_committee(
        path_to_xml, num_committors, return_core_wrapper=True)
    calc_committee = [committor.as_potential() for committor in xml_committee]

    if return_core_wrapper:
        return calc_committee, gap_wrapper
    else:
        return calc_committee

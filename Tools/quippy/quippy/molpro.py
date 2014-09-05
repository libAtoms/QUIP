#
#Molpro class - adapted from Castep class by James Kermode
#Alan Nichol
#

import sys, string, os, operator, itertools, logging, glob, re, subprocess
import numpy as np

from quippy.atoms import Atoms
from io import AtomsReaders, AtomsWriters, atoms_reader
from quippy.dictionary import Dictionary
from quippy.units import AU_FS, HARTREE, BOHR, BOLTZMANN_K, GPA, DEBYE
from quippy.periodictable import atomic_number
from quippy.atoms import make_lattice, get_lattice_params

from ordereddict import OrderedDict
from farray import *
from math import pi
import xml.dom.minidom as minidom
from HTMLParser import HTMLParser

#list of methods ( I think!?)
__all__ = ['MolproDatafile']

#will have to include units mapping

# this list is not complete and there is some crap in it as well - may have to rethink this.
valid_datafile_keywords = ['ANGSTROM','SYMMETRY','ORIENT','GEOMTYP','***', '--', 'ACCURACY', 'ACPF', 'ACTIVE', 'ADD', 'AIMS', 'ALTERN', 'ANGULAR', 'AOINT', 'AQCC', 'BASIS', 'BMAT', 'BRUECKNER', 'CANONICAL', 'CANORB', 'CASPROJ', 'CASVB', 'CCSD', 'CCSD(T)', 'CEPA', 'CHARGE', 'CHECK', 'CI', 'CI-PRO', 'CIGUESS', 'CIS', 'CISD', 'CIWEIGHTS', 'CLEAR', 'CLEARALL', 'CLOSED', 'COEFFS', 'COMPRESS', 'CON', 'CONFIG', 'CONICAL', 'COORD', 'COPT', 'CORE', 'CPMCSCF', 'CPP', 'CRIT', 'CUT', 'DATA', 'DDR', 'DELETE', 'DELOCAL', 'DELSTRUC', 'DEMC', 'DENSITY', 'DF-RKS', 'DF-UKS', 'DFTBLOCK', 'DFTDUMP', 'DFTFACTOR', 'DFTTHRESH', 'Difference gradients', 'DH', 'DIIS', 'DIP', 'DIP+', 'DIRECT', 'DM', 'DMA', 'DO', 'DONT', 'DUMMY', 'Q,X', ')', 'DUMP', 'ECP', 'ELSEIF', 'ENDDO', 'ENDIF', 'ENDZ', 'EOM', 'EOMPAR', 'EOMPRINT', 'ERASE', 'EXCHANGE', 'EXPEC', 'EXPEC2', 'EXTRA', 'FCI', 'FIELD', 'FIELD+', 'FILE', 'FIXORB', 'FIXSTRUC', 'FOCK', 'FORCE', 'FREEZE', 'FREQUENCIES', 'FROZEN', 'FULL', 'G1', 'GDIRECT', 'GENERAL', 'GEOMETRY', 'GEXPEC', 'GOTO', 'GPARAM', 'GPRINT', 'GRADTYP', 'GRID', 'GRIDPRINT', 'GRIDSAVE', 'GRIDSYM', 'GRIDTHRESH', 'GROUP', 'GTHRESH', 'GUESS', 'HESS', 'HESSELEM', 'HESSIAN', 'HF', ',', 'options', 'HF-SCF', 'IF', 'INACTIVE', 'INCLUDE', 'INDIVIDUAL', 'INIT', 'INSTANTON', 'INTOPT', 'IPOL', 'IPRINT', 'IRREPS', 'ITERATIONS', 'KS', 'KS-SCF', 'LABEL', 'LATTICE', 'LIMIT', 'LINEAR', 'LINESEARCH', 'LOCAL', 'LOCALI', 'LOCAO', 'LOCORB', 'LQUANT', 'MASS', 'MATROP', 'MAXDAV', 'MAXITER', 'MCSCF', 'MEMORY', 'MERGE', 'METHOD', 'MOLDEN', 'molpro', 'MOVE', 'MULTI', 'NACM', 'NATORB', 'NBO', 'NELEC', 'NOCASPROJ', 'NOCHECK', 'NOENEST', 'NOEXC', 'NOEXTRA', 'NOGPRINT', 'NOGRIDSAVE', 'NOGRIDSYM', 'NONLINEAR', 'NONUCLEAR', 'NOORDER', 'NOPAIR', 'NOSINGLE', 'NOSYMPROJ', 'NUMERICAL', 'NUMHES', 'OCC', 'OFFDIAG', 'OFFSET', 'OPEN', 'OPTG', 'OPTIM', 'OPTION', 'ORB', 'ORBIT', 'ORBITAL', 'ORBPERM', 'ORBPRINT', 'ORBREL', 'ORTH', 'ORTHCON', 'PAIR', 'PAIRS', 'PARAM', 'POLARIZABILITY', 'POLY', 'POP', 'POTENTIAL', 'PRINT', 'PROC', 'PROJECT', 'PROPERTY', 'PSPACE', 'PUNCH', 'PUT', 'QCI', 'QUAD', 'QUAD+', 'RADIAL', 'RADIUS', 'RANGEHYBRID', 'READ', 'READPUN', 'READVAR', 'REF', 'REFSTATE', 'REL', 'RELAX', 'RESTART', 'RESTRICT', 'RHF', 'RHF-SCF', 'RKS', 'RKS-SCF', 'ROOT', 'ROTATE', 'ROTATEA', 'RS2', 'RS2C', 'RS3', 'SADDLE', 'SAMC', 'SAVE', 'SCALE', 'SCHMIDT', 'SCORR', 'SELECT', 'SERVICE', 'SET', 'SHIFT', 'SHOW', 'SPECIAL', 'SPIN', 'SPINBASIS', 'START', 'STATE', 'STATUS', 'STEP', 'STRONG', 'STRUC', 'SURF', 'SYM', 'SYMELM', 'WF', ' card', 'SYMPROJ', 'TABLE', 'TEST', 'THERMO', 'THRESH', 'TRAN', 'TRAN2', 'TRANH', 'TRANS', 'TRNINT', 'TRUST', 'UHF', 'UHF-SCF', 'UKS', 'UKS-SCF', 'UNCOMPRESS', 'UPDATE', 'VARIABLE', 'VB', 'VBDUMP', 'VBWEIGHTS', 'VCI', 'VMP2', 'VORONOI', 'VSCF', 'WEIGHT', 'WF', 'WRITE', 'ZMAT', 'Program control:', '***', 'MEMORY', 'PUNCH', 'FILE', 'RESTART', 'INCLUDE', 'BASIS', 'GEOMETRY', 'ZMAT', 'PARALLEL', 'STATUS', 'PRINT', ',', 'GPRINT', 'THRESH', ',', 'GTHRESH', 'DIRECT', ',', 'GDIRECT', 'EXPEC', ',', 'GEXPEC', 'TEXT', 'EXIT', 'DO', 'ENDDO', 'IF', 'ELSEIF', 'ENDIF', 'IF block', 'GOTO', 'LABEL', 'DATA', 'DELETE', ', ', '      ERASE', 'MATROP', 'GRID', 'CUBE', 'CARTESIAN', 'SPHERICAL', 'USER', '--', 'Variables:', 'SET', 'SETI', 'SETA', 'CLEAR', 'CLEARALL', 'GETVAR', 'SHOW', 'TABLE', 'Wave function optimization:', 'INT', 'LSINT', 'SORT', 'CPP', 'HF', ', ', 'RHF', ', ', 'HF-SCF', ', or ', 'RHF-SCF', 'UHF', ' or ', 'UHF-SCF', 'DFT', 'KS', ', ', 'RKS', 'UKS', 'MULTI', ', ', 'MCSCF', ', or ', 'CASSCF', 'CASVB', 'CI', ', ', 'MRCI', ', or ', 'CI-PRO', 'CIPT2', 'ACPF', ', ', 'AQCC', 'CEPA', 'RS2', ', ', 'RS3', 'RS2C', 'MP2', 'MP3', 'MP4', 'CISD', 'CCSD', 'BCCD', 'QCI', ',', 'QCSID', 'UCCSD', 'RCCSD', 'FCI', ' or ', 'FULLCI', 'Local correlation methods:', 'LMP2', 'LMP3', 'LMP4', 'LCISD', 'LCCSD', 'Explicitly correlated methods:', 'DF-MP2-R12', 'DF-MP2-F12', 'DF-LMP2-R12', 'DF-LMP2-F12', 'Orbital manipulation:', 'LOCALI', 'MERGE', 'Properties and wavefunction analysis:', 'POP', 'DMA', 'PROPERTY', 'DIP', 'QUAD', 'LATTICE', 'Gradients and geometry optimization:', 'FORCES', 'OPTG', 'MIN', 'PUT', 'HESSIAN', 'FREQUENCY', 'MASS', 'DDR', 'QCISD, CCSD, LQCISD, LCCSD', ' can be appended by ', '(T)', '\nand then a perturbative correction for triple excitations will be computed (e.g., ', 'CCSD(T)', ').\n\n', 'HF', ', ', 'KS', ', ', 'MP2', ' and all local correlation methods can be prepended by ', 'DF-', ' to invoke density fitting.\n\n']


class KeyWordGetter(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.keywords = []
        self.is_keyword = False
        self.ignore_these=[',','\n',' or ']
        
    def handle_starttag(self,tag,attributes):
        if tag == 'tt' or tag == 'em': #HTMLParser converts tags to lower case
            self.is_keyword=True
        else:
            self.is_keyword=False
            
    def handle_data(self, data):
        if self.is_keyword == True:
            if data not in self.ignore_these:
                self.keywords.append(data)
            

class MolproDatafile(OrderedDict):
    """Class to wrap a molpro datafile"""

    def __init__(self, datafile=None,xml=None,atoms=None):
        OrderedDict.__init__(self)
        if datafile is not None:
            self.read(datafile)
        elif xml is not None:
            self.read_xml(xml)
        elif atoms is not None:
            self.update_from_atoms(atoms)

    def copy(self):
        new = MolproDatafile()
        new.update(self) #this is inherited form OrderedDict
        return new

    def parse_line(self, line, key=None):

        # split a line by the first non-alphanumeric character, then cats that character
        # to start of second item in list, if there is one.
        if key is not None:
            self[key].append(line)
        else:
            nonalpha = re.compile('[^a-zA-Z0-9()\-]')
            separator = re.findall(nonalpha,line)

            if len(separator) > 0:
                separator = separator[0]
            else:
                separator = ","
            fields = line.split(separator,1)
            key = fields[0].upper()
            if key not in self._keys:
                self[key] = []
                if len(fields) > 1 and fields[1] != '':
                    self[key].append(separator+fields[1])
            else:
                if key == "BASIS": # this warning should be passed through logging
                    print('WARNING: CHANGING BASIS DOES UNPREDICTABLE THINGS')
                n=2
                testkey = key
                while testkey in self._keys: #find a key like charge#3 if charge and charge#2 are already taken
                    testkey = key+"#%s" % n
                    n+=1
                key = testkey
                self[key] = []
                if len(fields) > 1:
                        self[key].append(separator+fields[1])

    def read(self, datafile):
        if operator.isMappingType(datafile):
            self.update(datafile)
            return
        
        # If given a file first make a list containing the lines
        if type(datafile) == type(''):
            data = open(datafile,'r')
            datalines = data.readlines()
        elif type(datafile) == type([]):
            datalines = datafile

        sub_lines = []
        current_key = None
        for compound_line in datalines:
            #molpro allows multiple statements on single line separated by ';' , so split these up
            lines = compound_line.split(';')
            for sub_line in lines:
                sub_lines.append(sub_line.strip())

        for line in sub_lines:
            # Skip comments and blank lines, NB molpro files may begin with  the line'***'
            if line.startswith('#') or line.startswith('!') or line.startswith('*') or line == '':
                continue

            # Check if any brackets in line
            open_bracket = re.search('{',line)
            close_bracket = re.search('}',line)

            if open_bracket and close_bracket:
                # superfluous brackets, no actual multi-line command or data
                line = line.replace("{","")
                line = line.replace("}","")
                open_bracket = False
                close_bracket = False
                self.parse_line(line.strip(),current_key)

            #check if command block starting
            elif open_bracket:
                if current_key == None:
                    line = re.sub('{',"",line).strip()
                    if line != "":
                        self.parse_line(line,current_key)
                        current_key = self._keys[-1]
                    else:
                        raise ValueError("Parse error in datafile: standalone open bracket")
                else:
                    raise ValueError("Parse error in datafile: nesting of curly {} brackets is illegal")

            #check if end of command block reached
            elif close_bracket:
                if current_key == None:
                    raise ValueError("Parse error in datafile:  check pairing of curly {} brackets")
                else:
                    line = re.sub('}',"",line).strip()
                    if line != "":
                        self.parse_line(line,current_key)
                        current_key = None
                    else:
                        current_key = None
                        continue
            # normal line - no brackets
            else : 
                self.parse_line(line,current_key)


    def read_from_molpro_output(self, molpro_output):
        """Read the input file from molpro output. Input should be filename, file-like object or list of lines"""
        if type(molpro_output) == type(''):
            molpro_output = open(molpro_output, 'r')
            molpro_output = molpro_output.readlines()
        elif hasattr(molpro_output, 'read'):
            molpro_output = molpro_output.readlines()

        # Remove newlines from end of each line in molpro_output
        molpro_output = [line.strip() for line in molpro_output]

        # Find the echo of the datafile in the molpro output
        try:
            datafile_start = molpro_output.index('default implementation of scratch files=df')
        except ValueError:
            raise ValueError('Unable to find echo of datafile in molpro output')

        datafile_lines = []
        i = datafile_start + 2 # skip over the blank line produced by molpro
        
        # If geometry is specified as a filename molpro will insert a comment which we should remove
        include = re.compile("Including file")
        end_of_input = re.compile("Variables initialized")
        datafile_ended = False
        while datafile_ended == False:
            line = molpro_output[i]
            i=i+1
            if re.search("Variables initialized",line):
                datafile_ended = True
            elif re.search("Including file",line):
                continue
            elif line.strip() == '':
                continue # skip arbitrary blank lines
            else:
                datafile_lines.append(line)
            
        self.read(datafile_lines)

    def read_xml(self, xml_output):
        """Read the input file from molpro output. Input should be filename, file-like object or list of lines"""
        if type(xml_output) == type(''):
            xml_output = open(xml_output, 'r')
            xml_output = xml_output.readlines()
        elif hasattr(xml_output, 'read'):
            xml_output = xml_output.readlines()

        # Remove newlines from end of each line in molpro_output
        xml_output = [line.strip() for line in xml_output]

        # Find the echo of the datafile in the molpro output
        try:
            datafile_start = xml_output.index('--><job>')
        except ValueError:
            raise ValueError('Unable to find echo of datafile in molpro output')

        datafile_lines = []
        i = datafile_start + 2 # skip over the comment indicator produced by molpro
        
        end_of_input = re.compile("Variables initialized")
        datafile_ended = False
        while datafile_ended == False:
            line = xml_output[i]
            i=i+1
            if re.search("Variables initialized",line):
                datafile_ended = True
            elif line.strip() == '':
                continue # skip arbitrary blank lines
            else:
                datafile_lines.append(line)
            
        self.read(datafile_lines)

    def write(self,datafile=sys.stdout):
        "Write molpro input file. datafile can be a filename or an open file"

        if type(datafile) == type(''):
            datafile = open(datafile,'w')

        for key, value in self.iteritems():
            #iteritems important here because order of lines matters
            #if have multiple instances of a command, say, 'hf' or 'charge'
            #the n occurrences of that keyword after the first will have #n appended
            # e.g. hf, hf#2, hf#3, etc.
            if re.search('#',key):
                shortkey = key.split('#')[0]
            else:
                shortkey = key
            if len(value) > 1:
                datafile.write(shortkey+'={\n')
                for line in value:
                    datafile.write(line+'\n')
                datafile.write('}\n')
            elif value != []:
                datafile.write('%s%s\n' % (shortkey, value[0]))
            else:
                datafile.write(shortkey+'\n')

    def to_atoms(self):
        #check if necessary input there & all makes sense
        if self.has_key('GEOMETRY'):
            block = self['geometry=']
            
        #now create
        atoms = Atoms(n=len(block))
        field_list = [line.strip() for line in block]
        field_list = [line.split() for line in block]
        
        #Way this currently works any labels will be lost
        #if want to use these for molecule specification
        #will have to do something more clever
        #label = re.compile('[0-9]')
        #field_list = map(label.split, field_list[0])
        elements = map(operator.itemgetter(0), field_list)

        #Look up elements by atomic number
        elements = [ not el.isdigit() and atomic_number(el) or el for el in elements ]
        
        #Transfer positions to Atoms object
        # Set the element and pos data
        atoms.set_atoms(elements) #Elements still needs to be defined, farray is a function
        atoms.pos[:,:] = farray([ [float(x) for x in row] \
                                  for row in [field[1:4] for field in field_list]]).T
        return atoms

    def update_from_atoms(self, at, geomfile=None):
        #As this stands the atomic positions aren't actually updated, which is silly
        #Need to decide whether there is any point in allowing geometry spec inside datafile
        for p in at.params:
            if p.upper() in valid_datafile_keywords:
                self[p] = at.params[p]
        self['geom'] = []

        if geomfile is not None:
            self['geom'].append("=%s" %geomfile)
        else:
            self['geometry'] = []
            for i in frange(at.n):
                self['geometry'].append(at.species[:,i].stripstrings() +' %f %f %f' % tuple(np.dot(at.g,at.pos[:,i])))


def read_xml_output(xmlfile,energy_from=None, extract_forces=False, extract_dipole=False, datafile=None, cluster=None):
    #parse an xml output file and return cluster with updated info
    # datafile tells which energies, forces to look for, cluster Atoms object which gets returned, this is echoed in the xml file so can be left out
    # If extract_forces is not given and the FORCE keyword is found in datafile, the default is to set extract_forces=True

    log = logging.getLogger('molpro_driver')
    
    if datafile is None:
        datafile=MolproDatafile(xml=xmlfile)
        if 'FORCE' in datafile:
            extract_forces=True

    energy_names = OrderedDict()
    energy_names['CCSD(T)-F12'] = ["total energy"]
    energy_names['CCSD(T)'] = ["total energy"]
    energy_names['MP2'] = ["total energy"]
    energy_names['RKS'] = ["Energy"]
    energy_names['RHF'] = ["Energy"]
    energy_names['HF'] = ["Energy"]
    #etc
    
    gradient_names = OrderedDict()
    gradient_names['CCSD(T)'] =[""]
    gradient_names['RKS'] =['RKS GRADIENT']
    gradient_names['MP2'] =['MP2 GRADIENT']

    all_methods=OrderedDict()
    all_methods['HF']=["RHF"]
    all_methods['MP2']=["MP2"]
    all_methods['RKS']=["RKS"]
    all_methods['CCSD(T)-F12']=["CCSD(T)-F12a","CCSD(T)-F12b"]
    all_methods['CCSD(T)']=["CCSD(T)"]

    if energy_from is None:
        log.critical("don't know which energy to extract, use keyword energy_from with options "+str([all_methods[k] for k in iter(all_methods)]).replace('[','').replace(']',''))

    #loop through datafile to look for methods.
    calcs=[] #holds the keys for getting correct method, energy_name, gradient_name
    data_keys_upper = [key.upper() for key in datafile._keys]
    for key in all_methods._keys:
       if key in data_keys_upper:
           calcs.append(key)
    dom = minidom.parse(xmlfile)    

    elements=[]
    position_matrix=[]
    cml = dom.documentElement.getElementsByTagName('cml:atomArray')

    for l in cml[0].childNodes:
        if l.nodeType== 1:
            element=l.attributes['elementType'].value.encode('ascii','ignore')
            elements.append(atomic_number(element))
            posx = l.attributes['x3'].value.encode('ascii','ignore')
            posy = l.attributes['y3'].value.encode('ascii','ignore')
            posz = l.attributes['z3'].value.encode('ascii','ignore')
            position_matrix.append([float(posx),float(posy),float(posz)])
    if cluster is None:
        cluster = Atoms(n=len(elements))
        cluster.set_atoms(elements)
        position_matrix=farray(position_matrix).T
        if not 'ANGSTROM' in datafile._keys and not 'angstrom' in datafile._keys:
            position_matrix = position_matrix * (1.0/0.529177249)
        cluster.pos[:,:]=position_matrix
        #note this leaves the lattice undefined

    #now look for each of these energies in xml file
    energy_found=False
    props = dom.documentElement.getElementsByTagName('property')
    for prop in props:
        prop_name = prop.attributes['name'].value.encode('ascii','ignore')
        prop_method = prop.attributes['method'].value.encode('ascii','ignore')
        for calc in calcs:
            if prop_name in energy_names[calc] and prop_method in all_methods[calc]:
                energy_param_name="_".join([prop_method,prop_name])
                energy_param_name=energy_param_name.replace(" ","_")
                #log.info("found "+energy_param_name)
dated routines for finding monomer pairs, triplets in Topology module
                energy_param=prop.attributes['value'].value.encode('ascii','ignore')
                my_energy=energy_param_name
                i_en=1
                while my_energy in cluster.params.iterkeys():
                    i_en+=1
                    my_energy='_'.join([energy_param_name,str(i_en)])
                cluster.params[my_energy] = float(energy_param)*27.21138386
                if prop_method == energy_from:
                    cluster.params['Energy']=float(energy_param)*27.21138386
                    energy_found=True
            elif extract_dipole and prop_name=='Dipole moment':
                dipole_param_name="_".join([prop_method,prop_name])
                dipole_param_name=dipole_param_name.replace(" ","_")
                log.info("found dipole moment: "+dipole_param_name)
                dipole_param=prop.attributes['value'].value.encode('ascii','ignore')
                cluster.params[dipole_param_name]=dipole_param

    if not energy_found:
        log.critical("couldn't find energy from "+energy_from+" prop method : "+prop_method)
                      
        
                
    # read gradients if requested
    if extract_forces:
        if not cluster.has_property('force'):
            cluster.add_property('force', 0.0, n_cols=3)

        grads = dom.documentElement.getElementsByTagName('gradient')
        force_matrix = grads[0].childNodes[0].data.split('\n')
        force_matrix = [str(i).split() for i in force_matrix]
        for i in force_matrix:
            try:
                force_matrix.remove([])
            except ValueError:
                break
        force_matrix = [[(-27.2113961/0.529177249)*float(j) for j in i] for i in force_matrix] # check this negative sign
       
        cluster.force[:] =farray(force_matrix).T

        if len(grads) != 1:
            for k in range(1,len(grads)):
                my_force='force%s'%str(k+1)
                force_matrix = grads[k].childNodes[0].data.split('\n')
                force_matrix = [str(i).split() for i in force_matrix]
                for i in force_matrix:
                    try:
                        force_matrix.remove([])
                    except ValueError:
                        break
                force_matrix = [[(-27.2113961/0.529177249)*float(j) for j in i] for i in force_matrix]
                cluster.add_property(my_force,farray(force_matrix).T)

    return cluster

def run_molpro(datafile, molpro, stem, test_mode=False):
    #Invoke molpro and return true if completed successfully

    log = logging.getLogger('molpro_driver')

    #write datafile
    datafile.write(stem)

    #check command line
    if not '%s' in molpro: molpro = molpro + ' %s'

    if test_mode:
        log.info('test mode: not running molpro')

    else:
        # Remove old output file 
        try:
            os.remove(stem+'.out')
        except:
            pass

        os.system(molpro % stem)

        #    error = subprocess()
    got_error=False

    return not got_error

def get_valid_keywords():
    #Download some pages from the user manual, parse the HTML, and return a list of valid molpro syntax
    #This because not all valid words are supported by the molpro help function

    if os.path.exists('index.html'):
        os.system('rm index.html')
    if os.path.exists('keywords.html'):
        os.system('rm keywords.html')
    subprocess.check_call("curl -L -o 'index.html' 'http://www.molpro.net/info/2012.1/doc/manual/node871.html'",shell=True)
    index = open('index.html','r')
    subprocess.check_call("curl -L -o 'keywords.html' 'http://www.molpro.net/info/2012.1/doc/manual/node37.html'",shell=True)
    keywords = open('keywords.html','r')
    htmlstring = index.read() + keywords.read()
    index.close()
    keywords.close()

    getter = KeyWordGetter()
    getter.feed(htmlstring)
    getter.close()
    return getter.keywords

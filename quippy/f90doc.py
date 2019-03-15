#!/usr/bin/python
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Originally based on:
# f90doc - automatic documentation generator for Fortran 90
# Copyright (C) 2004 Ian Rutt
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA

import re, getopt, sys, string, copy, os.path

major, minor = sys.version_info[0:2]

if (major, minor) < (2, 4):
    sys.stderr.write('Python 2.4 or later is needed to use this script\n')
    sys.exit(1)

if (major, minor) < (2, 5):
    all = lambda seq: not False in seq
    any = lambda seq: True in seq

# ++++++++++++++++++++
# global definitions
# ++++++++++++++++++++

sections = ['\section','\subsection*','\subparagraph']

# Define some regular expressions

module      = re.compile('^module',re.IGNORECASE)
module_end  = re.compile('^end\s*module',re.IGNORECASE)

program      = re.compile('^program',re.IGNORECASE)
program_end  = re.compile('^end\s*program',re.IGNORECASE)

type_re   = re.compile(r'^type\s+(?!\()',re.IGNORECASE)
type_end  = re.compile('^end\s*type',re.IGNORECASE)

types       = r'recursive|pure|double precision|elemental|(real(\(.*?\))?)|(complex(\(.*?\))?)|(integer(\(.*?\))?)|(logical)|(character(\(.*?\))?)|(type\s*\().*?(\))'
attribs     = r'allocatable|pointer|save|dimension\(.*?\)|parameter|target' # jrk33 added target
a_attribs   = r'allocatable|pointer|save|dimension\(.*?\)|intent\(.*?\)|optional|target'

types_re    = re.compile(types,re.IGNORECASE)

quoted      = re.compile('(\".*?\")|(\'.*?\')') # A quoted expression
comment     = re.compile('!.*')                 # A comment
whitespace  = re.compile(r'^\s*')               # Initial whitespace
c_ret       = re.compile(r'\r')

iface       = re.compile('^interface',re.IGNORECASE)
iface_end   = re.compile('^end\s*interface',re.IGNORECASE)

subt        = re.compile(r'^(recursive\s+)?subroutine',re.IGNORECASE)
subt_end    = re.compile(r'^end\s*subroutine\s*(\w*)',re.IGNORECASE)

recursive = re.compile('recursive',re.IGNORECASE)

funct       = re.compile('^(('+types+r')\s+)*function',re.IGNORECASE)
#funct       = re.compile('^function',re.IGNORECASE)
funct_end   = re.compile('^end\s*function\s*(\w*)',re.IGNORECASE)

prototype   = re.compile(r'^module procedure ([a-zA-Z0-9_,\s]*)')

contains = re.compile('^contains',re.IGNORECASE)

uses =  re.compile('^use\s+',re.IGNORECASE)
only =  re.compile('only\s*:\s*',re.IGNORECASE)

decl        =  re.compile('^('+types+r')\s*(,\s*('+attribs+r')\s*)*(::)?\s*\w+(\s*,\s*\w+)*',re.IGNORECASE)
d_colon     = re.compile('::')

attr_re     = re.compile('(,\s*('+attribs+r')\s*)+',re.IGNORECASE)
s_attrib_re = re.compile(attribs,re.IGNORECASE)


decl_a        = re.compile('^('+types+r')\s*(,\s*('+a_attribs+r')\s*)*(::)?\s*\w+(\s*,\s*\w+)*',re.IGNORECASE)
attr_re_a     = re.compile('(,\s*('+a_attribs+r')\s*)+',re.IGNORECASE)
s_attrib_re_a = re.compile(a_attribs,re.IGNORECASE)

cont_line   = re.compile('&')

fdoc_comm     = re.compile(r'^!\s?\*FD')
fdoc_comm_mid = re.compile(r'!\s?\*FD')
fdoc_mark     = re.compile('_FD\s?')
fdoc_rv_mark  = re.compile('_FDRV\s?')

result_re = re.compile(r'result\s*\((.*?)\)')

latex_ = re.compile(r'([_])')
latex_special_chars = re.compile(r'([%#])')

arg_split = re.compile(r'\s*(\w*)\s*(\(.+?\))?\s*(=\s*[\w\.]+\s*)?,?\s*')

size_re = re.compile(r'size\(([^,]+),([^\)]+)\)')
dimension_re = re.compile(r'^([-0-9.e]+)|((rank\(.*\))|(size\(.*\))|(len\(.*\))|(slen\(.*\)))$')

verbatim = False
displaymath = False

alnum = string.ascii_letters+string.digits+'_'

do_debug = True


string_lengths = {
    'key_len':256,
    'value_len':1024,
    'value_length':1023,
    'field_length':1023,
    'string_length':10240,
    'table_string_length':10,
    'default': 1024,
    'dict_char_a_len': 10
    }

valid_dim_re = re.compile(r'^(([-0-9.e]+)|(size\([_a-zA-Z0-9\+\-\*\/]*\))|(len\(.*\)))$')

#args_str_re = re.compile(r"""^\s*call param_register\([a-zA-Z][a-zA-Z0-9_]*,\s+(['"])([a-zA-Z][a-zA-Z0-9_]*)\1,\s+(['"])(.*?)\3,\s+.+?help_string=([\'"])(.+?)\5\)""")

args_str_re = re.compile(r"""^\s*call param_register\([a-zA-Z][a-zA-Z0-9_%]*\s*,\s*(['"])([a-zA-Z][a-zA-Z0-9_]*)\1\s*,\s*(['"]?)(.*?)\3\s*,\s*([a-zA-Z_][a-zA-Z0-9_%]*)\s*,.+?help_string=(['"])(.+?)\6\)""")

def debug(str):
    if do_debug:
        sys.stderr.write(str+'\n')
        

def print_line(str):

    global verbatim, displaymath

    if str == '':
        print
        return

    # Lines starting with '>' are to be printed verbatim
    if verbatim:
        if str[0] == '>':
            print str[1:]
            return
        else:
            verbatim = False
            #print_line(r'\end{verbatim}')
            #print_line(r'\end{sidebar}')
            print r'''\end{verbatim}
            \end{boxedminipage}

                    '''
    else:
        if str[0] == '>':
            print r'''
            
            \begin{boxedminipage}{\textwidth}
            \begin{verbatim}'''

#            print_line(r'\begin{sidebar}')
#            print_line(r'\begin{verbatim}')


            verbatim = True
            print str[1:]
            return
        else:
            pass

    if displaymath:
        if re.search(r'\\end{(displaymath|equation|eqnarray)\*?}',str):
            displaymath = False
    else:
        if re.search(r'\\begin{(displaymath|equation|eqnarray)\*?}',str):
            displaymath = True

    # Escape % and # everywhere
    s = latex_special_chars.sub(r'\\\1',str)

    if not displaymath and not verbatim:
        # Put code examples in single quotes in \texttt{} font
        s = re.sub(r"\\'", r"\\quote\\", s)
        s = re.sub(r"'(.*?)'",r'\\texttt{\1}',s)
        s = re.sub(r"\\quote\\", r"'", s)

        # Escape '_' only when not between $...$
        L = re.split(r'\$',s)
        L[::2] = [latex_.sub(r'\\\1',p) for p in L[::2]]
    
        print '$'.join(L)

    else:
        print s
        

def uniq(L, idfun=None):
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in L:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result


def combine_elements(elements):
    element_dict = {}
    func_args = []
    i = 0 # counter for appearance order of args
    for a in elements:
        if isinstance(a,C_subt) or isinstance(a,C_funct):
            func_args.append(a)
            continue
        i = i + 1
        element_dict[a.name] = (a,i)

    # Combine names with the same type, attributes and doc string
    rev_dict = {}
    for type, name in zip( [ x[0].type.lower() + str([y.lower for y in x[0].attributes]) + str(x[0].doc) \
                             for x in element_dict.values() ], element_dict.keys()):
        if rev_dict.has_key(type):
            rev_dict[type].append(element_dict[name])
        else:
            rev_dict[type] = [element_dict[name]]

    for k in rev_dict:
        names = [x[0].name for x in rev_dict[k]]
        a = rev_dict[k][0][0]
        names.sort(key=lambda x: element_dict[x][1])
        alist = []
        while names:
            n = 0
            length = 0
            while (length < 30 and n < len(names)):
                length = length + len(names[n])
                n = n +1
            ns = names[:n]
            del names[:n]
            b = copy.copy(a)
            b.name = ', '.join(ns)
            alist.append(b)

        rev_dict[k] = (alist, min([x[1] for x in rev_dict[k]]))

    # Sort by original appearance order of first name
    keys = rev_dict.keys()
    keys.sort(key=lambda x: rev_dict[x][1])

    return keys, rev_dict, func_args


# ++++++++++++++++++++
# CLASS DEFINITIONS for f90doc
# ++++++++++++++++++++

class C_prog:

    def __init__(self):
        self.name=''
        self.doc=[]
        self.subts=[]
        self.functs=[]
        self.uses=[]

    def __eq__(self, other):
        if other is None: return False
        return (self.name == other.name and
                self.doc == other.doc and
                self.subts == other.subts and
                self.functs == other.functs and
                self.uses == other.uses)

    def __ne__(self, other):
        return not self.__eq__(other)

    def latex(self,depth,fn='', short_doc=False):

        if self.doc:
            if self.doc[0].strip() == 'OMIT':
                return

            if self.doc[0].strip() == 'OMIT SHORT':
                if short_doc:
                    return
                else:
                    self.doc = self.doc[1:]

        
        print_line( r"\newpage")
        print_line(r'\index{general}{'+self.name+' program}')
        print_line( sections[depth]+r'[Program \texttt{')
        print_line( self.name+'}]')
        print_line( r"""{Program \texttt{""")
        print_line( self.name )
        if depth==0:
            print_line( r"""} in file """+fn+"""}""")
        else:
            print_line( r"}}")
        print_line( sections[depth+1]+"""{Purpose}""")
        for a in self.doc:
            #print_line( a)
            print_line(a)
        if self.uses!=[]:
            print_line( sections[depth+1]+r"{Uses}")
            u_temp=''
            for a in self.uses:
                u_temp=u_temp+a+', '
                if len(u_temp)>50:
                    print_line( r"\texttt{"+u_temp[:-2]+"}")
                    u_temp=''
                    print_line("\n")
            if u_temp!='':
                print_line( r"\texttt{"+u_temp[:-2]+"}")
        for a in self.subts:
            a.latex(depth+1,short_doc=short_doc)
        for a in self.functs:
            a.latex(depth+1,short_doc=short_doc)

#+++++++++++++++++++++++++++++++++++++++++++++

class C_module:

    def __init__(self):
        self.name=''
        self.types=[]
        self.elements=[]
        self.subts=[]
        self.functs=[]
        self.doc=[]
        self.uses=[]
        self.interfaces=[]

    def __eq__(self, other):
        if other is None: return False
        return (self.name == other.name and
                self.types == other.types and
                self.elements == other.elements and
                self.subts == other.subts and
                self.functs == other.functs and
                self.doc == other.doc and
                self.uses == other.uses and
                self.interfaces == other.interfaces)
        
    def __ne__(self, other):
        return not self.__eq__(other)

    def display(self):
        print 'module',self.name,self.doc
        for a in self.types:
            a.display()
        print '    module variables:'
        for a in self.elements:
            a.display()
        for a in self.subts:
            a.display()

    def latex(self,depth,fn='',short_doc=False):

        if self.doc:
            if self.doc[0].strip() == 'OMIT':
                return

            if self.doc[0].strip() == 'OMIT SHORT':
                if short_doc:
                    return
                else:
                    self.doc = self.doc[1:]
                
        
        print_line( r"\newpage")
        print_line(r'\index{general}{'+self.name+' module}')
        print_line( sections[depth]+r'[Module \texttt{')
        print_line( self.name+'}]')
        print_line( r"""{Module \texttt{""")
        print_line( self.name)
        if depth==0:
            print_line( r"""} in file """+fn+"""}""")
        else:
            print_line( r"}}")
        print_line( sections[depth+1]+"""{Purpose}""")
        for a in self.doc:
            print_line(a)
            #print a
        print_line( sections[depth+1]+r'{Usage}')
        print_line('>    use '+self.name)
        if self.uses!=[]:
            print_line( sections[depth+1]+r"{Uses}")
            u_temp=''
            for a in self.uses:
                u_temp=u_temp+a+', '
                if len(u_temp)>50:
                    print_line( r"\texttt{"+u_temp[:-2]+"}")
                    u_temp=''
                    print_line("\n")
            if u_temp!='':
                print_line( r"\texttt{"+u_temp[:-2]+"}")

        if self.elements!=[]:
            print_line( sections[depth+1]+r"""{Module variables}""")

            keys, rev_dict, func_args = combine_elements(self.elements)

            print_line( r"\begin{description}")

            for k in keys:
                for a in rev_dict[k][0]:
                    a.latex(short_doc=short_doc)

            print_line( r"\end{description}")
            
        for a in self.types:
            a.latex(short_doc=short_doc)

        for a in self.interfaces:
            a.latex(depth+1,short_doc=short_doc)
        for a in self.subts:
            a.latex(depth+1,short_doc=short_doc)
        for a in self.functs:
            a.latex(depth+1,short_doc=short_doc)

#+++++++++++++++++++++++++++++++++++++++++++++

class C_subt:

    def __init__(self):
        self.name=''
        self.arguments=[]
        self.doc=[]
        self.uses=[]
        self.recur=''
        self.args_str = {}

    def __eq__(self, other):
        if other is None: return False
        return (self.name == other.name and
                self.arguments == other.arguments and
                self.doc == other.doc and
                self.uses == other.uses and
                self.recur == other.recur)
        
    def __ne__(self, other):
        return not self.__eq__(other)

    
    def display(self):
        print '    subroutine',self.name,'(',
        for i in range(len(self.arguments)-1):
            print self.arguments[i].name,',',
        print self.arguments[len(self.arguments)-1].name,
        print ')',self.doc
        for a in self.arguments:
            a.display()

    def latex_compact(self, depth, short_doc=False):
        if self.doc:
            if self.doc[0].strip() == 'OMIT':
                return

            if self.doc[0].strip() == 'OMIT SHORT':
                if short_doc:
                    return
                else:
                    self.doc = self.doc[1:]

        if self.arguments!=[]:
            argl='('+','.join([x.name for x in self.arguments])+')'
        else:
            argl=''

        d_ent=r"Subroutine \texttt{"+self.name+argl+"}"

        print_line(r"\item["+d_ent+r"]\mbox{} \par\noindent")

        if self.arguments!=[]:

            print_line(r'\begin{description}')

            keys, rev_dict, func_args = combine_elements(self.arguments)

            for k in keys:
                for a in rev_dict[k][0]:
                    a.latex(short_doc=short_doc)


            print_line(r'\end{description}')


        for a in self.doc:
            #print a
            print_line(a)
        print_line('')


    def latex(self,depth,fn='',short_doc=False):

        if self.doc:
            if self.doc[0].strip() == 'OMIT':
                return

            if self.doc[0].strip() == 'OMIT SHORT':
                if short_doc:
                    return
                else:
                    self.doc = self.doc[1:]


        if self.arguments!=[]:
            argl='('
            for a in range(len(self.arguments)):
                arg = self.arguments[a]
                if isinstance(arg,C_decl) and 'optional' in arg.attributes:
                    if argl[-2:] == '],':
                        argl = argl[:-2]+','+arg.name.rstrip()+'],'
                    elif argl.rstrip()[-4:] == '], &':
                        argl = argl.rstrip()[:-4]+', &\n                        '+arg.name.rstrip()+'],'
                    elif argl[-1] == ',':
                        argl = argl[:-1]+'[,'+arg.name.rstrip()+'],'
                    else:
                        argl = argl+'['+arg.name.rstrip()+'],'
                else:
                    argl=argl+arg.name.rstrip()+','
                if (a+1)%4==0.0 and a+1!=len(self.arguments):
                    argl=argl+' &\n                        '
            argl=argl[:-1]+')'
        else:
            argl=''

        print_line(r'\index{general}{'+self.name+' subroutine}')
        
        if self.recur=='':            
            print_line( sections[depth]+r""" {Subroutine \texttt{"""+self.name)
        else:
            print_line( sections[depth]+r"""{Recursive subroutine \texttt{"""+self.name)
            
            #        if depth==0:
            #            print "} (in file "+latex_escape(fn)+")}"
            #        else:
        print_line("""}}""")
        print_line('>    call '+self.name+argl)
        
        for a in self.doc:
            print_line(a)

        if self.uses!=[]:
            print_line( sections[depth+1]+r"{Uses}")
            u_temp=''
            for a in self.uses:
                u_temp=u_temp+a+', '
                if len(u_temp)>50:
                    print_line( r"\texttt{"+u_temp[:-2]+"}")
                    u_temp=''
                    print_line("\n")
            if u_temp!='':
                print_line( r"\texttt{"+u_temp[:-2]+"}")

        if self.arguments!=[]:

            keys, rev_dict, func_args = combine_elements(self.arguments)
            
            print_line( r"\begin{description}")

            for k in keys:
                for a in rev_dict[k][0]:
                    a.latex(short_doc=short_doc)

            for f in func_args:
                f.latex_compact(depth,short_doc=short_doc)

            print_line( r"\end{description}")


        



#+++++++++++++++++++++++++++++++++++++++++++++

class C_funct:

    def __init__(self):
        self.name=''
        self.arguments=[]
        self.procedures=[]
        self.doc=[]
        self.uses=[]
        self.ret_val=None
        self.ret_val_doc=[]
        self.recur=''
        self.args_str = {}

    def __eq__(self, other):
        if other is None: return False
        return (self.name == other.name and
                self.arguments == other.arguments and
                self.procedures == other.procedures and
                self.doc == other.doc and
                self.uses == other.uses and
                self.ret_val == other.ret_val and
                self.ret_val_doc == other.ret_val_doc and
                self.recur == other.recur)
        
    def __ne__(self, other):
        return not self.__eq__(other)

    
    def display(self):
        print '   function',self.name,'(',
        for i in range(len(self.arguments)-1):
            print self.arguments[i].name,',',
        print self.arguments[len(self.arguments)-1].name,
        print ')',self.doc
        for a in self.arguments:
            a.display()

    def latex_compact(self, depth, short_doc=False):
        if self.doc:
            if self.doc[0].strip() == 'OMIT':
                return

            if self.doc[0].strip() == 'OMIT SHORT':
                if short_doc:
                    return
                else:
                    self.doc = self.doc[1:]

        if self.arguments!=[]:
            argl='('+','.join([x.name for x in self.arguments])+')'
        else:
            argl=''


        d_ent=r"Function \texttt{"+self.name+argl+"} --- "+self.ret_val.type

        for a in self.ret_val.attributes:
            d_ent=d_ent+", "+a

        print_line(r"\item["+d_ent+r"]\mbox{} \par\noindent")

        if self.arguments!=[]:

            keys, rev_dict, func_args = combine_elements(self.arguments)

            print_line(r'\begin{description}')
            for k in keys:
                for a in rev_dict[k][0]:
                    a.latex(short_doc=short_doc)
            print_line(r'\end{description}')
        
        for a in self.doc:
            #print a
            print_line(a)
        print_line('')


    def latex(self,depth,fn='',short_doc=False):

        if self.doc:
            if self.doc[0].strip() == 'OMIT':
                return

            if self.doc[0].strip() == 'OMIT SHORT':
                if short_doc:
                    return
                else:
                    self.doc = self.doc[1:]


#        print_line( r"""\begin{center}\rule{10cm}{0.5pt}\end{center}""")
#        print_line( r"""\rule{\textwidth}{0.5pt}""")


        print_line(r'\index{general}{'+self.name+' function}')

        if self.arguments!=[]:
            argl='('
            for a in range(len(self.arguments)):
                arg = self.arguments[a]
                if isinstance(arg,C_decl) and 'optional' in arg.attributes:
                    if argl[-2:] == '],':
                        argl = argl[:-2]+','+arg.name.rstrip()+'],'
                    elif argl.rstrip()[-4:] == '], &':
                        argl = argl.rstrip()[:-4]+', &\n                        '+arg.name.rstrip()+'],'
                    elif argl[-1] == ',':
                        argl = argl[:-1]+'[,'+arg.name.rstrip()+'],'
                    else:
                        argl = argl+'['+arg.name.rstrip()+'],'
                else:
                    argl=argl+arg.name.rstrip()+','
                if (a+1)%4==0.0 and a+1!=len(self.arguments):
                    argl=argl+' &\n                        '
            argl=argl[:-1]+')'
        else:
            argl=''
        
        if self.recur=='':
            print_line( sections[depth]+r"""{Function \texttt{"""+self.name)
        else:
            print_line( sections[depth]+r"""{Recursive function\texttt{"""+self.name)
        if depth==0:
            print_line("} (in file "+fn+")}")
        else:
            print_line("""}}""")
#        print_line( sections[depth+1]+r'{Usage}')
#        print_line(r'\begin{boxedminipage}{\textwidth}')
        ret_name = self.ret_val.name
        if ret_name.lower() == self.name.lower():
            ret_name = ret_name[0].lower()
        print_line('>    '+ret_name+' = '+self.name+argl)
#        print_line(r'\end{boxedminipage}'+'\n\n')
        for a in self.doc:
            #print a
            print_line(a)

        if self.uses!=[]:
            print_line(sections[depth+1]+r"{Uses}")
            u_temp=''
            for a in self.uses:
                u_temp=u_temp+a+', '
                if len(u_temp)>50:
                    print_line(r"\texttt{"+u_temp[:-2]+"}")
                    u_temp=''
                    print_line("\n")
            if u_temp!='':
                print_line(r"\texttt{"+u_temp[:-2]+"}")

        print_line(r"\begin{description}")

        if self.arguments!=[]:

            keys, rev_dict, func_args = combine_elements(self.arguments)

            for k in keys:
                for a in rev_dict[k][0]:
                    a.latex(short_doc=short_doc)

            for f in func_args:
                f.latex_compact(depth,short_doc=short_doc)

        
        #        print_line(sections[depth+1]+"{Return value --- ")


        print_line(r"\item[Return value --- ")

        self.ret_val.latex_rv()
        print_line(r"]\mbox{} \par\noindent")
        for a in self.ret_val_doc:
            print_line(a)

        print_line(r"\end{description}")






#+++++++++++++++++++++++++++++++++++++++++++++
        
class C_decl:

    def __init__(self):
        self.name=''
        self.type=''
        self.attributes=[]
        self.doc=[]
        self.value=''

    def __eq__(self, other):
        if other is None: return False
        return (self.name == other.name and
                self.type == other.type and
                self.attributes == other.attributes and
                self.doc == other.doc and
                self.value == other.value)
        
    def __ne__(self, other):
        return not self.__eq__(other)

    def display(self):
        print '        ',self.name,'\t',self.type,
        print self.attributes,
        if self.value!='':
            print 'value='+self.value,
        print self.doc

    def latex(self, short_doc=False):

        if self.doc:
            if self.doc[0].strip() == 'OMIT':
                return

            if self.doc[0].strip() == 'OMIT SHORT':
                if short_doc:
                    return
                else:
                    self.doc = self.doc[1:]


        if type(self.type) == type([]) and len(self.type) > 1:
            d_ent=r'\texttt{'+self.name+'} --- '
            

            for a in self.attributes:
                d_ent=d_ent+' '+a+', '

            if d_ent[-1] == ',':
                d_ent=d_ent[:-2] # remove trailing ','

            if (sum([len(t) for t in self.type])+len(self.attributes) < 30):
                print_line(r"\item["+d_ent+' \emph{or} '.join(self.type)+r"]\mbox{} \par\noindent")
            else:
                print_line(r"\item["+d_ent+r"]\mbox{} \par\noindent")
                print_line(r'\bfseries{'+' \emph{or} '.join(self.type)+r'} \par\noindent')

        else:
            if (type(self.type) == type([])):
                typename = self.type[0]
            else:
                typename = self.type
            d_ent=r"\texttt{"+self.name+"} --- "+typename

            for a in self.attributes:
                d_ent=d_ent+", "+a
                
            print_line(r"\item["+d_ent+r"]\mbox{} \par\noindent")


#        if self.value!='':
#            d_ent=d_ent+r", value = \texttt{"+latex_escape(self.value)+'}'

        for a in self.doc:
            #print a
            print_line(a)
        print_line('')
        
    def latex_rv(self):

        d_ent=self.type

        for a in self.attributes:
            d_ent=d_ent+", "+a
        if self.value!='':
            d_ent=d_ent+r", value = \texttt{"+self.value+'}'
        print_line(d_ent)
    
#+++++++++++++++++++++++++++++++++++++++++++++

class C_type:

    def __init__(self):
        self.name=''
        self.elements=[]
        self.doc=[]

    def __eq__(self, other):
        if other is None: return False
        return (self.name == other.name and
                self.elements == other.elements and
                self.doc == other.doc)
        
    def __ne__(self, other):
        return not self.__eq__(other)
        
    
    def display(self):
        print '    type',self.name,self.doc
        for a in self.elements:
            a.display()

    def latex(self, short_doc=False):

        if self.doc:
            if self.doc[0].strip() == 'OMIT':
                return

            if self.doc[0].strip() == 'OMIT SHORT':
                if short_doc:
                    return
                else:
                    self.doc = self.doc[1:]


        #        print_line(r"""\begin{center}\rule{10cm}{0.5pt}\end{center}""")
        #        print_line( r"""\rule{\textwidth}{0.5pt}""")


        print_line(r'\index{general}{'+self.name+' type}')
        
        print_line(r"""\subsection*{Type \texttt{"""+self.name+"""}}""")
        for a in self.doc:
            #print a
            print_line(a)
            
        print_line(r"""\subsubsection*{Elements}""")

        keys, rev_dict, func_args = combine_elements(self.elements)
        print_line(r"\begin{description}")

        for k in keys:
            for a in rev_dict[k][0]:
                a.latex(short_doc=short_doc)

        print_line(r"\end{description}")


#+++++++++++++++++++++++++++++++++++++++++++++

class C_interface:
    def __init__(self):
        self.name = ''
        self.procedures = []
        self.subts = []
        self.functs = []
        self.doc = []

    def __eq__(self, other):
        if other is None: return False
        return (self.name == other.name and
                self.procedures == other.procedures and
                self.subts == other.subts and
                self.functs == other.functs and
                self.doc == other.doc)
        
    def __ne__(self, other):
        return not self.__eq__(other)        

    def display(self):
        print '     interface', self.name, self.doc
        for a in self.elements:
            a.display()

    def latex(self, depth, short_doc=False):
        
        if self.doc:
            if self.doc[0].strip() == 'OMIT':
                return

            if self.doc[0].strip() == 'OMIT SHORT':
                if short_doc:
                    return
                else:
                    self.doc = self.doc[1:]


        
        #        print_line( r"""\rule{\textwidth}{0.5pt}""")


        print_line(r'\index{general}{'+self.name+' interface}')
        
        print_line(sections[depth]+r'{Interface \texttt{'+self.name+'}}')

        #        print_line(sections[depth+1]+r"""{Usage}""")

        is_sub = len(self.subts) != 0

        printed_args = []
        #        print_line(r'\begin{boxedminipage}{\textwidth}')
        for sub in self.functs+self.subts:

            if sub.arguments!=[]:
                argl='('
                for a in range(len(sub.arguments)):
                    arg = sub.arguments[a]
                    if isinstance(arg,C_decl) and 'optional' in arg.attributes:
                        if argl[-2:] == '],':
                            argl = argl[:-2]+','+arg.name.rstrip()+'],'
                        elif argl.rstrip()[-4:] == '], &':
                            argl = argl.rstrip()[:-4]+', &\n                        '+arg.name.rstrip()+'],'
                        elif argl[-1] == ',':
                            argl = argl[:-1]+'[,'+arg.name.rstrip()+'],'
                        else:
                            argl = argl+'['+arg.name.rstrip()+'],'
                    else:
                        argl=argl+arg.name.rstrip()+','
                    if (a+1)%4==0.0 and a+1!=len(sub.arguments):
                        argl=argl+' &\n                        '
                argl=argl[:-1]+')'
            else:
                argl=''

            if not is_sub and sub.ret_val.name != sub.name:
                hash_value = argl
            else:
                hash_value = argl

            if hash_value in printed_args:
                continue

            printed_args.append(hash_value)

            if not is_sub:
                ret_name = sub.ret_val.name
                if ret_name.lower() == self.name.lower() or ret_name.lower() == sub.name.lower():
                    ret_name = ret_name[0].lower()+str((self.functs+self.subts).index(sub)+1)
                print_line('>    '+ret_name+' = '+self.name+argl)
            else:
                print_line('>    call '+self.name+argl)
                #        print_line(r'\end{boxedminipage}'+'\n\n')


        for a in self.doc:
            print_line(a)

        for sub in self.functs+self.subts:
            for a in sub.doc:
                print_line(a)
            print_line('\n\n')

        got_args = (self.subts != [] and \
                   (sum([len(x.arguments) for x in self.subts])+\
                    sum([len(x.arguments) for x in self.functs]) != 0)) or self.functs != []
           

        func_args = []
        if got_args:
            print_line(r'\begin{description}')


            arg_dict = {}
            i = 0 # counter for appearance order of args
            for sub in self.functs+self.subts:
                for a in sub.arguments:
                    if isinstance(a,C_subt) or isinstance(a,C_funct):
                        func_args.append(a)
                        continue
                    i = i + 1
                    if arg_dict.has_key(a.name):
                        if a.type.lower()+str(sorted(map(string.lower,a.attributes))) in \
                           [x[0].type.lower()+str(sorted(map(string.lower, x[0].attributes))) for x in arg_dict[a.name]]:
                            pass # already got this name/type/attribute combo 
                        else:
                            arg_dict[a.name].append((a,i))
                            
                    else:
                        arg_dict[a.name] = [(a,i)]

            # Combine multiple types with the same name
            for name in arg_dict:
                types = [x[0].type for x in arg_dict[name]]
                types = uniq(types, string.lower)
                attr_lists = [x[0].attributes for x in arg_dict[name]]
                attributes = []

                contains_dimension = [ len([x for x in y if x.find('dimension') != -1]) != 0 for y in attr_lists ]
                
                for t in attr_lists:
                    attributes.extend(t)
                attributes = uniq(attributes, string.lower)

                dims = [x for x in attributes if x.find('dimension') != -1]
                attributes = [x for x in attributes if x.find('dimension') == -1]

                # If some attribute lists contains 'dimension' and some don't then
                # there are scalars in there as well.
                if True in contains_dimension and False in contains_dimension:
                    dims.insert(0, 'scalar')

                
                if (len(dims) != 0):
                    attributes.append(' \emph{or} '.join(dims))

                a = arg_dict[name][0][0]
                a.type = types #r' \emph{or} '.join(types)
                a.attributes = attributes
                arg_dict[name] = (a, arg_dict[name][0][1])


            # Combine names with the same type, attributes and doc string
            rev_dict = {}
            for type, name in zip( [ str([y.lower for y in x[0].type]) + \
                                     str([y.lower for y in x[0].attributes]) + str(x[0].doc) \
                                     for x in arg_dict.values() ], arg_dict.keys()):
                if rev_dict.has_key(type):
                    rev_dict[type].append(arg_dict[name])
                else:
                    rev_dict[type] = [arg_dict[name]]

            for k in rev_dict:
                names = [x[0].name for x in rev_dict[k]]
                a = rev_dict[k][0][0]
                names.sort(key=lambda x: arg_dict[x][1])

                # Split into pieces of max length 30 chars
                alist = []
                while names:
                    n = 0
                    length = 0
                    while (length < 30 and n < len(names)):
                        length = length + len(names[n])
                        n = n +1
                    ns = names[:n]
                    del names[:n]
                    b = copy.copy(a)
                    b.name = ', '.join(ns)
                    alist.append(b)

                rev_dict[k] = (alist, min([x[1] for x in rev_dict[k]]))

            # Sort by original appearance order of first name
            keys = rev_dict.keys()
            keys.sort(key=lambda x: rev_dict[x][1])

            for k in keys:
                for a in rev_dict[k][0]:
                    a.latex(short_doc=short_doc)
                            
            for f in func_args:
                f.latex_compact(depth,short_doc=short_doc)


        if self.functs != []:
            #            print_line(sections[depth+1]+"{Return value --- ")

            ret_types = [a.ret_val.type+str(a.ret_val.attributes) for a in self.functs]

            if len(filter(lambda x: x != self.functs[0].ret_val.type+str(self.functs[0].ret_val.attributes), \
                          ret_types)) == 0:
                
                print_line(r"\item[Return value --- ")
                self.functs[0].ret_val.latex_rv()
                print_line("]")
                for a in self.functs[0].ret_val_doc:
                    print_line(a)
            else:
                print_line(r"\item[Return values:]\mbox{} \par\noindent")
                print_line(r'\begin{description}')
                for f in self.functs:
                    shortname = f.ret_val.name[0].lower()+str(self.functs.index(f)+1)
                    print_line(r"\item[\texttt{"+shortname+"} --- ")
                    f.ret_val.latex_rv()
                    print_line(']')
                    for a in f.ret_val_doc:
                        print_line(a)
                print_line(r'\end{description}')



        if got_args:
            print_line(r"\end{description}")


class f90file:
    
    def __init__(self,fname):
        self.file=open(fname,'r')
        self.lines=self.file.readlines()
        self.file.close()
        self.dquotes=[]
        self.squotes=[]

    def readlines(self):
        lines = []
        line = self.next_line()
        while line is not None:
            lines.append(line)
            line = self.next_line()
        return lines

    def next_line(self):

        cline=''

        while (cline=='' and len(self.lines)!=0):
            cline=self.lines[0].strip()
            if cline.find('_FD')==1:
                break

            # jrk33 - join lines before removing delimiters

            # Join together continuation lines
            FD_index=cline.find('_FD')
            com2_index=cline.find('_COMMENT')
            if (FD_index==0 or com2_index==0):
                pass
            else:
                cont_index=cline.find('&')
                comm_index=cline.find('!')
                while (cont_index!=-1 and (comm_index==-1 or comm_index>cont_index)):
                    self.lines=[cline[:cont_index].strip()+self.lines[1]]+self.lines[2:]
                    cline=self.lines[0].strip()
                    cont_index=cline.find('&')


            # Remove quoted sections and store
            # jrk33 - swapped order of " and ' for nesting
#            s_quote,cline=remove_delimited(cline,"'","'")
#            d_quote,cline=remove_delimited(cline,'"','"')


#            self.dquotes=self.dquotes+d_quote
#            self.squotes=self.squotes+s_quote

            # split by '!', if necessary            
            comm_index=cline.find('!')
            if comm_index!=-1:
                self.lines=[cline[:comm_index],cline[comm_index:]]+self.lines[1:]
                cline=self.lines[0]
                cline=cline.strip()
                # jrk33 - changed comment mark from '!*FD' to '!%'
                if self.lines[1].find('!%')!=-1:
                    self.lines=[self.lines[0]]+['_FD'+self.lines[1][2:]]+self.lines[2:]
                else:
                    self.lines=[self.lines[0]]+['_COMMENT'+self.lines[1][1:]]+self.lines[2:]
                 
##            # Split by ';', if necessary
##            else:
##                self.lines=cline.split(';')+self.lines[1:]
##                cline=self.lines[0]
##                cline=cline.strip()
        
            self.lines=self.lines[1:]
#            cline,self.dquotes=recover_delimited(cline,'"','"',self.dquotes)
#            cline,self.squotes=recover_delimited(cline,"'","'",self.squotes)

#            for i in range(len(self.lines)):
                # jrk33 - swapped order of ' and "
#                if self.squotes!=[]:
#                    self.lines[i],self.squotes=recover_delimited(self.lines[i],"'","'",self.squotes)            
#                if self.dquotes!=[]:
#                    self.lines[i],self.dquotes=recover_delimited(self.lines[i],'"','"',self.dquotes)


        if cline=='':
            return None
        else:
#            debug(cline)
            return cline

def remove_delimited(line,d1,d2):

    bk=0
    temp_str=''
    undel_str=''
    delimited=[]
    
    for i in range(len(line)):
        if bk==1:
            if line[i]==d2:
                bk=0
                delimited.append(temp_str[:])
                temp_str=''
                undel_str=undel_str+line[i]
                continue
            temp_str=temp_str+line[i]
            continue
        if line[i]==d1:
            bk=1
        undel_str=undel_str+line[i]

    if bk==1:
        undel_str=undel_str+temp_str

    return delimited,undel_str

def recover_delimited(line,d1,d2,delimited):

    if delimited==[]:
        return line,[]

    i=0
    while i<len(line):
        if line[i]==d1:
            line=line[0:i+1]+delimited[0]+line[i+1:]
            i=i+len(delimited[0])+1
            delimited=delimited[1:]
        i=i+1

    return line,delimited

#----------------------------------------------

hold_doc = None

# ++++++++++++++++++++
# FUNCTION DEFINITIONS
# ++++++++++++++++++++



def usage():
    print r"""
f90doc: documentation generator for Fortran 90.
Usage: f90doc [-t <title>] [-a <author>] [-i <intro>] [-n] [-s] [-l] [-f] <filenames>
Options:
    -t: specify a title. 
    -a: specify an author.
    -i: specify intro file to be included
    -n: don't use an intro file
    -s: use short documentation format
    -l: write latex output
    
Notes:
    Output is as LaTeX, to standard output.
    A single LaTeX document is generated from all input
    files.

Markup:
    Prepare your document by adding comments after the things
    to which they refer (or on the same line), prefixed with !*FD
    For comments relating to function return values, use !*FDRV

Copyright (C) 2004 Ian Rutt. f90doc comes with ABSOLUTELY NO WARRANTY.
Modifications for Fortran 95 (c) 2006-2008 James Kermode.

This is free software, and you are welcome to redistribute it under
certain conditions. See LICENCE for details.
"""


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_uses(cline,file):

    if re.match(uses,cline)!=None:
        cline=uses.sub('',cline)
        cline=cline.strip()
        out=re.match(re.compile(r"\w+"),cline).group()
        cline=file.next_line()
        return [out,cline]
    else:
        return [None,cline]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_doc(cline,file):

    if cline and re.match(fdoc_mark,cline)!=None:
        out=fdoc_mark.sub('',cline)
        out=out.rstrip()
        cline=file.next_line()
        return [out,cline]
    else:
        return [None,cline]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_doc_rv(cline,file):

    cl=cline

    if cl is None:
        return [None, cl]
    
    if re.match(fdoc_rv_mark,cl)!=None:
        out=fdoc_rv_mark.sub('',cl)
        out=out.rstrip()
        cl=file.next_line()
        return [out,cl]
    else:
        return [None,cl]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_cont(cline,file):

    cl=cline

    if re.match(contains,cl)!=None:
        cl=file.next_line()
        return ['yes',cl]
    else:
        return [None,cl]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_program(cl,file):

    global hold_doc

    out=C_prog()
    cont=0

    if re.match(program,cl)!=None:
        # Get program name

        cl=program.sub('',cl)
        out.name=re.search(re.compile('\w+'),cl).group().strip()
        if out.name=='':
            out.name='<Unnamed>'

        # Get next line, and check each possibility in turn

        cl=file.next_line()

        while re.match(program_end,cl)==None:

            # contains statement
            check=check_cont(cl,file)
            if check[0]!=None:
                cont=1
                cl=check[1]
                continue

            if cont==0:
                
                # use statements
                check=check_uses(cl,file)
                if check[0]!=None:
                    out.uses.append(check[0])
                    cl=check[1]
                    continue
                
                # Doc comment
                check=check_doc(cl,file)
                if check[0]!=None:
                    out.doc.append(check[0])
                    cl=check[1]
                    continue
            else:

                # jrk33 - hold doc comment relating to next subrt or funct 
                check=check_doc(cl,file)
                if check[0]!=None:
                    if hold_doc == None:
                        hold_doc = [check[0]]
                    else:
                        hold_doc.append(check[0])
                    cl=check[1]
                    continue

                # Subroutine definition
                check=check_subt(cl,file)
                if check[0]!=None:
                    debug('    subroutine '+check[0].name)
                    out.subts.append(check[0])
                    cl=check[1]
                    continue
                
                # Function definition
                check=check_funct(cl,file)
                if check[0]!=None:
                    debug('    function '+check[0].name)
                    out.functs.append(check[0])
                    cl=check[1]
                    continue


                
            
            # If no joy, get next line
            cl=file.next_line()

        cl=file.next_line()
        
        return [out,cl]
    else:
        return [None,cl]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_module(cl,file):

    global hold_doc

    out=C_module()
    cont=0

    if re.match(module,cl)!=None: 

        # jrk33 - if we're holding a doc comment from before
        # subroutine definition, spit it out now
        if hold_doc is not None:
            for line in hold_doc:
                out.doc.append(line)
            hold_doc = None

        # Get module name
        cl=module.sub('',cl)
        out.name=re.search(re.compile('\w+'),cl).group()
        
        # Get next line, and check each possibility in turn

        cl=file.next_line()
        
        while re.match(module_end,cl)==None:

            # contains statement
            check=check_cont(cl,file)
            if check[0]!=None:
                cont=1
                cl=check[1]
                continue

            if cont==0:

                # use statements
                check=check_uses(cl,file)
                if check[0]!=None:
                    out.uses.append(check[0])
                    cl=check[1]
                    continue
                
                # Doc comment
                check=check_doc(cl,file)
                if check[0]!=None:
                    if hold_doc == None:
                        hold_doc = [check[0]]
                    else:
                        hold_doc.append(check[0])
                    cl=check[1]
                    continue

                # jrk33 - Interface definition
                check=check_interface(cl,file)
                if check[0] != None:
                    debug('    interface '+check[0].name)
                    out.interfaces.append(check[0])
                    cl=check[1]
                    continue
                
                # Type definition
                check=check_type(cl,file)
                if check[0]!=None:
                    debug('    type '+check[0].name)
                    out.types.append(check[0])
                    cl=check[1]
                    continue

                # Module variable
                check=check_decl(cl,file)
                if check[0]!=None:
                    for a in check[0]:
                        out.elements.append(a)
                        cl=check[1]
                    continue
            else:

                # jrk33 - hold doc comment relating to next subrt or funct 
                check=check_doc(cl,file)
                if check[0]!=None:
                    if hold_doc == None:
                        hold_doc = [check[0]]
                    else:
                        hold_doc.append(check[0])
                    cl=check[1]
                    continue
                
                # Subroutine definition
                check=check_subt(cl,file)
                if check[0]!=None:

                    debug('    subroutine '+check[0].name)
                    found = None
                    for i in out.interfaces:
                        if check[0].name.lower() in i.procedures:
                            found = i
                            break
                    
                    if found != None:
                        found.subts.append(check[0])
                    else:
                        out.subts.append(check[0])
                    cl=check[1]
                    continue
                
                # Function definition
                check=check_funct(cl,file)
                if check[0]!=None:

                    debug('    function '+check[0].name)
                    found = None
                    for i in out.interfaces:
                        if check[0].name.lower() in i.procedures:
                            found = i
                            break
                    
                    if found != None:
                        found.functs.append(check[0])
                    else:
                        out.functs.append(check[0])
                    cl=check[1]
                    continue

            
            # If no joy, get next line
            cl=file.next_line()

        cl=file.next_line()

        return [out,cl]
    else:
        return [None,cl]
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_subt(cl,file, grab_hold_doc=True):

    global hold_doc

    out=C_subt()
    
    if re.match(subt,cl)!=None:
        
        # Check if recursive

        if re.match(recursive,cl)!=None:
            out.recur='recursive'
        
        # Get subt name

        cl=subt.sub('',cl)
        out.name=re.search(re.compile('\w+'),cl).group()

        # Check to see if there are any arguments

        has_args=0
        if re.search(r'\(.+',cl)!=None:
            has_args=1

        # get argument list
            
        if has_args:

            cl=re.sub('\w+','',cl,count=1)
            argl=re.split('[\W]+',cl)
        
            del(argl[0])
            del(argl[len(argl)-1])

            while cl.strip() == '' or re.search('&',cl)!=None:
                cl=file.next_line()
                if cl.strip() == '': continue
                arglt=re.split('[\W]+',cl)
                del(arglt[len(arglt)-1])
                for a in arglt:
                    argl.append()

        else:
            argl=[]

        argl = map(string.lower, argl)

        # Get next line, and check each possibility in turn

        cl=file.next_line()

        while True:

            # Use statement
            ##check=check_uses(cl,file)
            ##if check[0]!=None:
            ##    out.uses.append(check[0])
            ##    cl=check[1]
            ##    continue

            # Doc comment
            check=check_doc(cl,file)
            if check[0]!=None:
                out.doc.append(check[0])
                cl=check[1]
                continue


            if has_args:
                # Argument
                check=check_arg(cl,file)
                if check[0]!=None:
                    for a in check[0]:
                        out.arguments.append(a)
                    cl=check[1]
                    continue

                # Interface section
                check=check_interface_decl(cl,file)
                if check[0] != None:
                    for a in check[0].procedures:
                        out.arguments.append(a)
                    cl = check[1]
                    continue

                # args_str arguments
                check=check_args_str(cl,file, out.args_str)
                if check[0]!=None:
                    cl=check[1]
                    continue



            m = subt_end.match(cl)

            if m == None:
                cl=file.next_line()
                continue
            elif m.group(1).lower() == out.name.lower() or m.group(1) == '':
                break

            debug('WARNING - found "end subroutine" line with mismatching label %r (expected %r)' % (m.group(1), out.name))
            # If no joy, get next line
            cl=file.next_line()



        # Select only those declarations that match entries
        # in argument list
        
        if has_args:
            #t_re_str='(^'
            ag_temp=[]
            #for a in argl:
            #    t_re_str=t_re_str+a+'$)|(^'
            #t_re_str=t_re_str[:-3]
            #t_re=re.compile(t_re_str,re.IGNORECASE)

            for i in out.arguments:
                if i.name.lower() in argl:
                    ag_temp.append(i)
                elif i.name.lower() in out.args_str:
                    out.args_str[i.name.lower()].type = i.type

            out.arguments=ag_temp
            out.arguments.sort(key=lambda x:argl.index(x.name.lower()))

        else:
            out.arguments=[]
            
        cl=file.next_line()

        # jrk33 - if we're holding a doc comment from before
        # subroutine definition, spit it out now
        if grab_hold_doc and hold_doc is not None:
            for line in hold_doc:
                out.doc.append(line)
            hold_doc = None


        return [out,cl]
    else:
        return [None,cl]
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_funct(cl,file,grab_hold_doc=True):

    global hold_doc

    out=C_funct()
    out.ret_val=C_decl()
    
    if re.match(funct,cl)!=None:

        # Check if recursive

        if re.search(recursive,cl)!=None:
            out.recur='recursive'
            cl=recursive.sub('',cl)

        # Get return type, if present

        cl=cl.strip()
        if re.match(types_re,cl)!=None:
            out.ret_val.type=re.match(types_re,cl).group()

        # jrk33 - Does function header specify alternate name of
        # return variable?
        ret_var = None
        if re.search(result_re,cl)!=None:
            ret_var = re.search(result_re,cl).group(1)
            cl = result_re.sub('',cl)

        # Get func name

        cl=funct.sub('',cl)
        out.name=re.search(re.compile('\w+'),cl).group()

        # Check to see if there are any arguments

        if re.search(r'\([^\)]+',cl)!=None:
            has_args=1
        else:
            has_args=0

        if has_args:
            # get argument list

            cl=re.sub('\w+','',cl,count=1)
            argl=re.split('[\W]+',cl)
        
            del(argl[0])
            del(argl[len(argl)-1])
        
            while cl.strip() == '' or re.search('&',cl)!=None:
                cl=file.next_line()
                if cl.strip() == '': 
                    continue
                arglt=re.split('[\W]+',cl)
                del(arglt[len(arglt)-1])
                for a in arglt:
                    argl.append(a.lower())
        else:
            argl=[]

        argl = map(string.lower, argl)

        # Get next line, and check each possibility in turn

        cl=file.next_line()

        while True:

            # Use statement
            ##check=check_uses(cl,file)
            ##if check[0]!=None:
            ##    out.uses.append(check[0])
            ##    cl=check[1]
            ##    continue

            # Doc comment - return value
            check=check_doc_rv(cl,file)
            if check[0]!=None:
                out.ret_val_doc.append(check[0])
                cl=check[1]
                continue
            
            # Doc comment
            check=check_doc(cl,file)
            if check[0]!=None:
                out.doc.append(check[0])
                cl=check[1]
                continue

            # Interface section
            check=check_interface_decl(cl,file)
            if check[0] != None:
                for a in check[0].procedures:
                    out.arguments.append(a)
                cl = check[1]
                continue

            # Argument
            check=check_arg(cl,file)
            if check[0]!=None:
                for a in check[0]:
                    out.arguments.append(a)
                    cl=check[1]
                continue

            # args_str arguments
            check=check_args_str(cl,file, out.args_str)
            if check[0]!=None:
                cl=check[1]
                continue
                

            m = re.match(funct_end,cl)

            if m == None:
                cl=file.next_line()
                continue
            
            elif m.group(1).lower() == out.name.lower() or m.group(1) == '':
                break

            cl = file.next_line()

        # Select only those declarations that match entries
        # in argument list

        ag_temp=[]

        #if has_args:
        #    t_re_str='(^'
        #    for a in argl:
        #        t_re_str=t_re_str+a+'$)|(^'
        #   t_re_str=t_re_str[:-3]
        #    t_re=re.compile(t_re_str,re.IGNORECASE)
            
        name_re=re.compile(out.name,re.IGNORECASE)
        
        for i in out.arguments:
            if has_args and i.name.lower() in argl and \
                   len([a for a in ag_temp if a.name.lower() == i.name.lower()]) == 0:
                ag_temp.append(i)
            if i.name.lower() in out.args_str:
                out.args_str[i.name.lower()].type = i.type
            if re.search(name_re,i.name)!=None:
                out.ret_val=i
            if ret_var != None and i.name.lower() == ret_var.lower():
                out.ret_val=i

        out.arguments=ag_temp
        out.arguments.sort(key=lambda x:argl.index(x.name.lower()))

        cl=file.next_line()

        # jrk33 - if we're holding a doc comment from before
        # subroutine definition, spit it out now
        if grab_hold_doc and hold_doc is not None:
            for line in hold_doc:
                out.doc.append(line)
            hold_doc = None

        return [out,cl]
    else:
        return [None,cl]
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_type(cl,file):

    global hold_doc

    out=C_type()

    if re.match(type_re,cl)!=None:
        
        # jrk33 - see if it's a global variable of this type.
        # if so, do nothing - it will be found by check_decl
        if decl.match(cl) != None:
            return [None,cl]

        if hold_doc != None:
            for line in hold_doc:
                out.doc.append(line)
            hold_doc = None


        # Get type name
        cl=type_re.sub('',cl)
        out.name=re.search(re.compile('\w+'),cl).group()

        # Get next line, and check each possibility in turn

        cl=file.next_line()

        while re.match(type_end,cl)==None:
            check=check_doc(cl,file)
            if check[0]!=None:
                out.doc.append(check[0])
                cl=check[1]
                continue

            check=check_decl(cl,file)
            if check[0]!=None:
                for a in check[0]:
                    out.elements.append(a)
                cl=check[1]
                continue

            cl=file.next_line()

        cl=file.next_line()

        return [out,cl]
    else:
        return [None,cl]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def check_interface(cl,file):

    global hold_doc
    
    out=C_interface()

    if re.match(iface,cl) != None:
        
        cl = iface.sub('',cl)
        out.name=cl.strip()

        #if out.name == '':
        #    return [None, cl]

        if hold_doc is not None:
            for line in hold_doc:
                out.doc.append(line)
            hold_doc = None

        cl = file.next_line()
        while re.match(iface_end, cl) == None:

            check = check_doc(cl, file)
            if check[0] != None:
                out.doc.append(check[0])
                cl = check[1]
                continue

            check = check_prototype(cl, file)
            if check[0] != None:
                for a in check[0]:
                    out.procedures.append(a)
                cl = check[1]
                continue

            cl=file.next_line()

        cl=file.next_line()

        return [out,cl]
    
    else:
        return [None,cl]


def check_interface_decl(cl,file):

    out = C_interface()

    if cl and re.match(iface,cl) != None:
        
        cl = file.next_line()
        while re.match(iface_end, cl) == None:

            # Subroutine declaration
            check=check_subt(cl,file,grab_hold_doc=False)
            if check[0]!=None:
                out.procedures.append(check[0])
                cl=check[1]
                continue
                
            # Function declaration
            check=check_funct(cl,file,grab_hold_doc=False)
            if check[0]!=None:
                out.procedures.append(check[0])
                cl=check[1]
                continue
            
            cl=file.next_line()

        cl=file.next_line()

        return [out,cl]
    
    else:
        return [None,cl]


def check_prototype(cl,file):

    m = prototype.match(cl)
    if m != None:
        out = map(string.strip, map(string.lower, m.group(1).split(',')))
        cl = file.next_line()
        return [out, cl]

    else:
        return [None, cl]
    
        

            

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++


def check_decl(cl,file):

    out=[]

    if re.match(decl,cl)!=None:
        tp=re.match(types_re,cl).group()
        atr=re.search(attr_re,cl)
        if atr!=None:
            atrl=s_attrib_re.findall(atr.group())
            for j in range(len(atrl)):
                atrl[j]=atrl[j].rstrip()
        else:
            atrl=[]
        m = re.search(d_colon,cl)
        if m is not None:
            names=cl[m.end():]
        else:
            names=types_re.sub('',cl)
        
        # old line - doesn't handle array constants
        # nl=re.split(r'\s*,\s*',names)
        nl=split_attribs(names)


        alist=[]
        for j in range(len(atrl)):
            alist.append(atrl[j])
            
        cl=file.next_line()
        check=check_doc(cl,file)

        dc=[]
        while check[0]!=None:
            # Doc comment
            dc.append(check[0])
            cl=check[1]
            check=check_doc(cl,file)

        cl=check[1]

        for i in range(len(nl)):
            nl[i]=nl[i].strip()
            nlv=re.split(r'\s*=\s*',nl[i])

            temp=C_decl()
            if len(nlv)==2:
                temp.value=nlv[1]

            names, sizes = splitnames(nlv[0])
            temp.name=names[0]
            temp.type=tp
            temp.doc=dc
            temp.attributes=alist[:]
            if sizes[0] != '':
                temp.attributes.append('dimension'+sizes[0])
            out.append(temp)

        return [out,cl]
    else:
        return [None,cl]
    
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def splitnames(names):
   nl = []
   sizes = []
   b = 0
   namestart = 0
   sizestart = 0
   name = ''
   size = ''
   for i, n in enumerate(names):
      if n == '(':
         b += 1
         size += '('
      elif n == ')':
         b -= 1
         size += ')'
      elif n == ',' and b == 0:
         nl.append(name)
         name = ''
         sizes.append(size)
         size = ''
      elif b == 0:
         name += n
      else:
         size += n

   nl.append(name)
   sizes.append(size)

   return nl,sizes


def check_arg(cl,file):

    out=[]

    if cl and re.match(decl_a,cl)!=None:
        tp=re.match(types_re,cl).group()
        m = re.search(d_colon,cl)
        if m is not None:
            atr_temp=cl[re.match(types_re,cl).end():m.start()]
            names=cl[m.end():]
        else:
            atr_temp= ''
            names=types_re.sub('',cl)

        atrl=split_attribs(atr_temp)

#        names=cl[re.search(d_colon,cl).end():]
##        nl=re.split(',',names)
##        for i in range(len(nl)):
##            nl[i]=nl[i].strip()


        # jrk33 - added code to cope with array declarations with
        # size after variable name, e.g. matrix(3,3) etc.

        # Remove values
        names = re.sub(r'=.*$','',names)

        nl, sizes = splitnames(names)

        alist=[]
        for j in range(len(atrl)):
            alist.append(atrl[j])

        cl=file.next_line()
        check=check_doc(cl,file)

        dc=[]

        while check[0]!=None:
            # Doc comment
            dc.append(check[0])
            cl=check[1]
            check=check_doc(cl,file)

        cl=check[1]

        for i in range(len(nl)):
            nl[i]=nl[i].strip()
            temp=C_decl()
            temp.name=nl[i]
            temp.doc=dc
            temp.type=tp

            temp.attributes=alist[:]

            # Append dimension if necessary
            if sizes[i] != '':
                temp.attributes.append('dimension'+sizes[i])

            out.append(temp)
            
        return [out,cl]
    else:
        return [None,cl]


def check_args_str(cl,file, dest):

    if cl is None:
        return [None,cl]
    m = args_str_re.match(cl)
    if cl and m is not None:
        print "CSG ", cl
        d = C_decl()
        d.name = m.group(2)
        d.value = m.group(4)
        d.var = m.group(5).lower()
        d.doc = [m.group(7)]
        dest[d.var] = d # store in funct.args_str or subt.args_str dictionary
        cl = file.next_line()
        return [True,cl]
    else:
        return [None,cl]


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def split_attribs(atr):

    atr=atr.strip()
    if re.match('[,]',atr)!=None:
        atr=atr[1:]
        atr=atr.strip()

    atrc=atr
    bk=0
    atrl=[]

    for i in range(len(atrc)):
        if atrc[i]=='(':
            bk=bk+1
            if bk==1:
                continue
        if atrc[i]==')':
            bk=bk-1
        if bk>0:
            atrc=atrc[:i]+'0'+atrc[i+1:]

    while re.search('[,]',atrc)!=None:
        atrl.append(atr[:re.search('[,]',atrc).start()]) # jrk33 changed [\s,] to [,]
        atr=atr[re.search('[,]',atrc).end():]
        atrc=atrc[re.search('[,]',atrc).end():]
        
    if atr!='':
        atrl.append(atr)

    return map(string.strip,atrl) # jrk33 added strip

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++

def main():

    # +++++++++++++++++++++
    # Check the argument list and retrieve
    # +++++++++++++++++++++

    all_args=sys.argv[1:]

    try:
        optlist, args = getopt.getopt(all_args, 'sa:t:hni:lfp', 'help')
    except getopt.GetoptError:
        # print help information and exit:
        print "f90doc: use -h or --help for help"
        sys.exit(2)

    doc_title='Title'
    doc_author='Authors'
    do_short_doc=False
    do_latex=False
    intro='intro'
    do_header = True

    for o,a in optlist:
        if o=='-t':
            doc_title=a
        if o=='-a':
            doc_author=a
        if o=='-s':
            do_short_doc=True
        if o=='-i':
            intro=a
        if o=='-n':
            intro=None
        if o in ('-h','--help'):
            usage()
            sys.exit()
        if o in '-l':
            do_latex = True
        if o in '-p':
            do_header = False

    if len(args) < 1:
        debug('You need to supply at one least argument!')
        sys.exit()


    programs, modules, functs, subts = read_files(args)

    if do_latex:
        write_latex(programs, modules, functs, subts,
                    doc_title, doc_author, do_short_doc, intro, do_header)

def read_files(args):

    global hold_doc

    programs=[]
    modules=[]
    functs=[]
    subts=[]

    
    for fn in args:

        fname=fn

        # Open the filename for reading

        debug('processing file '+fname)
        file=f90file(fname)

        # Get first line
        
        cline=file.next_line()
        
        while cline!=None:

            # programs
            check=check_program(cline,file)
            if check[0]!=None:
                debug('  program '+check[0].name)
                programs.append((check[0],fn))
                cline=check[1]
                continue

            # modules
            check=check_module(cline,file)
            if check[0]!=None:
                debug('  module '+check[0].name)
                modules.append((check[0],fn))
                cline=check[1]
                continue

            # jrk33 - hold doc comment relating to next module, subrt or funct 
            check=check_doc(cline,file)
            if check[0]!=None:
                if hold_doc == None:
                    hold_doc = [check[0]]
                else:
                    hold_doc.append(check[0])
                cline=check[1]
                continue

        
            # stand-alone subroutines
            check=check_subt(cline,file)
            if check[0]!=None:
                debug('  subroutine '+check[0].name)
                subts.append((check[0],fn))
                cline=check[1]
                continue
        
            # stand-alone functions
            check=check_funct(cline,file)
            if check[0]!=None:
                debug('  function '+check[0].name)
                functs.append((check[0],fn))
                cline=check[1]
                continue
       
            cline=file.next_line()

    return programs, modules, functs, subts



def write_latex(programs, modules, functs, subts, doc_title, doc_author, do_short_doc, intro, header):
    # Print start
    if os.path.exists('COPYRIGHT'):
        for line in open('COPYRIGHT').readlines():
            print '%'+line.strip()

    if header:
        print r"""
\documentclass[11pt]{article}
\textheight 10in
\topmargin -0.5in
\textwidth 6.5in
\oddsidemargin -0.2in
\parindent=0.3in
\pagestyle{headings}

%Set depth of contents page
\setcounter{tocdepth}{2}

\usepackage {makeidx, fancyhdr, boxedminipage, multind, colortbl, sverb}
\usepackage[dvips]{graphicx}
\pagestyle{fancy}

%\renewcommand{\sectionmark}[1]{\markboth{\thesection.\ #1}}
\renewcommand{\sectionmark}[1]{\markboth{}{#1}}
\fancyhf{}
\fancyhead[R]{\bfseries{\thepage}}
\fancyhead[L]{\bfseries{\rightmark}}
\renewcommand{\headrulewidth}{0.5pt}

\makeindex{general}

\begin{document}


\title{"""+doc_title+r"""}
\date{\today}
\author{"""+doc_author+r"""}
\maketitle

\thispagestyle{empty}

\tableofcontents

% Defines paragraphs
\setlength{\parskip}{5mm}
\setlength{\parindent}{0em}

\newpage
"""

    if intro is not None:
        print r'\include{'+intro+'}'


    for a in programs:
        a[0].latex(0,a[1],short_doc=do_short_doc)

    for a in modules:
        a[0].latex(0,a[1],short_doc=do_short_doc)

    if len(subts)+len(functs) != 0:

        print_line(r'\section{Miscellaneous Subroutines and Functions}')

        subts[0][0].latex(1,subts[0][1],short_doc=do_short_doc)

        for a in subts[1:]:
            a[0].latex(1,a[1],short_doc=do_short_doc)

        for a in functs:
            a[0].latex(1,a[1],short_doc=do_short_doc)


    if header:
        print r"""
\printindex{general}{Index}

\end{document}

"""


if __name__ == "__main__":
    main()

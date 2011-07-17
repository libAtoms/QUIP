"""
mockNDarray makes a lists of ndarrays (or lists)
            immitate an ndarray (without data copy)
"""

# Priithon copyright and licensing notes
# =====================================

# Priithon is a collection of other open-source projects and code written 
# at UCSF (contained in the 'Priithon' subfolder).

# Unless indicated otherwise, files in this project are covered by a BSD-type
# license, included below. (http://www.opensource.org/licenses/bsd-license.php)

# Priithon ships with a copy of FFTW,PYX,... which are provided under the GPL licencse.
# (http://www.opensource.org/licenses/gpl-license.php)
# This might require you to consider all of Priithon being under GPL!
# If you cannot accept this please remove those files.
# (IANAL - some comments are also at 
# http://www.scipy.org/mailinglists/mailman?fn=scipy-user/2003-March/001484.html)

# Individual authors are the holders of the copyright for their code and are
# listed in each file.

# "Sebastian Haase <haase@msg.ucsf.edu>" is the original author of Priithon.


# Priithon package license
# ------------------------

# Copyright (c) 2005 The Regents of the University of California

# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#   a. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer.
#   b. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#   c. Neither the name of the University of California, San Francisco nor 
#      the names of its contributors may be used to endorse or promote products 
#      derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
# THE POSSIBILITY OF SUCH DAMAGE.
            
__author__  = "Sebastian Haase <seb.haase@gmail.com>"
__license__ = "BSD license - see LICENSE file"

import numpy as N

class mockNDarray(object):
    def __init__(self, *arrs):
        def conv(a):
            if hasattr(a,'view'):
                return a.view()   # we use view(), so that we can fix it up
            if hasattr(a,'__len__'):
                return mockNDarray(*a) # recursively mockify lists
            if a is None:
                return mockNDarray()   # None makes an empty mockNDarray
            
            return N.array(a) # use standard conversion to ndarray (copy!) 
            # raise ValueError, "don't know how to mockify %s" %(a,)
        self._arrs = [conv(a) for a in arrs]
        self._mockAxisSet(0)

    def _mockAxisSet(self, i):
        """
        this is the workhorse function, that makes the internal state consistent
        sets:
           self._mockAxis
           self._ndim
           self._shape
        """
        self._mockAxis = i
        if len(self._arrs)==0:
            self._ndim     = 0
            self._shape    = ()
            return
        self._ndim     = 1+max((a.ndim for a in self._arrs))
        self._shape    = [1]*self._ndim

        #fix up shape of sub array so that they all have ndim=self._ndim-1
        #  fill shape by prepending 1s
        #  unless ndim was 0, then set shape to tuple of 0s
        nd1 = self._ndim-1
        for a in self._arrs:
            if a.ndim == nd1:
                continue
            if isinstance(a, mockNDarray):
                if a._ndim ==0:
                    a._ndim= nd1
                    a._shape = (0,)*nd1
            else:
                if a.ndim == 0:
                    a.shape = (0,)*nd1
                else:
                    a.shape = (1,)*(nd1-a.ndim)+a.shape

        # fixup the shape to reflext the "biggest" shape possible, like a "convex hull"
        iSub=0 # equal to i, except for i>_mockAxis: its one less
        for i in range(self._ndim):
            if i == self._mockAxis:
                self._shape[i] = len(self._arrs)
                continue
            for a in self._arrs:
                # OLD: the a.ndim>iSub check here means, that sub arrays may be "less dimensional" then `self` would imply if it was not "mock"
                # OLD:   if a.ndim <= self._ndim and a.ndim>iSub:
                if self._shape[i] < a.shape[iSub]:
                    self._shape[i] = a.shape[iSub]
            iSub+=1
        self._shape = tuple( self._shape )

    def _getshape(self):
        return self._shape
    def _setshape(self, s):
        # minimal "dummy" implementation
        #  useful for functions like U.mean2d, which want to set shape to -1,s[-2],s[-1]
        __setShapeErrMsg = "mockNDarray supports only trivial set_shape functionality"
        foundMinus=False
        if len(self._shape) != len(s):
            raise ValueError, __setShapeErrMsg
        for i in range(len(self._shape)):
            if s[i] == -1:
                if foundMinus:
                    raise ValueError, __setShapeErrMsg
                else:
                    foundMinus = True
            elif s[i] != self._shape[i]:
                    raise ValueError, __setShapeErrMsg

    shape = property( _getshape,_setshape )


    def _getndim(self):
        return self._ndim
    ndim = property( _getndim )

    def _getdtype(self):
        return self._ndim and min((a.dtype for a in self._arrs)) or None
    dtype = property( _getdtype )

    def __len__(self):
        return self._shape[0]

    def __getitem__(self, idx):
        import copy
        if isinstance(idx, int):
            if self._mockAxis == 0:
                return self._arrs[idx]
            else:
                s = copy.copy(self)
                s._arrs = [a[idx] for a in self._arrs]
                s._mockAxisSet( self._mockAxis-1 )
                return s
        elif isinstance(idx, tuple):
            if idx == ():
                return self
            if Ellipsis in idx:
                # expand Ellipsis [...] to make slice-handling easier ....
                dimsGiven = len(idx)-1
                for EllipsisIdx in range(len(idx)):
                    if idx[EllipsisIdx] is Ellipsis:
                        break
                idx = idx[:EllipsisIdx ] + (slice(None),)*(self._ndim-dimsGiven) + idx[EllipsisIdx+1:]

            if len(idx) <= self._mockAxis:
                mockIdx = slice(None)
                idxSkipMock = idx
            else:
                mockIdx = idx[self._mockAxis]
                idxSkipMock = idx[:self._mockAxis] + idx[self._mockAxis+1:]
            

            if isinstance(mockIdx, slice):
                s = copy.copy(self)
                
                s._arrs = [a[idxSkipMock] for a in self._arrs[mockIdx]]
                shiftMockAxisBecauseOfInt = sum((1 for i in idx[:self._mockAxis] if not isinstance(i, slice)))
                s._mockAxisSet( self._mockAxis-shiftMockAxisBecauseOfInt )          
                return s
            elif mockIdx is None:
                s = copy.copy(self)
                s._arrs = [a[None][idxSkipMock] for a in self._arrs]
                s._mockAxisSet( self._mockAxis+1 )
                idxSkipMock = (slice(None),)+idxSkipMock  # adjust idxSkipMock to keep new axis 
                return s[idxSkipMock]
                
            else: # mockIdx is "normal" int - CHECK
                # return non-mock ndarray, (or mockNDarray, if there are nested ones)
                return self._arrs[mockIdx][idxSkipMock]

        elif idx is Ellipsis:
            return self
        elif isinstance(idx, slice):
            #raise RuntimeError, "mockarray: slice indices not implemented yet"
            s = copy.copy(self)

            if self._mockAxis ==0:
                s._arrs = self._arrs[idx]
                #s._shape[0] = len(self._arrs)
            else:
                s._arrs = [a[idx] for a in self._arrs[mockIdx]]
                #s._shape[0] = len(self._arrs[0])
            #shiftMockAxisBecauseOfInt = sum((1 for i in idx[:self._mockAxis] if not isinstance(i, slice)))
            s._mockAxisSet( self._mockAxis )
            return s
        elif idx is None: # N.newaxis:
            s = copy.copy(self)
            s._arrs = [a[None] for a in self._arrs]
            s._mockAxisSet( self._mockAxis+1 )
            return s


        raise IndexError, "should not get here .... " 

    def transpose(self, *axes):
        if len(axes) == 1:
            axes = axes[0] # convert  transpose(self, axes) to  transpose(self, *axes)

        if len(axes) != self._ndim:
            raise ValueError, "axes don't match mockarray"

        for newMockAxis in range(self._ndim):
            if axes[newMockAxis] == self._mockAxis:
                break
        else:
            raise ValueError, "axes don't contain mockAxis"

        othersAxes = (ax<newMockAxis and ax or ax-1 for ax in axes[:newMockAxis] + axes[newMockAxis+1:])

        othersAxes = tuple(othersAxes)
        import copy
        s = copy.copy(self)
        s._mockAxisSet(newMockAxis)
        #s._shape = tuple(N.array(s._shape)[list(axes)])

        for i,a in enumerate(s._arrs):
            s._arrs[i] = a.transpose( *othersAxes )

        return s

    def view(self):
        from copy import copy
        return copy(self)

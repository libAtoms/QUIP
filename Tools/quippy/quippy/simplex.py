# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

"""
Direct mimimisation and simulated annealing using the downhill simplex
method of Nelder and Mead.

Code adapted from chapter 10 of Numerical Recipes (3rd edition).
"""

import numpy as np

class Converged(Exception):
    """Exception raised when convergence achieved."""
    pass

class NotConverged(Exception):
    """Exception raised when convergence fails."""
    pass

class DownhillSimplex(object):
    """
    Multidimensional minimisation using the downhill simplex method of
    Nelder and Mead.

    From Numerical Recipes, 3rd edition, section 10.5.
    """

    def __init__(self, func, x0=None, deltas=None, p=None,
                 args=None, ftol=1e-6, nmax=5000, tiny=1e-10):
        """
        Minimise *func(x)* using a DownhillSimplex.

        Either the initial simplex *p* or both the starting position
        *x0* and increments *deltas*. *deltas* can be a single number,
        in which case it is repeated for all the dimensions in
        *x0*. *args* can be used to pass extra arguments to *func*.

        Sets minimisation tolerance to *ftol*, as a fraction of the
        function value (default 1e-6) and maximum number of function
        evaluations to *nmax* (default 5000). *tiny* is estimate of
        machine precision (default 1e-10).
        """
        self.func = func
        if args is None:
            args = ()
        self.args = args
        self.ftol = ftol
        self.NMAX = nmax
        self.TINY = tiny
        self.ybest = np.finfo('d').max
        self.pbest = None
        self.restart(x0, deltas, p)


    def restart(self, x0=None, deltas=None, p=None):
        self.p = p
        if self.p is None:
            if x0 is None or deltas is None:
                raise ValueError('either "p" or both of "x0" and "deltas" must be present')

            try:
                len(deltas)
            except TypeError:
                deltas = [deltas]*len(x0)

            # construct initial simplex
            self.p = np.tile(x0, [len(x0)+1, 1])
            for i in range(1,len(x0)+1):
                self.p[i,i-1] += deltas[i-1]

        self.mpts, self.ndim = self.p.shape
        self.y = np.zeros(self.mpts)

        # initial function evaluations
        for i in range(self.mpts):
            self.y[i] = self.func(self.p[i,:], *self.args)


    def extrapolate(self, ihi, yhi, fac, temperature=None):
        """
        Extrapolate by a factor *fac* through the face of the simplex
        across from the high point *ihi*, try it, and replace the
        high point if the new point is better. Returns result of
        evaluating function at new point.
        """
        fac1 = (1.0 - fac)/self.ndim
        fac2 = fac1 - fac

        ptry = self.p.sum(axis=0)*fac1 - self.p[ihi,:]*fac2
        ytry = self.func(ptry, *self.args)

        # Save best ever value
        if ytry <= self.ybest:
            self.pbest = ptry
            self.ybest = ytry

        if temperature is not None:
            # Subtract thermal fluctuation to give the simplex Brownian motion
            ytry = ytry - (-temperature)*np.log(np.random.uniform())

        if ytry < yhi:
            self.y[ihi] = ytry
            self.p[ihi,:] = ptry

        return ytry

    def step(self, y, temperature=None):
        # Find highest (worst), next-highest and lowest (best) points
        order = y.argsort()
        ilo, ihi, inhi = order[0], order[-1], order[-2]

        # compute fractional range from highest to lowest and
        # return if done, putting best point and value in slot 0
        rtol = 2*(abs(y[ihi] - y[ilo])/
                  (abs(y[ihi]) + abs(y[ilo]) + self.TINY))

        #print 'ylo', y[ilo], 'yhi', y[ihi], 'rtol', rtol

        if rtol < self.ftol:
            self.y[0], self.y[ilo] = self.y[ilo], self.y[0]
            self.p[0,:], self.p[ilo,:] = self.p[ilo,:], self.p[0,:]
            self.fmin = self.y[0]
            raise Converged

        nfunc = 2
        # new iteration. Start by reflecting simplex from the high point.
        ytry = self.extrapolate(ihi, y[ihi], -1.0)
        if ytry <= y[ilo]:
            # new point is better than best point, so try additional
            # extrapolation by factor of 2.
            print 'reflection and expansion', ihi
            ytry = self.extrapolate(ihi, y[ihi], 2.0)
        elif ytry >= y[inhi]:
            # Reflected point is worse than second-highest point, so look
            # for an intermediate lower point, i.e. do a one-dimensional
            # contraction.
            print 'contraction', ihi
            ysave = y[ihi]
            ytry = self.extrapolate(ihi, y[ihi], 0.5)
            if ytry >= ysave:
                # Can't get rid of high point. Contract around the lowest point
                for i in range(self.mpts):
                    if i != ilo:
                        #print 'multiple contraction', i
                        self.p[i,:] = 0.5*(self.p[i,:] + self.p[ilo,:])
                        self.y[i] = self.func(self.p[i,:], *self.args)
                nfunc += self.ndim
        else:
            print 'reflection', ihi
            nfunc -= 1 # correct evaluation count

        return nfunc


    def minimise(self):
        """
        Run the downhill simplex minimiser.

        Returns coordinates at minimum as vector *x*. After minimisation,
        function value at minimum is available as *fmin* attribute,
        simplex as *p* attribute, and number of function evaluations as
        *nfunc* attribute.
        """

        self.nfunc = 0
        while True:
            try:
                self.nfunc += self.step(self.y.copy())
            except Converged:
                return self.p[0,:]
            if self.nfunc >= self.NMAX:
                raise NotConverged("NMAX exceeded")


    def anneal_step(self, iter, temperature):
        """
        Anneal for *iter* steps at temperature *temperature*.

        Returns True if converged to global minimum, False otherwise.
        """

        while True:
            # add positive logarithmically distributed random variable
            # proportional to temperature to function values at each vertex
            y_fluc = self.y + (-temperature)*np.log(np.random.uniform(size=self.mpts))

            try:
                iter -= self.step(y_fluc, temperature)
            except Converged:
                return True

            if iter < 0:
                return False


    def anneal(self, T, m, epsilon, verbose=True, graphics=False, pause=False):
        step = 0
        while not self.anneal_step(m, T):
            T = (1.0 - epsilon)*T
            step += m
            if graphics or verbose:
                ilo = self.y.argsort()[0]
                ihi = self.y.argsort()[-1]
            if graphics:
                if len(gca().patches) >= 1:
                    del gca().patches[-1:]
                if len(gca().lines) >= 1:
                    del gca().lines[-1]
                triplot(ds.p[:,0], ds.p[:,1], 'r.--')
                draw()
                if pause:
                    raw_input('continue...')
            if verbose:
                print 'step %d temperature %f x %s y %f' % (step, T, self.p[ilo,:], self.y[ilo])


if __name__ == '__main__':
    
    do_graphics = True
    if do_graphics:
        from pylab import plot, draw, contour, clf, triplot, gca

    def func(x):
        return np.cos(14.5*x[0]-0.3) + (x[1]+0.2)*x[1] + (x[0]+0.2)*x[0]

    if do_graphics:
        X, Y = np.mgrid[-3:3:100j,-3:3:100j]
        clf()
        contour(X, Y, func((X,Y)), 10)
        draw()

    ds = DownhillSimplex(func, x0=[1.,1.], deltas=.5, ftol=1e-5)
    ds.anneal(2, 1, 0.01, graphics=do_graphics, pause=False)

    # restart minimiser from best position found during annealing
    ds.ftol = 1e-6
    ds.restart(x0=ds.pbest, deltas=.1)
    print ds.minimise()

    # global minimum expected at ~[-0.195, -0.1]

# -*- coding: utf-8 -*-
'''
Module for defining the class related to the curve softener.
'''


# %%
class SmoothCurve:
    '''Creates an instance of an object that defines a curve smoother than the
    input through the :math:`k`-order B-Spline method for interpolation. ::

        SmoothCurve(x, y, k=3, n=300)

    Attributes:
        x (`tuple`, `list` or `numpy.ndarray`): abscisses of the curve to
            smooth.
        y (`tuple`, `list` or `numpy.ndarray`): ordinates of the curve to
            smooth. It must have the same length of ``x``.
        k (`int`): interpolation order.
        n (`int`): number of points of the returned smooth curve

    Note:
        The class ``SmoothCurve`` requires `numpy <http://www.numpy.org/>`_,
        `matplotlib <https://matplotlib.org/>`_. and
        `scipy <https://scipy.org/>`_.

        '''

    def __init__(self, x, y, k=3, n=300):
        '''
        SmoothCurve(x, y, k=3, n=300)
        '''
        from numpy import array
        self.x = array(x)
        self.y = array(y)
        self.k = k
        self.n = n
        # Smoothing the curve
        self.smooth()

    def smooth(self):
        '''Method to generate a smooth curve from the points input through the
        :math:`k`-order B-Spline method.

        Returns:
            (`numpy.ndarray`): :math:`\\left(2 \\times n \\right)` array\
                where :math:`n` is the number of nodes where the path has\
                crossed; the first row of the array contains the abscisses and\
                the second one contains the ordinates of the nodes into the\
                grid-graph.

        Examples:
        '''
        import numpy as np
        from scipy.interpolate import splev

        # Defining the knots vector, with k ending equal knots.
        length = len(self.x)
        t = np.linspace(0, 1, length-self.k+1, endpoint=True)
        t = np.append(np.zeros(self.k), t)
        t = np.append(t, np.ones(self.k))
        # Sequence of length 3 containing the knots, coefficients, and degree
        # of the spline to pass it as the tck argument to splev, the function
        # that will evaluate the curve.
        tck = [t, [self.x, self.y], self.k]
        # Required array of the values of the parameter.
        u = np.linspace(0, 1, self.n, endpoint=True)
        # evaluation
        smoothing = np.array(splev(u, tck))
        # Setting the attribute to the instanced object.
        setattr(self, 'smoothing', smoothing)
        return smoothing

    def plot(self):
        '''Method for generating a graphic of the ``smooth`` method output.
        It allows visually to compare the no smoothed line and its smooothed
        version.

        Returns:
            (`matplotlib.figure.Figure`): object with the matplotlib structure\
                of the plot. You might use it to save the figure for example.

        '''
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # plot original points
        ax.plot(self.x, self.y, 'k--', lw=0.5, marker='x',
                label='Original line')
        # plot smoothed line
        ax.plot(self.smoothing[0], self.smoothing[1], 'k', lw=1.5,
                label='Smoothed curve')
        ax.grid(True, ls='--', lw=0.5)
        ax.legend()
        ax.axis('equal')
        fig.tight_layout()
        return fig


# %%
'''
BSD 2 license.

Copyright (c) 2018, Universidad Nacional de Colombia, Exneyder Andres Montoya
    Araque and Ludger O. Suarez-Burgoa.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

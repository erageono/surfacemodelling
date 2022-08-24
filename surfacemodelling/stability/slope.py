# -*- coding: utf-8 -*-
'''
Module for defining the class related to the slope geometry.
'''


# %%
class AnthropicSlope:
    '''Creates an instance of an object that defines the geometrical frame
    of the slope to perform the analysis given the geometric properties of the
    slope. ::

        AnthropicSlope(slopeHeight, slopeDip, crownDist, toeDist, depth=None)

    The geometry of the slope is as follow:

        * It is a right slope, i.e. its face points to the right side.
        * Crown and toe planes are horizontal.
        * The face of the slope is continuous, ie, it has not berms.

    Attributes:
        slopeHeight (`int` or `float`): Height of the slope, ie, vertical
            length betwen crown and toe planes.
        slopeDip ((2, ) `tuple`, `list` or `numpy.ndarray`): Both horizontal
            and vertical components of the slope inclination given in that
            order.
        crownDist (`int` or `float`): Length of the horizontal plane in
            the crown of the slope.
        toeDist (`int` or `float`): Length of the horizontal plane in the
            toe of the slope.
        depth (`int` or `float` or `None`): Length of the segment beneath the
            slope toe. ``None`` is the default value; if is None, then the
            maximum depth is calculated.

    Note:
        The class ``AnthropicSlope`` requires `numpy <http://www.numpy.org/>`_
        and `matplotlib <https://matplotlib.org/>`_.

        '''

    def __init__(self, slopeHeight, slopeDip, crownDist, toeDist, depth=None):
        '''
        AnthropicSlope(slopeHeight, slopeDip, crownDist, toeDist, depth=None)
        '''
        from numpy import array

        self.slopeHeight = slopeHeight
        self.slopeDip = array(slopeDip)
        self.crownDist = crownDist
        self.toeDist = toeDist
        # Obtaining the maximum depth of the slope from the toe
        if depth is None:
            self.maxDepth()
        else:
            self.depth = depth
        # Defining the boundary coordinates
        self.defineBoundary()

    def maxDepth(self):
        '''Method to obtain the maximum depth of a slope where a circular
        slope failure analysis can be performed.

        The maximum depth is such that the biggest circle satisfished the
        following conditions:

            * It is tangent to the bottom.
            * crosses both the extreme points at the crown and toe.
            * It is orthogonal to the crown plane.

        Returns:
            (`int` or `float`): Maximum depth of the slope measured vertically\
                from the toe plane.

        '''
        import numpy as np

        # Origin of auxiliar coordinates at the begin of the slope-toe
        # Coordinates of the end of the slope-toe
        extremeToePointVec = np.array([self.toeDist, 0])
        # Horizontal distance of the slope face
        slopeDist = self.slopeHeight * self.slopeDip[0] / self.slopeDip[1]
        # Coordinates of the begin of the slope-crown
        extremeCrownPointVec = (-(slopeDist+self.crownDist), self.slopeHeight)
        # Distance between the two extreme points
        differenceVec = extremeToePointVec - extremeCrownPointVec
        distExtrPts = np.linalg.norm(differenceVec)
        # Radius of the largest circle
        maximumCircleRadius = distExtrPts/2 * distExtrPts/differenceVec[0]
        # Toe depth is the difference between maximum-circle radius and
        # the slope-height
        maxDepth = maximumCircleRadius - self.slopeHeight
        # Setting the attribute to the instanced object.
        setattr(self, 'depth', maxDepth)
        return maxDepth

    def defineBoundary(self):
        '''Method to obtain the coordinates of the boundary vertices of the
        slope and plot it if it is wanted.

        The origin of the coordinates is in the corner of the bottom with the
        back of the slope. The coordinates define a close polygon, ie, the
        first pair of coordinates is the same than the last one.

        Returns:
            (`numpy.ndarray`): Coordinates of the boundary vertices of the\
                slope.

        Examples:
        '''
        import numpy as np

        # Slope-face horizontal projection (horizontal distance)
        slopeDist = self.slopeHeight * self.slopeDip[0] / self.slopeDip[1]
        # Creating the contour
        coords = np.array(
                [[0, 0],
                 [0, (self.depth + self.slopeHeight)],
                 [self.crownDist, self.depth + self.slopeHeight],
                 [self.crownDist + slopeDist, self.depth],
                 [self.crownDist + slopeDist + self.toeDist, self.depth],
                 [self.crownDist + slopeDist + self.toeDist, 0],
                 [0, 0]]).T
        # Setting the attribute to the instanced object.
        setattr(self, 'coords', coords)

        return coords

    def plot(self):
        '''Method for generating a graphic of the slope boundary.

        Returns:
            (`matplotlib.figure.Figure`): object with the matplotlib structure\
                of the plot. You might use it to save the figure for example.

        Examples:
        '''
        from matplotlib import pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        # Plot body
        ax.plot(self.coords[0], self.coords[1], '-k')
        ax.plot(self.coords[0], self.coords[1], '.k')
        # Plot settings
        ax.grid(True, ls='--', lw=0.5)
        ax.set_aspect(1)
        fig.tight_layout()
        return fig


# %%
class NaturalSlope:
    '''Creates an instance of an object that defines the geometrical frame
    of the slope to perform the analysis given the terrain coordinates. ::

        NaturalSlope(terrainCoords, depth=None)

    The geometry of the slope is as follow:

        * It is a right slope, i.e. its face points to the right side.
        * The slope is defined with its surface's coordinates.
        * The surface is defined as a polyline such that each segment's slope\
            are always zero or negative.
        * The coordinates' order is such that the highest (and leftmost) point\
            is the first one, and the lowest (and rightmost) is the last one.

    Attributes:
        terrainCoords (`numpy.ndarray`): (2, n) array with the coordinates of
            the slope surface. It must not be a closed polyline, just the
            terrain segment.
        depth (`int` or `float` or `None`): Length of the segment beneath the
            slope toe. ``None`` is the default value; if is None, then the
            maximum depth is calculated.

    Note:
        The class ``NaturalSlope`` requires `numpy <http://www.numpy.org/>`_
        and `matplotlib <https://matplotlib.org/>`_.

        '''

    def __init__(self, terrainCoords, depth=None):
        '''
        NaturalSlope(terrainCoords, depth=None)
        '''
        self.terrainCoords = terrainCoords
        self.slopeHeight = terrainCoords[1, 0] - terrainCoords[1, -1]
        # Obtaining the maximum depth of the slope from the toe
        if depth is None:
            self.maxDepth()
        else:
            self.depth = depth
        # Defining the boundary coordinates
        self.defineBoundary()

    def maxDepth(self):
        '''Method to obtain the maximum depth of a slope where a circular
        slope failure analysis can be performed.

        The maximum depth is such that the biggest circle satisfished the
        following conditions:

            * It is tangent to the bottom.
            * crosses both the extreme points at the crown and toe.
            * It is orthogonal to the crown plane.

        Returns:
            (`int` or `float`): Maximum depth of the slope measured\
                vertically from the rightmost point of the surface coordinates.
        '''
        import numpy as np
        # Distance between the two extreme points
        differenceVec = self.terrainCoords[:, -1] - self.terrainCoords[:, 0]
        distExtrPts = np.linalg.norm(differenceVec)
        # Radius of the largest circle
        maximumCircleRadius = distExtrPts/2 * distExtrPts/differenceVec[0]
        # Toe depth is the difference between maximum-circle radius and
        # the slope-height
        maxDepth = maximumCircleRadius - self.slopeHeight
        # Setting the attribute to the instanced object.
        setattr(self, 'depth', maxDepth)
        return maxDepth

    def defineBoundary(self):
        '''Method to obtain the coordinates of the boundary vertices of the
        slope and plot it if it is wanted.

        The origin of the coordinates is in the corner of the bottom with the
        back of the slope. The coordinates define a closed polygon, ie, the
        first pair of coordinates is the same than the last one. That closed
        polygon is sorted clockwise.

        Returns:
            (`numpy.ndarray`): Coordinates of the boundary vertices of the\
                slope.

        '''
        import numpy as np

        # Obtaining the origin vector to move the surface
        originVec = [self.terrainCoords[0, 0],
                     self.terrainCoords[1, 0] - self.slopeHeight - self.depth]

        # Creating the contour
        coords = self.terrainCoords.T - np.array(originVec)
        extraCoords = np.array([[coords[-1, 0], 0], [0, 0]])
        coords = np.vstack(([0, 0], np.vstack((coords, extraCoords)))).T

        # Setting the attribute to the instanced object.
        setattr(self, 'coords', coords)
        return coords

    def plot(self):
        '''Method for generating a graphic of the slope boundary.

        Returns:
            (`matplotlib.figure.Figure`): object with the matplotlib structure\
                of the plot. You might use it to save the figure for example.

        '''
        from matplotlib import pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        # Plot body
        ax.plot(self.coords[0], self.coords[1], '-k')
        ax.plot(self.coords[0], self.coords[1], '.k')
        # Plot settings
        ax.grid(True, ls='--', lw=0.5)
        ax.set_aspect(1)
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

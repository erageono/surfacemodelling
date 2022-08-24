import numpy as np
import pyvista as pv


class Normals:
    @staticmethod
    def simple(surface):
        """Beregne normaler på hver punkt av et overflatte

        Simple calculation, prone to error with more complicated datasets

        Args:
            surface (pyvista surface): pyvista modell av overflatte

        Returns:
            surface.points: list av punkter
            normals: list av vectorer (3D) av normals
        """
        surface.compute_normals(inplace=True)
        normals = surface["Normals"]

        return surface.points, normals

    @staticmethod
    def alpha(surface, alpha=100):
        """Beregne normaler på et simplifisert overflatte hvor trekanter med kant lengde
        mere en et "alpha" lengde er ignorert

        Seems to be the most sensitive, particularly sensitive with lines of points
        Some parts of the interpolation seem unfounded

        Args:
            surface (pyvista surface): pyvista modell av overflatte

        Returns:
            surface.points: list av punkter
            normals: list av vectorer (3D) av normals
        """
        pts = surface.points
        pts = np.asarray(pts)
        points = pv.PolyData(pts)
        remade = points.delaunay_2d(alpha=100)
        return Normals.simple(remade)

    @staticmethod
    def decimated_smoothed(surface):
        """Beregne normaler på et simplifisert overflatte hvor
        overflatte er decimert og desimert og glattet

        Generates orientations based off of a geometrically visible regional trend, seems less
        sensitive than alpha method, however some parts of the interpolation seem unfounded

        Args:
            surface (pyvista surface): pyvista modell av overflatte

        Returns:
            surface.points: list av punkter
            normals: list av vectorer (3D) av normals
        """
        decimated = surface.decimate(
            target_reduction=0.5
        )  # redusere mesh trekanter til 50% av original
        smoothed = decimated.smooth(n_iter=1000, relaxation_factor=1)  # glatt desimert overflatte
        return Normals.simple(smoothed)

    def nearest_neighbour(self):
        # TODO:
        """
        from sklearn.neighbors import NearestNeighbors
        https://github.com/cgre-aachen/gempy/pull/548/files
        """
        pass

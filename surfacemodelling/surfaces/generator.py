import numpy as np
import pyvista as pv

from ..utilities import get_boundary_limits

                    
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


class Surface:
    @staticmethod
    def simple(points, name="name", save=False):
        """Enkelt triangulasjon metode


        Simple triangulation with pyvista's delaunay_2d function

        Args:
            points (list): [x,y,z]


        Returns:
            surface.points: list av punkter
            normals: list av vectorer (3D) av normals
        """

        points = np.asarray(points)
        if len(points) < 3:
            # TODO: plane that connects 2 points or is planar at 1 z coord
            print("Kan ikke lag overflatte med bare ", len(points), " punkt(er)!")
        elif len(points) == 3:
            verticies = points
            faces = np.hstack([3, 0, 1, 2])
            surf = pv.PolyData(verticies, faces)

        else:
            cloud = pv.PolyData(points)
            surf = cloud.delaunay_2d()
        if save:
            pv.save_meshio(name + ".ply", surf)
            print(str(name + ".ply"), "lagret.")
        return surf

    @staticmethod
    def universal_krige(
        points,
        bounds = 0,
        name="name",
        save=False,
        resolution=1,
        nlags=8,
        weight=True,
        variogram_model="spherical",
        drift_terms = ["regional_linear"],
        evaluate=False,
        plot_variogram=False,
        stddev = False,
        returngrid=False
    ):
        buffer=[100, 100, 10]
        from pykrige.uk import UniversalKriging
        import numpy as np

        if len(points)<5:
            evaluate = False
            print("Kriging med standardparameterer")
        if evaluate:
            from ..surfaces.parameter_tuning import evaluate_kriging
            bestparams = evaluate_kriging(points)
            variogram_model = bestparams["variogram_model"]
            nlags = bestparams["nlags"]
            weight = bestparams["weight"]
        points
        data = np.asarray(points) if isinstance(points,list) else points
        if bounds == 0:
            simplesurf = Surface.simple(points, save=False)
            bounds = get_boundary_limits(simplesurf, buffer)

        gridx = np.arange(bounds[0], bounds[1], resolution)
        gridy = np.arange(bounds[2], bounds[3], resolution)

        UK = UniversalKriging(
            data[:, 0],
            data[:, 1],
            data[:, 2],
            variogram_model=variogram_model,
            nlags=nlags,
            weight=weight,
            drift_terms=drift_terms,
            enable_plotting=plot_variogram,
            exact_values=True
        )
        print("Kriging")
        z, ss = UK.execute("grid", gridx, gridy)

        if returngrid==True:
            return z,ss

        z = np.squeeze(np.array(z))


        UK_xyz = []

        for i_y, row in enumerate(z):
            for i_x, zval in enumerate(row):
                UK_xyz.append([gridx[i_x], gridy[i_y], zval])

        UK_xyz = np.asarray(UK_xyz)

        prediction = pv.PolyData(UK_xyz)
        print("Meshing kriged grid")
        surf = prediction.delaunay_2d()
        print("Mesh created")

        if save == True:
            pv.save_meshio(name + ".ply", surf)
            print(str(name + ".ply"), "lagret.")
        if stddev == True:
            ss = np.array(ss).flatten()
            sd = [point**0.5 for point in ss]
            return surf, sd
        return surf

    @staticmethod
    def gempy(
        points,
        bounds = 0,
        buffer=[100, 100, 10],
        name="name",
        save=False,
        orientation_method="simple",
        plot_intermediate=False,
    ):
        import gempy as gp
        import vtk

        try:
            simplesurf = Surface.simple(points, save=False)
        except:
            print("too few points")
        if bounds == 0:
            simplesurf = Surface.simple(points, save=False)
            bounds = get_boundary_limits(simplesurf, buffer)
        model = gp.create_model(name)
        model = gp.init_data(model, extent=bounds, resolution=[50, 50, 50])

        gp.set_interpolator(model, theano_optimizer="fast_compile", verbose=[])

        model.set_default_surfaces()

        for point in points:
            model.add_surface_points(
                X=point[0], Y=point[1], Z=point[2], surface="surface1"
            )

        if orientation_method == "simple":
            if len(points) == 3:
                opt = points[0]
                a = [
                    points[0][0] - points[1][0],
                    points[0][1] - points[1][1],
                    points[0][2] - points[1][2],
                ]
                b = [
                    points[0][0] - points[2][0],
                    points[0][1] - points[2][1],
                    points[0][2] - points[2][2],
                ]
                orientation = np.cross(a, b).tolist()
                x = orientation[0]
                y = orientation[1]
                z = orientation[2]
                model.add_orientations(
                    X=opt[0],
                    Y=opt[1],
                    Z=opt[2],
                    surface="surface1",
                    pole_vector=(x, y, z),
                )
            else:
                orientation_points, orientations = Normals.simple(simplesurf)
        elif orientation_method == "alpha":
            orientation_points, orientations = Normals.alpha(simplesurf)
        elif orientation_method == "decimated_smoothed":
            orientation_points, orientations = Normals.decimated_smoothed(simplesurf)

        if len(points) != 3:
            for opt, ori in zip(orientation_points, orientations):
                x = ori[0]
                y = ori[1]
                z = -1 * ori[2]
                i = 0
                if z / (x + y) > 0.176:
                    i += 1
                    model.add_orientations(
                        X=opt[0], Y=opt[1], Z=opt[2], surface="surface1", pole_vector=(x, y, z),
                    )

        solution = gp.compute_model(model)
        vertices, simplices = gp.get_surfaces(solution)

        if plot_intermediate:
            gp.plot_3d(model, show_surfaces=True)
        for e, values in enumerate(vertices):
            # setup points and vertices
            Points = vtk.vtkPoints()
            Triangles = vtk.vtkCellArray()
            Triangle = vtk.vtkTriangle()
            for p in values:
                Points.InsertNextPoint(p)

            for i in simplices[e]:
                Triangle.GetPointIds().SetId(0, i[0])
                Triangle.GetPointIds().SetId(1, i[1])
                Triangle.GetPointIds().SetId(2, i[2])

                Triangles.InsertNextCell(Triangle)

            polydata = vtk.vtkPolyData()
            polydata.SetPoints(Points)
            polydata.SetPolys(Triangles)

            surf = pv.PolyData(polydata)

        if save:
            pv.save_meshio(name + ".ply", surf)
            print(str(name + ".ply"), "lagret.")
        return surf

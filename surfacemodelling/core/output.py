import numpy as np
import pyvista as pv
from simplification.cutil import simplify_coords

from ..utilities import distfromline, get_bounds_xy, get_boundary_limits

def get_segments(points):
    segments = []
    for i in range(len(points)-1):
        if len(points[i]) == 2:
            a = [round(points[i][0],0),round(points[i][1],0),-11000]
        else: 
            a = points[i]
        if len(points[i+1]) == 2:
            b = [round(points[i+1][0],0),round(points[i+1][1],0),9000]
        else: 
            b = points[i+1]
        segments.append([a,b])
    return segments

def get_segment_line(segment):
    """Get line and limits object used for cutting meshes
    Args:
        coords (list): [(x1,y1,z1),(x1,y1,z1)]

    Returns:
        line: PyVista PolyData (line) object
        limits: PyVista PolyData (box) object

    """
    a = segment[0]
    b = segment[1]
    a[-1] = 10000
    b[-1] = -10000
    pts = np.stack((a, b))
    line = pv.Spline(pts)
    bounds = line.bounds
    limits = pv.Box(bounds=bounds)
    return line, limits

def distance_from_segment(segment, point):
    ax = segment[0][0]
    ay = segment[0][1]
    az = point[2]

    bx = segment[1][0]
    by = segment[1][1]      
    bz = point[2]

    px = point[0]
    py = point[1]
    pz = point[2]

    a = np.asarray([ax,ay,az])
    b = np.asarray([bx,by,bz])
    p = np.asarray([px,py,pz])

    # normalized tangent vector
    d = np.divide(b - a, np.linalg.norm(b - a))

    # signed parallel distance components
    s = np.dot(a - p, d)

    # clamped parallel distance
    # signed parallel distance components
    t = np.dot(p - b, d)
    h = np.maximum.reduce([s, t, 0])
    # perpendicular distance component
    c = np.cross(p - a, d)

    return round(np.hypot(h, np.linalg.norm(c)),2), round(abs(s),2)

def get_section_XY(segments, surface,return_slc = False):
    """Get lists of relative X and  Y points of a slice of a surface
    Args:
        segments (list): list of segments of the cross section line
        surface (PolyData): PyVista polydata surface object
        return_slc (bool): setting to return slice PolyData object
    Returns:
        X (list): x position of points, starts at 0 
        Y (list): y position of points, corresponds to actual height.
    """
    #empty arrays for x and Y coords of line
    Xglobal = []
    Yglobal = []
    startpoint = 0
    if return_slc == True:
        slces = []
    for segment in segments:
        a = segment[0]
        b = segment[1]
        xlim = (abs(a[0]-b[0])**2+abs(a[1]-b[1])**2)**0.5
        xlim = round(xlim,2)
        line,limits = get_segment_line(segment)
        #cut the line of the surface and trim it to the length of the line
        slc = surface.slice_along_line(line)
        slc = slc.clip_surface(limits)
        if return_slc == True:
            slces.append(slc)
        slicepoints = slc.points
        #initialize local xy meshes
        X = []
        Y = []
        for pt in slicepoints:
            x = distfromline(a,b,pt)
            if x < 0:
                X.append(0+startpoint)
            elif x > xlim:
                X.append(xlim+startpoint)
            else:
                X.append(round(x,2)+startpoint)
        Y = slicepoints[:, 2]
        Y = [round(y, 2) for y in Y]
        XY = list(zip(X,Y))
        XY = list(set(XY))
        XY = sorted(XY)
        startx = XY[0][0]
        endx = XY[-1][0]
        cleaned = []
        for pt in XY:
            if pt in [XY[0], XY[-1]]:
                cleaned.append(pt)
            elif pt[0] not in [startx, endx]:
                cleaned.append(pt)
        X = []
        Y = []
        for pt in cleaned:
            X.append(pt[0])
            Y.append(pt[1])

        for x,y in zip(X,Y):
            Xglobal.append(x)
            Yglobal.append(y)
        startpoint += xlim

    if return_slc == True:
        return Xglobal,Yglobal, slces
    else:
        return Xglobal,Yglobal

def simplify_XY(X,Y,e=0.1):
    """Simplify coordinates
    Args:
        X (list): x position of points, starts at 0 
        Y (list): y position of points, corresponds to actual height.
        e (float, optional): permited deviation of simplified line 

    Returns:
        Xout (list): simplified x position of points, starts at 0 
        Yout (list): simplified y position of points, corresponds to actual height.
    """
    coords = [(x,y) for x,y in zip(X,Y)]
    Xout = []
    Yout = []

    simpcords = simplify_coords(coords, e)
    for pt in simpcords:
        Xout.append(pt[0])
        Yout.append(pt[1])
    return Xout, Yout

def calculate_raster_placement(x_y,referance_pixel,pixels_per_meter,size):
    dxf_x = x_y[0]-((referance_pixel[0])/pixels_per_meter)
    dxf_y = x_y[1]-((size[1]-referance_pixel[1])/pixels_per_meter)
    return (dxf_x,dxf_y)


class Section:
    @staticmethod
    def xy(points,surface):
        """Get  X and Y 2d coordinate points along surface from input points
        Args:
            points (list): [[x1,y1],[x2,y2]]
            surface: PyVista PolyData mesh

        Returns:
            X (list): X values in 2d coordinates (altså lengde fra startpunktet av snitt)
            Y (list): Y values in 2d coordinates (altså høyde)

        """
        segments = get_segments(points)
        X,Y = get_section_XY(segments,surface)
        return X,Y

    @staticmethod
    def simplified_dxf_trace(points,surface,e=0.1):
        """Get dxf along points cut from surface, simplifed to e
        Args:
            points (list): [[x1,y1],[x2,y2]]
            surface: PyVista PolyData mesh
            e (float or int): maksimal tillat avikk fra cuttt linje

        Returns:
            doc (): ezdxf doc object
        """
        import ezdxf as dxf
        from simplification.cutil import simplify_coords
        X,Y = Section.xy(points,surface)
        coords = [(x,y) for x,y in zip(X,Y)]
        simpcords = simplify_coords(coords, e)
        doc = dxf.new("R2000")
        msp = doc.modelspace()
        msp.add_polyline2d(simpcords)
        return doc

    @staticmethod
    def hoyde_data_csv(infile):

        import ezdxf as dxf
        from simplification.cutil import simplify_coords
        import csv

        inreader = csv.reader(infile, delimiter=";",)

        profile = []
        trace = []
        for i, row in enumerate(inreader):
            if i != 0:
                trace.append((float(row[0]),float(row[1])))
                profile.append((float(row[3]),float(row[2])))

        points = simplify_coords(trace, 0.1)

        simpprofile = simplify_coords(profile, 0.1)

        doc = dxf.new("R2000")
        msp = doc.modelspace()
        msp.add_polyline2d(simpprofile)

        return doc, points


    @staticmethod 
    def plot_position_placement(doc,points,positions,e=0.1):
        segments = get_segments(points)
        msp = doc.modelspace()
        for position in positions:
            y_pos = position["north"]
            x_pos = position["east"]
            z_pos = position["z"]
            #z_pos = round(float(position["ComputedGrid"]["Elevation"]["_text"]))
            point = (x_pos,y_pos,z_pos)
            section_placements = []
            total_distance = 0 
            for segment in segments:
                dist, x = distance_from_segment(segment,point)
                section_placements.append([dist,x+total_distance])
                if len(segments) > 1:
                    a = segment[0]
                    b = segment[1]
                    total_distance += round((abs(a[0]-b[0])**2+abs(a[1]-b[1])**2)**0.5,2)
            closest = 0
            if len(segments) > 1:
                distances = [placement[0] for placement in section_placements]
                closest = distances.index(min(distances))
                placement = section_placements[closest]
            
            else:
                placement = section_placements[0]

            print(placement)
            msp.add_text(position["name"]).set_pos((placement[1], z_pos+4), align='BOTTOM_CENTER')
            msp.add_text("Trukket "+str(placement[0])+" m").set_pos((placement[1], z_pos +1 ), align='BOTTOM_CENTER')
            msp.add_text(str(z_pos)+" m").set_pos((placement[1], z_pos +7 ), align='BOTTOM_CENTER')
            trace = [(placement[1],z_pos),(placement[1],z_pos-10)]
            msp.add_polyline2d(trace)
        return doc



    @staticmethod 
    def plot_totalsounding_rasters(doc,points,instructionpositions,max_distance = 999):
        segments = get_segments(points)
        msp = doc.modelspace()
        for position in instructionpositions:
            x_pos = position.tot.grid_coordinates[0]
            y_pos = position.tot.grid_coordinates[1]
            z_pos = position.tot.grid_coordinates[2]
            raster = position.tot.raster
            #z_pos = round(float(position["ComputedGrid"]["Elevation"]["_text"]))
            point = (x_pos,y_pos,z_pos)
            section_placements = []
            total_distance = 0 
            for segment in segments:
                dist, x = distance_from_segment(segment,point)
                section_placements.append([dist,x+total_distance])
                if len(segments) > 1:
                    a = segment[0]
                    b = segment[1]
                    total_distance += round((abs(a[0]-b[0])**2+abs(a[1]-b[1])**2)**0.5,2)
            closest = 0
            if len(segments) > 1:
                distances = [placement[0] for placement in section_placements]
                closest = distances.index(min(distances))
                placement = section_placements[closest]
            
            else:
                placement = section_placements[0]

            if placement[0] < max_distance:
                msp.add_text(position.name).set_pos((placement[1], z_pos+6), align='BOTTOM_CENTER')
                msp.add_text("Trukket "+str(placement[0])+" m").set_pos((placement[1], z_pos +3 ), align='BOTTOM_CENTER')
                msp.add_text(str(z_pos)+" m").set_pos((placement[1], z_pos +9 ), align='BOTTOM_CENTER')
                
                raster_placement = calculate_raster_placement(
                    (placement[1],z_pos),
                    raster.reference_pixel,
                    raster.pixels_per_meter,
                    raster.size,
                    
                    )
                
                raster_image_def = doc.add_image_def(
                    filename=raster.raster_location, 
                    size_in_pixel=(raster.size[0], raster.size[1])
                    )
                msp.add_image(
                    insert=raster_placement, 
                    size_in_units=(raster.size[0]/raster.pixels_per_meter, raster.size[1]/raster.pixels_per_meter), 
                    image_def=raster_image_def, 
                    rotation=0
                    )
        return doc

    @staticmethod
    def stability(points,surface,material=False,resolution=20,search=False):
        """Calculate stability between points cut from surface
        Args:
            points (list): [[x1,y1],[x2,y2]]
            surface: PyVista PolyData mesh
            e (float or int): maksimal tillat avikk fra cuttt linje

        Returns:
            stability (float): stability calculated for section, default to lower bound shansep material
            fig (matplotlib figure): matplotlib figure object
        """
        from ..stability.main import calculate_stability
        X,Y = Section.xy(points,surface)
        X.remove(X[0])
        X.remove(X[-1])
        Y.remove(Y[0])
        Y.remove(Y[-1])
        X,Y = simplify_XY(X,Y,e=0.5)
        if Y[-1]>Y[0]:
            print(X)
            xmax =max(X)
            newX = [round(-x+xmax,2) for x in X]
            X = newX
            print(X)
        slopecoords = [X,Y]
        if material==False:
            material = {"cohesion":0, "frictAngle":0,"unitWeight":19,"wtUnitWeight":19,"adp_factors":{"Aa":1.0,"Ad":0.65,"Ap":0.33},"shansep_factors":{"OCR":1.0,"a":0.25,"m":0.65}}
        stability,fig = calculate_stability(slopecoords,material,resolution=resolution,check_all=search)
        return stability,fig

    @staticmethod
    def three_d_view(points,surface,bounds=False):
        """Calculate stability between points cut from surface
        Args:
            points (list): [[x1,y1],[x2,y2]]
            surface: PyVista PolyData mesh

        Returns:
            section_three_d: pyvista plotter object
        """

        #slices
        segments = get_segments(points)
        X,Y,slces = get_section_XY(segments,surface,return_slc=True)

        #Background
        if bounds == False:
            outergeom = surface.extract_geometry()
            outerpoints = outergeom.points
            x = outerpoints[:, 0]
            y = outerpoints[:, 1]
            z = outerpoints[:, 2]
            bounds = [
                min(x),
                max(x),
                min(y),
                max(y),
                min(z),
                max(z),
            ]
        # texture = Texture.topo_map(bounds)
        # surface.texture_map_to_plane(inplace=True)


        pv.set_plot_theme('document')
        section_three_d=pv.Plotter()
        for slc in slces:
            spline = pv.PolyData()
            spline.points = slc.points
            cells = np.full((len(slc.points)-1, 3), 2, dtype=np.int_)
            cells[:, 1] = np.arange(0, len(slc.points)-1, dtype=np.int_)
            cells[:, 2] = np.arange(1, len(slc.points), dtype=np.int_)
            spline.lines = cells
            spline = spline.tube(radius=0.5)
            section_three_d.add_mesh(spline, color="red")

        # section_three_d.add_mesh(surface,texture=texture,opacity=0.5)
        section_three_d.show_grid()

        return section_three_d

    @staticmethod
    def html_snippet(plotter):
        from ipywidgets.embed import embed_snippet, embed_minimal_html
        from .pv_pythreejs import convert_plotter
        pythreejs_renderer = convert_plotter(plotter)
        embed_minimal_html("test.html",views=[pythreejs_renderer])
        #return embed_snippet(views=[pythreejs_renderer])

class InterpolatedSurface():
    def plot_points(points):
        return [pv.Sphere(radius=0.5,center=point) for point in points]
            
    def universial_krige(points,variogram_model="power",plot_variogram=True):
        from ..surfaces.generator import Surface
        return Surface.universal_krige(
            points,
            variogram_model=variogram_model,
            evaluate=False,
            plot_variogram=plot_variogram,
        )
    def surface_dxf(surface):
        import ezdxf 
        doc = ezdxf.new('R2000')
        msp = doc.modelspace()
        mesh = msp.add_mesh()
        mesh.dxf.subdivision_levels = 0
        with mesh.edit_data() as mesh_data:
            mesh_data.vertices = surface.points
            mesh_data.faces = surface.faces.reshape(-1, 4)[:,1:]
        return doc

    def plot_surface(
        points,
        plotter=False,
        plot_points=True,
        color = "b",
        variogram_model="power",
        plot_variogram=True,
        dxf = False
    ):
        """Standard method to plot surface generated from points
        Args:
            points (list): [[x1,y1,z1],[x2,y2,z2]]
            plotter: PyVista plotter
            plot_points: bool
            color: color (valgfri)
            dxf = False

        Returns:
            plotter: PyVista plotter
        """
        if not plotter:
            plotter = pv.Plotter()
        surface = InterpolatedSurface.universial_krige(
            points,
            variogram_model=variogram_model,
            plot_variogram=plot_variogram
            )
        plotter.add_mesh(surface,color=color,opacity=0.5)
        if plot_points:
            points = InterpolatedSurface.plot_points(points)
            for point in points:
                plotter.add_mesh(point,color=color)
        if dxf == True:
            return plotter, InterpolatedSurface.surface_dxf(surface)
        return plotter
    


        



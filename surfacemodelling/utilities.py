import numpy as np

def get_boundary_limits(surface, buffer, from_xyz=False, from_xy = False):
    if isinstance(surface,dict):
        kofnames = surface.keys()
        points = []
        for kofn in kofnames:
            for point in surface[kofn]:
                points.append(point)
        surface = points
        from_xyz = True

    elif from_xy == True:
        points = []
        x=[]
        y=[]
        z=[0]
        for point in surface:
            x.append(point[0])
            y.append(point[1])
            if from_xyz == True:
                z.append(point[-1])

    else:
        outergeom = surface.extract_geometry()
        outergeom.plot()
        outerpoints = outergeom.points
        x = outerpoints[:, 0]
        y = outerpoints[:, 1]
        z = outerpoints[:, 2]
    print(min(x))
    bounds = [
        min(x) - buffer[0],
        max(x) + buffer[0],
        min(y) - buffer[1],
        max(y) + buffer[1],
        min(z) - buffer[2],
        max(z) + buffer[2],
    ]
    print(bounds)
    return bounds

def get_bounds_xy(points,buffer):
    X = []
    Y = []
    for point in points:
        X.append(point[0])
        Y.append(point[-1])
    bounds = [
        min(X) - buffer[0],
        max(X) + buffer[0],
        min(Y) - buffer[1],
        max(Y) + buffer[1],
        -buffer[2],
        buffer[2],
    ]
    return bounds
    

def clean_to_xyz(points):
    xyz = []
    for row in points:
        xyz.append([row[1], row[2], row[3]])
    return xyz

def clip(item,surface,above=True,plot=False):
    
    import pyvista as pv
    #if type(surface)
    clipper = pv.Plotter()
    clipper.add_mesh(item, color="r", opacity = 0.5)
    clipper.add_mesh(surface, opacity = 0.5)
    
    clipped = item.clip_surface(surface,invert=above)
    try:
        clipper.add_mesh(clipped,color="g")
    except:
        print("Clipping Fail")
    if plot == True:
        clipper.show()

    return clipped

def distfromline(a,b,pt):
    ax = a[0]
    ay = a[1]
    bx = b[0]
    by = b[1]
    if len(pt) > 3:
        az = pt[3]
        bz = pt[3]
        px = pt[1]
        py = pt[2]
        pz = pt[3]
    elif len(pt) == 3:
        az = pt[2]
        bz = pt[2]
        px = pt[0]
        py = pt[1]
        pz = pt[2]
    a = np.asarray([ax, ay, az])
    b = np.asarray([bx, by, bz])
    p = np.asarray([px, py, pz])

    # normalized tangent vector
    d = np.divide(b - a, np.linalg.norm(b - a))

    # signed parallel distance components
    s = np.dot(a - p, d)

    if len(pt) > 3:
        # clamped parallel distance
        # signed parallel distance components
        t = np.dot(p - b, d)
        h = np.maximum.reduce([s, t, 0])
        # perpendicular distance component
        c = np.cross(p - a, d)

        return np.hypot(h, np.linalg.norm(c)), -s

    elif len(pt) == 3:
        return -s


def simplify_XY(X,Y,e=0.1):
    coords = []
    for x,y in zip(X,Y):
        coords.append((x,y))
    Xout = []
    Yout = []

    simpcords = simplify_coords(coords, e)
    for pt in simpcords:
        Xout.append(pt[0])
        Yout.append(pt[1])
    return X,Y

def find_svg_dim(svg_string):
    SVG_SCALE = 23.4145
    
    import ast
    #Find size
    width_index_start = svg_string.find("width=")
    height_index = svg_string.find("height=")
    
    width = int(float(svg_string[width_index_start+7:height_index-2]))
    
    
    height = int(svg_string[height_index+8:height_index+8+3])
    
    size = [width,height]
    #find referance_pixel
    
    translate_index_start = svg_string.find("translate")
    translate_index_end = svg_string.find(")",translate_index_start)
    translate_tuple = ast.literal_eval(svg_string[translate_index_start+9:translate_index_end+1])
    
    referance_pixel = [int(round(translate_tuple[0],0)),int(round(translate_tuple[1],0))]
    

    return size, referance_pixel, SVG_SCALE
    
ERA_colour_list = ["b","g","r","c","m","y","k","w"]

from numpy import array,arange

from .slope import NaturalSlope
from .slipsurface import CircularSurface
from .slices import MaterialParameters, Slices
from .slopestabl import SlopeStabl

def round_base(x, base):
    return base * round(x/base)

def find_possible_distances(terrain_list,resolution=5):
    
    possible_distances = []
    X = terrain_list[0]
    Y = terrain_list[1]
    n = int(X[-1]/resolution)
    for i in range(n):
        x_i = i*resolution
        y_i = Y[X.index(min(X, key=lambda x:abs(x-x_i)))]
        for i_i in range(n-i):
            x_i_i = x_i + i_i*resolution
            y_i_i = Y[X.index(min(X, key=lambda x:abs(x-x_i_i)))]
            if y_i-y_i_i>5:
                possible_distances.append([x_i,x_i_i])
    possible_distances = [t for t in {tuple(i) for i in possible_distances}]
    possible_distances.sort()
    return possible_distances

def find_possible_radii(distances,resolution):
    min_radius = ((distances[1]-distances[0])/2)+(distances[1]-distances[0])/(resolution/2)+1
    max_radius = min_radius*3
    return arange(min_radius,max_radius,resolution)

def find_minimum_fs(terrain_list,material,resolution=5,check_all=True):
    min_fs = 10
    last_fs = min_fs
    terrainCoords = array(terrain_list)
    slope = NaturalSlope(terrainCoords,depth=50)
    calc_mat = MaterialParameters(
        cohesion=material["cohesion"],
        frictAngle=material["frictAngle"],
        unitWeight=material["unitWeight"],
        wtUnitWeight=material["wtUnitWeight"],
        adp_factors=material["adp_factors"],
        shansep_factors=material["shansep_factors"],
        )

    all_calculations = {}
    possible_distances = find_possible_distances(terrain_list,resolution=resolution)
    if check_all==False:
        possible_distances = [((terrain_list[0][0]),(terrain_list[0][-1])),((terrain_list[0][0])+5,(terrain_list[0][-1]-5))]
    print("Calculating possible failures between:")
    print(possible_distances)
    for distances in possible_distances:
        radii = find_possible_radii(distances,resolution)
        for radius in radii:
            
            surface = CircularSurface(slopeCoords=slope.coords,dist1=distances[0], dist2=distances[1],radius=radius)
            slices = Slices(material=calc_mat, slipSurfCoords=surface.coords,slopeCoords=slope.coords, numSlices=10,watertabCoords=None, bim=None)
            try:
                stabAnalysis = SlopeStabl(slices, seedFS=1, Kh=0,interSlcFunc=1)
                fs = stabAnalysis.__dict__["FS"]["fs"]
                if last_fs < fs:
                    break
                last_fs = fs
                if fs < min_fs and fs >0:
                    min_fs = fs
                    #moved here to try and speed it up
                    all_calculations[fs] = {"obj": stabAnalysis}
            except:
                fs = False
                stabAnalysis = False


    return all_calculations,min_fs


def calculate_stability(terrain_list,material,resolution=20,check_all=False):
    if terrain_list[-1][0]<terrain_list[-1][-1]:
        X = terrain_list[0]
        Y = terrain_list[1]
        X.reverse()
        Y.reverse()
        terrain_list=[X,Y]
    all_calculations,min_fs=find_minimum_fs(terrain_list,material,resolution=resolution,check_all=check_all)
    print("min_fs")
    print(min_fs)
    import matplotlib.pyplot as plt
    
    try:
    
        fig = plt.figure()
        X = terrain_list[0]
        Y  = terrain_list[-1]
        ax = fig.add_subplot(111)
        ax.axis('equal')
        ax.plot(X, Y, 'b-')
        for slice_ in all_calculations[min_fs]["obj"].slices.slices:  # Plotting each slice
            ax.plot(slice_.coords, ':r', lw=0.5)
            print(slice_.coords)
        ax.plot(all_calculations[min_fs]["obj"].slices.slipSurfCoords, '-r')
    except: 
        pass
    return min_fs, fig
                                          
    
    
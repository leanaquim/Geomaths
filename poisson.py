import numpy as np
import matplotlib.pyplot as plt
from mesh import Mesh


def premiers_test_poisson() :
    # Reading the original mesh
    slice = Mesh("./shell/slice.obj")
    horiz0 = Mesh("./shell/horizon0.obj")

    # Extracting the vertices of the horizon we want to move
    z = slice.V[0][2]
    print('z', z)

    f = []
    for point in horiz0.V:
        if point[2] == z:
            f.append(point[1])
    print('f', f)

    n = len(f)
    print('n', n)

    # Defining the matrices
        # Target function
    g = np.array([f[0]] * n)
    print('g', g)

        # 
    b = np.array([[g[i]-g[i-1]] for i in range(1, n)])
    b[0] = g[0]
    b[-1] = g[0]
    print('b', b)

    result = f

    A = np.matrix(np.zeros((n,n)))
    #A[0, 0] = result[1]
    #A[-1, -1] = - result[n-2]
    for i in range(1, n-2):
        #A[i, i] = result[i + 1]
        #A[i, i-1] = - result[i]
        A[i,0] = 1
    # np.fill_diagonal(A,      1)
    # np.fill_diagonal(A[1:], -1)

    # result = [g[0]] + (np.linalg.inv(A.T*A)*(A.T)*b).T.tolist()[0] + [g[-1]]
    result = A.dot(f)
    result = [result[0, i] for i in range(56)]

    # Neanderthal version: 
    # for _ in range(10) :
    #     for i in range(1 , n-1) :
    #         result[i] = (f[i-1] + f[i+1] + (2*g[i]-g[i-1]-g[i+1]) )/2.

    print(A)
    print(result)

    t = np.linspace(0, 1, n)
    plt.plot(t, f, linewidth=5, label="f_init")
    plt.plot(t, g, label="g")
    plt.plot(t, result, label="result")
    plt.legend()
    plt.show()


def extract_border(mesh):
    borders = []
    for v in range(mesh.nverts) :
        if mesh.on_border(v) :
            borders.append(v)
            
    return borders


def flatten_horizon(mesh, horizon):
    modified_mesh = mesh
    list_borders = extract_border(mesh)
    nvert = mesh.nverts
    ncorn = mesh.ncorners
    nborders = len(list_borders)
    nhoriz = len(horizon)

    list_y = np.zeros(ncorn + nhoriz + nborders)
    list_y = [[list_y[i]] for i in range(len(list_y))]
    
    # Defining the 'transition' matrix
    A = np.matrix(np.zeros((ncorn + nhoriz + nborders, nvert)))
        
        # Constraints on horizon vertices
    y0 = mesh.org(horizon[0]) # id of the first point of the horizon
    
        # Laplace
    for i in range(ncorn):
        A[i,mesh.org(i)] = -1
        A[i,mesh.dst(i)] = 1

        # horizon
    for e in range(nhoriz):
        i = mesh.org(horizon[e])
        if y0==i: continue
        A[ncorn + e, y0] = 1*100
        A[ncorn + e, i] = -1*100
        list_y[ncorn + nhoriz + e][0] = 0*100

        #  fix borders
    for i in range(len(list_borders)):
        A[ncorn + len(horizon) + i, list_borders[i]] = 1*10
        list_y[ncorn + nhoriz + i][0] = mesh.V[list_borders[i]][1]*10

    # compute
    result = np.empty(nvert)
    result = (np.linalg.inv(A.T*A)*(A.T)*list_y).T.tolist()[0]

    # edit mesh
    for i in range(nvert) :
        modified_mesh.V[i,1] = result[i]

    # for i in list_borders:
    #     modified_mesh.V[i,2] = .1

    modified_mesh.print_to_file('chevron/modified_mesh.obj')


def flatten_horizons(model, mesh, horizon_list):
    modified_mesh = mesh
    list_borders = extract_border(mesh)
    nvert = mesh.nverts
    ncorn = mesh.ncorners
    nborders = len(list_borders)
    nb_horizons = len(horizon_list)
    nhoriz_global = sum(len(hor) for hor in horizon_list)

    list_y = np.zeros(ncorn + nhoriz_global + nborders)
    list_y = [[list_y[i]] for i in range(len(list_y))]
    
    # Defining the 'transition' matrix
    A = np.matrix(np.zeros((ncorn + nhoriz_global + nborders, nvert)))
    
        # Laplace
    for i in range(ncorn):
        A[i,mesh.org(i)] = -1
        A[i,mesh.dst(i)] = 1

        #  fix borders
    for i in range(len(list_borders)):
        A[ncorn + i, list_borders[i]] = 1*100
        list_y[ncorn + i][0] = mesh.V[list_borders[i]][1]*100
        
        # horizon
    nb_lines_hor = 0
    for hor in range(nb_horizons):
        nhoriz = len(horizon_list[hor])
        y0 = mesh.org(horizon_list[hor][0]) # Constraints on horizon vertices : id of the first point of the horizon
        for e in range(nhoriz):
            i = mesh.org(horizon_list[hor][e])
            if y0==i: continue
            A[ncorn + nborders + nb_lines_hor + e, y0] = 1*100
            A[ncorn + nborders + nb_lines_hor + e, i] = -1*100
            list_y[ncorn + nborders + nb_lines_hor + e][0] = 0*100
        nb_lines_hor += nhoriz

    # compute
    result = np.empty(nvert)
    result = (np.linalg.inv(A.T*A)*(A.T)*list_y).T.tolist()[0]

    # edit mesh
    for i in range(nvert) :
        modified_mesh.V[i,1] = result[i]

    for i in list_borders:
        modified_mesh.V[i,2] = .1

    modified_mesh.print_to_file(model + '/slice_horizons.obj')


def flatten_fault(mesh):
    return False

if __name__ == '__main__' :
    # mesh = Mesh("./chevron/slice.obj")
    # horizon = [258, 267, 273, 285, 294, 307, 317, 816, 825, 837, 849, 861, 870, 882, 3975, 4044, 4056,
    # 4065, 4071, 4074, 4173, 4206, 4272, 4284, 4308, 4317, 4368, 4377, 4452, 4458, 4479, 4485, 4494, 4503,
    # 4512, 4521, 4527, 4548, 4569, 4584, 4587, 4606, 4620, 4647, 4656, 4662, 4677, 4692, 4716, 4734, 4740,
    # 4752, 4791, 4797, 4833, 4851, 4854, 4866, 4879, 4890, 4920, 4932, 4980, 5124, 5130, 5136, 5148, 5154,
    # 5181, 5190, 5523, 5541, 5560, 5583, 6117, 6126, 6144, 6147, 6237, 6279, 6297, 6321, 6369, 6474, 6582,
    # 6600, 6606, 6645, 6672, 6687, 6729, 6768, 6780, 6786, 6837, 6849, 6885, 6909, 6927, 6951, 6957, 6969,
    # 7017, 7035, 7059, 7074, 7095, 7101, 7110, 7122, 7140, 7152, 7182, 7200, 7218, 7230, 7242, 7266, 7281,
    # 7293, 7308, 7329, 7359, 7377, 7404, 7413, 7428, 7437, 7563, 7575, 9084, 9621, 10002, 10692]

    # horizon1 = [90, 93, 152, 159, 164, 170, 174, 200, 207, 217, 1206, 1242, 1477, 1512, 1539, 1599, 1608,
    # 1653, 1659, 1671, 1695, 1707, 1716, 1726, 1737, 1746, 1749, 1752, 1767, 1770, 1773, 1779, 1794, 1800,
    # 1819, 1828, 1830, 1836, 1843, 1845, 1848, 1851, 1860, 1872, 1881, 1896, 1902, 1917, 1923, 1929, 1938,
    # 1941, 1948, 1953, 1962, 1965, 1977, 1986, 1995, 2004, 2013, 2022, 2031, 2040, 2049, 2058, 2067, 2076,
    # 2091, 3963, 3981, 3999, 4014, 4029, 4038, 4062, 4080, 4098, 4116, 4128, 4143, 4158, 4170, 4188, 4203,
    # 4218, 4233, 4245, 4257, 4287, 4299, 4326, 4341, 4362, 4389, 4401, 4467, 4902, 4953, 5010, 5019, 5022,
    # 5031, 5040, 5052, 5058, 5262, 5268, 5274, 5283, 5358, 5382, 5391, 5403, 5412, 5610, 5625, 5640, 5655,
    # 5667, 5676, 5685, 5694, 5703, 5712, 5721, 5730, 5742]

    
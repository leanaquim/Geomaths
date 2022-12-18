import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import lsmr
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


def flatten_fault(model, mesh):
    return False

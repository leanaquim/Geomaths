import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import lsmr
from mesh import Mesh


def extract_border(mesh):
    borders = []
    for v in range(mesh.nverts) :
        if mesh.on_border(v) :
            borders.append(v)
            
    return borders

def extract_faults(mesh, is_fault):
    faults = []
    for v in range(mesh.nverts) :
        if is_fault[v] :
            faults.append(mesh.org(v))
            faults.append(mesh.dst(v))
            
    return faults

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


def flatten_fault(model, mesh, is_fault):
    list_faults = extract_faults(mesh, is_fault)
    modified_mesh = mesh

    for i in list_faults:
        modified_mesh.V[i,2] = .1

    modified_mesh.print_to_file(model + '/slice_faults.obj')
    
    return False

import numpy as np
from scipy.sparse.linalg import lsmr
import scipy.sparse

def extract_border(mesh):
    borders = []
    for v in range(mesh.nverts) :
        if mesh.on_border(v) :
            borders.append(v)
    return borders

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
    A = scipy.sparse.lil_matrix((ncorn + nhoriz_global + nborders, nvert))
    
        # Laplace
    for i in range(ncorn):
        A[i,mesh.org(i)] = -1
        A[i,mesh.dst(i)] = 1
        list_y[i][0] = mesh.V[mesh.dst(i)][1] - mesh.V[mesh.org(i)][1]

        
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
    A = A.tocsr()
    result = lsmr(A, list_y, atol=1e-20, btol=1e-20)[0]

    # Accentuate the borders for visualisation
    for i in list_borders:
        modified_mesh.V[i,2] = .05

    # edit mesh
    for i in range(nvert):
        modified_mesh.V[i,1] = result[i]
    modified_mesh.print_to_file(model + '/slice_horizons.obj')
    
    return modified_mesh

def flatten_fault(model, mesh, is_fault, fault_opposite):
    modified_mesh = mesh
    list_borders = extract_border(mesh)
    nvert = mesh.nverts
    ncorn = mesh.ncorners
    nborders = len(list_borders)
    
    list_x = np.zeros(2 * ncorn + nborders)
    list_x = [[list_x[i]] for i in range(len(list_x))]

    
    # Defining the 'transition' matrix
    A = scipy.sparse.lil_matrix((2 * ncorn + nborders, nvert))
    
        # Laplace
    for i in range(ncorn):
        A[i,mesh.org(i)] = -1
        A[i,mesh.dst(i)] = 1
        list_x[i][0] = mesh.V[mesh.dst(i)][0] - mesh.V[mesh.org(i)][0]


        # verticalize faults
    for i in range(ncorn) :
        if is_fault[i] :
            # Half-edge i = fault
            A[i,mesh.org(i)] = -50
            A[i,mesh.dst(i)] = 50
        
        # Join the two sides of the faults
    for i in range(ncorn) :
        if is_fault[i] :
            A[ncorn + nborders + i, mesh.org(i)] = -1*1
            A[ncorn + nborders + i, mesh.dst(fault_opposite[i])] = 1*1

    # compute
    A = A.tocsr()
    result = lsmr(A, list_x, atol=1e-20, btol=1e-20)[0]

    # Accentuate the borders for visualization
    for i in list_borders:
        modified_mesh.V[i,2] = .05

    # edit mesh
    for i in range(nvert) :
        modified_mesh.V[i,0] = result[i]
    modified_mesh.print_to_file(model + '/slice_faults.obj')
    
    return modified_mesh

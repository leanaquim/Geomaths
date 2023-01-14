import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import lsmr
import scipy.sparse
from mesh import Mesh


def extract_border(mesh):
    borders = []
    for v in range(mesh.nverts) :
        if mesh.on_border(v) :
            borders.append(v)
            
    return borders

def extract_faults(mesh, is_fault):
    nvert = mesh.nverts
    nb_faults = 0
    fault_list = [] # list of corners in each fault
    tag = [-1 for i in range(mesh.ncorners)]

    for half_edge in range(nvert): 
        if is_fault[half_edge] :
            
            if tag[mesh.org(half_edge)] != -1  or tag[mesh.dst(half_edge)] != -1 :
                # if the halfedge links two faults : fusionner failles
                if tag[mesh.org(half_edge)] != -1 and tag[mesh.dst(half_edge)] != -1 :
                    nb_faults += 1 # create a new fault constituted by the two previous ones

                    fault_number1 = tag[mesh.org(half_edge)]
                    fault_number2 = tag[mesh.dst(half_edge)]
                    
                    for corner in fault_list[fault_number1 - 1] :
                        tag[corner] = nb_faults

                    for corner in fault_list[fault_number2 - 1] :
                        tag[corner] = nb_faults

                    fault_list.append([])
                    fault_list[nb_faults - 1] = fault_list[fault_number1 - 1] + fault_list[fault_number2 - 1]
                    fault_list[fault_number1 - 1] = []
                    fault_list[fault_number2 - 1] = []

                # fault already known -> add the corner to the fault
                elif tag[mesh.org(half_edge)] != -1 :
                    fault_number = tag[mesh.org(half_edge)]
                    fault_list[fault_number - 1].append(mesh.dst(half_edge))

                elif tag[mesh.dst(half_edge)] != -1 :
                    fault_number = tag[mesh.dst(half_edge)]
                    fault_list[fault_number - 1].append(mesh.org(half_edge))

            else :
                nb_faults += 1
                fault_list.append([])
            
            fault_list[nb_faults - 1].append(mesh.org(half_edge))

            tag[mesh.org(half_edge)] = nb_faults # tag the two corners linked by the current halfedge
            tag[mesh.dst(half_edge)] = nb_faults

    fault_list = [fault_list[i] for i in range (len(fault_list)) if fault_list[i] != []]
    print("Nb de failles trouvées : ", len(fault_list))
    return fault_list

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
    # A = np.matrix(np.zeros((ncorn + nhoriz + nborders, nvert)))
    A = scipy.sparse.lil_matrix((ncorn + nhoriz + nborders, nvert))
        
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
        A[ncorn + len(horizon) + i, list_borders[i]] = 1*0
        list_y[ncorn + nhoriz + i][0] = mesh.V[list_borders[i]][1]*0

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
    # A = np.matrix(np.zeros((ncorn + nhoriz_global + nborders, nvert)))
    A = scipy.sparse.lil_matrix((ncorn + nhoriz_global + nborders, nvert))
    
        # Laplace
    for i in range(ncorn):
        A[i,mesh.org(i)] = -1
        A[i,mesh.dst(i)] = 1
        list_y[i][0] = mesh.V[mesh.dst(i)][0] - mesh.V[mesh.org(i)][0]

        #  fix borders
    for i in range(len(list_borders)):
        A[ncorn + i, list_borders[i]] = 1*1
        list_y[ncorn + i][0] = mesh.V[list_borders[i]][1]*1
        
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
    # result = np.empty(nvert)
    # result = (np.linalg.inv(A.T*A)*(A.T)*list_y).T.tolist()[0]
    A = A.tocsr()
    result = lsmr(A, list_y, atol=1e-20, btol=1e-20)[0]

    # edit mesh
    for i in range(nvert):
        modified_mesh.V[i,1] = result[i]

    # for i in list_borders:
    #     modified_mesh.V[i,2] = .1

    modified_mesh.print_to_file(model + '/slice_horizons.obj')
    
    return modified_mesh


def flatten_fault(model, mesh, is_fault, fault_opposite):
    modified_mesh = mesh
    list_borders = extract_border(mesh)
    nvert = mesh.nverts
    ncorn = mesh.ncorners
    nborders = len(list_borders)
    
    list_y = np.zeros(ncorn + nborders)
    list_y = [[list_y[i]] for i in range(len(list_y))]
    
    # Defining the 'transition' matrix
    A = np.matrix(np.zeros((ncorn + nborders, nvert)))
    
        # Laplace
    for i in range(ncorn):
        A[i,mesh.org(i)] = -1
        A[i,mesh.dst(i)] = 1

    list_borders = extract_border(mesh)
    list_faults = extract_faults(mesh, is_fault)
    
    nvert = mesh.nverts
    ncorn = mesh.ncorners
    nborders = len(list_borders)
    # nfaults = len(list_faults)                         #############
    
    list_x = np.zeros(2 * ncorn + nborders)
    list_x = [[list_x[i]] for i in range(len(list_x))]
    
    # Defining the 'transition' matrix
    # A = np.matrix(np.zeros((2 * ncorn + nborders, nvert)))
    A = scipy.sparse.lil_matrix((2 * ncorn + nborders, nvert))
    
        # Laplace
    for i in range(ncorn):
        A[i,mesh.org(i)] = -1
        A[i,mesh.dst(i)] = 1
        list_x[i][0] = mesh.V[mesh.dst(i)][0] - mesh.V[mesh.org(i)][0]

        #  fix borders
    for i in range(nborders):
        if not is_fault[i] :
            A[ncorn + i, list_borders[i]] = 1*1
            list_x[ncorn + i][0] = mesh.V[list_borders[i]][0]*1

        # verticalize faults
    for i in range(ncorn) :
        if is_fault[i] :
            # Half-edge i = fault
            A[i,mesh.org(i)] = -50
            A[i,mesh.dst(i)] = 50
        
        # Join the two sides of the faults
    for i in range(ncorn) :
        if is_fault[i] :
            A[ncorn + nborders + i, mesh.org(i)] = -1*2
            A[ncorn + nborders + i, mesh.dst(fault_opposite[i])] = 1*2

    # compute
    # result = np.empty(nvert)
    # result = (np.linalg.inv(A.T*A)*(A.T)*list_x).T.tolist()[0]
    A = A.tocsr()
    result = lsmr(A, list_x, atol=1e-20, btol=1e-20)[0]

    # edit mesh
    for i in range(nvert) :
        modified_mesh.V[i,0] = result[i]

    # for i in range(len(list_faults)):
    #     for j in list_faults[i] :
    #         modified_mesh.V[j,2] = i * 0.01

    modified_mesh.print_to_file(model + '/slice_faults.obj')
    
    return modified_mesh

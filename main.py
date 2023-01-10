import numpy as np 
from mesh import Mesh
import matplotlib.pyplot as plt
import poisson

from chevron_attributes import horizon_id_chevron
from ifp1_attributes import horizon_id_ifp1, is_fault_ifp1, fault_opposite_ifp1
from ifp2_attributes import horizon_id_ifp2, is_fault_ifp2
from shell_attributes import horizon_id_shell


# Stores the indexes of the vertices of faults/horizons
def extract_horizons(mesh, horizon_id):
    """
    @param horizon_id list of horizons computed in the attributes file
    
    Extracts the list of horizons of the model.
    For each horizon found, it stores the index of all half-edges constituing the horizons (line number in the slice.obj file)
    """
    nb_vertices = mesh.ncorners
    nb_horizons = max(horizon_id) + 1
    horizon_list = [[] for _ in range(nb_horizons)]
    for hor in range(nb_horizons):
        for vertex in range(nb_vertices):
            if horizon_id[vertex] == hor :
                horizon_list[hor].append(vertex)
    return horizon_list


def extract_faults(mesh, is_fault):
    """

    """
    nvert = mesh.nverts
    nb_faults = 0
    fault_list = [] # list of corners in each fault
    tag = [-1 for i in range(mesh.ncorners)]

    for half_edge in range(nvert): 
        if is_fault[half_edge] :

            if tag[mesh.org(half_edge)] != -1  or tag[mesh.dst(half_edge)] != -1 :
                # fault already known -> add the corner to the fault
                if tag[mesh.org(half_edge)] != -1 :
                    fault_number = tag[mesh.org(half_edge)]
                    fault_list[fault_number - 1].append(mesh.dst(half_edge))

                if tag[mesh.dst(half_edge)] != -1 :
                    fault_number = tag[mesh.dst(half_edge)]
                    fault_list[fault_number - 1].append(mesh.org(half_edge))

                # if the halfedge links two faults : fusionner failles
                if tag[mesh.org(half_edge)] != -1 and tag[mesh.dst(half_edge)] != -1 :
                    nb_faults += 1 # create a new fault constituted by the two previous ones

                    fault_number1 = tag[mesh.org(half_edge)]
                    fault_number2 = tag[mesh.dst(half_edge)]
                    
                    fault_list.append([])
                    fault_list[nb_faults - 1] = fault_list[fault_number1 - 1] + fault_list[fault_number2 - 1]
                    fault_list[fault_number1 - 1] = []
                    fault_list[fault_number2 - 1] = []
                
            else :
                nb_faults += 1
                fault_list.append([])
            
            fault_list[nb_faults - 1].append(half_edge)

            tag[mesh.org(half_edge)] = nb_faults # tag the two corners linked by the current halfedge
            tag[mesh.dst(half_edge)] = nb_faults

    fault_list = [fault_list[i] for i in range (len(fault_list)) if fault_list[i] != []]

    return fault_list


# def count_nb_faults(mesh, is_fault):
#     nb_faults = 0
#     dict = {}
#     ncorn = mesh.ncorners
#     for i in range(ncorn):
#         if is_fault[i]:
#             if nb_faults not in dict.keys():
#                 dict[nb_faults] = []
#                 nb_faults += 1
#             if mesh.dst(i) in dict[nb_faults]:
#                 dict[nb_faults].append(mesh.org(i))
#                 dict[nb_faults].append(mesh.dst(i))
#     return False


def rec_count(mesh, halfedge, nb_faults, is_fault, tag):
    """
    Counts the number of faults present in the mesh
    """
    next_halfedge = mesh.opposite(mesh.prev(mesh.opposite(halfedge)))
    
    if is_fault[next_halfedge]:
        tag[next_halfedge] = nb_faults
        rec_count(mesh, next_halfedge, nb_faults, is_fault, tag)
    
    else :
        tag[next_halfedge] = -1



def count_nb_faults(mesh, is_fault):
    nb_faults = 0
    ncorn = mesh.ncorners
    tag = np.zeros(ncorn)
    
    for halfedge in range(ncorn):
        if tag[halfedge] == 0: # new halfedge 
            if is_fault[halfedge]: # new fault
                nb_faults += 1
                
                tag[halfedge] = nb_faults
                tag[mesh.opposite(halfedge)] = nb_faults

                rec_count(mesh, halfedge, nb_faults, is_fault, tag)
            else :
                tag[halfedge] = -1
    print("Tag : ", tag)
    return nb_faults


def flatten(model, mesh, horizon_id, is_fault, fault_opposite):
    horizon_list = extract_horizons(mesh, horizon_id)
    mesh = poisson.flatten_horizons(model, mesh, horizon_list)
    poisson.flatten_fault(model, mesh, is_fault, fault_opposite)


if __name__ == '__main__' :
    
    # flatten('chevron', Mesh("chevron/slice.obj"), horizon_id_chevron)
    flatten('ifp1', Mesh("ifp1/slice.obj"), horizon_id_ifp1, is_fault_ifp1, fault_opposite_ifp1)
    # flatten('ifp2', Mesh("ifp2/slice.obj"), horizon_id_ifp2)
    # flatten('shell', Mesh("shell/slice.obj"), horizon_id_shell)

    #print(extract_faults(Mesh("ifp2/slice.obj"), is_fault_ifp2))
    # fault_list = extract_faults(Mesh("ifp2/slice.obj"), is_fault_ifp2)
    # print([len(fault_list[i]) for i in range(len(fault_list))])
    # Probleme : 216 fois la meme liste...

    # faults_ifp2 = extract_faults(Mesh('ifp2/slice.obj'), is_fault_ifp2)
    # print(len(faults_ifp2))
    # poisson.flatten_fault('ifp1', Mesh('ifp1/slice_horizons.obj'), is_fault_ifp1)
    # print(count_nb_faults(Mesh('ifp2/slice.obj'), is_fault_ifp2))
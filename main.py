import numpy as np 
from mesh import Mesh
import matplotlib.pyplot as plt
import poisson

from chevron_attributes import horizon_id_chevron
from ifp1_attributes import horizon_id_ifp1
from ifp2_attributes import horizon_id_ifp2, is_fault_ifp2
from shell_attributes import horizon_id_shell


# Stores the indexes of the vertices of faults/horizons
def extract_horizons(mesh, horizon_id):
    """
    \param horizon_id list of horizons computed in the attributes file
    
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
    ncorn = mesh.ncorners
    nb_faults = count_nb_faults(mesh, is_fault)
    fault_list = [[] for _ in range(nb_faults)]
    for fault in range(nb_faults):
        for vertex in range(ncorn):
            if is_fault[vertex]:
                fault_list[fault].append(vertex)
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
    next_halfedge = mesh.opposite(mesh.prev(mesh.opposite(halfedge)))
    if is_fault[next_halfedge]:
        tag[next_halfedge] = nb_faults
        rec_count(mesh, next_halfedge, nb_faults, tag)
    else :
        tag[next_halfedge] = -1


def count_nb_faults(mesh, is_fault):
    nb_faults = 0
    ncorn = mesh.ncorners
    tag = np.zeros(ncorn)
    for halfedge in range(ncorn):
        if tag[halfedge] == 0:
            if is_fault[halfedge]:
                nb_faults += 1
                tag[halfedge] = nb_faults
                rec_count(mesh, halfedge, nb_faults, is_fault, tag)
            else :
                tag[halfedge] = -1
    return nb_faults


def flatten(model, mesh, horizon_id):
    horizon_list = extract_horizons(mesh, horizon_id)
    poisson.flatten_horizons(model, mesh, horizon_list)


if __name__ == '__main__' :
    
    # flatten('chevron', Mesh("chevron/slice.obj"), horizon_id_chevron)
    # flatten('ifp1', Mesh("ifp1/slice.obj"), horizon_id_ifp1)
    # flatten('ifp2', Mesh("ifp2/slice.obj"), horizon_id_ifp2)
    # flatten('shell', Mesh("shell/slice.obj"), horizon_id_shell)

    #print(extract_faults(Mesh("ifp2/slice.obj"), is_fault_ifp2))
    fault_list = extract_faults(Mesh("ifp2/slice.obj"), is_fault_ifp2)
    print([len(fault_list[i]) for i in range(len(fault_list))])
    # Probleme : 216 fois la meme liste...
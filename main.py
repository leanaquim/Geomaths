import numpy as np 
from mesh import Mesh
import matplotlib.pyplot as plt
import poisson

from chevron_attributes import horizon_id_chevron
from ifp1_attributes import horizon_id_ifp1
from ifp2_attributes import horizon_id_ifp2
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


def main(model, mesh, horizon_id):
    horizon_list = extract_horizons(mesh, horizon_id)
    poisson.flatten_horizons(model, mesh, horizon_list)


if __name__ == '__main__' :
    
    main('chevron', Mesh("chevron/slice.obj"), horizon_id_chevron)
    main('ifp1', Mesh("ifp1/slice.obj"), horizon_id_ifp1)
    main('ifp2', Mesh("ifp2/slice.obj"), horizon_id_ifp2)
    main('shell', Mesh("shell/slice.obj"), horizon_id_shell)

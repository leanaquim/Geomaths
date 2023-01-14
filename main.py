import numpy as np 
from mesh import Mesh
import matplotlib.pyplot as plt
import lsmr

import chevron.attributes as chevron
import ifp1.attributes as ifp1
import ifp2.attributes as ifp2
import shell.attributes as shell

# Stores the indexes of the vertices of horizons
def extract_horizons(mesh, horizon_id):
    nb_vertices = mesh.ncorners
    nb_horizons = max(horizon_id) + 1
    horizon_list = [[] for _ in range(nb_horizons)]
    for hor in range(nb_horizons):
        for vertex in range(nb_vertices):
            if horizon_id[vertex] == hor :
                horizon_list[hor].append(vertex)
    return horizon_list

# Main function
def flatten(model, mesh, horizon_id, is_fault, fault_opposite):
    horizon_list = extract_horizons(mesh, horizon_id)
    mesh = lsmr.flatten_horizons(model, mesh, horizon_list)
    lsmr.flatten_fault(model, mesh, is_fault, fault_opposite)


if __name__ == '__main__' :
    
    flatten('chevron', Mesh("chevron/slice.obj"), chevron.horizon_id, chevron.is_fault, chevron.fault_opposite)
    flatten('ifp1', Mesh("ifp1/slice.obj"), ifp1.horizon_id, ifp1.is_fault, ifp1.fault_opposite)
    flatten('ifp2', Mesh("ifp2/slice.obj"), ifp2.horizon_id, ifp2.is_fault, ifp2.fault_opposite)
    flatten('shell', Mesh("shell/slice.obj"), shell.horizon_id, shell.is_fault, shell.fault_opposite)

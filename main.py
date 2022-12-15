import numpy as np 
from mesh import Mesh
import matplotlib.pyplot as plt
import poisson

model = 'ifp1'
mesh = Mesh(model+"/slice.obj")

from ifp1_attributes import horizon_id, is_fault, fault_opposite

nb_vertices = mesh.ncorners

# Stores the indexes of the vertices of faults/horizons
def extract_horizons(horizon_id):
    """
    \param horizon_id list of horizons computed in the attributes file
    
    Extracts the list of horizons of the model.
    For each horizon found, it stores the index of all half-edges constituing the horizons (line number in the slice.obj file)
    """
    nb_horizons = max(horizon_id) + 1
    horizon_list = [[] for _ in range(nb_horizons)]
    for hor in range(nb_horizons):
        for vertex in range(nb_vertices):
            if horizon_id[vertex] == hor :
                horizon_list[hor].append(vertex)
    return horizon_list

def extract_faults(is_fault):
    #TODO
    return False

def main(mesh, horizon_id):
    horizon_list = extract_horizons(horizon_id)
    poisson.flatten_horizons(mesh, horizon_list)

if __name__ == '__main__' :
    h0 = extract_horizons(horizon_id)[0] # half-edge indexes
    list_v = [mesh.org(i) for i in h0]
    list_y = [mesh.V[i,1] for i in list_v]
    list_x = [mesh.V[i,0] for i in list_v]

    # print(extract_horizons(horizon_id)[1])

    # plt.scatter(list_x, list_y)
    # plt.ylim(0,1)
    # plt.show()
    
    main(mesh, horizon_id)
    
    main('ifp1', mesh, horizon_id)

    



        
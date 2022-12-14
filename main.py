import numpy as np 
from mesh import Mesh

model = 'chevron'
mesh = Mesh(model+"/slice.obj")

from chevron_attributes import horizon_id, is_fault, fault_opposite

nb_vertices = mesh.ncorners

# Stores the indexes of the vertices of faults/horizons
def extract_horizons(horizon_id):
    """
    \param horizon_id list of horizons computed in the attributes file
    
    Extracts the list of horizons of the model.
    For each horizon found, it stores the index of all half-edges constituing the horizons (line number in the slice.obj file)
    """
    nb_horizons = max(horizon_id)
    horizon_list = [[] for _ in range(nb_horizons)]
    for hor in range(nb_horizons):
        for vertex in range(nb_vertices):
            if horizon_id[vertex] == hor :
                horizon_list[hor].append(vertex)
    return horizon_list

def extract_faults(is_fault):
    #TODO
    return False

if __name__ == '__main__' :
    print(extract_horizons(horizon_id)[0])



        
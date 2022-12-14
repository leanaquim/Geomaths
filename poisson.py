import numpy as np
import matplotlib.pyplot as plt
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

def flatten_horizon(mesh, horizon):
    modified_mesh = mesh
    n = mesh.ncorners
    
    # Defining the 'transition' matrix
    A = np.matrix((n,n))
        # Constraints on horizon vertices
    for i in horizon :
        A[i,0] = 1 
        # Constraints on the rest of the mesh (Laplace ?)
    # TODO
    
    # Computing the new mesh


    modified_mesh.print_to_file('chevron/modified_mesh.obj')


if __name__ == '__main__' :
    mesh = Mesh("./chevron/slice.obj")
    horizon = [258, 267, 273, 285, 294, 307, 317, 816, 825, 837, 849, 861, 870, 882, 3975, 4044, 4056,
    4065, 4071, 4074, 4173, 4206, 4272, 4284, 4308, 4317, 4368, 4377, 4452, 4458, 4479, 4485, 4494, 4503,
    4512, 4521, 4527, 4548, 4569, 4584, 4587, 4606, 4620, 4647, 4656, 4662, 4677, 4692, 4716, 4734, 4740,
    4752, 4791, 4797, 4833, 4851, 4854, 4866, 4879, 4890, 4920, 4932, 4980, 5124, 5130, 5136, 5148, 5154,
    5181, 5190, 5523, 5541, 5560, 5583, 6117, 6126, 6144, 6147, 6237, 6279, 6297, 6321, 6369, 6474, 6582,
    6600, 6606, 6645, 6672, 6687, 6729, 6768, 6780, 6786, 6837, 6849, 6885, 6909, 6927, 6951, 6957, 6969,
    7017, 7035, 7059, 7074, 7095, 7101, 7110, 7122, 7140, 7152, 7182, 7200, 7218, 7230, 7242, 7266, 7281,
    7293, 7308, 7329, 7359, 7377, 7404, 7413, 7428, 7437, 7563, 7575, 9084, 9621, 10002, 10692]

    flatten_horizon(mesh, horizon)
    
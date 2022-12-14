import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import lsmr
from mesh import Mesh

slice = Mesh("./shell/slice.obj")
horiz0 = Mesh("./shell/horizon0.obj")

# Target horizon coordinate
z = slice.V[0][2]
print('z', z)

# List of x coordinates of the horizon 
f = []
for point in horiz0.V:
    if point[2] == z:
        f.append(point[1])
print('f', f)

# Number of points
n = len(f)
print('n', n)

# Target function (constant value)
g = np.array([f[0]] * n)
print('g', g)

# b = np.array([[g[i]-g[i-1]] for i in range(1, n)])
# b[0] += g[0]
# b[-1] += g[0]
# # b = np.array([[1] for _ in range(1, n)])
# print('b', b)

phi = np.array([[1], [1]])

A = np.matrix(np.ones((n, 2)))
for i in range(n):
    A[i, 0] = f[i]
print("A", A)

# result = f

# A = np.matrix(np.zeros((n-1,n-2)))
# A[0, 0] = result[1]
# A[-1, -1] = - result[n-2]
# for i in range(1, n-2):
#     A[i, i] = result[i + 1]
#     A[i, i-1] = - result[i]
# print("A", A)

result = np.linalg.inv(A.T * A) * (A.T) * g

# result = [g[0]] + (np.linalg.inv(A.T*A)*(A.T)*phi).T.tolist()[0] + [g[-1]]
print(len(result))

# Neanderthal version: 
# for _ in range(10) :
#     for i in range(1 , n-1) :
#         result[i] = (f[i-1] + f[i+1] + (2*g[i]-g[i-1]-g[i+1]) )/2.

# least squares from scipy:
# result = np.concatenate(([g[0]], lsmr(A, b)[0], [g[-1]]))

t = np.linspace(0, 1, n)
plt.plot(t, f, linewidth=3, label="f_init")
plt.plot(t, g, linewidth=3, label="g")
plt.plot(t, result, label="result")
plt.legend()
plt.show()

slice_2 = slice
horiz0_2 = horiz0

i = 0
for point in horiz0_2.V:
    if point[2] == z and i < 56:
        point[1] = result[i]
        slice_2.V[1] = result[i]
        i += 1

horiz0_2.print_to_file("horizon0_2.obj")
slice_2.print_to_file("slice_2.obj")

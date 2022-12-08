import numpy as np
import matplotlib.pyplot as plt
from mesh import Mesh

slice = Mesh("./shell/slice.obj")
horiz0 = Mesh("./shell/horizon0.obj")

z = slice.V[0][2]
print('z', z)

f = []
for point in horiz0.V:
    if point[2] == z:
        f.append(point[1])
print('f', f)

n = len(f)
print('n', n)

g = np.array([f[0]] * n)
print('g', g)

b = np.array([[g[i]-g[i-1]] for i in range(1, n)])
b[0] = g[0]
b[-1] = g[0]
print('b', b)

result = f

A = np.matrix(np.zeros((n-1,n-2)))
A[0, 0] = result[1]
A[-1, -1] = - result[n-2]
for i in range(1, n-2):
    A[i, i] = result[i + 1]
    A[i, i-1] = - result[i]
# np.fill_diagonal(A,      1)
# np.fill_diagonal(A[1:], -1)

result = [g[0]] + (np.linalg.inv(A.T*A)*(A.T)*b).T.tolist()[0] + [g[-1]]

# Neanderthal version: 
# for _ in range(10) :
#     for i in range(1 , n-1) :
#         result[i] = (f[i-1] + f[i+1] + (2*g[i]-g[i-1]-g[i+1]) )/2.

t = np.linspace(0, 1, n)
plt.plot(t, f, linewidth=5, label="f_init")
plt.plot(t, g, label="g")
plt.plot(t, result, label="result")
plt.legend()
plt.show()

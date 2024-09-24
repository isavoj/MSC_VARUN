import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve

# Define parameters
N = 100              # Number of Grid poitns
u0_t, u0_b = 0, 0    # Dirichlet boundary values on the top and bottom edges
u0_l, u0_r = 0, 0    # Dirichlet boundary values on the left and right edges
dx = 1.0 / (N + 1)   # Grid spacing

# 1d Diff matix
A_1d = (sp.eye(N, k=-1) + sp.eye(N, k=1) - 2 * sp.eye(N)) / (dx**2)

#2D diff matrix using Kronecker products
A = sp.kron(sp.eye(N), A_1d) + sp.kron(A_1d, sp.eye(N))


def source_function(x, y):
    # Source term of poission (Charge distribution)
    return np.sin(np.pi * x) * np.sin(np.pi * y)

# Generate grid points
x = np.linspace(0, 1, N + 2)  # N+2 includes the boundary points
y = np.linspace(0, 1, N + 2)
X, Y = np.meshgrid(x, y)

# Create the source term on the interior grid points - no boundaries that is
f = source_function(X[1:-1, 1:-1], Y[1:-1, 1:-1])


b = f.flatten()
b = b.reshape((N, N))

# ADD boundary conditions
b[0, :] -= u0_b / dx**2
b[-1, :] -= u0_t / dx**2
b[:, 0] -= u0_l / dx**2
b[:, -1] -= u0_r / dx**2

# Reshape b to a vector for the matrix system A*u = b
b = b.flatten()

# Solve the linear system A*u = b
u = spsolve(A, b)

# Reshape the solution vector back to a 2D grid of interior points
u = u.reshape((N, N))

# Combine the interior solution with Dirichlet boundary values
U = np.vstack([
    np.ones((1, N + 2)) * u0_b,
    np.hstack([np.ones((N, 1)) * u0_l, u, np.ones((N, 1)) * u0_r]),
    np.ones((1, N + 2)) * u0_t
])

fig, ax = plt.subplots(figsize=(8, 6))

X_plot, Y_plot = np.meshgrid(np.linspace(0, 1, N + 2), np.linspace(0, 1, N + 2))

c = ax.pcolor(X_plot, Y_plot, U, cmap='RdBu_r', shading='auto')
cb = plt.colorbar(c, ax=ax)

def exact_solution(x, y):
    return -(1 / (2 * np.pi ** 2)) * np.sin(np.pi * x) * np.sin(np.pi * y)


exact = np.array([[exact_solution(xi, yi) for yi in y] for xi in x])

plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.contourf(X_plot, Y_plot, U, cmap='viridis')
plt.colorbar(label='Numerical Solution')
plt.title('Finite Differences')
plt.xlabel('x')
plt.ylabel('y')

plt.subplot(1, 2, 2)
plt.contourf(x, y, exact, cmap='viridis')
plt.colorbar(label='Exact Solution')
plt.title('Exact Solution')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.show()


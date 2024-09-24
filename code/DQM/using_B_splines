import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline
from scipy.linalg import solve

# Define grid parameters
n = 100  # Grid Points
x = np.linspace(0, 1, n)
y = np.linspace(0, 1, n)
h = x[1] - x[0]  # Step size

X, Y = np.meshgrid(x, y)

f = np.sin(np.pi * X) * np.sin(np.pi * Y)
f_flat = f.flatten()  # Flatten for the solver

degree = 3
##### x
bspl_x = make_interp_spline(x, np.ones_like(x), k=degree) 
t = bspl_x.t  # Knot vector
c = np.eye(n)  # Coefficients 
splines_x = [BSpline(t, c[i], degree) for i in range(n)]  # B-spline basis for x

##### y
bspl_y = make_interp_spline(y, np.ones_like(y), k=degree)  
t_y = bspl_y.t  # Knot vector
splines_y = [BSpline(t_y, c[i], degree) for i in range(n)]  

D2_x = np.zeros((n, n))
D2_y = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        D2_x[i, j] = splines_x[j](x[i], nu=2)  # Second derivative in x
        D2_y[i, j] = splines_y[j](y[i], nu=2)  # second derivative in y

# Combine into a full 2D Laplacian matrix with Kronecker product
Laplacian = np.kron(np.eye(n), D2_x) + np.kron(D2_y, np.eye(n))

# Set Dirichlet conditions: u = 0 on the boundary
boundary_indices = np.unique(
    np.concatenate([
        np.arange(0, n),
        np.arange(n * (n - 1), n**2),
        np.arange(0, n**2, n),
        np.arange(n - 1, n**2, n)
    ])
)

# Modify the laplacian fo r the boundaries
for index in boundary_indices:
    Laplacian[index, :] = 0
    Laplacian[index, index] = 1
    f_flat[index] = 0

# Solve the linear system: Laplacian * U = f
u_approx = solve(Laplacian, f_flat)

#For  visalization purposes
u_approx = u_approx.reshape((n, n))

u_exact = -np.sin(np.pi * X) * np.sin(np.pi * Y) / (2 * np.pi ** 2)  # Example of exact solution

# Plotting the results
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.contourf(X, Y, u_approx, levels=50, cmap='viridis')
plt.colorbar()
plt.title('B-Spline DQM')

plt.subplot(1, 2, 2)
plt.contourf(X, Y, u_exact, levels=50, cmap='viridis')
plt.colorbar()
plt.title('Exact Solution')

plt.show()

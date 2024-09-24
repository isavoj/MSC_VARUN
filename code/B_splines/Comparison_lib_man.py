import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline, make_interp_spline

# This script compares the B-spline basis function and its first and second derivatives
# between a manually implemented method and the SciPy library to validate the manual implementation.

def b_spline_basis(x, i, k, t):
    """
    Recursive definition of the B-spline basis function B_{i,k}(x).

    Parameters:
        x (float): The point at which to evaluate the basis function.
        i (int): The index of the basis function.
        k (int): The degree of the basis function.
        t (list): The knot vector.

    Returns:
        float: The value of the B-spline basis function at x.
    """
    if k == 0:
        return 1.0 if t[i] <= x < t[i + 1] else 0.0

    left = 0 if (t[i + k] == t[i]) else (x - t[i]) / (t[i + k] - t[i]) * b_spline_basis(x, i, k - 1, t)
    right = 0 if (t[i + k + 1] == t[i + 1]) else (t[i + k + 1] - x) / (t[i + k + 1] - t[i + 1]) * b_spline_basis(x, i + 1, k - 1, t)
    return left + right

def first_derivative_b_spline(x, i, k, t):
    """
    Compute the first derivative of the B-spline basis function B'_{i,k}(x).

    Parameters:
        x (float): The point at which to evaluate the derivative.
        i (int): The index of the basis function.
        k (int): The degree of the basis function.
        t (list): The knot vector.

    Returns:
        float: The value of the first derivative of the B-spline basis function at x.
    """
    if k == 0:
        return 0.0

    left = 0 if (t[i + k] == t[i]) else k / (t[i + k] - t[i]) * b_spline_basis(x, i, k - 1, t)
    right = 0 if (t[i + k + 1] == t[i + 1]) else k / (t[i + k + 1] - t[i + 1]) * b_spline_basis(x, i + 1, k - 1, t)
    return left - right

def second_derivative_b_spline(x, i, k, t):
    """
    Compute the second derivative of the B-spline basis function B''_{i,k}(x).

    Parameters:
        x (float): The point at which to evaluate the second derivative.
        i (int): The index of the basis function.
        k (int): The degree of the basis function.
        t (list): The knot vector.

    Returns:
        float: The value of the second derivative of the B-spline basis function at x.
    """
    if k == 0:
        return 0.0

    left = 0 if (t[i + k] == t[i]) else k / (t[i + k] - t[i]) * first_derivative_b_spline(x, i, k - 1, t)
    right = 0 if (t[i + k + 1] == t[i + 1]) else k / (t[i + k + 1] - t[i + 1]) * first_derivative_b_spline(x, i + 1, k - 1, t)
    return left - right

###################################Parameters#####################################################

degree = 3  # Cubic spline
n = 20
x = np.linspace(0, 1, n)
x_vals = np.linspace(0.1, 0.4, 100)  # Points to evaluate
i = 5  # Basis function index within the active interval


###################################SciPy Computation#####################################################

bspl_x = make_interp_spline(x, np.ones_like(x), k=degree)  # Dummy coefficients just to get knots
t = bspl_x.t  # Knot vector from SciPy
spline = BSpline(t, np.eye(len(t) - degree - 1)[i], degree)  # B-spline basis for the given index i

# Compute SciPy B-spline basis function, first derivative, and second derivative
b_spline_vals_scipy = [spline(x_val) for x_val in x_vals]
first_derivative_vals_scipy = [spline(x_val, nu=1) for x_val in x_vals]
second_derivative_vals_scipy = [spline(x_val, nu=2) for x_val in x_vals]

######################################## MANUAL Computation################################################

b_spline_vals_manual = [b_spline_basis(x_val, i, degree, t) for x_val in x_vals] # B-spline basis for the given index i

# Compute Manual B-spline basis function, first derivative, and second derivative
first_derivative_vals_manual = [first_derivative_b_spline(x_val, i, degree, t) for x_val in x_vals]
second_derivative_vals_manual = [second_derivative_b_spline(x_val, i, degree, t) for x_val in x_vals]



plt.figure(figsize=(14, 8))

plt.subplot(3, 1, 1)
plt.plot(x_vals, b_spline_vals_manual, label='Manual Basis Function', linestyle='--', linewidth=2.5)
plt.plot(x_vals, b_spline_vals_scipy, label='SciPy Basis Function', linestyle='-', linewidth=1.5)
plt.title('B-spline Basis Function Comparison')
plt.xlabel('x')
plt.ylabel('Value')
plt.legend()
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(x_vals, first_derivative_vals_manual, label='Manual First Derivative', linestyle='--', linewidth=2.5)
plt.plot(x_vals, first_derivative_vals_scipy, label='SciPy First Derivative', linestyle='-', linewidth=1.5)
plt.title('First Derivative Comparison')
plt.xlabel('x')
plt.ylabel('Value')
plt.legend()
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(x_vals, second_derivative_vals_manual, label='Manual Second Derivative', linestyle='--', linewidth=2.5)
plt.plot(x_vals, second_derivative_vals_scipy, label='SciPy Second Derivative', linestyle='-', linewidth=1.5)
plt.title('Second Derivative Comparison')
plt.xlabel('x')
plt.ylabel('Value')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()

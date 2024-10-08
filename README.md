# MSC_VARUN
---

:warning: **Warning:** This project is currently under development. Additional information and updates will be added continuously. Stay tuned!

---
### Literature and Code for Varun's Master Thesis Project

Welcome to the repository for Varun's Master Thesis. This repository contains important literature and code related to Varun's research on solving Partial Differential Equations (PDEs) with a focus on B-splines. The purpose of this work is to compare his findings with the results of Isabella's network.

## Contents

- **Literature**: The `literature` folder contains various resources, papers, and notes that focus on Differential Quadrature and related numerical methods used throughout the project. This collection was essential in developing the code and understanding the methodologies.

- **Code**: The code in this repository is developed based on the literature reviewed. However, these are  just very simple examples i added.
The DQ -code was inspired by the paper [*Solving 2D-Poisson equation using modified cubic B-spline differential quadrature method*](https://www.sciencedirect.com/science/article/pii/S2090447917301521). Additionally, I added finite difference methods.
Both these cases solves the potential for an example source term `sin(pi * x) * sin(pi * y)` and compares the numerical solution with the exact solution. The results of these comparisons can be seen in the figures below.

> ### Some info about the DQ code 
> 
> We consider the well-known two-dimensional Poisson equation:
> 
> ```math
> \frac{\partial^2 z}{\partial x^2} + \frac{\partial^2 z}{\partial y^2} = f(x, y) = sin(pi * x) * sin(pi * y), \quad x, y \in [a, b]
> ```
> 
> under Dirichlet and mixed boundary conditions. The Differential Quadrature (DQ) method approximates the spatial partial derivatives in the PDE by a sum of weights multiplied by the function values at discrete knots over the domain \([a, b]\). The first and second-> > > order partial derivatives with respect to \(x\) over domain \([a, b]\) are defined as:
> 
> ```math
> U_x(x_i) = \sum_{j=1}^{N} W^{(1)}_{i,j} U(x_j), \quad \text{for } i = 1, \ldots, N
> ```
> 
> 
> ```math
> U_{xx}(x_i) = \sum_{j=1}^{N} W^{(2)}_{i,j} U(x_j), \quad \text{for } i = 1, \ldots, N
> ```
> 
> where $`W^{(1)}_{i,j}`$ and $`W^{(2)}_{i,j}`$ are the weighting coefficients of the first and second-order partial derivatives, respectively. 
> In this project, we use the modified cubic-B-spline method along with Shu's method for obtaining the weighting coefficients.
> 
> 
> Especially we have the B-spline basis functions $`\psi_k(x_i)`$  and its derivatives in the following relaion:
> 
> ```math
> \psi'_k(x_i) = \sum_{j=1}^{N} W^{(1)}_{i,j} \psi_k(x_j) \quad \text{for } i = 1, 2, \ldots, N; \, k = 1, 2, \ldots, N
> ```
> 
> To get the weighting coefficients of the second-order derivatives:
> 
> ```math
> W^{(2)}_{i,j} = 2W^{(1)}_{i,j} \cdot W^{(1)}_{i,i} - \frac{1}{x_i - x_j}, \quad \text{for } i, j = 1, \ldots, N; i \neq j
> ```
> 
> and
> 
> ```math
> W^{(2)}_{i,i} = - \sum_{j=1, j \neq i}^{N} W^{(2)}_{i,j}
> ```
> 
> 
> The differential quadrature method is applied which leads to a set of DQ algebraic equations as follows:
> 
> ```math
> \sum_{k=1}^{N} W^{(2)}_{i,k} Z_{k,j} + \sum_{k=1}^{M} W^{(2)}_{j,k} Z_{i,k} = f(x, y)
> ```
> 
> which can be rewritten as:
> 
> ```math
> \sum_{k=2}^{N-1} W^{(2)}_{i,k} Z_{k,j} + \sum_{k=2}^{M-1} W^{(2)}_{j,k} Z_{i,k} = f(x_i, y_j) - W_{i,1} Z_{1,j} - W_{i,N} Z_{N,j} - W_{j,1} Z_{i,1} - W_{j,M} Z_{i,M}
> ```
> 
> And then solved for. 
> 
> **_NOTE:_** 
> - In the provided DQ code, the second derivatives are directly obtained from the B-splines. Although the method differs slightly from the above formulation, the underlying principles remain consistent, and the equations are used to solve for the potential.
> 
> - The code that demonstrates these methods will be added soon.



## Recommended Reading about PDEs

For a comprehensive summary of solving PDEs in the electrostatic case, I highly recommend the following resource:

- **[Numerical Integration of Partial Differential Equations: Stationary Problems, Elliptic PDE](https://www.mps.mpg.de/phd/numerical-integration-partial-differential-equations-stationary-problems-elliptic-pde)**  
  This material provides an excellent overview of different methods used to solve PDEs, especially in the context of electrostatics. It covers various numerical techniques and is a great starting point for understanding the theoretical background.
Additionally, this GitHub repository contains course materials for **Fast Methods for Partial Differential and Integral Equations**, which includes a wealth of information valuable to this project. Be sure to explore the "readings" folder for further insights.

- **[MIT 18.336 - Fast Methods for Partial Differential and Integral Equations](https://github.com/mitmath/18336/tree/master)**  
  This resource originally hosted course material that covered fast computational methods for solving PDEs and integral equations, including valuable insights on spectral methods and other efficient techniques.

## Recommended Reading on B-Splines

For a review of B-splines and their properties, I recommend the materials provided by the University of Oslo. The following chapters were particularly useful:

- **[Chapter 1: Introduction to B-Splines](https://www.uio.no/studier/emner/matnat/ifi/nedlagte-emner/INF-MAT5340/v07/undervisningsmateriale/kap1.pdf)**: This chapter provides a foundational understanding of B-splines, including basic definitions and concepts.
- **[Chapter 3: Properties of B-Splines](https://www.uio.no/studier/emner/matnat/ifi/nedlagte-emner/INF-MAT5340/v07/undervisningsmateriale/kap3.pdf)**: This chapter covers the detailed properties of B-splines, which were extensively used in this project.

These resources offer comprehensive insights into B-splines, which were necessary in the development and understanding of the code in this repository.


## Additional Numerical Methods for PDE-Solving

For further exploration of numerical methods that can be discussed in the thesis, I introduce some documents here:

- **[Implementing FFTs in Practice](https://www.csd.uwo.ca/~mmorenom/CS433-CS9624/Resources/Implementing_FFTs_in_Practice.pdf)**: This document offers valuable insights into the practical aspects of FFTs.


### Results

Below are the results comparing the numerical solution with the exact solution:

**Figure 1: Using B-Spline DQM (code/DQM/using_B_splines.py)**

![B-splines](Figures/B_spline_DQ.png)

**Figure 2: Using Finite-Differences (code/finite_differences/simple_solution.py)**

![Finite-Differences](Figures/fin_dir.png)

**Figure 3: Comparison between Library SciPy and Manual implementation (code/B_splines/comparison_lib_man.py)**

![Comparison](Figures/comp_lib_man.png)


By solving, for example the following problem using DQ-Spline method
![solution](Figures/SOL_41.png)

We can compare it with the solution that KAN-neural network learns. 
Knowing the solution of the function, we can also fix the activation functions.
![easier](Figures/easier.png)








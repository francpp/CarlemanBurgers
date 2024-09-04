# PACS Project - Inhomogeneous, Viscous Burgers Equation Solver

### Author: Francesco Pettenon

## Project Overview

This project focuses on solving the inhomogeneous, viscous Burgers equation, potentially including linear damping. The solution is achieved using a direct application of the Carleman linearization method combined with Euler's method. The results from this solver are compared with solutions generated by MATLAB's inbuilt solvers, offering a benchmark for accuracy and efficiency.

The current implementation is based on the methodology described in the paper ["Efficient quantum algorithm for dissipative nonlinear differential equations"](https://arxiv.org/abs/2011.03185) by Jin-Peng Liu, Herman Øie Kolden, Hari K. Krovi, Nuno F. Loureiro, Konstantina Trivisa, and Andrew M. Childs.

## Problem Statement

We focus on solving an inhomogeneous, viscous Burgers equation:

```math
\partial_t u + u \partial_x u = \nu \partial_x^2 u + f
```

After space discretization, we describe the problem by the following \(n\)-dimensional quadratic ordinary differential equation (ODE):

```math
\frac{du}{dt} = F_2 u^{\otimes 2}+ F_1 u + F_0(t), \quad u(0) = u_{\text{in}}.
```

Here

```math
u = \begin{bmatrix} u_1, \dots, u_n \end{bmatrix}^T \in \mathbb{R}^n, \quad u^{\otimes 2} = \begin{bmatrix} u_1^2, u_1 u_2, \dots, u_1 u_n, u_2 u_1, \dots, u_n u_{n-1}, u_n^2 \end{bmatrix}^T \in \mathbb{R}^{n^2},
```

where each \(u_j = u_j(t)\) is a function of \(t\) on the interval \([0, T]\).

This problem is then linearized and solved using the **Carleman method**. We compare the solutions with those obtained with classical numerical methods, such as:

- Euler Method
- RK45
- Direct PDE solvers

## Usage Instructions

To set up and run the solver, follow these steps:

1. **Clone Repositories**

First, clone the necessary repositories:

```bash
git clone https://github.com/francpp/CarlemanBurgers.git
git clone --recursive git@github.com:pacs-course/pacs-examples.git
```

Ensure that `pacs-examples` is cloned **outside** of your project directory (i.e., it should not reside inside the `CarlemanBurgers` folder).

2. **Modify Makefile**

After cloning, modify the `PACS_ROOT` variable in the `Makefile` inside the project folder (`CarlemanBurgers`) to point to the correct directory where `pacs-examples` is located. This will ensure that all the necessary paths are correctly configured.

3. **Load Necessary Modules**

Load the required modules before proceeding:

```bash
module load gcc-glibc/11
module load eigen
module load boost
```

4. **Install Required Packages**

Install additional required tools:

```bash
sudo apt-get install gnuplot doxygen
```

These will be used for plotting and generating documentation.

5. **Installation**

Navigate to the project directory (`CarlemanBurgers`) and run the following command to install dependencies:

```bash
make install
```

6. **Compilation**

Compile the code by running:

```bash
make
```

This will generate the necessary executable files for running simulations.

7. **Set the Library Path**

Set the library path to include the `lib` directory:

```bash
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(pwd)/lib
```

8. **Run the Simulation**

Execute the main simulation:

```bash
./main_executable
```

9. **Generate Documentation**

To generate documentation, run the following command:

```bash
make docs
```

10. **Clean the Build**

To clean up the build files, use:

```bash
make clean
```

---

By following these steps, you will successfully set up, compile, and run the solver, along with generating any necessary documentation.


## Folder Structure

The project directory is organized as follows:

- **build/**: Stores the object (`.o`) files generated during the compilation process. These files are intermediate steps in building the final executable.

- **data/**: Contains input data files and parameter settings for the simulations. The primary files to be edited are `parameters.pot` and `parameters.json`. These files control various aspects of the simulation, including initial conditions, solver settings, and physical parameters.

- **include/**: Header files used throughout the project. This directory contains declarations for functions, classes, and templates that are implemented in the source files.

- **lib/**: This directory contains the compiled libraries and object files required by the project. These are generated during the build process and are used to link different components of the codebase efficiently.

- **src/**: Contains the source code of the project, further divided into subdirectories:
  - **discretization/**: Code related to the discretization of the PDEs and ODEs.
  - **error analysis/**: Modules for analyzing and calculating the error in the numerical solutions.
  - **initial conditions/**: Functions and classes for setting up the initial conditions of the simulations.
  - **params/**: Handling and processing of parameter files.
  - **plots/**: Generate plots with solution and errors.
  - **solvers/**: Implementation of the solvers, including Carleman linearization and Euler's method.
  - **utils/**: Utility functions and helpers used across the project.
  - **mainSimulation.hpp** and **mainSimulation.cpp**: The main simulation implementation files that drive the core functionality.

- **main.cpp**: The entry point of the program where the main function resides. It orchestrates the flow of the program by initializing parameters, invoking solvers, and handling output.

- **Makefile**: Automates the build process, including installation, compilation, and cleaning of the project.

- **output/**: Stores the output plots and data generated during the simulation. These files are created automatically as the simulation runs and are used to visualize results and analyze the behavior of the system.

- **README.md**: This file, providing an overview, usage instructions, and a description of the project's structure.

This structure is designed to ensure that the project is modular, easy to navigate, and scalable for future developments.

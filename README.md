# Usador

[![C++](https://img.shields.io/badge/C++-%2300599C.svg?logo=c%2B%2B&logoColor=white)](https://cplusplus.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21904035.svg)](https://doi.org/10.5281/zenodo.21904035)

* [PRESENTATION](#presentation)
    - [Matrix diffusion equation](#matrix-diffusion-equation)
    - [Boundary conditions](#boundary-conditions)
    - [Derived observables](#derived-observables)
    - [Solution algorithm](#solution-algorithm)
    - [Acknowledgements](#acknowledgements)
* [INSTALLATION AND USAGE](#installation-and-usage)
    - [Dependencies](#dependencies)
    - [Usage](#usage)
* [REFERENCES](#references)

## PRESENTATION

Usador is a C++ 2017 program to solve the matrix diffusion equation known as the *Usadel equation*[^1] describing the coherent propagation of a wave in a two-dimensional disordered medium.
The solution of this equation provides the distribution of singular values of the transmission matrix associated with the propagation of a wave between two edges of the medium[^2].
It also provides the disorder-averaged intensity profile of transmission eigenstates (also known as transmission eigenchannels[^3][^4]).
The name is an acronym for *"Usadel equation Solver for Arbitrary DisOrdered Regions"*.

### Matrix diffusion equation

The matrix diffusion equation reads

<p>$$ \mathbf{\nabla}_{\mathbf{r}}\cdot\mathbf{\mathsf{J}}(\mathbf{r}) = -\frac{[\mathsf{\sigma}_3, \mathsf{Q}]}{2\ell_{\rm a}} ,\qquad  \mathbf{\mathsf{J}}(\mathbf{r}) = -\frac{\ell_{\rm t}}{d} \mathsf{Q}\mathbf{\nabla}_{\mathbf{r}}\mathsf{Q} $$</p>

where 

- $d$ the number of spatial dimensions (limited to 2 in the program by design), 
- $\ell_{\rm t}$ is the transport [mean free path](https://en.wikipedia.org/wiki/Mean_free_path), 
- $\ell_{\rm a}$ is the ballistic [absorption length](https://en.wikipedia.org/wiki/Attenuation_length), which is infinite in the absence of absorption,
- $[\mathsf{A}, \mathsf{B}] = \mathsf{A}\mathsf{B} - \mathsf{B}\mathsf{A}$ denotes the [matrix commutator](https://en.wikipedia.org/wiki/Commutator),
- $\mathsf{\sigma}_3$ is the third [Pauli matrix](https://en.wikipedia.org/wiki/Pauli_matrices),
- $\mathsf{Q}(\mathbf{r})$ is a 2-by-2 complex matrix field obeying the constraints $\mathsf{Q}(\mathbf{r})^2=\mathsf{1}$ and $\mathrm{Tr}\mathsf{Q}(\mathbf{r})=0$ at each point of space.
- $\mathbf{\mathsf{J}}(\mathbf{r})$ is a $d$-vector of 2-by-2 complex matrices, which will be referred to as the matrix current.

The solution of this equation for $\mathsf{Q}(\mathbf{r})$ provides the distribution of transmission eigenvalues, the intensity profile of transmission eigenstates, and other observables related to transmission eigenstates (see below).
This equation closely resembles the standard [diffusion equation](https://en.wikipedia.org/wiki/Diffusion_equation) and this is not by chance: It is based on the same fundamental assumption that is the smallness of the mean free path compared to the system size ($\ell_{\rm t}\ll L$).
Apart from the matrix nature of this equation, this equation is distinguished by its *nonlinearity* in $\mathsf{Q}$ which makes its solution much richer than that of the standard diffusion equation.
Indeed, it can describe coherent effects (i.e., effects dependent on the phase of the wave) despite the large amount of scatterings.
A very similar equation is the *Usadel equation*[^1] which governs the dynamics of the electron Green's function in dirty (strongly scattering) superconductors.

### Boundary conditions

The matrix diffusion equation is accompanied by boundary conditions which depend on the transmission eigenvalue $T$ (between 0 and 1).
If the edge of the disordered medium is perfectly reflecting, then the normal component of the matrix flux is zero, or in other words,

<p>$$ \mathbf{n}\cdot\mathbf{\nabla}_{\mathbf{r}}\mathsf{Q}(\mathbf{r}) = 0 $$</p>

where $\mathbf{n}$ is the outward normal to the medium.
If the edge is connected to an input or an output duct which are controlled so as to achieve transmission $T$, then the boundary conditions read

<p>$$ \mathsf{Q}(\mathbf{r}_{\rm in}  + z_0\mathbf{n}) = \begin{pmatrix}1 & \tfrac{-2\mathrm{i}}{\sqrt{T}}\\ 0 & -1\end{pmatrix} $$</p>
<p>$$ \mathsf{Q}(\mathbf{r}_{\rm out} + z_0\mathbf{n}) = \begin{pmatrix}1 & 0\\ \tfrac{-2\mathrm{i}}{\sqrt{T}} & -1\end{pmatrix} $$</p>

where $\mathbf{r}_{\rm in}$ and $\mathbf{r}_{\rm out}$ denote positions on the input or output edges, respectively, and $z_0=\mu\ell_{\rm t}$ is the diffusive extrapolation length, with $\mu=\frac{\pi}{4}$ in two dimensions.
If, on the contrary, the duct is leaky or uncontrolled, leaving the wave escape freely (hence acting as an absorber), then the boundary condition is given by

<p>$$ \mathsf{Q}(\mathbf{r} + z_0\mathbf{n}) = \begin{pmatrix}1 & 0\\ 0 & -1\end{pmatrix} $$</p>

for all $\mathbf{r}$ on the leaky edge.

### Derived observables

The matrix diffusion equation predicts several observables related to transmission eigenstates.
The first one is the intensity profile of transmission eigenstates defined by

<p>$$ I_{T}(\mathbf{r}) = \frac{1}{N_{\rm e}\rho(T)} \left\langle \sum_{e=1}^{N_{\rm e}} |\psi_{T_e}(\mathbf{r})|^2 \delta(T-T_e) \right\rangle $$</p>

where

- $\langle\cdots\rangle$ denotes the average over the disorder, 
- $T_e$ is the $e$-th transmission eigenvalue[^2], i.e., the eigenvalue of $\mathsf{t}^\dagger \mathsf{t}$, where $\mathsf{t}$ is the transmission matrix which relates the wavefront in the output duct to the wavefront in the input duct,
- $N_{\rm e}=\min(N_{\rm in},N_{\rm out})$ is the number of transmission eigenchannels, $N_{\rm in}$ and $N_{\rm out}$ being the number of waveguide eigenmodes in the input duct and the output duct, respectively, and
- $\psi_{T_e}(\mathbf{r})$ is the transmission eigenstate (or transmission eigenchannel[^3][^4]) normalized to unit incident intensity on average.

The prediction of the matrix diffusion equation for this observable is

<p>$$ I_{T}(\mathbf{r}) = \frac{S_{\rm in}}{\pi S_{\rm e}\rho(T)\sqrt{T}} \Re Q_{12}(\mathbf{r}) $$</p>

where $Q_{12}(\mathbf{r})$ is the upper-right component of the matrix intensity $\mathsf{Q}(\mathbf{r})$.
Note that the matrix diffusion equation also predicts the current density associated with transmission eigenstates by replacing $I_{T}(\mathbf{r})$ by $\mathbf{J}_{T}(\mathbf{r})$ and $Q_{12}(\mathbf{r})$ by $\mathbf{J}_{12}(\mathbf{r})$ in the formula hereabove.
In this case, $\mathbf{J}_{12}(\mathbf{r})$ represents the upper-right component of the matrix current $\mathbf{\mathsf{J}}(\mathbf{r})$ in the matrix diffusion equation.

In addition, the distribution of transmission eigenvalues defined by[^2]

<p>$$ \rho(T) = \frac{1}{N_{\rm e}} \sum_{e=1}^{N_{\rm e}} \langle\delta(T-T_{\rm e})\rangle $$</p>

also turns out to be captured by the matrix diffusion equation according to[^5][^6]

<p>$$ \rho(T) = \frac{d\mu}{\pi S_{\rm e}T^{\frac{3}{2}}} \Re \int_{\mathcal{S}_{\rm out}} \mathrm{d}{\mathbf{y}}\cdot\mathbf{J}_{12}(\mathbf{r}) $$</p>

where the integral is over the output surface $\mathcal{S}_{\rm out}$ with the surface element $\mathrm{d}\mathbf{y}$ pointing in the outward direction with respect to the disordered medium.
The resulting distribution $\rho(T)$ is normalized by $\int_{0^+}^1 \mathrm{T}\,\rho(T)=1$, excluding possible transmission eigenvalues exactly equal to zero.

### Solution algorithm

The matrix diffusion equation given above can be solved analytically in the absence of absorption or leaky edges (that is when all the degrees of freedom of the wave can be controlled).
However, it cannot be solved in the general case of an absorbing medium or with losses at the edges, for which a numerical solver is required.
The main difficulty encountered in solving the matrix diffusion equation is its nonlinearity with respect to $\mathsf{Q}(\mathbf{r})$, which arises from the normalization condition $\mathsf{Q}(\mathbf{r})^2 = \mathsf{1}$.
This normalization reduces the number of unknowns from four complex functions in the 2-by-2 matrix $\mathsf{Q}(\mathbf{r})$ to only two.
In order to take into account this constraint, it is appropriate to consider the following parameterization

<p>$$ \mathsf{Q}(\mathbf{r}) = \sin\varphi(\mathbf{r})\cos\vartheta(\mathbf{r}) \mathsf{\sigma}_1 - \sin\vartheta(\mathbf{r}) \mathsf{\sigma}_2 + \cos\varphi(\mathbf{r}) \cos\vartheta(\mathbf{r}) \mathsf{\sigma}_3 $$</p>

where $\vartheta(\mathbf{r})$ and $\varphi(\mathbf{r})$ are two complex functions which can be geometrically interpreted as angles over the spherical manifold $\mathsf{Q}(\mathbf{r})^2 = \mathsf{1}$, and $\mathsf{\sigma}_1,\mathsf{\sigma}_2,\mathsf{\sigma}_3$ are the three [Pauli matrices](https://en.wikipedia.org/wiki/Pauli_matrices).
This equation is accompanied by boundary conditions at input and output edges, which read

<p>$$ \vartheta(\mathbf{r}_{\rm out} + z_0\mathbf{n}) = -\vartheta(\mathbf{r}_{\rm in} + z_0\mathbf{n}) = \frac{\pi}{2} + \mathrm{i}\mathrm{arccosh}(\frac{1}{\sqrt{T}}) $$</p>

<p>$$ \varphi(\mathbf{r}_{\rm out} + z_0\mathbf{n}) = \varphi(\mathbf{r}_{\rm in} + z_0\mathbf{n}) = \frac{\pi}{2} - \mathrm{i}\mathrm{arccosh}(\frac{1}{\sqrt{1-T}}) , $$</p>

and by boundary conditions at leaky edges, which read

<p>$$ \vartheta(\mathbf{r} + z_0\mathbf{n}) = \varphi(\mathbf{r} + z_0\mathbf{n}) = 0 $$</p>

where $\mathbf{r}$ lies on the leaky edge.
In this parametrization, the matrix diffusion equation becomes a system of two nonlinear equations for the two unknown functions $\vartheta(\mathbf{r})$ and $\varphi(\mathbf{r})$:

<p>$$ \nabla^2_{\mathbf{r}}\vartheta + \sin\vartheta\cos\vartheta(\mathbf{\nabla}_{\mathbf{r}}\varphi)^2 = \frac{d}{\ell_{\rm a}\ell_{\rm t}} \cos\varphi\sin\vartheta $$</p>

<p>$$ \cos\vartheta\nabla^2_{\mathbf{r}}\varphi - 2\sin\vartheta\mathbf{\nabla}_{\mathbf{r}}\vartheta\cdot\mathbf{\nabla}_{\mathbf{r}}\varphi = \frac{d}{\ell_{\rm a}\ell_{\rm t}} \sin\varphi $$</p>

In the program, these equations are discretized over a square lattice, and the [Newton-Raphson algorithm](https://en.wikipedia.org/wiki/Newton's_method) is then used to iteratively solve the system.
This algorithm is stabilized using linear backtracking[^7].
The Jacobian matrix is computed used first-order [finite difference](https://en.wikipedia.org/wiki/Finite_difference) only between pairs of neighboring points in order to preserve the sparsity of the Jacobian matrix.
An appropriate initial guess for this algorithm is

<p>$$ \vartheta^{(0)}(\mathbf{r}) = 0, \qquad
\varphi^{(0)}(\mathbf{r}) = \frac{\pi}{2} - \mathrm{i} \mathrm{arccosh}\left( \frac{1}{\sqrt{1-T}} \right) $$</p>

which has the particularity of falling midway between the input and output boundary conditions given above.

### Acknowledgements

This program has been written by David Gaspard ([Institut Langevin](https://ror.org/00kr24y60), [ESPCI Paris](https://ror.org/03zx86w41), [PSL University](https://ror.org/013cjyk83), [CNRS](https://ror.org/02feahw73)) mainly in July 2025.
This research has been supported by the [ANR](https://ror.org/00rbzpz17) project MARS_light under reference [ANR-19-CE30-0026](https://anr.fr/Project-ANR-19-CE30-0026), by the program "Investissements d'Avenir" launched by the French Government.
It also received support from a grant of the [Simons Foundation](https://ror.org/01cmst727) (No. 1027116).

## INSTALLATION AND USAGE 

The source files can be downloaded using the [`git clone`](https://git-scm.com/docs/git-clone) command.
To compile the program, call the [`make`](https://en.wikipedia.org/wiki/Make_(software)) utility in the root directory:
```shell
make all
```
This should generate an executable.

### Dependencies

The program requires a [C++](https://en.wikipedia.org/wiki/C++) compiler, such as from the [GNU Compiler Collection](https://en.wikipedia.org/wiki/GNU_Compiler_Collection), and the libraries [OpenBLAS](https://en.wikipedia.org/wiki/OpenBLAS), [UMFPACK](https://en.wikipedia.org/wiki/UMFPACK), [MUMPS](https://en.wikipedia.org/wiki/MUMPS_(software)), and [libpng](https://en.wikipedia.org/wiki/Libpng).
It also calls [Python 3](https://en.wikipedia.org/wiki/Python_(programming_language)) with the [NumPy](https://en.wikipedia.org/wiki/NumPy), [Matplotlib](https://en.wikipedia.org/wiki/Matplotlib), and [csv](https://docs.python.org/3/library/csv.html) modules, but also the command `pdflatex` from [TeX Live](https://en.wikipedia.org/wiki/TeX_Live) with the [PGFPlots](https://ctan.org/pkg/pgfplots) package to process graphical outputs automatically.

### Usage

The recommended way to launch the program is through the [`run`](run) script:
```shell
./run
```
This script defines the environment variables and calls the executable.
The executable is built from the [`src/Main.cpp`](src/Main.cpp) file which contains the computation instructions in [C++](https://en.wikipedia.org/wiki/C++).

#### Initialization

The first step is to define the mesh over which the wave equation will be solved.
This can be done by instantiating a `SquareMesh` object as follows
```cpp
SquareMesh mesh("model/image.png");
```
where `model/image.png` can be replaced by the file path of any image in the [`model/`](model) folder or elsewhere on the host machine.
This model image completely defines the mesh structure and uses the following color code:

- White (`#FFFFFF`): The pixel is void and should be skipped from the mesh.
- Black (`#000000`): The pixel is added to the mesh. The default boundary condition between black and white pixels are Dirichlet boundary conditions $\psi(\mathbf{r})=0$.
- Red (`#FF0000`): Black pixels surrounding a red pixel are flagged as belonging to an input duct for the construction of the transmission matrix.
  Red pixels are not added to the mesh.
- Blue (`#0000FF`): Black pixels surrounding a blue pixel are flagged as belonging to an output duct for the construction of the transmission matrix. 
  Blue pixels are not added to the mesh.
- Green (`#00FF00`): Black pixels surrounding a green pixels are flagged as belonging to a duct which is neither input nor output and thus ignored from the transmission matrix.
  Green pixels are not added to the mesh.

Note that red (`#FF0000`), blue (`#0000FF`), and green (`#00FF00`) pixels assume free-escape boundary conditions presented before.
It is worth noting that `SquareMesh` objects can also be generated procedurally using the methods found in the file [`src/SquareMesh.hpp`](src/SquareMesh.hpp).
This can be helpful if the boundaries are complicated or change randomly from one simulation to the other.

A `UsadelSystem` object (containing the solver of the Usadel equation and computation methods of physical observables) is then instantiated by the command:
```cpp
UsadelSystem usys(sysname, mesh, holscat, holabso, tval);
```
where

- `sysname` is a string containing the full name of the system, and which will be used for file output.
- `mesh` is the `SquareMesh` object initialized from a model image (see above).
- `holscat` is the ratio of the mesh step over the scattering [mean free path](https://en.wikipedia.org/wiki/Mean_free_path).
- `holabso` is the ratio of the mesh step over the ballistic [absorption length](https://en.wikipedia.org/wiki/Attenuation_length).
- `tval` is the initial value of the tranmsission eigenvalue between 0 and 1.

#### Computations

The `UsadelSystem` object is the main computational object of the program.
Computations can be performed by calling the methods of `UsadelSystem` (see also [`src/UsadelSystem.hpp`](src/UsadelSystem.hpp)).
A typical example of computation of the matrix field $\mathsf{Q}(\mathbf{r})$ for a given value of transmission eigenvalue could be
```cpp
usys.setTransmission(tval);  // Reset the transmission eigenvalue.
usys.initConstant();  // Initialize the Q field to a constant.
int niter = usys.solveNewton(maxit, nsub, toldf, tolr, verbose);  // Solves the Usadel equation with the Newton-Raphson method.
```
where

- `tval` is the desired transmission eigenvalue $T$ between 0 and 1.
- `maxit` is the maximum number of iterations of the Newton-Raphson solver, typically between 20 and 100.
- `nsub` is the maximum number of substep used for backtracking line search, typically between 20 and 50 for double precision.
- `toldf` is tolerance over the relative displacement prescribed by the Newton-Raphson step, typically $10^{-7}$.
- `tolr` is tolerance over the norm of the residual compared to the norm of the initial residual, typically $10^{-10}$.
- `verbose` is the verbosity level of the Newton-Raphson solver. 0 is no output, 1 displays each iteration to stdout.
- `niter` is the number of iterations returned by the Newton-Raphson solver to meet the convergence criterion.

The resulting data can then be saved to a CSV file and plotted with:
```cpp
const std::string filepath = "path/to/file";  // File path without extension (CSV by default).
usys.savePlot(filepath);
```
On output, the CSV file contains the position of each point of the lattice and the values of $\vartheta(\mathbf{r})$, $\varphi(\mathbf{r})$, each component of $\mathsf{Q}(\mathbf{r})$, and several observables such as $I_{T}(\mathbf{r})$, the intensity of the transmission eigenstate at the given value of $T$.
This command also calls the Python script [`plot/plot_map.py`](plot/plot_map.py) in order to plot the field.
More precisely, this script generates a [PGF/TikZ](https://en.wikipedia.org/wiki/PGF/TikZ) file and a PNG file, and compiles the final PDF file using LaTeX.
By default, this script plots the intensity of the transmission eigenstate at the given value of $T$.

After solving the equation with the `solveNewton()` method, the transmission eigenvalue density $\rho(T)$ at the given value of $T$ can also be computed by:
```cpp
double rho = usys.getRho();
```
This command is called recursively in the following function in order to produce the full transmission eigenvalue distribution $\rho(T)$ (see [`src/Main.cpp`](src/Main.cpp)):
```cpp
computeDistributionSerial(usys, tmin, tmax, ntval);
```
where

- `tmin` and `tmax` are the minimum and maximum number of transmission eigenvalues, typically $10^{-2}$ and $1$.
- `ntval` is the desired number of samples of transmission eigenvalues, typically $300$.

This function iterates through the transmission eigenvalues ​​and, at each step, uses the solution of the previous step to initialize the Newton-Raphson algorithm.
It stops automatically as soon as the algorithm fails to converge, saves the results in a CSV file, and calls the Python script [`plot/plot_distrib.py`](plot/plot_distrib.py) to plot the data using LaTeX and [PGF/TikZ](https://en.wikipedia.org/wiki/PGF/TikZ).

## REFERENCES

[^1]: K. D. Usadel, *Generalized Diffusion Equation for Superconducting Alloys*, [Phys. Rev. Lett. **25**, 507-509 (1970)](https://doi.org/10.1103/PhysRevLett.25.507).
[^2]: C. W. J. Beenakker, *Random-matrix theory of quantum transport*, [Rev. Mod. Phys. **69**, 731-808 (1997)](https://doi.org/10.1103/RevModPhys.69.731).
[^3]: W. Choi, A. P. Mosk, Q.-H. Park, and W. Choi, *Transmission eigenchannels in a disordered medium*, [Phys. Rev. B **83**, 134207 (2011)](https://doi.org/10.1103/PhysRevB.83.134207).
[^4]: M. Davy, Z. Shi, J. Park, C. Tian, and A. Z. Genack, *Universal structure of transmission eigenchannels inside opaque media*, [Nat. Commun. **6**, 6893 (2015)](https://doi.org/10.1038/ncomms7893).
[^5]: D. Gaspard and A. Goetschy, *Radiant Field Theory: A Transport Approach to Shaped Wave Transmission through Disordered Media*, [Phys. Rev. Lett. **135**, 033804 (2025)](https://doi.org/10.1103/g3kd-sg4x).
[^6]: D. Gaspard and A. Goetschy, *Transmission eigenvalue distribution in disordered media from radiant field theory*, [Phys. Rev. Res. **7**, 033071 (2025)](https://doi.org/10.1103/djhy-16mh).
[^7]: W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, [*Numerical Recipes: The Art of Scientific Computing*](https://www.cambridge.org/9780521880688), (Cambridge University Press, 2007), 3rd ed., sec. 9.7.1.

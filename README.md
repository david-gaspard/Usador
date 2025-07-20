# Usador

[![C++](https://img.shields.io/badge/C++-%2300599C.svg?logo=c%2B%2B&logoColor=white)](https://cplusplus.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)

* [PRESENTATION](#presentation)
    - [Usadel equation](#usadel-equation)
    - [Transmission eigenvalue distribution](#transmission-eigenvalue-distribution)
    - [Profile of transmission eigenchannels](#profile-of-transmission-eigenchannels)
    - [Acknowledgements](#acknowledgements)
* [USAGE AND OPTIONS](#usage-and-options)
* [INTERNAL IMPLEMENTATION](#internal-implementation)
* [REFERENCES](#references)

## PRESENTATION

Usador is a C++ 2017 program to solve the matrix diffusion equation known as the *Usadel equation* [[1](#1)] describing the coherent propagation of a wave in a two-dimensional disordered medium.
The name is an acronym for *"Usadel equation Solver for Arbitrary DisOrdered Regions"*.
The solution of this equation provides the distribution of eigenvalues of the transmission matrix associated with the propagation of a wave between two surfaces located inside or at the edge of the medium.
It also provides the disorder-averaged intensity profile of transmission eigenchannels (aka transmission eigenstates).

### Usadel equation

The matrix diffusion equation reads

<p>$$ \nabla\cdot(\frac{\ell_{\mathrm{s}}}{d} \mathsf{Q}\nabla\mathsf{Q}) = \frac{1}{2\ell_{\mathrm{a}}} [\mathsf{\Lambda}_3, \mathsf{Q}] $$</p>

where $d$ the number of spatial dimensions (limited to 2 in the program by design), $\ell_{\mathrm{s}}$ is the scattering [mean free path](https://en.wikipedia.org/wiki/Mean_free_path), $\ell_{\mathrm{a}}$ is the ballistic [absorption length](https://en.wikipedia.org/wiki/Attenuation_length).
In this equation, $[\mathsf{A}, \mathsf{B}] = \mathsf{A}\mathsf{B} - \mathsf{B}\mathsf{A}$ denotes the [matrix commutator](https://en.wikipedia.org/wiki/Commutator), and $\mathsf{\Lambda}_3$ is the third [Pauli matrix](https://en.wikipedia.org/wiki/Pauli_matrices).
The field $\mathsf{Q}(\mathbf{r})$, a 2-by-2 complex matrix, is the main unknown of the equation and its solution provides the distribution of transmission eigenvalues or the intensity profile of transmission eigenchannels.
This equation closely resembles the standard [diffusion equation](https://en.wikipedia.org/wiki/Diffusion_equation) and this is not by chance: It is based on the same fundamental assumption that is the smallness of the mean free path compared to the system size ($\ell_{\mathrm{s}}\ll L$).
Apart from the matrix nature of this equation, this equation is distinguished by its nonlinearity in $\mathsf{Q}$ which makes its solution much richer than that of the standard diffusion equation.
Indeed, it can describe coherent effects---effects dependent on the phase of the wave---despite the large amount of scatterings.

### Transmission eigenvalue distribution

Without heterogeneities, the solution of this equation is trivial: $\mathsf{Q}(\mathbf{r}) = \mathsf{\Lambda}_3$.
Things change when looking at specific observables because then the equation must be supplemented by nontrivial boundary conditions at the edge or in the bulk of the medium.
In particular, if we consider the transmission eigenvalue distribution from one surface to another, we need to impose the following boundary conditions on these surfaces:

* $\mathsf{Q}_{\mathrm{out}} = \mathrm{e}^{\mathrm{i}\gamma_{\mathrm{a}}\mathsf{\Lambda}_+} \mathsf{Q}_{\mathrm{in}} \mathrm{e}^{-\mathrm{i}\gamma_{\mathrm{a}}\mathsf{\Lambda}_+}$ for the input surface.
    
* $\mathsf{Q}_{\mathrm{out}} = \mathrm{e}^{\mathrm{i}\gamma_{\mathrm{b}}\mathsf{\Lambda}_-} \mathsf{Q}_{\mathrm{in}} \mathrm{e}^{-\mathrm{i}\gamma_{\mathrm{b}}\mathsf{\Lambda}_-}$ for the output surface.

Here, $\mathsf{Q}_{\mathrm{in}}$ and $\mathsf{Q}_{\mathrm{out}}$ denote the input and output fields which are located on either side of the surface, $\gamma_{\mathrm{a}}$ and $\gamma_{\mathrm{b}}$ are the contact parameters which are related to the transmission probability $T$ by $\gamma=\gamma_{\mathrm{a}}\gamma_{\mathrm{b}} = T^{-1} + \mathsf{i}0^+$, and $\mathsf{\Lambda}_\pm = (\mathsf{\Lambda}_1 \pm \mathrm{i}\mathsf{\Lambda}_2)/2$ are the ladder Pauli matrices.

The transmission eigenvalue distribution is then given by

<p>$$ \rho(T) = ... $$</p>

<span style="color:red">TODO: Calculation of $\rho(T)$ remains to be implemented...................................</span>

### Profile of transmission eigenchannels

<span style="color:red">TODO: Computation of $I_{\mathrm{a},T}(\mathbf{r})$ remains to be implemented...................................</span>

### Acknowledgements

This program was developed by David Gaspard ([Institut Langevin](https://ror.org/00kr24y60), [ESPCI Paris](https://ror.org/03zx86w41), [PSL University](https://ror.org/013cjyk83), [CNRS](https://ror.org/02feahw73)) mainly in July 2025.
This research has been supported by the [ANR](https://ror.org/00rbzpz17) project MARS_light under reference [ANR-19-CE30-0026](https://anr.fr/Project-ANR-19-CE30-0026), by the program "Investissements d'Avenir" launched by the French Government.
It also received support from a grant of the [Simons Foundation](https://ror.org/01cmst727) (No. 1027116).

## USAGE AND OPTIONS

<span style="color:red">TODO: Write a few words on the commands...................................</span>

## INTERNAL IMPLEMENTATION

The Usadel equation is solved iteratively using the [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton's_method) which reads in summary:

> Repeat $s=1,\ldots,s_{\mathrm{max}}$ until $\|\mathbf{p}_s\|$ and $\|\mathbf{r}_s\|$ are small enough:
>     $\mathsf{M}_s\cdot\mathbf{p}_s = -\mathbf{r}_s$
>     $\mathbf{q}_{s+1} = \mathbf{q}_s + \tau_s\mathbf{p}_s$

where $\mathbf{q}_s$ is a vector representation of the $\mathsf{Q}(\mathbf{r})$ field on each point of the lattice, $\mathbf{r}_s$ is the residual vector which must be cancelled out on the solution, $\mathsf{M}_s = \partial\mathbf{r}/\partial\mathbf{q}$ is the Jacobian obtained by differentiating the residual $\mathbf{r}_s$ with respect to the unknown vector $\mathbf{q}_s$.
The parameter $\tau_s\in[0,1]$ is a line search parameter ensuring the residual norm $\|\mathbf{r}_s\|$ decreases at each iteration, a requirement which is not guaranteed by the strict Newton-Raphson algorithm ($\tau_s=1$).

<span style="color:red">TODO: Write a few words of explanations on the algorithm, the representation of $\mathsf{Q}$ as a vector form, the angle parameters $(\theta,\eta)$, the expression of the Usadel residual equation (equation that must vanish). Also explain what does the Newton-Raphson method (the sparse Jacobian), and add a comment on the initial ansatz (constant or random).</span>

## REFERENCES

<a id="1">[1]</a>
Klaus D. Usadel,
*Generalized Diffusion Equation for Superconducting Alloys*,
[Phys. Rev. Lett. **25**, 507-509 (1970)](https://doi.org/10.1103/PhysRevLett.25.507)


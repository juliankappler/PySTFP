# PySTFP: Python module for the short-time Fokker-Planck propagator

## About

This module contains sympy-based symbolic tools for the short-time propagator for one-dimensional Fokker-Planck dynamics. This python module accompanies the paper Ref. <a href="#ref_1">[1]</a>.

## Introduction

For a given drift $a(x)$ and diffusivity $D(x)$ we consider the Fokker-Planck equation

$$\partial_t P = - \partial_x \left( a P\right) + \partial_x^2 \left(D P\right),$$

where $P \equiv P(x,t \mid x_0,t_0)$ is the transition density, or propagator, to find a particle that started at time $t_0$ at location $x_0$ at a later time $t$ at a location $x$.

We consider the Fokker-Planck equation on a spatial domain $x \in \mathbb{R} = (-\infty,\infty)$, with boundary conditions $P(x,t) \rightarrow 0$ as $|x| \rightarrow \infty$,
and the initial conditions $P(x,t_0) = \delta (x-x_0)$,  where $\delta$ denotes the Dirac-delta distribution.

In this module, we implement an approximate short-time propagator, in both the normalization-preserving representation

$$P_K(x,t \mid x_0 ,t_0) = \frac{1}{\sqrt{ 4 D(x_0) \pi \Delta t}} \exp\left[ - \frac{ \Delta x^2}{4 D(x_0) \Delta t} \right] \times \left[ 1 +\sum_{k=1}^K \sqrt{\Delta t}^k \mathcal{Q}_k(\Delta x,\Delta t,x_0)\right],$$

and the positivity-preserving representation

$$P_K(x,t \mid x_0 ,t_0) = \frac{1}{\sqrt{ 4 D(x_0) \pi \Delta t}} \exp\left[ - \frac{ \Delta x^2}{4 D(x_0) \Delta t} + \sum_{k=1}^K \sqrt{\Delta t}^k \hat{\mathcal{Q}}_k(\Delta x, \Delta t,x_0)\right].$$

In the formulas above, $\Delta x \equiv x - x_0$, $\Delta t = t - t_0$, and $K \in \mathbb{N}_0$ is the truncation order for the approximate propagator. The explicit form of the coefficients $\mathcal{Q}_k$, $\hat{\mathcal{Q}}_k$, are a main result of Ref. <a href="#ref_1">[1]</a>. We have precalculated them up to $K = 8$ (corresponding to an accuracy $\Delta t^4$ for the propagator), but include code to calculate also higher-order coefficients.

In PySTFP we furthermore implement the explicit perturbation expansions for the first four moments $\langle \Delta x^n \rangle$ ($n = 0, 1, 2, 3, 4$),
as well as the medium entropy production rate, the total entropy production rate,
and the Gibbs entropy. See Ref. <a href="#ref_1">[1]</a> and the examples below for more details.


## Examples

In the subfolder [examples/](examples/) we provide several notebooks that illustrate the use of PySTFP, 
as well as allow to re-derive and extend the symbolic perturbative expressions used.

Here is a list of the notebooks:

**Subfolder [examples/derivations](examples/derivations)**

These notebooks implement the symbolic derivations of all the perturbative results from Ref. <a href="#ref_1">[1]</a>:

* [perturbative propagators and moments.ipynb](examples/derivations/perturbative%20propagators%20and%20moments.ipynb): 
Symbolically evaluate the perturbation series for both the normalization-preserving and positivity-preserving propagator. Additionally, use the normalization-preserving propagator to calculate some moments $\langle \Delta x^k \rangle$.
* [entropy production rates.ipynb](examples/derivations/entropy%20production%20rates.ipynb): Symbolically calculate the medium entropy production rate, the total entropy production rate, and the Gibbs entropy.


**Subfolder [examples/](examples/)**

These notebooks showcase basic functionality of PySTFP:

* [plot propagator.ipynb](examples/plot%20propagator.ipynb): Plot the short-time propagator for a given drift and diffusivity profile.
* [propagator in physical units.ipynb](examples/propagator%20in%20physical%20units.ipynb): Access the symbolic propagator stored in PySTFP, which are expressed in dimensionless variables, and rewrite them in physical dimensions.
* [entropy production in physical units.ipynb](examples/entropy%20production%20in%20physical%20units.ipynb): Access the symbolic entropy production (rates) stored in PySTFP, which are expressed in dimensionless variables, and rewrite them in physical dimensions.


**Subfolder [examples/numerical example](examples/numerical%20example)**

These notebooks allow to reproduce the plots of Ref. <a href="#ref_1">[1]</a>:

* [analytical solution.ipynb](examples/numerical%20example/analytical%20solution.ipynb): Evaluate the drift, diffusivity, and propagator exactly.
* [Fig1 diffusivity and drift.ipynb](examples/numerical%20example/Fig1%20diffusivity%20and%20drift.ipynb): Plot the drift and diffusivity.
* [Fig2 exact vs approximate propagator.ipynb](examples/numerical%20example/Fig2%20exact%20vs%20approximate%20propagator.ipynb): Compare the exact and perturbative propagators.
* [Fig3 exact vs approximate Kramers-Moyal coefficients.ipynb](examples/numerical%20example/Fig3%20exact%20vs%20approximate%20Kramers-Moyal%20coefficients.ipynb): Compare the exact and perturbative finite-time Kramers-Moyal coefficients.
* [Fig4 exact vs approximate medium entropy production.ipynb](examples/numerical%20example/Fig4%20exact%20vs%20approximate%20medium%20entropy%20production.ipynb): Compare the exact and perturbative medium entropy production rate.
* [Fig5 exact vs second order propagators.ipynb](examples/numerical%20example/Fig5%20exact%20vs%20second%20order%20propagators.ipynb): Compare exact and various second-order perturbative propagators.


## Installation

PySTFP requires sympy and numpy. To install these requirements as well as PySTFP, you can run the following commands.

```bash
>> git clone https://github.com/juliankappler/PySTFP.git .
>> cd PySTFP
>> pip install -r requirements.txt
>> pip install .
```

## References

<a id="ref_1">[1] **Short-time Fokker-Planck propagator beyond the Gaussian approximation**. J. Kappler. arXiv: add link once paper has been uploaded to the arXiv.</a>

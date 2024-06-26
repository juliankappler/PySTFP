# PySTFP: Python module for the short-time Fokker-Planck propagator

## About

This module contains sympy-based symbolic tools for the short-time propagator for one-dimensional Fokker-Planck dynamics. This python module accompanies the paper Ref. <a href="#ref_1">[1]</a>.

## Introduction

For a given drift $a(x)$ and diffusivity $D(x)$ we consider the Fokker-Planck equation

$$\partial_t P = - \partial_x \left( a P\right) + \partial_x^2 \left(D P\right),$$

where $P \equiv P(x,t \mid x_0,t_0)$ is the transition density, or propagator, to find a particle that started at time $t_0$ at location $x_0$ at a later time $t$ at a location $x$. This means we consider the initial conditions $P(x,t_0\mid x_0, t_0) = \delta (x-x_0)$,  where $\delta$ denotes the Dirac-delta distribution. As spatial domain we consider $x \in \mathbb{R} = (-\infty,\infty)$, with boundary conditions $P(x,t\mid x_0, t_0) \rightarrow 0$ as $|x| \rightarrow \infty$.

In this module, we implement an approximate short-time propagator $P_K(x,t\mid x_0,t_0)$, in both the normalization-preserving representation

$$P_K(x,t \mid x_0 ,t_0) = \frac{1}{\sqrt{ 4 \pi D(x_0)  \Delta t}} \exp\left[ - \frac{ \Delta x^2}{4 D(x_0) \Delta t} \right] \times \left[ 1 +\sum_{k=1}^K \sqrt{\Delta t}^k \mathcal{Q}_k\left(\tilde{x},x_0\right)\right],$$

and the positivity-preserving representation

$$P_K(x,t \mid x_0 ,t_0) = \frac{1}{\sqrt{ 4\pi D(x_0)  \Delta t}} \exp\left[ - \frac{ \Delta x^2}{4 D(x_0) \Delta t} + \sum_{k=1}^K \sqrt{\Delta t}^k \hat{\mathcal{Q}}_k\left(\tilde{x},x_0\right)\right].$$

Here, $\Delta x \equiv x - x_0$, $\Delta t = t - t_0$, are the spatial and temporal increments, $\tilde{x} = \Delta x/ \sqrt{ 2 D(x_0) \Delta t}$ is the rescaled spatial increment, and $K \in \mathbb{N}_0$ is the truncation order for the approximate propagator. The explicit form of the coefficients $\mathcal{Q}_k$, $\hat{\mathcal{Q}}_k$, are a main result of Ref. <a href="#ref_1">[1]</a>. The module contains the precalculated coefficients up to $K = 8$ (corresponding to an accuracy $\Delta t^4$ for the propagator), but the module also includes code to calculate the coefficients to even higher order. Note that for a typical spatial increment $\Delta x$, for which the support of the propagator is non-negligible, we have $\tilde{x} \lesssim 1$; this is the reason why the power-series coefficients $\mathcal{Q}_k$, $\hat{\mathcal{Q}}_k$ have $\tilde{x}$ as an argument, see Ref. <a href="#ref_1">[1]</a> for more details.

In PySTFP we furthermore implement the explicit perturbation expansions for the spatial moments $\langle \Delta x^n \rangle$ ($n = 0, 1, 2, 3, 4$),
as well as the medium entropy production rate, the total entropy production rate,
and the Gibbs entropy. See Ref. <a href="#ref_1">[1]</a> and the examples below for more details.


## Examples

In the subfolder [examples/](examples/) we provide several notebooks that illustrate the use of PySTFP, 
as well as allow to re-derive and extend the symbolic perturbative expressions used.

Here is a list of the notebooks:


**Subfolder [examples/](examples/)**

These notebooks showcase basic functionality of PySTFP:

* [plot propagator.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/plot%20propagator.ipynb): Plot the short-time propagator for a given drift and diffusivity profile.
* [propagator in physical units.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/propagator%20in%20physical%20units.ipynb): Access the symbolic propagator stored in PySTFP, which are expressed in dimensionless variables, and rewrite them in physical dimensions.
* [entropy production in physical units.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/entropy%20production%20in%20physical%20units.ipynb): Access the symbolic entropy production (rates) stored in PySTFP, which are expressed in dimensionless variables, and rewrite them in physical dimensions.



**Subfolder [examples/derivations](examples/derivations)**

These notebooks implement the symbolic derivations of all the perturbative results from Ref. <a href="#ref_1">[1]</a>:

* [perturbative propagators and moments.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/derivations/perturbative%20propagators%20and%20moments.ipynb): 
Symbolically evaluate the perturbation series for both the normalization-preserving and positivity-preserving propagator. Additionally, use the normalization-preserving propagator to calculate some moments $\langle \Delta x^k \rangle$.
* [entropy production rates.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/derivations/entropy%20production%20rates.ipynb): Symbolically calculate the medium entropy production rate, the total entropy production rate, and the Gibbs entropy.


**Subfolder [examples/numerical example](examples/numerical%20example)**

These notebooks allow to reproduce the plots of Ref. <a href="#ref_1">[1]</a>:

* [analytical solution.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/numerical%20example/analytical%20solution.ipynb): Evaluate the drift, diffusivity, and propagator exactly.
* [Fig1 diffusivity and drift.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/numerical%20example/Fig1%20diffusivity%20and%20drift.ipynb): Plot the drift and diffusivity.
* [Fig2 exact vs approximate propagator.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/numerical%20example/Fig2%20exact%20vs%20approximate%20propagator.ipynb): Compare the exact and perturbative propagators.
* [Fig3 exact vs approximate Kramers-Moyal coefficients.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/numerical%20example/Fig3%20exact%20vs%20approximate%20Kramers-Moyal%20coefficients.ipynb): Compare the exact and perturbative finite-time Kramers-Moyal coefficients.
* [Fig4 exact vs approximate medium entropy production.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/numerical%20example/Fig4%20exact%20vs%20approximate%20medium%20entropy%20production.ipynb): Compare the exact and perturbative medium entropy production rate.
* [Fig5 exact vs second order propagators.ipynb](https://github.com/juliankappler/PySTFP/blob/main/examples/numerical%20example/Fig5%20exact%20vs%20second%20order%20propagators.ipynb): Compare exact and various second-order perturbative propagators.


## Installation

PySTFP requires sympy and numpy, the example notebooks furthermore use matplotlib and h5py. To install these requirements as well as PySTFP, you can run the following commands.

```bash
>> git clone https://github.com/juliankappler/PySTFP.git .
>> cd PySTFP
>> pip install -r requirements.txt
>> pip install .
```

## References

<a id="ref_1">[1] **Short-time Fokker-Planck propagator beyond the Gaussian approximation**. J. Kappler. arXiv: [2405.18381](http://arxiv.org/abs/2405.18381).</a>

## Acknowledgements

We acknowledge funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 101068745.

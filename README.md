# Riemann-solver
This project come form owner's project reports during studying Advanced Computational Fluid Dynamics (AAE6201-20241-A) in the Hong Kong Polytechnic University.

1-D FDM Riemann solver for sod shock tube and 2-D FVM Riemann solver for bow shock. 

Numerical method for 1-D Riemann solver will be introduced briefly. For more information please refer to ```./doc```.

## Mathematical Derivation

We consider the 1D Euler equations governing the flow of an ideal gas:
![](./img/Project1.png)

1. Conservation of mass:

$$
\begin{equation}
       \frac{\partial \rho}{\partial t} + \frac{\partial (\rho u)}{\partial x} = 0 
\end{equation}
$$

3. Conservation of momentum:

$$
\begin{equation}
      \frac{\partial (\rho u)}{\partial t} + \frac{\partial (\rho u^2 + p)}{\partial x} = 0 
\end{equation}
$$

5. Conservation of energy:

$$
\begin{equation}
      \frac{\partial (\rho E)}{\partial t} + \frac{\partial (\rho u E + p)}{\partial x} = 0
\end{equation}
$$

where $\rho$ is the density, $u$ is the velocity, $p$ is the pressure, and $E$ is the total energy per unit volume.
## Requirements
1. Compare the numerical solution with the exact solution at different time instants (do not let the wave arrive at the boundaries).
2. Use different schemes for space discretization (around 100 grid points).
3. Use first-order difference formula for time discretization.
4. Use S-W or L-F flux vector splitting for original flux and characteristic flux, and compare their difference.

## Exact solution
According to the knowledge of gas dynamics, three types of waves may appear in the Sod shock tube: 
1. shock wave. After passing through the shock wave, the density, velocity, and pressure of the fluid all experience sudden changes, satisfying the Rankine-Hugoniot (R-H) relation; 
2. contact discontinuity. After passing through the contact discontinuity, only the density of the fluid changes suddenly, while the velocity and pressure remain unchanged;
3. expansion wave or rarefaction wave. It is an entropy wave with continuous and smooth internal physical quantities, and the physical quantities at the head and tail are continuous but the derivatives are discontinuous (weak discontinuity), with the Riemann invariants remaining invariant. Considering the general case, there are five possibilities of combination waves in the tube. According to the conservation of mass flux, momentum flux, and energy flux, by taking a control volume moving with the shock wave and sufficiently small in thickness, equations can be written and solved for analysis following different scenarios.
In this project we only consider a circumstance that expansion waves and shock wave occur at two end end of the shock tube and  propagate to different direction.
The relationships can be categorized and written as follows for the shock wave at right and expansion waves at left :
In Region 1 and 3

$$
\begin{equation}
    p^{*}/\left(\rho^{*L}\right)^{\gamma}=p_{1}/\left(\rho_{1}\right)^{\gamma}\\u_{1}+\frac{2c_{1}}{\gamma-1}=u^{*}+\frac{2c^{L}}{\gamma-1}
\end{equation}
$$

in which 

$$
\begin{equation}
    c^{L}=\sqrt{\gamma p^{*}/\rho^{*L}}
\end{equation}
$$

In Region 2 and 4

$$
\begin{equation}
    \begin{aligned}&\rho_{2}\left(u_{2}-Z_{2}\right)=\rho^{*R}\left(u^{*}-Z_{2}\right)\\&\rho_{2}u_{2}\left(u_{2}-Z_{2}\right)+p_{2}=\rho^{*R}u^{*}\left(u^{*}-Z_{2}\right)+p^{*}\\&E_{2}\left(u_{2}-Z_{2}\right)+u_{2}p_{2}=E^{*R}\left(u^{*}-Z_{2}\right)+p^{*}u^{*}\end{aligned}
\end{equation}
$$

In this project, analytic solutions was calculated using Ho Taihsiang's code.
## Numerical simulation
### Basic equations
The original equations can be written as follow:

$$
\begin{equation}
    \frac{\partial \mathbf{U}}{\partial t}+\frac{\partial \mathbf{F}}{\partial x}=0
\end{equation}
$$

in which 

$$
\begin{equation}
   \mathbf{U}=\begin{bmatrix}\rho\\\rho u\\\rho E\end{bmatrix},
   \mathbf{F}=\begin{bmatrix}\rho u\\\rho u^2+p\\\rho uE+pu\end{bmatrix}
\end{equation}
$$

set 

$$
A=\frac{\partial F}{\partial U}
$$

we have 

$$
\begin{equation}
    \frac{\partial\mathbf{U}}{\partial t}+A\frac{\partial\mathbf{U}}{\partial x}=0
\end{equation}
$$

in which 

$$
\begin{equation}
    A=\begin{bmatrix}0&1&0\\\\\left(\frac{\gamma-3}2\right)u^2&(3-\gamma)u&\gamma-1\\\\(\gamma-1)u^3-\gamma uE&\gamma E-\frac{3(\gamma-1)}2u^2&\gamma u\end{bmatrix}
\end{equation}
$$

it should be noticed that $E$ can be written as

$$
\begin{equation}
    E=\frac{p}{(\gamma-1)\rho}+\frac{u^2}{2}
\end{equation}
$$

assuming that the fluid is ideal gas.
$A$ is diagonalizable, 

$$
\begin{equation}
    \mathbf{F}=R\Lambda^+R^{-1} \mathbf{U}+R\Lambda^-R^{-1} \mathbf{U}
\end{equation}
$$

$$
\begin{equation}
    \frac{\partial\mathbf{U}}{\partial t}+\frac{\partial\mathbf{F}^+}{\partial x}+\frac{\partial\mathbf{F}^-}{\partial x}=0
\end{equation}
$$

where

$$
\begin{equation}
    \boldsymbol{F}^+=R\Lambda^+R^{-1} \mathbf{U},\quad\boldsymbol{F}^-=R\Lambda^-R^{-1} \mathbf{U}
\end{equation}
$$

the right eigenvectors

$$
\begin{equation}
    \boldsymbol{R}=\begin{bmatrix}1&1&1\\u-a&u&u+a\\H-ua&u^2/2&H+ua\end{bmatrix}
\end{equation}
$$

### Second order upwind scheme
Consider a general scheme in conservative form
$$
\begin{equation}
\frac{\partial \mathbf{U}}{\partial t}=-\mathbf{A}\frac{\mathbf{U}_{i+\frac{1}{2}}^n-\mathbf{U}_{i-\frac{1}{2}}^n}{\Delta x}
\end{equation}
$$
Consider the second order upwind scheme with limiter
$$
\begin{equation}
  \mathbf{U}_{i+1/2}^n=\mathbf{U}_i^n+\frac12\varphi(r_i)(\mathbf{U}_i^n-\mathbf{U}_{i-1}^n),\quad r_i=\frac{\mathbf{U}_{i+1}-\mathbf{U}_i}{\mathbf{U}_i-\mathbf{U}_{i-1}}
\end{equation}
$$
set r=1 if $\mathbf{U}_i-\mathbf{U}_{i-1}=0$.

$$
\begin{equation}
    \frac{\partial \mathbf{U}}{\partial t}=-\frac{\mathbf{A}}{\Delta x}\left[\left(1+\frac{1}{2}\varphi(r_i)-\frac{1}{2}\frac{\varphi(r_{i-1})}{r_{i-1}}\right)(\mathbf{U}_i^n-\mathbf{U}_{i-1}^n)\right]
\end{equation}
$$
substitute equation 14 into equation 8-10 we have
$$
\begin{equation}
    \frac{\partial \mathbf{U}}{\partial t}=-\frac{\mathbf{A}^{+}}{\Delta x}\left[\left(1+\frac{1}{2}\varphi(r_i)-\frac{1}{2}\frac{\varphi(r_{i-1})}{r_{i-1}}\right)(\mathbf{U}_i^n-\mathbf{U}_{i-1}^n)\right]+\frac{\mathbf{A}^{-}}{\Delta x}\left[\left(1+\frac{1}{2}\varphi(r_i)-\frac{1}{2}\frac{\varphi(r_{i+1})}{r_{i+1}}\right)(\mathbf{U}_i^n-\mathbf{U}_{i+1}^n)\right]
\end{equation}
$$

Which is
$$
\begin{equation}
\begin{aligned}
      \mathbf{U}_{i}^{n+1}=-\mathbf{A}^{+}\frac{\Delta t}{\Delta x}\left[\left(1+\frac{1}{2}\varphi(r_i)-\frac{1}{2}\frac{\varphi(r_{i-1})}{r_{i-1}}\right)(\mathbf{U}_i^n-\mathbf{U}_{i-1}^n)\right] \\
      +\mathbf{A}^{-}\frac{\Delta t}{\Delta x}\left[\left(1+\frac{1}{2}\varphi(r_i)-\frac{1}{2}\frac{\varphi(r_{i+1})}{r_{i+1}}\right)(\mathbf{U}_i^n-\mathbf{U}_{i+1}^n)\right]+\mathbf{U}_{i}^n
\end{aligned}
\end{equation}
$$
### Van Leer limiter
Use Van Leer’s limiter, 
$$
\begin{equation}
  \varphi(r)=\frac{r+|r|}{1+|r|}
\end{equation}
$$

### Lax-Wendroff scheme
Discretize equation 9 with Lax-Wendroff scheme
$$
\begin{equation}
\mathbf{U}_{i}^{n+1}=\mathbf{U}_{i}^{n}-\mathbf{A}\frac{\Delta t}{2 \Delta x}(\mathbf{U}_{i+1}^{n}-\mathbf{U}_{i-1}^{n})+\mathbf{A}^2\frac{\Delta t^2}{2\Delta x^2}(\mathbf{U}_{i+1}^{n}-2\mathbf{U}_{i}^{n}+\mathbf{U}_{i-1}^{n})
\end{equation}
$$
### Steger-Warming flux vector splitting
As mentioned before, the Jacobian matrix A for Euler equations can be diagonalized.
$$
\begin{equation}
  \mathbf{A}=\mathbf{P} \Lambda \mathbf{P}^{-1}
\end{equation}
$$
with the left eigenvectors
$$
\begin{equation}
    \boldsymbol{R}^{-1}=\frac{(\gamma-1)}{2a^2}\begin{bmatrix}H+\frac{a(u-a)}{\gamma-1}&-u-\frac{a}{\gamma-1}&1\\\frac{4}{\gamma-1}a^2-2H&2u&-2\\H-\frac{a(u-a)}{\gamma-1}&-u+\frac{a}{\gamma-1}&1\end{bmatrix}
\end{equation}
$$
split the eigenvalues as
$$
\begin{equation}
    \boldsymbol{\Lambda}=\boldsymbol{\Lambda}^++\boldsymbol{\Lambda}^-\quad\boldsymbol{\Lambda}^+=\begin{bmatrix}\lambda_1^+&&\\&\ddots&\\&&\lambda_m^+\end{bmatrix}\quad\boldsymbol{\Lambda}^-=\begin{bmatrix}\lambda_1^-&&\\&\ddots&\\&&\lambda_m^-\end{bmatrix}
\end{equation}
$$
use Steger-Warming scheme to evaluate $\lambda^+$ and $\lambda^-$
$$
\begin{equation}
    \lambda_i^+=\frac12(\lambda_i+|\lambda_i|),\lambda_i^-=\frac12(\lambda_i-|\lambda_i|)
\end{equation}
$$
### Lax-Friedrichs flux vector splitting
Lax-Friedrichs flux vector splitting the positive and negative Jacobian matrix are created following
$$
\begin{equation}
    A^+=\frac12(A+\lambda_{max}I)
\end{equation}
$$
as well as 
$$
\begin{equation}
    A^-=\frac12(A-\lambda_{max}I)
\end{equation}
$$
For a local L-F splitting, $\lambda_{max}$ is evaluated at each point.





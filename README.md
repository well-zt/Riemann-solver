# Riemann-solver
This project come form owner's project reports during studying Advanced Computational Fluid Dynamics (AAE6201-20241-A) in the Hong Kong Polytechnic University.

This project includes two main parts, 1-D FDM Riemann solver for sod shock tube and 2-D FVM Riemann solver for bow shock. 

## 1-D sod shock tube

We consider the 1D Euler equations governing the flow of an ideal gas:

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

![](./img/Project1.jpg)

## Simulation result
![](./img/project1_space.png)

## 2-D bow shock
For more information please refer to ```./doc```.

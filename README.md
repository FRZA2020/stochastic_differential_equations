# stochastic_differential_equations


This R script is an implementation of the implicit Euler method for a linear 2 dimensional differential equation. 
It uses the sde package for simulation a Brownian motion.
The method is described in: P. Kloeden, Numerical Solution of Stochastic Differential Equations (Chapter 12)
The plots show the solution of the differential equation. component11 is the first, component21 the second component of the solution calculated via the implicit euler (3 samples). EXPMA is the exact solution and should be identical to the first solution, because we have reused its sample of the Brownian Motion in the exact solution. The fifth plot is the difference between the exact solution and the first solution and should be close to zero.

# FLDHM
Supplementary material for the paper "Functional Linear Dynamics for Human Movement" by Arlo Hook, Edward Gunning, Giles Hooker, Mark Watsford, Paul Pen Wu, John Warmenhoven.

This repository contains code and sythetic data to estimate a single joint model of sagital ankle kinematics and a multi-joint model of the saggital hip, knee, and ankle kinematics during walking gait using functional linear ordinary differential equations. Results produced by this code will not mimic the exact results found in the paper and should instead be used as a resource for practitioners wishing to implement this modelling technique.

## Synthetic Data Generation

Data used in this repository is synthetic data generated from the multi-joint model estimated in the paper. This is done by taking the estimated parameters and forming the coupled dynamic system which is then solved using initial conditions for each variable sampled from a normal distribution with mean and variance of the inital conditions of the real data. Additionally an error term is added by generating smooth errors from the empirical error covariance obtained from the residuals of solutions to the ODE using the true initial conditions of the data.

# Thermal Analysis of a Multi-Core Processor

This MATLAB script simulates the transient heat conduction in a 2D chip containing two heat-generating cores. The heat source switches from one core to the other after 0.05 seconds, demonstrating thermal crosstalk and the evolution of temperature distribution over time. The simulation uses an explicit finite difference method (Euler time integration) with Dirichlet boundary conditions.
Features
2D heat diffusion solved on a uniform grid using finite differences.
Two separate core regions with time‑dependent heat generation.
Constant temperature (Dirichlet) boundary conditions at the chip edges.
Real‑time visualisation – the script updates two plots every 20 time steps:
    A 2D heat map (temperature colormap, axes in cm).
    A 3D surface plot of the temperature distribution.
 Easily adjustable parameters – chip size, grid resolution, material properties, heat flux, and simulation time can be modified at the top of the script.

# Lyman-Alpha-Scatter
A Monte Carlo model that simulates the "scattering" effect of Lyman Alpha photons in nebulae of neutral hydrogen using the Matplot library in Python. 

This code consists of two files, a 2-dimensional one and a 3-dimensional one. The 3-dimensional model has a very long run time and will not produce files simulating more than fifty particles cumulatively.

Both work in the same way. A larger function called "simulate" wraps around the parameters and commands. "Simulate" has one parameter: the amount of images you want to produce, which is defined outside of the function by the variable "frames." A while loop calls the function "frame" amount of times to create iterative images built on the simulated data. Data is plotted with the quiver function, and every coordinate and component of velocity is held in its own array within a larger array. These values change per condition, velocity/position and collisions. Collisions between Ly-a and hydrogen particles are simulated using a for loop that re-calculates the velocity and position of a Lyman alpha photon that collides with a hydrogen particle. This section can be commented out of the 3d simulation to yield actual results. 

The 2d simulation uses polar coordinates, while the 3d simulation uses spherical coordinates. 

There is an external object component that prints background photons onto the image. It does not interact with the other particles, and its location and size can be changed. 

The size of the box of the simulation, the number of photons and atoms and the temperature can be changed.

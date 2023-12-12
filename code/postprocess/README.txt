Post-Processing Scripts
-------------------------------------

LoadData.m 
----------
The latest version of /src/rigid2D.m outputs options and parameters as a .mat file for each simulation.
It also outputs separate data files for particle velocities and density profiles.

It the flag for tracers is turned on, then tracer files will be written.

This script can be used for calculate fluid pressure, stress tensor, and velocity. 
The script for calculating Yukawa solutions is included as well.

MidArcLen.m
-----------
For a biayer structure, this script calculates the midplane configurations, midplane arc length, midplane enclosed area, and the reduced area by using fft. The test data file is included.


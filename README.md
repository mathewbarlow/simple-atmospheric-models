# simple-atmospheric-models
python code for simple models of the atmosphere and ocean

inertial_oscillations_1.0.py numerically calculates trajectories for inertial oscillations, given a starting location and wind.  The analytic solution for constant f is included for comparison and validation. The code solves the equation set:  
<br> du/dt = fv
<br> dv/dt = -fu

<img src="output-figures-animation/traj.png" width="200" height="200">

inertial_oscillations_anim_1.0.py does the same calculation, but makes an animated gif of the trajectory

<img src="output-figures-animation/inert.gif" width="200" height="200">

If you find any errors, please contact me at Mathew_Barlow@uml.edu

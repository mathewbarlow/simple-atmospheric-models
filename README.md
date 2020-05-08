# simple-atmospheric-models
Code for simple models of the atmosphere (and ocean)

This is a repository for code I have used to make class examples; some animations using the output of the codes is shown further below. Unfortunately, several of the programs need to be completely re-coded for clarity and none of them have been carefully tested, so please use with caution. My goal is to eventually produce well-tested Python versions of the all programs here but I haven't figured out a way to get funding for this effort, so it remains a very slowly moving hobby project. If you find any errors, please contact me at Mathew_Barlow@uml.edu

<b>inertial_oscillations_1.1.py</b>
</br>calculates trajectories for inertial oscillations

<b>baro_1d_1.0.f</b>
</br>solves the linearized barotropic vorticity equation in one spatial dimension (longitude)

</br>
</br>
<b>More Information</b>
</br>
<b>inertial_oscillations_1.1.py</b> numerically calculates trajectories for inertial oscillations, given a starting location and wind.  The analytic solution for constant f is included for comparison and validation. The code solves the equation set:  <br> du/dt = fv
<br> dv/dt = -fu

<img src="output-figures-animation/traj.png" width="200" height="200">

<b>inertial_oscillations_anim_1.1.py</b> does the same calculation, but makes an animated gif of the trajectory

<img src="output-figures-animation/inert.gif" width="200" height="200">



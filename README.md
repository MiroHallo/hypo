# Bayesian earthquake hypocenter location from wave arrival times
Robust 3D Bayesian hypocenter localization from body P- and S-wave
arrival times in a layered 1D medium including full Uncertainty Quantification
***************************************

This toolbox provides a robust framework for 3D earthquake hypocenter localization 
and rigorous Uncertainty Quantification (UQ) in 1D layered media. This toolbox utilizes 
body P- and S-wave arrival times to assess hypocenter locations through a grid search 
method. A key feature of the implementation is the estimation of arrival time 
uncertainties via Monte Carlo simulations and ray tracing within randomly perturbed 
velocity models. Following the probabilistic inverse theory of Tarantola (2005), these 
uncertainties are rigorously evaluated within a Bayesian framework to determine the 
final posterior Probability Density Function (PDF). This approach ensures that the 
solution provides not only the most likely coordinates but also an assessment of the 
location uncertainty, making it an ideal tool for both seismic research and industrial 
applications.

1 METHODOLOGY
===================

  Tarantola, A. (2005, Chapter 7.1). Inverse Problem Theory and Methods 
for Model Parameter Estimation, Society for Industrial and Applied 
Mathematics, Philadelphia, USA.

  Hallo, M., Opršal, I., Asano, K., Gallovič, F. (2019). Seismotectonics
of the 2018 Northern Osaka M6.1 earthquake and its aftershocks: joint
movements on strike-slip and reverse faults in inland Japan, Earth,
Planets and Space, 71:34. [https://doi.org/10.1186/s40623-019-1016-8](https://doi.org/10.1186/s40623-019-1016-8)

2 TECHNICAL IMPLEMENTATION
===================

Bayesian Inference, Grid Search (3D), Ray Tracing (1D), Monte Carlo simulations, 
Uncertainty Quantification, Cross-Platform (Windows, Linux)

The official software version is archived on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19343031.svg)](https://doi.org/10.5281/zenodo.19343031)

3 PACKAGE CONTENT
===================

  1. `crustal.dat` - Input 1D velocity model of Earth's crust
  2. `a1_Pert_Model.m` - Perturb velocity model and perform Monte Carlo simulations to get uncertainties of P and S-wave arrival times
  3. `a2_LatLon2xyz.m` - Converts coordinates of stations from WGS84 Lat/Lon into local coordinates Easting/Northing
  4. `a3_hypo_homo.m` - Compute earthquake hypocenter location by grid search method in the homogeneous medium
  5. `a4_hypo_layer.m` - Compute earthquake hypocenter location by grid search method in the 1D layered medium

4 REQUIREMENTS
===================
  
  MATLAB: Version R2016b, Codes do not require any additional Matlab Toolboxes

5 COPYRIGHT
===================

Copyright (C) 2017-2019  Miroslav Hallo

This program is published under the GNU General Public License (GNU GPL).

This program is free software: you can modify it and/or redistribute it
or any derivative version under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version.

This code is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
and don't remove their names from the code.

You should have received copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

# Earthquake hypocenter location
Bayesian earthquake hypocenter location from direct P- and S-wave
arrival times in a layered 1D medium
***************************************

  Open-source Matlab functions for assessment of the earthquake
hypocenter location by grid search method in the homogeneous or 1D
layered medium. The posterior probability density function (PDF) of the
earthquake hypocenter location is following theory by Tarantola (2005).
Hence, the solution contains also information about hypocenter location
uncertainty.

1 METHODOLOGY
===================

  Hallo, M., Opršal, I., Asano, K., Gallovič, F. (2019). Seismotectonics
of the 2018 Northern Osaka M6.1 earthquake and its aftershocks: joint
movements on strike-slip and reverse faults in inland Japan, Earth,
Planets and Space, 71:34.

2 PACKAGE CONTENT
===================

  a) "crustal.dat" - 1D velocity model of Earth's crust
  
  b) "a1_Pert_Model.m" - Matlab function that perturb velocity model and
  perform Monte-Carlo simulations to get uncertainties of P and S wave
  arrival times
  
  c) "a2_LatLon2xyz.m" - Matlab function that converts coordinates of
  stations from WGS84 Lat/Lon into local coordinates Easting/Northing
  
  d) "a3_hypo_homo.m" - Matlab function that compute earthquake
  hypocenter location by grid search method in the homogeneous medium
  
  e) "a4_hypo_layer.m" - Matlab function that compute earthquake
  hypocenter location by grid search method in the 1D layered medium

3 REQUIREMENTS
===================

  MATLAB: Codes do not require any additional Matlab Toolboxes.

4 COPYRIGHT
===================

Copyright (C) 2017,2018  Miroslav Hallo

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

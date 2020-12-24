# Master-Equation-for-Quantum-Nonlinear-Mixing-of-Thermal-Photons

A doubly resonant nonlinear nanophotonic system, when heated, can emit thermal radiation beyond the known limit imposed by Kirchhoff's laws.
This finding of surpassing the blackbody limit on thermal emission was recently published in our work:
https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-28-2-2045&id=425964

This repository contains a script to reproduce the subfigures in fig3 of the paper. 
You can simply run the script 'produce_fig3' in matlab to see how nonlinear mixing of photons at first harmonic (w1) can be used to enhance the thermal emission at its second harmonic (w2=2*w1).  

The above script uses a function 'thermalshg' which computes various steady state quantities in presence of nonlinear mixing. 
One can also use this script to validate the consistency of this master equation approach (from quantum optics) with thermodynamic laws by showing:
1. At thermal equilibrium, there is no heat flow between the system and the environment.
2. Under nonequilibrium, heat always flows from high temperature body to low temperature body. 

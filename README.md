# MC Simulation Programm   

This is a Monte Carlo Simulation Programm (Fortran implementation) to compute the 
  - internal energy (with Ewald summation) 
  - radial distribution function 
  - relative permittivity 
  - 
for fluid mixtures of dipolar and charged hard spheres 

Note: 
  - The coupling parameter λ linearly incorporates the electrostatic pair interactions (within the acceptance criterion) to compute a λ-scaled canonical ensemble average of the internal energy (necessary for thermodynamic integration, see papers) 

Written by: M. Theiss 

To-Do's:
  - comments: work in progress

The implementation aligns with the following literature. 

[1] Theiss, M. & Gross, J., Dipolar Hard Spheres: Comprehensive Data from Monte Carlo Simulations, Journal of Chemical & Engineering Data, 64 (2),
827-832 (2019)

[2] Theiss, M. & Gross, J., Nonprimitive Model Electrolyte Solutions: Comprehensive Data from Monte Carlo Simulations, Journal of Chemical &
Engineering Data, 65 (2), 634-639 (2020)

## execute programm

thermodynamic and simulation input parameters in INPUT.inp file

1. compile with "make"
2. run with "./main INPUT.inp"

About
-----

This is a modified version of the ABC Quantum Reactive Scattering Program
originally developed by the group of D. E. Manolopoulos [1].



Modifications
-------------

1) Limit of 1000 Mb scratch files removed (setrec routine).

2) All local arrays are dynamically allocated.

3) The 'nmax' parameter moved to ABC's input file (driver routine).

4) Some Fortran 77 features promoted to Fortran 90 (comments, continuation line,
source code extensions etc).

5) Added a new routine named pes() (pes.f90) for user's defined potential
energy surfaces (potsub routine).

6) Added wrappers to invoke linear algebra operations from the Magma library
(hybrid CPU+GPU), if desired.



References
----------

[1] D. Skouteris et al., Computer Physics Communications 133 128–135 (2000).

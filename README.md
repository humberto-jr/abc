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



How to build
-----

First compile the external user defined PES into object files *.o (-c flag).
Then, ABC is built for many scenarios as follows:

1. Intel ifort + Intel MKL (default):
`make`

Option (1) requires that the env variable MKLROOT is properly defined. Check the
content of `$MKLROOT` first.

If instead GNU gfortran is intended,

2. GNU gfortran + Intel MKL:
`make FC=gfortran`

Make sure both the PES and ABC are compiled with the same kind of compiler in
order to avoid missing libraries during the link.

3. Intel ifort + Intel MKL + MAGMA:
`make use_magma LINEAR_ALGEBRA=MAGMA`

If MAGMA is not installed in standard locations, you can provide the path to its
components;

4. Intel ifort + Intel MKL + MAGMA:
`make use_magma LINEAR_ALGEBRA=MAGMA MAGMAROOT=[path/to/magma-2.5.2] CUDAROOT=[path/to/cuda/10.1]`

Wrappers for MAGMA are written in C and are built with GNU gcc. If Intel icc is
needed, then

5. Intel ifort + Intel icc + Intel MKL + MAGMA:
`make use_magma CC=icc LINEAR_ALGEBRA=MAGMA MAGMAROOT=[path/to/magma-2.5.2] CUDAROOT=[path/to/cuda/10.1]`

6. GNU gfortran + Intel icc + Intel MKL + MAGMA:
`make use_magma FC=gfortran CC=icc LINEAR_ALGEBRA=MAGMA MAGMAROOT=[path/to/magma-2.5.2] CUDAROOT=[path/to/cuda/10.1]`



References
----------

[1] D. Skouteris et al., Computer Physics Communications 133 128–135 (2000).

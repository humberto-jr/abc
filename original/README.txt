README file for the CCP6 Quantum Reactive Scattering Program ABC
----------------------------------------------------------------


1. To extract the ABC distribution from abc.tar.gz, execute the
   following commands:

   gunzip abc.tar.gz
   tar -xvf abc.tar



2. This should create a directory abc in your present working
   directory that contains the following files:

   abc.f     The main ABC program
   fun.f     Additional special function subroutines
   lin.f     Linear algebra subroutines
   pot.f     Potential energy surface subroutines

   calpha    An example compile script for Compaq computers
   csgi      An example compile script for SGI computers

   BW.3p     3-body parameters for the BW Cl+H2 potential
   SW.2p     2-body parameters for the SW F+H2 potential
   SW.3p     3-body parameters for the SW potential

   fhd.d     An example data set for the F+HD reaction [1]
   fhd.out   ... and the corresponding output file

   clhd.d    An example data set for the Cl+HD reaction [2]
   clhd.out  ... and the corresponding output file

   hd2.d     An example data set for the H+D2 reaction [3] 
   hd2.out   ... and the corresponding output file



3. In order to compile the ABC program on a SGI workstation
   with an f90 compiler and the -lscs library (Silicon Graphics
   -Cray Scientific Library), type 

   cd abc
   csgi

   On a Compaq alpha workstation with an f90 compiler and the 
   -ldxml library (Digital Extended Maths Library) type

   cd abc
   calpha

   These commands should generate the executable program abc.x. 



4. If you do not have either of the above computers, but you
   do have a unix workstation that has an f90 compiler and
   BLAS and LAPACK libraries, you could try

   f90 abc.f fun.f lin.f pot.f -o abc.x -llapack -lblas.

   If you don't have a Fortran 90 compiler, but you do have
   a Fortran 77 compiler that supports the automatic array
   extension to the f77 standard (such as the GNU/Linux g77
   compiler), you could try

   f77 abc.f fun.f lin.f pot.f -o abc.x -llapack -lblas.

   If you are lucky, one or other of these commands will work
   correctly and generate an acceptable executable file abc.x.



5. To ensure that the program is working properly, try

   time abc.x < fhd.d > fhd.new  
   diff fhd.new fhd.out

   time abc.x < clhd.d > clhd.new
   diff clhd.new clhd.out

   time abc.x < hd2.d > hd2.new
   diff hd2.new hd2.out

   There shouldn't be any differences. 

   For reference, the (user) times taken to run these three
   examples on a Compaq DS20 workstation are approximately
   30 mins, 4 hours 15 mins, and 1 hour 15 mins, respectively.

 

6. More information about how to use the ABC program once it
   is working correctly is given in Ref. [4]. 



References
----------

[1] R.T. Skodje et al., J. Chem. Phys. 112 (2000) 4536.

[2] D. Skouteris et al., Science 286 (1999) 1713.

[3] M.P. de Miranda et al., J. Chem. Phys. 108 (1998) 3142. 
   
[4] D. Skouteris et al., ABC: A Quantum Reactive Scattering Program, 
    Computer Physics Communications (submitted)

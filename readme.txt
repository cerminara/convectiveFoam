Solver for the incompressible Bounssinesq equation, based on pimpleFoam and
buoyantBoussinesqPimpleFoam. See www.openfoam.org.

In this solver, the buoyancy term is written explicitly in his original
formulation.

The test case is analogous to the convective plume
described in Whitehead et al. 2013, doi:10.1017/jfm.2013.330.

/* ---- *\
Equations:

div(U) = 0              ---- (PISO algorithm)
ddt(U) + div( U x U ) = -grad(p) + div(stress) - beta(T-TRef) g
ddt(T) + div(U T) = div(diffusion)

where beta is the thermal expansion.

In this new version the gravity is a function of the vertical coordinate, in
particular:
g = g0 ( 1 + gamma (g0 . x)/|g0| )
where g0 and gamma are given in the transportProperties dictionary.
\* ---- */


How to use convectiveFoam?

- clone the source files in your preferred directory:
    git clone https://github.com/cerminara/convectiveFoam.git

- Be sure the environment variables of OpenFoam 4.x have been
  loaded in your terminal window

- enter in the convectiveFoam folder and install it running 
  the script:
    ./Allwmake

- try it with the singlePlume case, in your tutorial folder, by
  using the following commands:
    cd $WM_PROJECT_USER_DIR/tutorials/convectiveFoam/singlePlume/Pr1_Ra1e6_limLin
    ./Allrun

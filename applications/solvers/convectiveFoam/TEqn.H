{
    //alphat = turbulence->nut()/Prt;
    //alphat.correctBoundaryConditions();

    //volScalarField alphaEff("alphaEff", turbulence->nu()/Pr);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alpha, T)
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);
}

// Solve the Momentum equation


MRF.correctBoundaryVelocity(U);

dimensionedScalar Dt("Dt", dimTime, runTime.deltaTValue());
tmp<fvVectorMatrix> tUEqn
(
//    fvm::ddt(U) + fvm::div(phi, U)
//  + MRF.DDt(U)+ 
    //turbulence->divDevReff(U)
//    turbulence->nu()*fvm::laplacian(U)

    (fvm::Sp(scalar(1)/Dt, U) - U/Dt)
  - fvm::laplacian(turbulence->nu(), U)
  - fvc::div(turbulence->nu()*dev(fvc::grad(U)().T()))
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    //for(int ii=0; ii<10; ii++)
    {
        solve
        (
            UEqn
         ==
           -fvc::grad(p)
          + fvc::reconstruct
            (
                -fvc::interpolate(beta*(T-TRef))*(g & mesh.Sf())
            )
        );
    }

    fvOptions.correct(U);
}

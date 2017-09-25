//volScalarField momRes
//(
//    (
//       -fvc::laplacian(turbulence->nu(), U)
//      - fvc::div(turbulence->nu()*dev(fvc::grad(U)().T()))
//      + fvc::grad(p)
//      - fvc::reconstruct
//        (
//            -fvc::interpolate(beta*(T-TRef))*(g & mesh.Sf())
//        )
//    ).component(2)
//);
//
//volScalarField momB
//(
//    (
//        fvc::grad(p)
//      - fvc::reconstruct
//        (
//            -fvc::interpolate(beta*(T-TRef))*(g & mesh.Sf())
//        )
//    ).component(2)
//);

fvVectorMatrix UEqn1
(
   -fvm::laplacian(turbulence->nu(), U)
  - fvc::div(turbulence->nu()*dev(fvc::grad(U)().T()))
 ==
   -fvc::grad(p)
  + fvc::reconstruct
    (
        -fvc::interpolate(beta*(T-TRef))*(g & mesh.Sf())
    )
);

//scalar maxRes = max(momRes.internalField()).value();
//scalar meanB = sum(momB.internalField()*mesh.V()).value()/sum(mesh.V()).value();
//scalar res = solve(UEqn1).max().initialResidual();
//vectorField res = UEqn1.residual();
volVectorField URef
(
    IOobject
    (
        "URef",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("URef", dimensionSet(0,1,-1,0,0,0,0), gAverage(U))
);
vectorField wA = (UEqn1 & U);
vectorField pA = (UEqn1 & URef);

scalar normFactorX = gSum(mag(wA.component(0) - pA.component(0)) + mag(UEqn1.source().component(0) - pA.component(0))); 
scalar residualX = gSum(mag(UEqn1.source().component(0) - wA.component(0)))/(normFactorX + SMALL); 
scalar normFactorY = gSum(mag(wA.component(1) - pA.component(1)) + mag(UEqn1.source().component(1) - pA.component(1))); 
scalar residualY = gSum(mag(UEqn1.source().component(1) - wA.component(1)))/(normFactorY + SMALL); 
scalar normFactorZ = gSum(mag(wA.component(2) - pA.component(2)) + mag(UEqn1.source().component(2) - pA.component(2))); 
scalar residualZ = gSum(mag(UEqn1.source().component(2) - wA.component(2)))/(normFactorZ + SMALL); 

Info
  << "res: "
  << "(" << residualX << ", " << residualY << ", " << residualZ << ")"
  << endl;


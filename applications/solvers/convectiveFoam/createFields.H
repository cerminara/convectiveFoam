Info<< "Reading thermophysical properties\n" << endl;

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Thermal diffusivity [m2/s]
dimensionedScalar alpha
(
    "alpha",
    dimLength*dimLength/dimTime,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Gravity
dimensionedVector g0("g0", dimAcceleration, laminarTransport);
dimensionedScalar gamma("gamma", dimless/dimLength, laminarTransport);

//Info << g0 << "   " << gamma << endl;
//volVectorField g("g", g0*(scalar(1) + (g0 & mesh.C())/mag(g0)*gamma));
surfaceVectorField gf("gf", g0*(scalar(1) + (g0 & mesh.Cf())/mag(g0)*gamma));
Info << "min(g) = " << gMin(gf) << ",    max(g) = " << gMax(gf) << endl;

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// kinematic turbulent thermal thermal conductivity m2/s
//Info<< "Reading field alphat\n" << endl;
//volScalarField alphat
//(
//    IOobject
//    (
//        "alphat",
//        runTime.timeName(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh
//);


//#include "readGravitationalAcceleration.H"


#include "createMRF.H"

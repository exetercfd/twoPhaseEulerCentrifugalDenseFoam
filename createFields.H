#include "readGravitationalAcceleration.H"
#include "readhRef.H"
 
Info<< "Creating twoPhaseSystem\n" << endl;
 
twoPhaseSystem fluid(mesh, g);
 
phaseModel& phase1 = fluid.phase1();
phaseModel& phase2 = fluid.phase2();
 
volScalarField& alpha1 = phase1;
volVectorField& U1 = phase1.U();
volVectorField& U2 = phase2.U();
 
volScalarField& p = phase1.thermo().p();
 
dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);
 
#include "gh.H"
 
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
 
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.U()
);
 
 
// APR ADDED thetaMax and omega
 
Info<< "Reading gravityProperties\n" << endl;
 
IOdictionary gravityProperties
(
    IOobject
    (
        "gravityProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
 
dimensionedScalar coneAngle
(
    gravityProperties.lookupOrDefault("coneAngle", 0.0)
);
 
dimensionedScalar omega
(
    gravityProperties.lookupOrDefault("omega", 0.0)
);
 
dimensionedVector omega_vector
(
    gravityProperties.lookupOrDefault("omega_vector", vector::zero)
);

const scalar PI = constant::mathematical::pi; 

Info << "phase 1 is\n" << phase1.name() << endl;

Info << "coneAngle in radians\n" << coneAngle*PI/180.0 << endl;

Info << "coneAngle in degrees\n" << coneAngle << endl;
 
Info << "omega\n" << omega << endl;
 
Info << "omega_vector\n" << omega_vector << endl;
 
volVectorField acceleration
(
    IOobject
    (
        "acceleration",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("acceleration", dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0))
);
 
volVectorField centrifugalAcceleration
(
    IOobject
    (
        "centrifugalAcceleration",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("centrifugalAcceleration", dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0))
);

// APR ADDED thetaMax and omega
 
 
 
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());
 
Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);
 
 
Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K1(IOobject::groupName("K", phase1.name()), 0.5*magSqr(U1));
volScalarField K2(IOobject::groupName("K", phase2.name()), 0.5*magSqr(U2));
 
#include "createMRF.H"

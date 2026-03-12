Info<< "Reading thermal properties\n" << endl;
IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoCp(thermalProperties.lookup("rhoCp"));//rho*Cp used to convert volumetric heating into the temperature source term
dimensionedScalar kf(thermalProperties.lookup("kf"));//Thermal Conductivity of fluid
dimensionedScalar ks(thermalProperties.lookup("ks"));//Thermal Conductivity of solid
dimensionedScalar heatSource(thermalProperties.lookup("heatSource"));//uniform volumetric heat generation (e.g. W/m^3)
dimensionedScalar thermalSource("thermalSource", heatSource/rhoCp);//temperature-equation source term: q'''/(rho*Cp)
dimensionedScalar unitThermalSource("unitThermalSource", thermalSource.dimensions(), 1.0);//unit source used for adjoint forcing
dimensionedScalar b2("b2", dimensionSet(0,0,-1,2-Pnorm,0,0,0),1);
Switch DDHS(thermalProperties.lookup("DDHS"));//design dependent heat source
volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    (ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoCp
);
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
Info<< "Reading field Tb\n" << endl;
   volScalarField Tb
   (
       IOobject
       (
           "Tb",
           runTime.timeName(),
           mesh,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh
   ); 
volScalarField T2(T);
volScalarField fsensMeanT
(
    IOobject
    (
        "fsensMeanT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField fsenshMeanT(fsensMeanT);

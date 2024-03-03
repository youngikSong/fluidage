/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "timeTransport.H"
#include "surfaceFields.H"
#include "fvmDdt.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvcFlux.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "incompressibleMomentumTransportModel.H"
#include "compressibleMomentumTransportModel.H"

#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "interfaceCompression.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(timeTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        timeTransport,
        dictionary
    );
}
}


template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::timeTransport::diffusivityType,
    3
>::names[] =
{
    "none",
    "constant",
    "viscosity"
};

const Foam::NamedEnum
<
    Foam::functionObjects::timeTransport::diffusivityType,
    3
> Foam::functionObjects::timeTransport::diffusivityTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::functionObjects::timeTransport::D() const
{
    const word Dname("D" + fieldName_);
    
    //do you have diffusivityType::constant stored?
    if (diffusivity_ == diffusivityType::constant) 
    {
        return volScalarField::New
        (
            Dname,
            mesh_,
            dimensionedScalar(Dname, dimViscosity, D_)
        );
    }
    //if no diffusivity data stored, search for momentumTransportModel file
    // alphat is kinematic turbulent thermal conductivity, alphal maybe for laminar? idk
    else 
    {
        const momentumTransportModel& turbulence =
            mesh_.lookupType<momentumTransportModel>();

        return alphal_*turbulence.nu() + alphat_*turbulence.nut();
    } 
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// & as reference, see https://www.learncpp.com/cpp-tutorial/lvalue-references/
Foam::functionObjects::timeTransport::timeTransport 
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
: 
  //member initializer list, to call base class constructors
  //and initialize data members before the body of the constructor executes
  //uses constructor below
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(dict.lookupOrDefault<word>("field", "s")), //fieldname s
    diffusivity_(diffusivityType::none),
    D_(0),
    nCorr_(0),
    s_
    (
        IOobject
        (
            fieldName_,
            time_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    MULES_(false),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    ),
    sRestart_(false)
{
    read(dict);

    if (mesh_.solution().solversDict().found(fieldName_))
    {
        const dictionary& controls = mesh_.solution().solverDict(fieldName_);

        if (controls.found("nSubCycles"))
        {
            MULES_ = true;

            typeIOobject<surfaceScalarField> sPhiHeader
            (
                IOobject::groupName("sPhi", s_.group()),
                runTime.name(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            );

            sRestart_ = sPhiHeader.headerOk();

            if (sRestart_)
            {
                Info << "Restarting s" << endl;
            }

            const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(phiName_);

            // Scalar's volumetric flux
            // this is what phi means https://openfoamwiki.net/index.php/Main_FAQ#What_is_the_field_phi_that_the_solver_is_writing
            tsPhi_ = new surfaceScalarField 
            (
                sPhiHeader,
                phi*fvc::interpolate(s_)
            );

            if (controls.lookupOrDefault<Switch>("MULESCorr", false))
            {
                mesh_.schemes().setFluxRequired(fieldName_);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::timeTransport::~timeTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::timeTransport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    schemesField_ = dict.lookupOrDefault<word>("schemesField", fieldName_); 
    // find schemesfield by name fieldname, which is s

    diffusivity_ = diffusivityTypeNames_.read(dict.lookup("diffusivity"));

    switch(diffusivity_)
    {
        case diffusivityType::none:
            break;

        case diffusivityType::constant:
            dict.lookup("D") >> D_;
            break;

        case diffusivityType::viscosity:
            dict.lookup("alphal") >> alphal_;
            dict.lookup("alphat") >> alphat_;
            break;
    }

    dict.readIfPresent("nCorr", nCorr_);

    return true;
}


Foam::wordList Foam::functionObjects::timeTransport::fields() const
{
    return wordList{phiName_};
}


bool Foam::functionObjects::timeTransport::execute()
{
    Info<< type() << " write:" << endl;

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    const word divScheme("div(phi," + schemesField_ + ")"); //div(phi,s), divergence of phi*s

    // Set under-relaxation coeff
    scalar relaxCoeff = 0.0;
    if (mesh_.solution().relaxEquation(schemesField_))
    {
        relaxCoeff = mesh_.solution().equationRelaxationFactor(schemesField_);
    }

    const Foam::fvModels& fvModels(Foam::fvModels::New(mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(mesh_)
    );
    //for incompressible maybe?
    if (phi.dimensions() == dimVolume/dimTime)
    {
        if (MULES_)
        {
            subCycleMULES();

            fvConstraints.constrain(s_);
        }
        else
        {
            for (int i=0; i<=nCorr_; i++)
            {
                fvScalarMatrix sEqn
                (
                    fvm::ddt(s_)
                  + fvm::div(phi, s_, divScheme)
                  - dimensionedScalar(1) //fvModels.source(s_)
                );
                //source will be 1 for fluid-age, so can use default code

                if (diffusivity_ != diffusivityType::none)
                {
                    sEqn -= fvm::laplacian(D(), s_);
                }

                sEqn.relax(relaxCoeff);

                fvConstraints.constrain(sEqn);

                sEqn.solve(schemesField_);

                fvConstraints.constrain(s_);
            }
        }
    }
    //for compressible maybe
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        for (int i=0; i<=nCorr_; i++)
        {
            fvScalarMatrix sEqn
            (
                fvm::ddt(rho, s_)
              + fvm::div(phi, s_, divScheme)
              - rho //fvModels.source(rho, s_)
            );
            //same here with upper side

            if (diffusivity_ != diffusivityType::none)
            {
                sEqn -= fvm::laplacian(rho*D(), s_);
            }

            sEqn.relax(relaxCoeff);

            fvConstraints.constrain(sEqn);

            sEqn.solve(schemesField_);

            fvConstraints.constrain(s_);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Incompatible dimensions for phi: " << phi.dimensions() << nl
            << "Dimensions should be " << dimMass/dimTime << " or "
            << dimVolume/dimTime << exit(FatalError);
    }

    Info<< endl;

    return true;
}


void Foam::functionObjects::timeTransport::subCycleMULES()
{
    const dictionary& controls = mesh_.solution().solverDict(fieldName_);
    const label nSubCycles(controls.lookup<label>("nSubCycles"));
    const bool LTS = fv::localEulerDdt::enabled(mesh_);

    if (nSubCycles > 1)
    {
        tmp<volScalarField> trSubDeltaT;

        if (LTS)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh_, nSubCycles);
        }

        for
        (
            subCycle<volScalarField> sSubCycle(s_, nSubCycles);
            !(++sSubCycle).end();
        )
        {
            solveMULES();
        }
    }
    else
    {
        solveMULES();
    }


    // Apply the diffusivity term separately to allow implicit solution
    // and boundedness of the explicit advection
    if (diffusivity_ != diffusivityType::none)
    {
        fvScalarMatrix sEqn
        (
            fvm::ddt(s_) - fvc::ddt(s_)
          - fvm::laplacian(D(), s_)
        );

        sEqn.solve(controls.subDict("diffusivity"));

        Info<< fieldName_ << " volume fraction = "
            << s_.weightedAverage(mesh_.V()).value()
            << "  Min(" << fieldName_ << ") = " << min(s_).value()
            << "  Max(" << fieldName_ << ") = " << max(s_).value()
            << endl;
    }
}


void Foam::functionObjects::timeTransport::solveMULES()
{
    const dictionary& controls = mesh_.solution().solverDict(fieldName_);
    const label nCorr(controls.lookup<label>("nCorr"));
    const label nSubCycles(controls.lookup<label>("nSubCycles"));
    const bool MULESCorr(controls.lookupOrDefault<Switch>("MULESCorr", false));

    // Apply the compression correction from the previous iteration
    // Improves efficiency for steady-simulations but can only be applied
    // once the s field is reasonably steady, i.e. fully developed
    const bool applyPrevCorr
    (
        controls.lookupOrDefault<Switch>("applyPrevCorr", false)
    );

    const bool LTS = fv::localEulerDdt::enabled(mesh_);

    const word divScheme("div(phi," + schemesField_ + ")");

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);

    surfaceScalarField& sPhi_ = tsPhi_.ref();

    const word sScheme(mesh_.schemes().div(divScheme)[1].wordToken());

    // If a compressive convection scheme is used
    // the interface normal must be cached
    tmp<surfaceScalarField> nHatf;

    if (compressionSchemes.found(sScheme))
    {
        const surfaceVectorField gradsf(fvc::interpolate(fvc::grad(s_)));

        nHatf = new surfaceScalarField
        (
            IOobject
            (
                "nHatf",
                s_.time().name(),
                mesh_
            ),
            gradsf/(mag(gradsf) + deltaN_) & mesh_.Sf()
        );
    }

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    {
        tmp<fv::ddtScheme<scalar>> tddtS
        (
            fv::ddtScheme<scalar>::New
            (
                mesh_,
                mesh_.schemes().ddt("ddt(s)")
            )
        );
        const fv::ddtScheme<scalar>& ddtS = tddtS();

        if
        (
            isType<fv::EulerDdtScheme<scalar>>(ddtS)
         || isType<fv::localEulerDdtScheme<scalar>>(ddtS)
        )
        {
            ocCoeff = 0;
        }
        else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtS))
        {
            if (nSubCycles > 1)
            {
                FatalErrorInFunction
                    << "Sub-cycling is not supported "
                       "with the CrankNicolson ddt scheme"
                    << exit(FatalError);
            }

            if
            (
                sRestart_
             || time_.timeIndex() > time_.startTimeIndex() + 1
            )
            {
                ocCoeff =
                    refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtS)
                   .ocCoeff();
            }
        }
        else
        {
            FatalErrorInFunction
                << "Only Euler and CrankNicolson ddt schemes are supported"
                << exit(FatalError);
        }
    }

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    tmp<surfaceScalarField> phiCN(phi);

    // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = surfaceScalarField::New
        (
            "phiCN",
            cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime()
        );
    }

    if (MULESCorr)
    {
        fvScalarMatrix sEqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh_).fvmDdt(s_)
              : fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(s_)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh_,
                phiCN,
                upwind<scalar>(mesh_, phiCN)
            ).fvmDiv(phiCN, s_)
        );

        sEqn.solve();

        Info<< fieldName_ << " volume fraction = "
            << s_.weightedAverage(mesh_.Vsc()).value()
            << "  Min(" << fieldName_ << ") = " << min(s_).value()
            << "  Max(" << fieldName_ << ") = " << max(s_).value()
            << endl;

        tmp<surfaceScalarField> tsPhiUD(sEqn.flux());
        sPhi_ = tsPhiUD();

        if (applyPrevCorr && tsPhiCorr0_.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                geometricOneField(),
                s_,
                sPhi_,
                tsPhiCorr0_.ref(),
                oneField(),
                zeroField()
            );

            sPhi_ += tsPhiCorr0_();
        }

        // Cache the upwind-flux
        tsPhiCorr0_ = tsPhiUD;
    }

    for (int sCorr=0; sCorr<nCorr; sCorr++)
    {
        // Split operator
        tmp<surfaceScalarField> tsPhiUn
        (
            fvc::flux
            (
                phiCN(),
                (cnCoeff*s_ + (1.0 - cnCoeff)*s_.oldTime())(),
                mesh_.schemes().div(divScheme)
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tsPhiCorr(tsPhiUn() - sPhi_);
            volScalarField s0("s0", s_);

            MULES::correct
            (
                geometricOneField(),
                s_,
                tsPhiUn(),
                tsPhiCorr.ref(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (sCorr == 0)
            {
                sPhi_ += tsPhiCorr();
            }
            else
            {
                s_ = 0.5*s_ + 0.5*s0;
                sPhi_ += 0.5*tsPhiCorr();
            }
        }
        else
        {
            sPhi_ = tsPhiUn;

            MULES::explicitSolve
            (
                geometricOneField(),
                s_,
                phiCN,
                sPhi_,
                oneField(),
                zeroField()
            );
        }
    }

    if (applyPrevCorr && MULESCorr)
    {
        tsPhiCorr0_ = sPhi_ - tsPhiCorr0_;
        tsPhiCorr0_.ref().rename("sPhiCorr0");
    }
    else
    {
        tsPhiCorr0_.clear();
    }

    if
    (
        word(mesh_.schemes().ddt("ddt(s)"))
     == fv::CrankNicolsonDdtScheme<scalar>::typeName
    )
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step s flux
            sPhi_ = (sPhi_ - (1.0 - cnCoeff)*sPhi_.oldTime())/cnCoeff;
        }
    }

    Info<< fieldName_ << "volume fraction = "
        << s_.weightedAverage(mesh_.Vsc()).value()
        << "  Min(" << fieldName_ << ") = " << min(s_).value()
        << "  Max(" << fieldName_ << ") = " << max(s_).value()
        << endl;
}


bool Foam::functionObjects::timeTransport::write()
{
    return true;
}


// ************************************************************************* //

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

#include "avgageTransport.H"
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
    defineTypeNameAndDebug(avgageTransport, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        avgageTransport,
        dictionary
    );
}
}


template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::avgageTransport::diffusivityType,
    3
>::names[] =
{
    "none",
    "constant",
    "viscosity"
};

const Foam::NamedEnum
<
    Foam::functionObjects::avgageTransport::diffusivityType,
    3
> Foam::functionObjects::avgageTransport::diffusivityTypeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::functionObjects::avgageTransport::D() const
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
Foam::functionObjects::avgageTransport::avgageTransport 
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
    fieldName_(dict.lookupOrDefault<word>("field", "cphi")), //fieldname cphi
    diffusivity_(diffusivityType::none),
    D_(0),
    nCorr_(0),
    cphi_
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
    cphiRestart_(false)
{
    read(dict);

    if (mesh_.solution().solversDict().found(fieldName_))
    {
        const dictionary& controls = mesh_.solution().solverDict(fieldName_);

        if (controls.found("nSubCycles"))
        {
            MULES_ = true;

            typeIOobject<surfaceScalarField> cphiPhiHeader
            (
                IOobject::groupName("cphiPhi", cphi_.group()),
                runTime.name(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            );

            cphiRestart_ = cphiPhiHeader.headerOk();

            if (cphiRestart_)
            {
                Info << "Restarting cphi" << endl;
            }

            const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(phiName_);

            // Scalar's volumetric flux
            // this is what phi means https://openfoamwiki.net/index.php/Main_FAQ#What_is_the_field_phi_that_the_solver_is_writing
            tcphiPhi_ = new surfaceScalarField 
            (
                cphiPhiHeader,
                phi*fvc::interpolate(cphi_)
            );

            if (controls.lookupOrDefault<Switch>("MULESCorr", false))
            {
                mesh_.schemes().setFluxRequired(fieldName_);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::avgageTransport::~avgageTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::avgageTransport::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");
    rhoName_ = dict.lookupOrDefault<word>("rho", "rho");
    zetaName_ = dict.lookupOrDefault<word>("zeta", "zeta");
    schemesField_ = dict.lookupOrDefault<word>("schemesField", fieldName_); 
    // find schemesfield by name fieldname, which is cphi

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


Foam::wordList Foam::functionObjects::avgageTransport::fields() const
{
    return wordList{phiName_};
}


bool Foam::functionObjects::avgageTransport::execute()
{
    Info<< type() << " write:" << endl;

    const surfaceScalarField& phi =
        mesh_.lookupObject<surfaceScalarField>(phiName_);
    const volScalarField& zeta =
        mesh_.lookupObject<volScalarField>(zetaName_);

    const word divScheme("div(phi," + schemesField_ + ")"); //div(phi,cphi), divergence of phi*cphi

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

            fvConstraints.constrain(cphi_);
        }
        else
        {
            for (int i=0; i<=nCorr_; i++)
            {
                fvScalarMatrix cphiEqn
                (
                    fvm::ddt(cphi_)
                  + fvm::div(phi, cphi_, divScheme)
                 ==
                    zeta
                );
                //source will be 1 for fluid-age, so can use default code

                if (diffusivity_ != diffusivityType::none)
                {
                    cphiEqn -= fvm::laplacian(D(), cphi_);
                }

                cphiEqn.relax(relaxCoeff);

                fvConstraints.constrain(cphiEqn);

                cphiEqn.solve(schemesField_);

                fvConstraints.constrain(cphi_);
            }
        }
    }
    //for compressible maybe
    else if (phi.dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);
        const volScalarField& zeta =
            mesh_.lookupObject<volScalarField>(zetaName_);

        for (int i=0; i<=nCorr_; i++)
        {
            fvScalarMatrix cphiEqn
            (
                fvm::ddt(rho, cphi_)
              + fvm::div(phi, cphi_, divScheme)
             ==
                zeta
            );
            //same here with upper side

            if (diffusivity_ != diffusivityType::none)
            {
                cphiEqn -= fvm::laplacian(rho*D(), cphi_);
            }

            cphiEqn.relax(relaxCoeff);

            fvConstraints.constrain(cphiEqn);

            cphiEqn.solve(schemesField_);

            fvConstraints.constrain(cphi_);
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


void Foam::functionObjects::avgageTransport::subCycleMULES()
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
            subCycle<volScalarField> cphiSubCycle(cphi_, nSubCycles);
            !(++cphiSubCycle).end();
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
        fvScalarMatrix cphiEqn
        (
            fvm::ddt(cphi_) - fvc::ddt(cphi_)
          - fvm::laplacian(D(), cphi_)
        );

        cphiEqn.solve(controls.subDict("diffusivity"));

        Info<< fieldName_ << " volume fraction = "
            << cphi_.weightedAverage(mesh_.V()).value()
            << "  Min(" << fieldName_ << ") = " << min(cphi_).value()
            << "  Max(" << fieldName_ << ") = " << max(cphi_).value()
            << endl;
    }
}


void Foam::functionObjects::avgageTransport::solveMULES()
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

    surfaceScalarField& cphiPhi_ = tcphiPhi_.ref();

    const word cphiScheme(mesh_.schemes().div(divScheme)[1].wordToken());

    // If a compressive convection scheme is used
    // the interface normal must be cached
    tmp<surfaceScalarField> nHatf;

    if (compressionSchemes.found(cphiScheme))
    {
        const surfaceVectorField gradsf(fvc::interpolate(fvc::grad(cphi_)));

        nHatf = new surfaceScalarField
        (
            IOobject
            (
                "nHatf",
                cphi_.time().name(),
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
                mesh_.schemes().ddt("ddt(cphi)")
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
                cphiRestart_
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
        fvScalarMatrix cphiEqn
        (
            (
                LTS
              ? fv::localEulerDdtScheme<scalar>(mesh_).fvmDdt(cphi_)
              : fv::EulerDdtScheme<scalar>(mesh_).fvmDdt(cphi_)
            )
          + fv::gaussConvectionScheme<scalar>
            (
                mesh_,
                phiCN,
                upwind<scalar>(mesh_, phiCN)
            ).fvmDiv(phiCN, cphi_)
        );

        cphiEqn.solve();

        Info<< fieldName_ << " volume fraction = "
            << cphi_.weightedAverage(mesh_.Vsc()).value()
            << "  Min(" << fieldName_ << ") = " << min(cphi_).value()
            << "  Max(" << fieldName_ << ") = " << max(cphi_).value()
            << endl;

        tmp<surfaceScalarField> tcphiPhiUD(cphiEqn.flux());
        cphiPhi_ = tcphiPhiUD();

        if (applyPrevCorr && tcphiPhiCorr0_.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct
            (
                geometricOneField(),
                cphi_,
                cphiPhi_,
                tcphiPhiCorr0_.ref(),
                oneField(),
                zeroField()
            );

            cphiPhi_ += tcphiPhiCorr0_();
        }

        // Cache the upwind-flux
        tcphiPhiCorr0_ = tcphiPhiUD;
    }

    for (int cphiCorr=0; cphiCorr<nCorr; cphiCorr++)
    {
        // Split operator
        tmp<surfaceScalarField> tcphiPhiUn
        (
            fvc::flux
            (
                phiCN(),
                (cnCoeff*cphi_ + (1.0 - cnCoeff)*cphi_.oldTime())(),
                mesh_.schemes().div(divScheme)
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> tcphiPhiCorr(tcphiPhiUn() - cphiPhi_);
            volScalarField cphi0("cphi0", cphi_);

            MULES::correct
            (
                geometricOneField(),
                cphi_,
                tcphiPhiUn(),
                tcphiPhiCorr.ref(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (cphiCorr == 0)
            {
                cphiPhi_ += tcphiPhiCorr();
            }
            else
            {
                cphi_ = 0.5*cphi_ + 0.5*cphi0;
                cphiPhi_ += 0.5*tcphiPhiCorr();
            }
        }
        else
        {
            cphiPhi_ = tcphiPhiUn;

            MULES::explicitSolve
            (
                geometricOneField(),
                cphi_,
                phiCN,
                cphiPhi_,
                oneField(),
                zeroField()
            );
        }
    }

    if (applyPrevCorr && MULESCorr)
    {
        tcphiPhiCorr0_ = cphiPhi_ - tcphiPhiCorr0_;
        tcphiPhiCorr0_.ref().rename("cphiPhiCorr0");
    }
    else
    {
        tcphiPhiCorr0_.clear();
    }

    if
    (
        word(mesh_.schemes().ddt("ddt(cphi)"))
     == fv::CrankNicolsonDdtScheme<scalar>::typeName
    )
    {
        if (ocCoeff > 0)
        {
            // Calculate the end-of-time-step cphi flux
            cphiPhi_ = (cphiPhi_ - (1.0 - cnCoeff)*cphiPhi_.oldTime())/cnCoeff;
        }
    }

    Info<< fieldName_ << "volume fraction = "
        << cphi_.weightedAverage(mesh_.Vsc()).value()
        << "  Min(" << fieldName_ << ") = " << min(cphi_).value()
        << "  Max(" << fieldName_ << ") = " << max(cphi_).value()
        << endl;
}


bool Foam::functionObjects::avgageTransport::write()
{
    return true;
}


// ************************************************************************* //

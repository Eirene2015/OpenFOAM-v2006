/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "heatTransferCoeffModel.H"
#include "fvMesh.H"
#include "fluidThermo.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatTransferCoeffModel, 0);
    defineRunTimeSelectionTable(heatTransferCoeffModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::FieldField<Foam::Field, Foam::scalar>>
Foam::heatTransferCoeffModel::q() const
{
    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);
    const volScalarField::Boundary& Tbf = T.boundaryField();

    auto tq = tmp<FieldField<Field, scalar>>::New(Tbf.size());
    auto& q = tq.ref();

    forAll(q, patchi)
    {
        q.set(patchi, new Field<scalar>(Tbf[patchi].size(), Zero));
    }

    typedef compressible::turbulenceModel cmpTurbModel;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel& turb =
            mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        const volScalarField& he = turb.transport().he();
        const volScalarField::Boundary& hebf = he.boundaryField();

        const volScalarField alphaEff(turb.alphaEff());
        const volScalarField::Boundary& alphaEffbf = alphaEff.boundaryField();

        for (const label patchi : patchSet_)
        {
            q[patchi] = alphaEffbf[patchi]*hebf[patchi].snGrad();
        }
    }
    else if (mesh_.foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const fluidThermo& thermo =
            mesh_.lookupObject<fluidThermo>(fluidThermo::dictName);

        const volScalarField& he = thermo.he();
        const volScalarField::Boundary& hebf = he.boundaryField();

        const volScalarField& alpha(thermo.alpha());
        const volScalarField::Boundary& alphabf = alpha.boundaryField();

        for (const label patchi : patchSet_)
        {
            q[patchi] = alphabf[patchi]*hebf[patchi].snGrad();
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find a valid thermo model to evaluate q" << nl
            << exit(FatalError);
    }

    // Add radiative heat flux contribution if present

    const volScalarField* qrPtr = mesh_.cfindObject<volScalarField>(qrName_);

    if (qrPtr)
    {
        const volScalarField::Boundary& qrbf = qrPtr->boundaryField();

        for (const label patchi : patchSet_)
        {
            q[patchi] += qrbf[patchi];
        }
    }

    return tq;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferCoeffModel::heatTransferCoeffModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& TName
)
:
    mesh_(mesh),
    patchSet_(),
    TName_(TName),
    qrName_("qr")
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::heatTransferCoeffModel::read(const dictionary& dict)
{
    patchSet_ = mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches"));

    dict.readIfPresent("qr", qrName_);

    return true;
}


bool Foam::heatTransferCoeffModel::calc
(
    volScalarField& result,
    const FieldField<Field, scalar>& q
)
{
    htc(result, q);

    return true;
}


// ************************************************************************* //

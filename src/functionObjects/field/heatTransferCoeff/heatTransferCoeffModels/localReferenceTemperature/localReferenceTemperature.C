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

#include "localReferenceTemperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{
    defineTypeNameAndDebug(localReferenceTemperature, 0);
    addToRunTimeSelectionTable
    (
        heatTransferCoeffModel,
        localReferenceTemperature,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferCoeffModels::localReferenceTemperature::
localReferenceTemperature
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& TName
)
:
    heatTransferCoeffModel(dict, mesh, TName)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::heatTransferCoeffModels::localReferenceTemperature::read
(
    const dictionary& dict
)
{
    return heatTransferCoeffModel::read(dict);
}


void Foam::heatTransferCoeffModels::localReferenceTemperature::htc
(
    volScalarField& htc,
    const FieldField<Field, scalar>& q
)
{
    ///const FieldField<Field, scalar> qBf(q());
    const volScalarField& T = mesh_.lookupObject<volScalarField>(TName_);
    const volScalarField::Boundary& Tbf = T.boundaryField();
    const scalar eps = ROOTVSMALL;

    volScalarField::Boundary& htcBf = htc.boundaryFieldRef();

    for (const label patchi : patchSet_)
    {
        const scalarField Tc(Tbf[patchi].patchInternalField());
        htcBf[patchi] = q[patchi]/(Tc - Tbf[patchi] + eps);
    }
}


// ************************************************************************* //

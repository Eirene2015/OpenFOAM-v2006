/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

Description
    Finite area edge-based patch fields

\*---------------------------------------------------------------------------*/

#include "faePatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeFaePatchField(faePatchTypeField)                                  \
                                                                              \
defineNamedTemplateTypeNameAndDebug(faePatchTypeField, 0);                    \
template<>                                                                    \
int                                                                           \
faePatchTypeField::disallowGenericFaePatchField                               \
(                                                                             \
    debug::debugSwitch("disallowGenericFaPatchField", 0)                      \
);                                                                            \
defineTemplateRunTimeSelectionTable(faePatchTypeField, patch);                \
defineTemplateRunTimeSelectionTable(faePatchTypeField, patchMapper);          \
defineTemplateRunTimeSelectionTable(faePatchTypeField, dictionary);

makeFaePatchField(faePatchScalarField)
makeFaePatchField(faePatchVectorField)
makeFaePatchField(faePatchSphericalTensorField)
makeFaePatchField(faePatchSymmTensorField)
makeFaePatchField(faePatchTensorField)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

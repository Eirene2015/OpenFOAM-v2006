/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR YEAR AUTHOR,AFFILIATION
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

#include "codedPatchFunction1Template.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"

//{{{ begin codeInclude
${codeInclude}
//}}} end codeInclude


namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = ${SHA1sum}
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void ${typeName}_${SHA1sum}(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}


namespace PatchFunction1Types
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//makePatchFunction1(${typeName}PatchFunction1${FieldType});
defineTypeNameAndDebug
(
    ${typeName}PatchFunction1${FieldType},
    0
);
PatchFunction1<${TemplateType}>::adddictionaryConstructorToTable
    <${typeName}PatchFunction1${FieldType}>
    add${typeName}PatchFunction1${FieldType}ConstructorToTable_;



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}PatchFunction1${FieldType}::
${typeName}PatchFunction1${FieldType}
(
    const polyPatch& pp,
    const word& type,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<${TemplateType}>(pp, entryName, dict, faceValues)
{
    if (${verbose:-false})
    {
        printMessage("Construct ${typeName} from components");
    }
}


${typeName}PatchFunction1${FieldType}::
${typeName}PatchFunction1${FieldType}
(
    const ${typeName}PatchFunction1${FieldType}& rhs
)
:
    PatchFunction1<${TemplateType}>(rhs)
{}


${typeName}PatchFunction1${FieldType}::
${typeName}PatchFunction1${FieldType}
(
    const ${typeName}PatchFunction1${FieldType}& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<${TemplateType}>(rhs, pp)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<Field<${TemplateType}>>
${typeName}PatchFunction1${FieldType}::value
(
    const scalar x
) const
{
//{{{ begin code
    ${code}
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam


// ************************************************************************* //

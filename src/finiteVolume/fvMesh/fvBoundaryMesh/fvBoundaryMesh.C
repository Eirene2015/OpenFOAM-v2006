/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "fvBoundaryMesh.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvBoundaryMesh::addPatches(const polyBoundaryMesh& basicBdry)
{
    setSize(basicBdry.size());

    // Set boundary patches
    fvPatchList& Patches = *this;

    forAll(Patches, patchi)
    {
        Patches.set(patchi, fvPatch::New(basicBdry[patchi], *this));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvBoundaryMesh::fvBoundaryMesh
(
    const fvMesh& m
)
:
    fvPatchList(),
    mesh_(m)
{}


Foam::fvBoundaryMesh::fvBoundaryMesh
(
    const fvMesh& m,
    const polyBoundaryMesh& basicBdry
)
:
    fvPatchList(basicBdry.size()),
    mesh_(m)
{
    addPatches(basicBdry);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fvBoundaryMesh::indices
(
    const keyType& key,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(key, useGroups);
}


Foam::label Foam::fvBoundaryMesh::findPatchID(const word& patchName) const
{
    const fvPatchList& patches = *this;

    forAll(patches, patchi)
    {
        if (patches[patchi].name() == patchName)
        {
            return patchi;
        }
    }

    // Not found, return -1
    return -1;
}


void Foam::fvBoundaryMesh::movePoints()
{
    fvPatchList& patches = *this;

    for (fvPatch& p : patches)
    {
        p.initMovePoints();
    }

    for (fvPatch& p : patches)
    {
        p.movePoints();
    }
}


Foam::lduInterfacePtrsList Foam::fvBoundaryMesh::interfaces() const
{
    lduInterfacePtrsList interfaces(size());

    forAll(interfaces, patchi)
    {
        if (isA<lduInterface>(this->operator[](patchi)))
        {
            interfaces.set
            (
                patchi,
               &refCast<const lduInterface>(this->operator[](patchi))
            );
        }
    }

    return interfaces;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fvBoundaryMesh::readUpdate(const polyBoundaryMesh& basicBdry)
{
    clear();
    addPatches(basicBdry);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::fvPatch& Foam::fvBoundaryMesh::operator[]
(
    const word& patchName
) const
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << abort(FatalError);
    }

    return operator[](patchi);
}


Foam::fvPatch& Foam::fvBoundaryMesh::operator[]
(
    const word& patchName
)
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << abort(FatalError);
    }

    return operator[](patchi);
}


// ************************************************************************* //

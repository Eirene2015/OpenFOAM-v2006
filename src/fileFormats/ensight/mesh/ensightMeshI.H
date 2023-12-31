/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::polyMesh& Foam::ensightMesh::mesh() const
{
    return mesh_;
}


inline const Foam::ensightMesh::options& Foam::ensightMesh::option() const
{
    return *options_;
}


inline const Foam::Map<Foam::ensightCells>&
Foam::ensightMesh::cellZoneParts() const
{
    return cellZoneParts_;
}


inline const Foam::Map<Foam::ensightFaces>&
Foam::ensightMesh::faceZoneParts() const
{
    return faceZoneParts_;
}


inline const Foam::Map<Foam::ensightFaces>&
Foam::ensightMesh::boundaryParts() const
{
    return boundaryParts_;
}


inline bool Foam::ensightMesh::expire()
{
    clear();

    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


inline bool Foam::ensightMesh::needsUpdate() const
{
    return needsUpdate_;
}


inline bool Foam::ensightMesh::empty() const
{
    return
    (
        cellZoneParts_.empty()
     && faceZoneParts_.empty()
     && boundaryParts_.empty()
    );
}


inline Foam::label Foam::ensightMesh::size() const
{
    return
    (
        cellZoneParts_.size()
      + faceZoneParts_.size()
      + boundaryParts_.size()
    );
}


inline void Foam::ensightMesh::write
(
    autoPtr<ensightGeoFile>& os,
    bool parallel
) const
{
    write(os.ref(), parallel);
}


// ************************************************************************* //

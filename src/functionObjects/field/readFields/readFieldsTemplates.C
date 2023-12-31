/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "readFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::readFields::loadField(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal IntVolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    /// typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    if (foundObject<VolFieldType>(fieldName))
    {
        DebugInfo
            << "readFields : " << VolFieldType::typeName
            << " " << fieldName << " already in database"
            << endl;
    }
    else if (foundObject<IntVolFieldType>(fieldName))
    {
        DebugInfo
            << "readFields : " << IntVolFieldType::typeName
            << " " << fieldName << " already in database"
            << endl;
    }
    else if (foundObject<SurfaceFieldType>(fieldName))
    {
        DebugInfo
            << "readFields: " << SurfaceFieldType::typeName
            << " " << fieldName << " already exists in database"
            << " already in database" << endl;
    }
    /// else if (foundObject<SurfFieldType>(fieldName))
    /// {
    ///     DebugInfo
    ///         << "readFields: " << SurfFieldType::typeName
    ///         << " " << fieldName << " already exists in database"
    ///         << " already in database" << endl;
    /// }
    else
    {
        IOobject fieldHeader
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (fieldHeader.typeHeaderOk<VolFieldType>(true, true, false))
        {
            // Store field on mesh database
            Log << "    Reading " << fieldName << endl;
            auto* fldPtr(new VolFieldType(fieldHeader, mesh_));
            mesh_.objectRegistry::store(fldPtr);
            return true;
        }
        else if (fieldHeader.typeHeaderOk<IntVolFieldType>(true, true, false))
        {
            // Store field on mesh database
            Log << "    Reading " << fieldName << endl;
            auto* fldPtr(new IntVolFieldType(fieldHeader, mesh_));
            mesh_.objectRegistry::store(fldPtr);
            return true;
        }
        else if (fieldHeader.typeHeaderOk<SurfaceFieldType>(true, true, false))
        {
            // Store field on mesh database
            Log << "    Reading " << fieldName << endl;
            auto* fldPtr(new SurfaceFieldType(fieldHeader, mesh_));
            mesh_.objectRegistry::store(fldPtr);
            return true;
        }
        /// else if (fieldHeader.typeHeaderOk<SurfFieldType>(true, true, false))
        /// {
        ///     const surfMesh* surfptr = isA<surfMesh>(obr());
        ///     if (surfptr)
        ///     {
        ///         const surfMesh& s = surfptr;
        ///
        ///         // Store field on surfMesh database
        ///         Log << "    Reading " << fieldName << endl;
        ///         auto* fldPtr(new SurfFieldType(fieldHeader, s));
        ///         s.store(fldPtr);
        ///         return true;
        ///     }
        /// }
    }

    return false;
}


// ************************************************************************* //

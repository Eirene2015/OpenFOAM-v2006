/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

Class
    Foam::areaWrite

Description
    Write finite area mesh/fields to standard output formats.

    Example of function object specification:

    \verbatim
    surfaces
    {
        type    areaWrite;
        libs    (utilityFunctionObjects);

        // Write at same frequency as fields
        writeControl    outputTime;
        writeInterval   1;

        // Fields to be sampled
        fields          (p U);

        // Output surface format
        surfaceFormat   vtk;

        formatOptions
        {
            vtk
            {
                precision  10;
            }
        }
    }
    \endverbatim

    Entries:
    \table
        Property | Description                              | Required | Default
        type    | Type-name: \c areaWrite                   | yes |
        region  | name for a single region                  | no  | region0
        area    | select a single area                      | no  |
        areas   | wordRe list of multiple areas             | no  |
        fields  | wordRe list of fields                     | yes |
        surfaceFormat | output surface format               | yes |
        formatOptions | dictionary of format options        | no  |
    \endtable

SourceFiles
    areaWrite.C
    areaWriteTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef areaWrite_H
#define areaWrite_H

#include "fvMeshFunctionObject.H"
#include "polyMesh.H"
#include "areaFieldsFwd.H"
#include "surfaceWriter.H"
#include "HashPtrTable.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;
class polySurface;

/*---------------------------------------------------------------------------*\
                          Class areaWrite Declaration
\*---------------------------------------------------------------------------*/

class areaWrite
:
    public functionObjects::fvMeshFunctionObject
{
    // Static Data Members

        //- Tolerance for merging points (fraction of mesh bounding box)
        static scalar mergeTol_;


    // Private Data

        //- Load fields from files (not from objectRegistry)
        const bool loadFromFiles_;

        //- Output verbosity
        bool verbose_;

        //- Output path
        fileName outputPath_;


    // Read from dictionary

        //- Names of areas to select
        wordRes selectAreas_;

        //- Names of fields to write
        wordRes fieldSelection_;

        //- Pointers to the requested mesh regions
        HashTable<const faMesh*> meshes_;

        //- Hold intermediate surfaces.
        //  The faMesh has an indirect face list but we require real ones.
        autoPtr<objectRegistry> surfaces_;


    // Output control

        //- Surface writers (one per surface)
        HashPtrTable<surfaceWriter> writers_;


    // Private Member Functions

        //- Write fieldName on surface and on outputDir path
        template<class Type>
        void writeSurface
        (
            surfaceWriter& writer,
            const Field<Type>& values,
            const word& fieldName
        );

        //- Write all applicable fields
        template<class GeoField>
        void performAction
        (
            surfaceWriter& writer,
            const faMesh& areaMesh,
            const IOobjectList& objects
        );

        //- Mark intermediate surfaces and writers as needing an update.
        void expire();

        //- No copy construct
        areaWrite(const areaWrite&) = delete;

        //- No copy assignment
        void operator=(const areaWrite&) = delete;


public:

    //- Runtime type information
    TypeName("areaWrite");


    // Constructors

        //- Construct from Time and dictionary
        areaWrite
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct for given objectRegistry and dictionary
        //  allow the possibility to load fields from files
        areaWrite
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~areaWrite() = default;


    // Member Functions

        //- Set verbosity level
        void verbose(const bool verbosity = true);

        //- Read the areaWrite dictionary
        virtual bool read(const dictionary& dict);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Sample and write
        virtual bool write();

        //- Update for changes of mesh - expires the surfaces
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion - expires the surfaces
        virtual void movePoints(const polyMesh& mesh);

        //- Update for changes of mesh due to readUpdate - expires the surfaces
        virtual void readUpdate(const polyMesh::readUpdateState state);

        //- Get merge tolerance
        static scalar mergeTol();

        //- Set merge tolerance and return old value
        static scalar mergeTol(const scalar tol);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "areaWriteTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

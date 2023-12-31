/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::ensightMesh

Description
    Encapsulation of volume meshes for writing in ensight format.
    It manages cellZones, facesZone, patches.

    When cellZones are present (and not disabled), the cells are grouped
    in parts according to the zone.
    Any remaining \em unzoned cells are placed into the "internalMesh" part,
    which is always part 0. If cellZones are missing or disabled,
    all cells are placed into the "internalMesh" part.

    If one or more cellZones are explicitly requested, all other cells
    (including any unzoned cells) are ignored.

    The converted patch faces are restricted by the volume mesh coverage.
    Except when the entire internal mesh has been explicitly suppressed.

Note
    The internal data management uses a Map for cellZones, faceZones and
    patches. The internalMesh is treated as cellZone with a special index.

    Since the patches are subsetted by the internal mesh coverage,
    they are treated as indirect patches rather than regular poly patches.

SourceFiles
    ensightMesh.C
    ensightMeshI.H
    ensightMeshOptions.C

\*---------------------------------------------------------------------------*/

#ifndef ensightMesh_H
#define ensightMesh_H

#include "Map.H"
#include "ensightCells.H"
#include "ensightFaces.H"
#include "wordRes.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class ensightGeoFile;
class ensightMesh;

/*---------------------------------------------------------------------------*\
                         Class ensightMesh Declaration
\*---------------------------------------------------------------------------*/

class ensightMesh
{
public:

    // Forward Declarations
    class options;

    //- The zone-id for internal mesh or unzoned cells.
    static const label internalZone;


private:

    // Private Data

        //- Writer options
        const std::unique_ptr<options> options_;

        //- Reference to the OpenFOAM mesh
        const polyMesh& mesh_;

        //- Volume elements per cellZone, lookup by zone index.
        //  The zone -1 is reserved for internal mesh (unzoned cells)
        Map<ensightCells> cellZoneParts_;

        //- Face elements per faceZone, lookup by zone index.
        Map<ensightFaces> faceZoneParts_;

        //- Face elements per selected patch, lookup by patch index
        Map<ensightFaces> boundaryParts_;

        //- Track if it needs an update
        mutable bool needsUpdate_;


    // Private Member Functions

        //- Clear all storage
        void clear();

        //- Enforce consistent index/part numbering
        void renumber();

        //- No copy construct
        ensightMesh(const ensightMesh&) = delete;

        //- No copy assignment
        void operator=(const ensightMesh&) = delete;


public:

    // Constructors

        //- Construct from mesh with all default options
        explicit ensightMesh(const polyMesh& mesh);

        //- Construct from components
        ensightMesh(const polyMesh& mesh, const options& opts);


    // Member Functions

    // Access

        //- Reference to the underlying polyMesh
        inline const polyMesh& mesh() const;

        //- Reference to the writer/mesh options
        inline const ensightMesh::options& option() const;

        //- Face elements per selected patch, lookup by patch index
        //  Process in sorted order.
        //  May require special treatment for zone -1 (internal).
        inline const Map<ensightCells>& cellZoneParts() const;

        //- Face elements per faceZone, lookup by zone index.
        //  Process in sorted order.
        inline const Map<ensightFaces>& faceZoneParts() const;

        //- Face elements per selected patch, lookup by patch index
        //  Process in sorted order.
        inline const Map<ensightFaces>& boundaryParts() const;


    // Sizing Information

        //- Any parts?
        inline bool empty() const;

        //- Number of parts
        inline label size() const;


    // Other

        //- Does the content need an update?
        inline bool needsUpdate() const;

        //- Mark as needing an update.
        //  May also free up unneeded data.
        //  Return false if already marked as expired.
        inline bool expire();

        //- Update for new mesh
        void correct();


    // Output

        //- Write geometry to file. Normally in parallel
        inline void write
        (
            autoPtr<ensightGeoFile>& os,
            bool parallel = Pstream::parRun()
        ) const;

        //- Write geometry to file. Normally in parallel
        void write
        (
            ensightGeoFile& os,
            bool parallel = Pstream::parRun()
        ) const;
};


/*---------------------------------------------------------------------------*\
                    Class ensightMesh::options Declaration
\*---------------------------------------------------------------------------*/

//- Configuration options for the ensightMesh
class ensightMesh::options
{
    // Private Data

        //- Create in 'expired' mode
        bool lazy_;

        //- Use the internal mesh
        bool internal_;

        //- Use the boundary mesh
        bool boundary_;

        //- Handle cellZones (if internal_ is true)
        bool cellZones_;

        //- Selected patches only
        wordRes patchInclude_;

        //- Deselected patches
        wordRes patchExclude_;

        //- Selected cellZones
        wordRes cellZoneInclude_;

        //- Selected faceZones
        wordRes faceZoneInclude_;


public:

    // Constructors

        //- Default construct. Non-lazy with internal/boundary/cellZones.
        options();


    // Member Functions

    // Access

        //- Lazy creation? (ie, ensightMesh starts as needsUpdate)
        bool lazy() const;

        //- Using internal?
        bool useInternalMesh() const;

        //- Using boundary?
        bool useBoundaryMesh() const;

        //- Using faceZones?
        bool useFaceZones() const;

        //- Using cellZones?
        bool useCellZones() const;

        //- Selection of patches. Empty if unspecified.
        const wordRes& patchSelection() const;

        //- Selection of black listed patches. Empty if unspecified.
        const wordRes& patchExclude() const;

        //- Selection of faceZones. Empty if unspecified.
        const wordRes& faceZoneSelection() const;

        //- Selection of faceZones. Empty if unspecified.
        const wordRes& cellZoneSelection() const;


    // Edit

        //- Reset to defaults
        void reset();

        //- Lazy creation - ensightMesh starts as needsUpdate
        void lazy(bool beLazy);

        //- Alter the useBoundaryMesh state
        void useInternalMesh(bool on);

        //- Alter the useBoundaryMesh state
        void useBoundaryMesh(bool on);

        //- Alter the useCellZones state
        void useCellZones(bool on);

        //- Define patch selection matcher
        void patchSelection(const UList<wordRe>& patterns);

        //- Define patch selection matcher
        void patchSelection(List<wordRe>&& patterns);

        //- Define patch selection blacklist
        void patchExclude(const UList<wordRe>& patterns);

        //- Define patch selection blacklist
        void patchExclude(List<wordRe>&& patterns);

        //- Define faceZone selection matcher
        void faceZoneSelection(const UList<wordRe>& patterns);

        //- Define faceZone selection matcher
        void faceZoneSelection(List<wordRe>&& patterns);

        //- Define cellZone selection matcher
        void cellZoneSelection(const UList<wordRe>& patterns);

        //- Define cellZone selection matcher
        void cellZoneSelection(List<wordRe>&& patterns);


    // Output

        //- Report values
        void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensightMeshI.H"

#endif

// ************************************************************************* //

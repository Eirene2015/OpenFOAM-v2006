/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::isoSurfaceTopo

Description
    Marching tet iso surface algorithm with optional filtering to keep only
    points originating from mesh edges.

SourceFiles
    isoSurfaceTopo.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurfaceTopo_H
#define isoSurfaceTopo_H

#include "labelPair.H"
#include "bitSet.H"
#include "edgeList.H"
#include "isoSurfaceBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class tetMatcher;

/*---------------------------------------------------------------------------*\
                       Class isoSurfaceTopo Declaration
\*---------------------------------------------------------------------------*/

class isoSurfaceTopo
:
    public isoSurfaceBase
{
    // Private Data

        enum cellCutType
        {
            NOTCUT,     //!< Not cut
            SPHERE,     //!< All edges to cell centre cut
            CUT         //!< Normal cut
        };


        //- Reference to mesh
        const polyMesh& mesh_;

        const scalarField& cVals_;

        const scalarField& pVals_;

        //- Optional cells to ignore
        const bitSet& ignoreCells_;

        //- Optional boundary faces to ignore. Used to exclude cyclicACMI
        //  (since duplicate faces)
        bitSet ignoreBoundaryFaces_;

        //- Corrected version of tetBasePtIs
        labelList tetBasePtIs_;

        //- Per point: originating mesh vertex/cc. See encoding above
        edgeList pointToVerts_;

        //- For every point the originating face in mesh
        labelList pointToFace_;


    // Private Member Functions

        scalar minTetQ
        (
            const label facei,
            const label faceBasePtI
        ) const;

        void fixTetBasePtIs();

        //- Determine whether cell is cut
        cellCutType calcCutType
        (
            const bool isTet,
            const label
        ) const;

        //- Determine for all mesh whether cell is cut
        //  \return number of cells cut
        label calcCutTypes
        (
            tetMatcher& tet,
            List<cellCutType>& cellCutTypes
        );

        //- Generate single point on edge
        label generatePoint
        (
            const label facei,
            const bool edgeIsDiag,
            const edge& vertices,

            DynamicList<edge>& pointToVerts,
            DynamicList<label>& pointToFace,
            DynamicList<bool>& pointFromDiag,
            EdgeMap<label>& vertsToPoint
        ) const;

        //- Generate triangles from tet
        void generateTriPoints
        (
            const label facei,
            const FixedList<scalar, 4>& s,
            const FixedList<point, 4>& p,
            const FixedList<label, 4>& pIndex,
            const FixedList<bool, 6>& edgeIsDiag,

            DynamicList<edge>& pointToVerts,
            DynamicList<label>& pointToFace,
            DynamicList<bool>& pointFromDiag,

            EdgeMap<label>& vertsToPoint,
            DynamicList<label>& verts
        ) const;

        //- Generate triangles from cell
        void generateTriPoints
        (
            const label celli,
            const bool isTet,

            DynamicList<edge>& pointToVerts,
            DynamicList<label>& pointToFace,
            DynamicList<bool>& pointFromDiag,

            EdgeMap<label>& vertsToPoint,
            DynamicList<label>& verts,
            DynamicList<label>& faceLabels
        ) const;


        // Simplification

            void triangulateOutside
            (
                const bool filterDiag,
                const primitivePatch& pp,
                const boolList& pointFromDiag,
                const labelList& pointToFace,
                const label cellID,

                DynamicList<face>& compactFaces,
                DynamicList<label>& compactCellIDs
            ) const;

            MeshStorage removeInsidePoints
            (
                const bool filterDiag,
                const MeshStorage& s,
                const boolList& pointFromDiag,
                const labelList& pointToFace,
                const labelList& start,              // Per cell:starting tri
                DynamicList<label>& pointCompactMap, // Per point the original
                DynamicList<label>& compactCellIDs   // Per face the cellID
            ) const;


public:

    //- Filtering type
    using isoSurfaceBase::filterType;


    //- Runtime type information
    TypeName("isoSurfaceTopo");


    // Constructors

        //- Construct from dictionary
        isoSurfaceTopo
        (
            const polyMesh& mesh,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso,
            const filterType filter = filterType::DIAGCELL,
            const boundBox& bounds = boundBox::invertedBox,
            const bitSet& ignoreCells = bitSet()
        );


    // Member Functions

        //- For every point originating face (pyramid) in mesh
        const labelList& pointToFace() const
        {
            return pointToFace_;
        }

        //- Per point: originating mesh vertex/cc. See encoding above
        const edgeList& pointToVerts() const
        {
            return pointToVerts_;
        }

        //- Interpolates cCoords,pCoords.
        template<class Type>
        tmp<Field<Type>> interpolate
        (
            const Field<Type>& cCoords,
            const Field<Type>& pCoords
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "isoSurfaceTopoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

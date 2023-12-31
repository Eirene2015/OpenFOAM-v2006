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

Class
    Foam::faMesh

Description
    Finite area mesh.  Used for 2-D non-Euclidian finite area method.

SourceFiles
    faMesh.C
    faMeshDemandDrivenData.C

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef faMesh_H
#define faMesh_H

#include "GeoMesh.H"
#include "MeshObject.H"
#include "polyMesh.H"
#include "lduMesh.H"
#include "faBoundaryMesh.H"
#include "edgeList.H"
#include "faceList.H"
#include "primitiveFieldsFwd.H"
#include "DimensionedField.H"
#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "indirectPrimitivePatch.H"
#include "edgeInterpolation.H"
#include "labelIOList.H"
#include "FieldFields.H"
#include "faGlobalMeshData.H"
#include "faSchemes.H"
#include "faSolution.H"
#include "data.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class faMeshLduAddressing;
class faMeshMapper;

/*---------------------------------------------------------------------------*\
                           Class faMesh Declaration
\*---------------------------------------------------------------------------*/

class faMesh
:
    public GeoMesh<polyMesh>,
    public MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>,
    public lduMesh,
    public edgeInterpolation,
    public faSchemes,
    public faSolution,
    public data
{
    // Private data

        //- Face labels
        labelIOList faceLabels_;

        //- Boundary mesh
        faBoundaryMesh boundary_;


        // Primitive mesh data

            //- Edges, addressing into local point list
            edgeList edges_;

            //- Edge owner
            labelList edgeOwner_;

            //- Edge neighbour
            labelList edgeNeighbour_;


        // Primitive size data

            //- Number of points
            mutable label nPoints_;

            //- Number of edges
            mutable label nEdges_;

            //- Number of internal edges
            mutable label nInternalEdges_;

            //- Number of faces
            mutable label nFaces_;


        // Communication support

            //- Communicator used for parallel communication
            label comm_;


    // Demand-driven data

        //- Primitive patch
        mutable indirectPrimitivePatch* patchPtr_;

        //- Ldu addressing data
        mutable faMeshLduAddressing* lduPtr_;

        //- Current time index for motion
        //  Note.  The whole mechanism will be replaced once the
        //  dimensionedField is created and the dimensionedField
        //  will take care of the old-time levels.
        mutable label curTimeIndex_;

        //- Face areas
        mutable DimensionedField<scalar, areaMesh>* SPtr_;

        //- Face areas old time level
        mutable DimensionedField<scalar, areaMesh>* S0Ptr_;

        //- Face areas old-old time level
        mutable DimensionedField<scalar, areaMesh>* S00Ptr_;

        //- Patch starts in the edge list
        mutable labelList* patchStartsPtr_;

        //- Edge length vectors
        mutable edgeVectorField* LePtr_;

        //- Mag edge length vectors
        mutable edgeScalarField* magLePtr_;

        //- Face centres
        mutable areaVectorField* centresPtr_;

        //- Edge centres
        mutable edgeVectorField* edgeCentresPtr_;

        //- Face area normals
        mutable areaVectorField* faceAreaNormalsPtr_;

        //- Edge area normals
        mutable edgeVectorField* edgeAreaNormalsPtr_;

        //- Edge area normals
        mutable vectorField* pointAreaNormalsPtr_;

        //- Face curvatures
        mutable areaScalarField* faceCurvaturesPtr_;

        //- Edge transformation tensors
        mutable FieldField<Field, tensor>* edgeTransformTensorsPtr_;

        //- Whether point normals must be corrected for a patch
        mutable boolList* correctPatchPointNormalsPtr_;


        // Other mesh-related data

            //- Parallel info
            mutable faGlobalMeshData* globalMeshDataPtr_;


    // Static Private Data

        //- Use quadrics fit
        static const int quadricsFit_;


    // Private Member Functions

        //- No copy construct
        faMesh(const faMesh&) = delete;

        //- No copy assignment
        void operator=(const faMesh&) = delete;


        //- Set primitive mesh data
        void setPrimitiveMeshData();


        // Private member functions to calculate demand driven data

            //- Calculate ldu addressing
            void calcLduAddressing() const;

            //- Calculate patch starts in the edge list
            void calcPatchStarts() const;

            //- Calculate edge lengths
            void calcLe() const;

            //- Calculate mag edge lengths
            void calcMagLe() const;

            //- Calculate face centres
            void calcAreaCentres() const;

            //- Calculate edge centres
            void calcEdgeCentres() const;

            //- Calculate face areas
            void calcS() const;

            //- Calculate face area normals
            void calcFaceAreaNormals() const;

            //- Calculate edge area normals
            void calcEdgeAreaNormals() const;

            //- Calculate point area normals
            void calcPointAreaNormals() const;

            //- Calculate point area normals by quadrics fit
            void calcPointAreaNormalsByQuadricsFit() const;

            //- Calculate face curvatures
            void calcFaceCurvatures() const;

            //- Calculate edge transformation tensors
            void calcEdgeTransformTensors() const;

            //- Clear geometry but not the face areas
            void clearGeomNotAreas() const;

            //- Clear geometry
            void clearGeom() const;

            //- Clear addressing
            void clearAddressing() const;

            //- Clear demand-driven data
            void clearOut() const;

public:

    // Public typedefs

        typedef faMesh Mesh;
        typedef faBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("faMesh");


    //- Return the mesh sub-directory name (usually "faMesh")
    static word meshSubDir;


    // Constructors

        //- Construct from polyMesh
        explicit faMesh(const polyMesh& m);

        //- Construct from components without boundary.
        //  Boundary is added using addFaPatches() member function
        faMesh
        (
            const polyMesh& m,
            const labelList& faceLabels
        );

        //- Construct from finite area mesh definition file
        faMesh
        (
            const polyMesh& m,
            const fileName& defFile
        );

        //- Construct from polyPatch
        faMesh
        (
            const polyMesh& m,
            const label polyPatchID
        );


    //- Destructor
    virtual ~faMesh();


    // Member Functions

        // Helpers

            //- Add boundary patches. Constructor helper
            void addFaPatches(const List<faPatch*> &);


        // Database

            //- Return access to polyMesh
            const polyMesh& mesh() const
            {
                return
                    MeshObject
                    <
                        polyMesh,
                        Foam::UpdateableMeshObject,
                        faMesh
                    >::mesh();
            }

            //- Return the local mesh directory (dbDir()/meshSubDir)
            fileName meshDir() const;

            //- Return reference to time
            const Time& time() const;

            //- Return the current instance directory for points
            //  Used in the construction of geometric mesh data dependent
            //  on points
            const fileName& pointsInstance() const;

            //- Return the current instance directory for faces
            const fileName& facesInstance() const;


            // Mesh size parameters

                inline label nPoints() const
                {
                    return nPoints_;
                }

                inline label nEdges() const
                {
                    return nEdges_;
                }

                inline label nInternalEdges() const
                {
                    return nInternalEdges_;
                }

                inline label nFaces() const
                {
                    return nFaces_;
                }


            // Primitive mesh data

                //- Return mesh points
                const pointField& points() const;

                //- Return edges
                const edgeList& edges() const;

                //- Return faces
                const faceList& faces() const;

                //- Edge owner addressing
                inline const labelList& edgeOwner() const
                {
                    return edgeOwner_;
                }

                //- Edge neighbour addressing
                inline const labelList& edgeNeighbour() const
                {
                    return edgeNeighbour_;
                }


            // Communication support

                //- Return communicator used for parallel communication
                label comm() const;

                //- Return communicator used for parallel communication
                label& comm();


        // Access

            //- Return true if thisDb() is a valid DB
            virtual bool hasDb() const;

            //- Return reference to the mesh database
            virtual const objectRegistry& thisDb() const;

            //- Name function is needed to disambiguate those inherited
            //  from base classes
            const word& name() const
            {
                return thisDb().name();
            }

            //- Return constant reference to boundary mesh
            const faBoundaryMesh& boundary() const;

            //- Return faMesh face labels
            const labelList& faceLabels() const
            {
                return faceLabels_;
            }


            //- Return parallel info
            const faGlobalMeshData& globalData() const;

            //- Return ldu addressing
            virtual const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return boundary().interfaces();
            }

            //- Internal face owner
            const labelUList& owner() const
            {
                return lduAddr().lowerAddr();
            }

            //- Internal face neighbour
            const labelUList& neighbour() const
            {
                return lduAddr().upperAddr();
            }

            //- Return true if given edge label is internal to the mesh
            inline bool isInternalEdge(const label edgeIndex) const
            {
                return edgeIndex < nInternalEdges();
            }


        // Mesh motion and morphing

            //- Is mesh moving
            bool moving() const
            {
                return mesh().moving();
            }

            //- Update after mesh motion
            virtual bool movePoints();

            //- Update after topo change
            virtual void updateMesh(const mapPolyMesh&);


        // Mapping

            //- Map all fields in time using given map.
            virtual void mapFields(const faMeshMapper& mapper) const;

            //- Map face areas in time using given map.
            virtual void mapOldAreas(const faMeshMapper& mapper) const;


        // Demand-driven data

            //- Return constant reference to primitive patch
            const indirectPrimitivePatch& patch() const;

            //- Return reference to primitive patch
            indirectPrimitivePatch& patch();

            //- Return patch starts
            const labelList& patchStarts() const;

            //- Return edge length vectors
            const edgeVectorField& Le() const;

            //- Return edge length magnitudes
            const edgeScalarField& magLe() const;

            //- Return face centres as areaVectorField
            const areaVectorField& areaCentres() const;

            //- Return edge centres as edgeVectorField
            const edgeVectorField& edgeCentres() const;

            //- Return face areas
            const DimensionedField<scalar, areaMesh>& S() const;

            //- Return old-time face areas
            const DimensionedField<scalar, areaMesh>& S0() const;

            //- Return old-old-time face areas
            const DimensionedField<scalar, areaMesh>& S00() const;

            //- Return face area normals
            const areaVectorField& faceAreaNormals() const;

            //- Return edge area normals
            const edgeVectorField& edgeAreaNormals() const;

            //- Return point area normals
            const vectorField& pointAreaNormals() const;

            //- Return face curvatures
            const areaScalarField& faceCurvatures() const;

            //- Return edge transformation tensors
            const FieldField<Field, tensor>& edgeTransformTensors() const;

            //- Return internal point labels
            labelList internalPoints() const;

            //- Return boundary point labels
            labelList boundaryPoints() const;

            //- Return edge length correction
            tmp<edgeScalarField> edgeLengthCorrection() const;

            //- Whether point normals should be corrected for a patch
            bool correctPatchPointNormals(const label patchID) const;

            //- Set whether point normals should be corrected for a patch
            boolList& correctPatchPointNormals() const;

        //- Write mesh
        virtual bool write(const bool valid = true) const;


    // Member Operators

        bool operator!=(const faMesh& m) const;

        bool operator==(const faMesh& m) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faPatchFaMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

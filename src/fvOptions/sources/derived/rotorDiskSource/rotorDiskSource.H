/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::fv::rotorDiskSource

Group
    grpFvOptionsSources

Description
    Rotor disk source

    Cell based momemtum source which approximates the mean effects of
    rotor forces on a cylindrical region within the domain.

Usage
    Example usage:
    \verbatim
        fields          (U);    // names of fields on which to apply source
        nBlades         3;      // number of blades
        tipEffect       0.96;   // normalised radius above which lift = 0

        inletFlowType   local;  // inlet flow type specification

        geometryMode    auto;   // geometry specification

        refDirection    (-1 0 0); // reference direction
                                  // - used as reference for psi angle

        trimModel       fixed;  // fixed || targetForce

        flapCoeffs
        {
            beta0           0;  // coning angle [deg]
            beta1c          0;  // lateral flapping coeff (cos coeff)
            beta2s          0;  // longitudinal flapping coeff (sin coeff)
        }

        blade
        {
            // see bladeModel.H for documentation
        }

        profiles
        {
            profile1
            {
                type    lookup; // lookup || series
                ...
                // see lookupProfile.H or seriesProfile.H for documentation
            }
            profile2
            {
                ...
            }
        }
    \endverbatim

    Where:
    Valid options for the \c geometryMode entry include:
    - auto          : determine rotor coordinate system from cells
    - specified     : specified coordinate system

    Valid options for the \c inletFlowType entry include:
    - fixed         : specified velocity
    - local         : use local flow conditions
    - surfaceNormal : specified normal velocity (positive towards rotor)

See also
    Foam::bladeModel
    Foam::lookupProfile
    Foam::seriesProfile

SourceFiles
    rotorDiskSource.C
    rotorDiskSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef rotorDiskSource_H
#define rotorDiskSource_H

#include "cellSetOption.H"
#include "cylindricalCS.H"
#include "Enum.H"
#include "bladeModel.H"
#include "profileModelList.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class trimModel;

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class rotorDiskSource Declaration
\*---------------------------------------------------------------------------*/

class rotorDiskSource
:
    public cellSetOption
{
public:

    enum geometryModeType
    {
        gmAuto,
        gmSpecified
    };
    static const Enum<geometryModeType> geometryModeTypeNames_;

    enum inletFlowType
    {
        ifFixed,
        ifSurfaceNormal,
        ifLocal
    };
    static const Enum<inletFlowType> inletFlowTypeNames_;


protected:

    // Helper structures to encapsulate flap and trim data
    // Note: all input in degrees (converted to radians internally)

        struct flapData
        {
            scalar beta0;   // coning angle
            scalar beta1c;  // lateral flapping coeff (cos coeff)
            scalar beta2s;  // longitudinal flapping coeff (sin coeff)
        };


    // Protected data

        //- Reference density for incompressible case
        scalar rhoRef_;

        //- Rotational speed [rad/s]
        //  Positive anti-clockwise when looking along -ve lift direction
        scalar omega_;

        //- Number of blades
        label nBlades_;

        //- Inlet flow type
        inletFlowType inletFlow_;

        //- Inlet velocity for specified inflow
        vector inletVelocity_;

        //- Tip effect [0-1]
        //  Ratio of blade radius beyond which lift=0
        scalar tipEffect_;

        //- Blade flap coefficients [rad/s]
        flapData flap_;

        //- Cell centre positions in local rotor frame
        //  (Cylindrical r-theta-z)
        List<point> x_;

        //- Rotation tensor for flap angle
        List<tensor> Rcone_;

        //- Area [m2]
        List<scalar> area_;

        //- Rotor local cylindrical coordinate system (r-theta-z)
        coordSystem::cylindrical coordSys_;

        //- Cached rotation tensors for cylindrical coordinates
        autoPtr<tensorField> Rcyl_;

        //- Maximum radius
        scalar rMax_;

        //- Trim model
        autoPtr<trimModel> trim_;

        //- Blade data
        bladeModel blade_;

        //- Profile data
        profileModelList profiles_;


    // Protected Member Functions

        //- Check data
        void checkData();

        //- Set the face areas per cell, and optionally correct the rotor axis
        void setFaceArea(vector& axis, const bool correct);

        //- Create the coordinate system
        void createCoordinateSystem();

        //- Construct geometry
        void constructGeometry();

        //- Return the inlet flow field
        tmp<vectorField> inflowVelocity(const volVectorField& U) const;

        //- Helper function to write rotor values
        template<class Type>
        void writeField
        (
            const word& name,
            const List<Type>& values,
            const bool writeNow = false
        ) const;


public:

    //- Runtime type information
    TypeName("rotorDisk");


    // Constructors

        //- Construct from components
        rotorDiskSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~rotorDiskSource();


    // Member Functions

        // Access

            //- Return the reference density for incompressible case
            inline scalar rhoRef() const;

            //- Return the rotational speed [rad/s]
            //  Positive anti-clockwise when looking along -ve lift direction
            inline scalar omega() const;

            //- Return the cell centre positions in local rotor frame
            //  (Cylindrical r-theta-z)
            inline const List<point>& x() const;

            //- Return the rotor coordinate system (r-theta-z)
            inline const coordSystem::cylindrical& coordSys() const;


        // Evaluation

            //- Calculate forces
            template<class RhoFieldType>
            void calculate
            (
                const RhoFieldType& rho,
                const vectorField& U,
                const scalarField& thetag,
                vectorField& force,
                const bool divideVolume = true,
                const bool output = true
            ) const;


        // Source term addition

            //- Source term to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Source term to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rotorDiskSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "rotorDiskSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::projectCurveEdge

Description
    Defines the edge from the projection onto a surface (single surface)
    or intersection of two surfaces.

SourceFiles
    projectCurveEdge.C

\*---------------------------------------------------------------------------*/

#ifndef projectCurveEdge_H
#define projectCurveEdge_H

#include "blockEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class pointConstraint;

/*---------------------------------------------------------------------------*\
                      Class projectCurveEdge Declaration
\*---------------------------------------------------------------------------*/

class projectCurveEdge
:
    public blockEdge
{
    // Private Data

        const searchableSurfaces& geometry_;

        //- The indices of surfaces onto which the points are projected
        labelList surfaces_;


    // Private Member Functions

        //- No copy construct
        projectCurveEdge(const projectCurveEdge&) = delete;

        //- No copy assignment
        void operator=(const projectCurveEdge&) = delete;


public:

    //- Runtime type information
    TypeName("projectCurve");


    // Constructors

        //- Construct from Istream setting pointsList
        projectCurveEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField& points,
            Istream& is
        );


    //- Destructor
    virtual ~projectCurveEdge() = default;


    // Member Functions

        //- Return the point positions corresponding to the curve parameters
        //  0 <= lambda <= 1
        virtual point position(const scalar) const
        {
            NotImplemented;
            return point::max;
        }

        //- Return the point positions corresponding to the curve parameters
        //  0 <= lambda <= 1
        virtual tmp<pointField> position(const scalarList&) const;

        //- Return the length of the curve
        virtual scalar length() const
        {
            NotImplemented;
            return 1;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

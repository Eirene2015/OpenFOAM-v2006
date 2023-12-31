/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 Bernhard Gschaider
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
    Foam::expressions::patchExprFieldBase

Description
    Base class for managing patches with expressions.
    The expected input supports values, gradients and mixed conditions

Usage
    \table
        Property     | Description                          | Required | Default
        valueExpr    | expression for fixed value           | no  | 0
        gradientExpr | expression for patch normal gradient | no  | 0
        fractionExpr | expression for value fraction weight | no  | 1
    \endtable

SourceFiles
    patchExprFieldBase.C

\*---------------------------------------------------------------------------*/

#ifndef expressions_patchExprFieldBase_H
#define expressions_patchExprFieldBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "exprString.H"

namespace Foam
{

// Forward Declarations
class facePointPatch;
class fvPatch;
class polyPatch;

namespace expressions
{

/*---------------------------------------------------------------------------*\
                     Class patchExprFieldBase Declaration
\*---------------------------------------------------------------------------*/

class patchExprFieldBase
{
protected:

    // Protected Data

        bool debug_;
        bool allowGradient_;

        //- Slightly dodgy concept here
        bool evalOnConstruct_;

        // The expressions
        expressions::exprString valueExpr_;
        expressions::exprString gradExpr_;
        expressions::exprString fracExpr_;


public:

    // Constructors

        //- Default construct
        patchExprFieldBase();

        //- Construct with specified gradient handling
        explicit patchExprFieldBase(bool allowGradient);

        //- Construct from dictionary
        explicit patchExprFieldBase
        (
            const dictionary& dict,
            bool allowGradient = false,
            bool isPointVal = false
        );

        //- Copy constructor
        patchExprFieldBase(const patchExprFieldBase& rhs);


    // Member Functions

        //- Write
        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::codedPoints0MotionSolver

Group

Description
    Provides a general interface to enable dynamic code compilation of
    mesh motion solvers.

    The entries are:
    \plaintable
       codeInclude | include files
       codeOptions | include paths; inserted into EXE_INC in Make/options
       codeLibs    | link line; inserted into LIB_LIBS in Make/options
       localCode   | c++; local static functions;
       code        | c++; upon motionSolver::curPoints();
    \endplaintable

    Note that the dynamically generated motionSolver is an points0MotionSolver,
    i.e. it has access to the points0 (points from constant/polyMesh).

Usage
    Example of simplistic specification:
    \verbatim
    difference
    {
        type coded;
        // Name of on-the-fly generated motion solver
        name myMotion;
        code
        #{
            const pointField& p0 = points0();
            const tensor rotT =
                rotationTensor(vector(1, 0, 0), vector(0, 1, 0));
            return transform(rotT, p0);
        #};
    }
    \endverbatim

See also
    Foam::motionSolver
    Foam::codedBase

SourceFiles
    codedPoints0MotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef codedPoints0MotionSolver_H
#define codedPoints0MotionSolver_H

#include "motionSolver.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                  Class codedPoints0MotionSolver Declaration
\*---------------------------------------------------------------------------*/

class codedPoints0MotionSolver
:
    public motionSolver,
    public codedBase
{
protected:

    // Protected Data

        //- Name of redirected motion solver
        word name_;

        //- Underlying motionSolver
        mutable autoPtr<motionSolver> redirectMotionSolverPtr_;


    // Protected Member Functions

        //- Get the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Return a description (type + name) for the output
        virtual string description() const;

        //- Clear any redirected objects
        virtual void clearRedirect() const;

        // The dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;


        //- No copy construct
        codedPoints0MotionSolver(const codedPoints0MotionSolver&) = delete;

        //- No copy assignment
        void operator=(const codedPoints0MotionSolver&) = delete;


public:

    // Static Data Members

        //- Name of the C code template to be used
        static constexpr const char* const codeTemplateC
            = "codedPoints0MotionSolverTemplate.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
            = "codedPoints0MotionSolverTemplate.H";


    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from mesh and dictionary
        codedPoints0MotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~codedPoints0MotionSolver() = default;


    // Member Functions

        //- Dynamically compiled motionSolver
        motionSolver& redirectMotionSolver() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField& fld);

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

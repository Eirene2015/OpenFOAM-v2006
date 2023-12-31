/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::PatchFunction1Types::CodedField

Description
    PatchFunction1 with the code supplied by an on-the-fly compiled C++
    expression.

Usage
    Example:
    \verbatim
    <patchName>
    {
        type            uniformFixedValue;
        uniformValue
        {
            type            coded;
            // Explictly supply name of generated PatchFunction1. Only needed
            // if entryname ('uniformValue') would clash with existing
            // runtime selection tables.
            name            myExpression;

            code
            #{
                const polyPatch& pp = this->patch();
                Pout<< "** Patch size:" << pp.size() << endl;
                return tmp<vectorField>::New(pp.size(), vector(1, 0, 0))
            #};
        }

        //codeInclude
        //#{
        //    #include "fvCFD.H"
        //#};

        //codeOptions
        //#{
        //    -I$(LIB_SRC)/finiteVolume/lnInclude
        //#};
    }
    \endverbatim

See also
    Foam::dynamicCode
    Foam::codedFixedValue
    Foam::functionEntries::codeStream

SourceFiles
    CodedField.C

\*---------------------------------------------------------------------------*/

#ifndef PatchFunction1Types_CodedField_H
#define PatchFunction1Types_CodedField_H

#include "PatchFunction1.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                     Class CodedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CodedField
:
    public PatchFunction1<Type>,
    protected codedBase
{
    // Private Data

        //- Dictionary contents for the function
        const dictionary dict_;

        const word name_;

        mutable autoPtr<PatchFunction1<Type>> redirectFunctionPtr_;


    // Private Member Functions

        //- Get the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        // Return a description (type + name) for the output
        virtual string description() const;

        // Clear the ptr to the redirected object
        virtual void clearRedirect() const;

        //- Get reference to the underlying Function1
        const PatchFunction1<Type>& redirectFunction() const;

        // Get the (sub)dictionary to initialize the codeContext
        virtual const dictionary& codeDict(const dictionary& fullDict) const;

        // Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;

        //- No copy assignment
        void operator=(const CodedField<Type>&) = delete;


public:

    // Static Data Members

        //- Name of the C code template to be used
        static constexpr const char* const codeTemplateC
            = "codedPatchFunction1Template.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
            = "codedPatchFunction1Template.H";


    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct from entry name and dictionary
        CodedField
        (
            const polyPatch& pp,
            const word& type,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Copy construct
        explicit CodedField(const CodedField<Type>& rhs);

        //- Copy construct, setting patch
        explicit CodedField
        (
            const CodedField<Type>& rhs,
            const polyPatch& pp
        );

        //- Construct and return a clone
        virtual tmp<PatchFunction1<Type>> clone() const
        {
            return tmp<PatchFunction1<Type>>
            (
                new CodedField<Type>(*this)
            );
        }

        //- Construct and return a clone setting patch
        virtual tmp<PatchFunction1<Type>> clone(const polyPatch& pp) const
        {
            return tmp<PatchFunction1<Type>>
            (
                new CodedField<Type>(*this, pp)
            );
        }


    //- Destructor
    virtual ~CodedField() = default;


    // Member Functions

        // Evaluation

            //- Return CodedField value
            virtual tmp<Field<Type>> value(const scalar x) const;

            //- Is value constant (i.e. independent of x)
            virtual inline bool constant() const
            {
                return false;
            }

            //- Is value uniform (i.e. independent of coordinate)
            virtual inline bool uniform() const
            {
                return false;
            }

            //- Integrate between two values
            virtual tmp<Field<Type>> integrate
            (
                const scalar x1,
                const scalar x2
            ) const;


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const FieldMapper& mapper);

            //- Reverse map the given PatchFunction1 onto this PatchFunction1
            virtual void rmap
            (
                const PatchFunction1<Type>& pf1,
                const labelList& addr
            );


        // I-O

            //- Write in dictionary format
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CodedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

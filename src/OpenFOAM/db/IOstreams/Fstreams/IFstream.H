/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::IFstream

Description
    Input from file stream, using an ISstream

SourceFiles
    IFstream.C

\*---------------------------------------------------------------------------*/

#ifndef IFstream_H
#define IFstream_H

#include "ISstream.H"
#include "fileName.H"
#include "className.H"
#include <fstream>
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace Detail
{

/*---------------------------------------------------------------------------*\
                  Class Detail::IFstreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- A std::istream with the ability to handle compressed files
class IFstreamAllocator
{
protected:

    // Member Data

        //- The allocated stream pointer (ifstream or igzstream).
        std::unique_ptr<std::istream> allocatedPtr_;

        //- The detected compression type
        IOstream::compressionType detectedCompression_;


    // Constructors

        //- Construct from pathname
        IFstreamAllocator(const fileName& pathname);
};

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                          Class IFstream Declaration
\*---------------------------------------------------------------------------*/

class IFstream
:
    public Detail::IFstreamAllocator,
    public ISstream
{
public:

    //- Declare type-name (with debug switch)
    ClassName("IFstream");


    // Constructors

        //- Construct from pathname
        explicit IFstream
        (
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Construct from pathname, format (version)
        IFstream
        (
            const fileName& pathname,
            streamFormat fmt,
            versionNumber ver = currentVersion
        )
        :
            IFstream(pathname, IOstreamOption(fmt, ver))
        {}


    //- Destructor
    ~IFstream() = default;


    // Member Functions

    // Access

        //- Read/write access to the name of the stream
        using ISstream::name;


    // STL stream

        //- Access to underlying std::istream
        virtual std::istream& stdStream();

        //- Const access to underlying std::istream
        virtual const std::istream& stdStream() const;

        //- Rewind the stream so that it may be read again
        virtual void rewind();


    // Print

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const;


    // Member Operators

        //- Return a non-const reference to const IFstream
        //  Needed for read-constructors where the stream argument is temporary:
        //  e.g. thing thisThing(IFstream("fileName")());
        IFstream& operator()() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

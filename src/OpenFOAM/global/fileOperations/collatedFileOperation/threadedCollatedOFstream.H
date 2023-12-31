/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::threadedCollatedOFstream

Description
    Master-only drop-in replacement for OFstream.

SourceFiles
    threadedCollatedOFstream.C

\*---------------------------------------------------------------------------*/

#ifndef threadedCollatedOFstream_H
#define threadedCollatedOFstream_H

#include "StringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstreamCollator;

/*---------------------------------------------------------------------------*\
                  Class threadedCollatedOFstream Declaration
\*---------------------------------------------------------------------------*/

class threadedCollatedOFstream
:
    public OStringStream
{
    // Private data

        OFstreamCollator& writer_;

        const fileName pathName_;

        const IOstream::compressionType compression_;

        const bool useThread_;


public:

    // Constructors

        //- Construct and set stream status
        threadedCollatedOFstream
        (
            OFstreamCollator& writer,
            const fileName& pathname,
            IOstreamOption streamOpt = IOstreamOption(),
            const bool useThread = true
        );

        //- Construct and set stream status
        threadedCollatedOFstream
        (
            OFstreamCollator& writer,
            const fileName& pathname,
            streamFormat fmt,
            versionNumber ver = currentVersion,
            compressionType comp = compressionType::UNCOMPRESSED,
            const bool useThread = true
        )
        :
            threadedCollatedOFstream
            (
                writer,
                pathname,
                IOstreamOption(fmt, ver, comp),
                useThread
            )
        {}


    //- Destructor
    ~threadedCollatedOFstream();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

\*---------------------------------------------------------------------------*/

#include "ensightOutput.H"
#include "ensightPTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<typename> class FieldContainer, class Type>
void Foam::ensightOutput::Detail::copyComponent
(
    scalarField& res,
    const FieldContainer<Type>& input,
    const direction cmpt
)
{
    res.resize(input.size());

    auto iter = res.begin();

    for (const Type& val : input)
    {
        *iter = component(val, cmpt);
        ++iter;
    }
}


template<template<typename> class FieldContainer>
bool Foam::ensightOutput::Detail::writeCoordinates
(
    ensightGeoFile& os,
    const label partId,
    const word& partName,
    const label nPoints,
    const FieldContainer<Foam::point>& fld,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    const label nSlaves = (parallel ? Pstream::nProcs() : 0);

    // Using manual copyComponent(...) instead of fld.component() to support
    // indirect lists etc.

    scalarField send(fld.size());

    if (Pstream::master())
    {
        // Serial output, or parallel (master)

        os.beginPart(partId, partName);
        os.beginCoordinates(nPoints);

        for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
        {
            copyComponent(send, fld, cmpt);
            os.writeList(send);

            for (int slave=1; slave < nSlaves; ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                scalarField recv(fromSlave);
                os.writeList(recv);
            }
        }
    }
    else if (nSlaves)
    {
        // Parallel (slaves)

        for (direction cmpt=0; cmpt < point::nComponents; ++cmpt)
        {
            copyComponent(send, fld, cmpt);

            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster << send;
        }
    }

    return true;
}


template<template<typename> class FieldContainer, class Type>
bool Foam::ensightOutput::Detail::writeFieldComponents
(
    ensightFile& os,
    const char* key,
    const FieldContainer<Type>& fld,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    const label nSlaves = (parallel ? Pstream::nProcs() : 0);

    // Preliminary checks
    {
        bool hasField = !fld.empty();

        if (parallel)
        {
            reduce(hasField, orOp<bool>());
        }

        // No field
        if (!hasField) return false;
    }


    // Using manual copyComponent(...) instead of fld.component() to support
    // indirect lists etc.

    scalarField send(fld.size());

    if (Pstream::master())
    {
        // Serial output, or parallel (master)

        os.writeKeyword(key);

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const direction cmpt = ensightPTraits<Type>::componentOrder[d];

            copyComponent(send, fld, cmpt);
            os.writeList(send);

            for (label slave=1; slave < nSlaves; ++slave)
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);
                scalarField recv(fromSlave);
                os.writeList(recv);
            }
        }
    }
    else if (nSlaves)
    {
        // Parallel (slaves)

        for (direction d=0; d < pTraits<Type>::nComponents; ++d)
        {
            const direction cmpt = ensightPTraits<Type>::componentOrder[d];

            copyComponent(send, fld, cmpt);

            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            toMaster << send;
        }
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::Detail::writeFaceSubField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Preliminary checks: total() contains pre-reduced information
    {
        // No geometry
        if (parallel ? !part.total() : !part.size()) return false;

        bool hasField = !fld.empty();

        if (parallel)
        {
            reduce(hasField, orOp<bool>());
        }

        // No field
        if (!hasField) return false;
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    labelList localAddr;

    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        ensightOutput::Detail::writeFieldComponents
        (
            os,
            ensightFaces::key(etype),
            SubField<Type>(fld, part.range(etype)),
            parallel
        );
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::writeField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightCells& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Preliminary checks: total() contains pre-reduced information
    {
        // No geometry
        if (parallel ? !part.total() : !part.size()) return false;

        bool hasField = !fld.empty();

        if (parallel)
        {
            reduce(hasField, orOp<bool>());
        }

        // No field
        if (!hasField) return false;
    }

    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightCells::nTypes; ++typei)
    {
        const auto etype = ensightCells::elemType(typei);

        ensightOutput::Detail::writeFieldComponents
        (
            os,
            ensightCells::key(etype),
            UIndirectList<Type>(fld, part.cellIds(etype)),
            parallel
        );
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::writeField
(
    ensightFile& os,
    const Field<Type>& fld,
    const ensightFaces& part,
    bool parallel
)
{
    parallel = parallel && Pstream::parRun();

    // Preliminary checks: total() contains pre-reduced information
    {
        // No geometry
        if (parallel ? !part.total() : !part.size()) return false;

        bool hasField = !fld.empty();

        if (parallel)
        {
            reduce(hasField, orOp<bool>());
        }

        // No field
        if (!hasField) return false;
    }


    if (Pstream::master())
    {
        os.beginPart(part.index());
    }

    for (int typei=0; typei < ensightFaces::nTypes; ++typei)
    {
        const auto etype = ensightFaces::elemType(typei);

        ensightOutput::Detail::writeFieldComponents
        (
            os,
            ensightFaces::key(etype),
            UIndirectList<Type>(fld, part.faceIds(etype)),
            parallel
        );
    }

    return true;
}


// ************************************************************************* //

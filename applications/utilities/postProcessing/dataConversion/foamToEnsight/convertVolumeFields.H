/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

Description
    Code chunk for converting volume and dimensioned fields
    included by foamToEnsight.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Cell field data output
{
    Info<< "Write volume field (";

    writeAllVolFields(ensCase, ensMesh, objects, nearCellValue);
    writeAllDimFields(ensCase, ensMesh, objects);

    Info<< " )" << nl;

    // PointData
    // - only construct pointMesh on request (it constructs edge addressing)
    if (doPointValues)
    {
        Info<< "Write point field (";
        writeAllPointFields(ensCase, ensMesh, objects);
        Info<< " )" << nl;
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
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

#include "objectivePartialVolume.H"
#include "createZeroField.H"
#include "IOmanip.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectivePartialVolume, 1);
addToRunTimeSelectionTable
(
    objectiveIncompressible,
    objectivePartialVolume,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectivePartialVolume::objectivePartialVolume
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objectiveIncompressible(mesh, dict, adjointSolverName, primalSolverName),
    initVol_(Zero),
    objectivePatches_
    (
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        ).sortedToc()
    )
{
    // Read target volume if present. Else use the current one as a target
    if (!dict.readIfPresent("initialVolume", initVol_))
    {
        const scalar oneThird(1.0/3.0);
        for (const label patchi : objectivePatches_)
        {
            const fvPatch& patch = mesh_.boundary()[patchi];
            initVol_ -= oneThird*gSum(patch.Sf() & patch.Cf());
        }
    }
    // Allocate boundary field pointers
    bdxdbDirectMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    bdSdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar objectivePartialVolume::J()
{
    J_ = Zero;
    const scalar oneThird(1.0/3.0);
    for (const label patchi : objectivePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        J_ -= oneThird*gSum(patch.Sf() & patch.Cf());
    }
    J_ -= initVol_;
    J_ /= initVol_;
    return J_;
}


void objectivePartialVolume::update_dxdbDirectMultiplier()
{
    const scalar oneThird(1.0/3.0);
    for (const label patchi : objectivePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        tmp<vectorField> tnf = patch.nf();
        const vectorField& nf = tnf();
        bdxdbDirectMultPtr_()[patchi] = -oneThird*nf/initVol_;
    }
}


void objectivePartialVolume::update_dSdbMultiplier()
{
    const scalar oneThird(1.0/3.0);
    for (const label patchi : objectivePatches_)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        bdSdbMultPtr_()[patchi] = -oneThird*patch.Cf()/initVol_;
    }
}


bool objectivePartialVolume::write(const bool valid) const
{
    if (Pstream::master())
    {
        // File is opened only upon invocation of the write function
        // in order to avoid various instantiations of the same objective
        // opening the same file
        unsigned int width = IOstream::defaultPrecision() + 6;
        if (objFunctionFilePtr_.empty())
        {
            setObjectiveFilePtr();
            objFunctionFilePtr_() << setw(4)     << "#"               << " ";
            objFunctionFilePtr_() << setw(width) << "(V - VInit)/VInit" << " ";
            objFunctionFilePtr_() << setw(width) << "VInit" << endl;
        }

        objFunctionFilePtr_() << setw(4)     << mesh_.time().value() << " ";
        objFunctionFilePtr_() << setw(width) << J_ << " ";
        objFunctionFilePtr_() << setw(width) << initVol_ << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// ************************************************************************* //

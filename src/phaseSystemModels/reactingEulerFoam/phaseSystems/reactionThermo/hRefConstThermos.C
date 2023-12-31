/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "makeReactionThermo.H"
#include "makeThermo.H"

#include "rhoReactionThermo.H"
#include "heRhoThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "rPolynomial.H"
#include "perfectFluid.H"
#include "rhoConst.H"

#include "sensibleEnthalpy.H"

#include "hRefConstThermo.H"
#include "eRefConstThermo.H"

#include "constTransport.H"
#include "H2_gTransport.H"
#include "H2_lTransport.H"

#include "pureMixture.H"
#include "multiComponentMixture.H"

#include "thermoPhysicsTypes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Thermo type typedefs:

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectGas<specie>
        >,
        sensibleEnthalpy
    >
> constRefGasHThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectFluid<specie>
        >,
        sensibleEnthalpy
    >
> constRefFluidHThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            rPolynomial<specie>
        >,
        sensibleEnthalpy
    >
> constRefrPolFluidHThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        eRefConstThermo
        <
            perfectGas<specie>
        >,
        sensibleInternalEnergy
    >
> constRefGasEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        eRefConstThermo
        <
            perfectFluid<specie>
        >,
        sensibleInternalEnergy
    >
> constRefFluidEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        eRefConstThermo
        <
            rPolynomial<specie>
        >,
        sensibleInternalEnergy
    >
> constRefrPolFluidEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        eRefConstThermo
        <
            rhoConst<specie>
        >,
        sensibleInternalEnergy
    >
> constRefRhoConstEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            rhoConst<specie>
        >,
        sensibleEnthalpy
    >
> constRefRhoConstHThermoPhysics;


// pureMixture, sensibleEnthalpy:

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    rPolynomial,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectFluid,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    rhoConst,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_gTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_gTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    rPolynomial,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_gTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectFluid,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_gTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    rhoConst,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_lTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_lTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    rPolynomial,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_lTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectFluid,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_lTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    rhoConst,
    specie
);

// pureMixture, sensibleInternalEnergy:

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    rPolynomial,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    perfectFluid,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    rhoConst,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_gTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_gTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    rPolynomial,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_gTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    perfectFluid,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_gTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    rhoConst,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_lTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_lTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    rPolynomial,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_lTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    perfectFluid,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    H2_lTransport,
    sensibleInternalEnergy,
    eRefConstThermo,
    rhoConst,
    specie
);

// multiComponentMixture, sensibleInternalEnergy:

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefGasEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefFluidEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefrPolFluidEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefRhoConstEThermoPhysics
);


// multiComponentMixture, sensibleEnthalpy:

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefRhoConstHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefFluidHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefrPolFluidHThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefGasHThermoPhysics
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

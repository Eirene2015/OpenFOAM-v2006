/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::H2_gTransport

Group
    grpSpecieTransport

Description
    properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    H2_gTransportI.H
    H2_gTransport.C

\*---------------------------------------------------------------------------*/

#ifndef H2_gTransport_H
#define H2_gTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo> class H2_gTransport;

template<class Thermo>
inline H2_gTransport<Thermo> operator+
(
    const H2_gTransport<Thermo>&,
    const H2_gTransport<Thermo>&
);

template<class Thermo>
inline H2_gTransport<Thermo> operator*
(
    const scalar,
    const H2_gTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const H2_gTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class H2_gTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class H2_gTransport
:
    public Thermo
{
    // Private Data

        //- Reciprocal Prandtl Number [-]
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline H2_gTransport
        (
            const Thermo& t,
            const scalar Pr
        );


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct as named copy
        inline H2_gTransport(const word&, const H2_gTransport&);

        //- Construct from dictionary
        explicit H2_gTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<H2_gTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<H2_gTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "H2_g<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Heat capacity at constant pressure [J/(kmol K)]
        inline scalar cp(const scalar p, const scalar T) const;

        //- Heat capacity at constant volume [J/(kmol K)]
        inline scalar cv(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const H2_gTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend H2_gTransport operator+ <Thermo>
        (
            const H2_gTransport&,
            const H2_gTransport&
        );

        friend H2_gTransport operator* <Thermo>
        (
            const scalar,
            const H2_gTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const H2_gTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "H2_gTransportI.H"

#ifdef NoRepository
    #include "H2_gTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

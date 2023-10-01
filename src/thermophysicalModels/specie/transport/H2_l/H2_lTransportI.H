/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
//- #include <cmath> // Include the C++ math library for pow

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::H2_lTransport<Thermo>::H2_lTransport
(
    const Thermo& t,
    const scalar Pr
)
:
    Thermo(t),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::H2_lTransport<Thermo>::H2_lTransport
(
    const word& name,
    const H2_lTransport& ct
)
:
    Thermo(name, ct),
    rPr_(ct.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::H2_lTransport<Thermo>>
Foam::H2_lTransport<Thermo>::clone() const
{
    return autoPtr<H2_lTransport<Thermo>>::New(*this);
}


template<class Thermo>
inline Foam::autoPtr<Foam::H2_lTransport<Thermo>>
Foam::H2_lTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<H2_lTransport<Thermo>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::H2_lTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{

    //- Fitted Curve For P = 80000 Pa: 276.3824 * Foam::pow(T - 65.4049, 3) + -14.3984 * Foam::pow(T - 65.4049, 2) + 0.5491 * (T - 65.4049) + -0.0154
    //- Fitted Curve For P = 1000000 Pa: 92.5502 * Foam::pow(T - 93.7469, 3) + -2.8942 * Foam::pow(T - 93.7469, 2) + 0.0540 * (T - 93.7469) + -0.0013
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar muValue = 0.0;

    if (p <= P1)
    {
        muValue = ( 276.3824 * Foam::pow(T - 65.4049, 3) - 14.3984 * Foam::pow(T - 65.4049, 2) + 0.5491 * (T - 65.4049) - 0.0154 ) * 1e-6;
    }
    else if (p <= P2)
    {
        scalar muP1 = ( 276.3824 * Foam::pow(T - 65.4049, 3) - 14.3984 * Foam::pow(T - 65.4049, 2) + 0.5491 * (T - 65.4049) - 0.0154) * 1e-6;
        scalar muP2 = ( 92.5502 * Foam::pow(T - 93.7469, 3) - 2.8942 * Foam::pow(T - 93.7469, 2) + 0.0540 * (T - 93.7469) - 0.0013 ) * 1e-6;
        scalar muWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        muValue = muP1 + muWeight * (muP2 - muP1);
    }
    else
    {
        muValue = ( 92.5502 * Foam::pow(T - 93.7469, 3) - 2.8942 * Foam::pow(T - 93.7469, 2) + 0.0540 * (T - 93.7469) - 0.0013 ) * 1e-6;
    }

    return muValue;
}

template<class Thermo>
inline Foam::scalar Foam::H2_lTransport<Thermo>::cp
(
    const scalar p,
    const scalar T
) const
{
    //- TODO: i have to check to see what is the p value in openfoam. (is it pressure or prgh?)

    //- Fitted Curve P = 800000 Pa: 2134.4985 * Foam::pow(T - 90.3987, 3) + 12.2402 * Foam::pow(T - 90.3987, 2) + 0.4811 * Foam::pow(T - 90.3987) + 0.0182
    //- Fitted Curve P = 1000000 Pa: 2175.0394 * Foam::pow(T - 93.7469, 3) + 15.5442 * Foam::pow(T - 93.7469, 2) + 0.6391 * Foam::pow(T - 93.7469) + 0.0276
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar cpValue = 0.0;

    if (p <= P1)
    {
        cpValue = 2134.4985 * Foam::pow(T - 90.3987, 3) + 12.2402 * Foam::pow(T - 90.3987, 2) + 0.4811 * (T - 90.3987) + 0.0182;
    }
    else if (p <= P2)
    {
        scalar cpP1 = 2134.4985 * Foam::pow(T - 90.3987, 3) + 12.2402 * Foam::pow(T - 90.3987, 2) + 0.4811 * (T - 90.3987) + 0.0182;
        scalar cpP2 = 2175.0394 * Foam::pow(T - 93.7469, 3) + 15.5442 * Foam::pow(T - 93.7469, 2) + 0.6391 * (T - 93.7469) + 0.0276;
        scalar cpWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        cpValue = cpP1 + cpWeight * (cpP2 - cpP1);
    }
    else
    {
        cpValue = 2175.0394 * Foam::pow(T - 93.7469, 3) + 15.5442 * Foam::pow(T - 93.7469, 2) + 0.6391 * (T - 93.7469) + 0.0276;
    }

    return cpValue;
}



template<class Thermo>
inline Foam::scalar Foam::H2_lTransport<Thermo>::cv
(
    const scalar p,
    const scalar T
) const
{
    //- TODO: i have to check to see what is the p value in openfoam. (is it pressure or prgh?)


    //- Fitted Curve P = 80000 Pa: 1706.9617 * pow(T - 65.4049, 3) + 1.8709 * pow(T - 65.4049, 2) + 0.0727 * (T - 65.4049) + 0.0031
    //- Fitted Curve P = 1000000 Pa: 1878.2529 * pow(T - 93.7469, 3) + 15.5442 * pow(T - 93.7469, 2) + 0.6391 * (T - 93.7469) + 0.0276
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar cpValue = 0.0;

    if (p <= P1)
    {
        cpValue = 1706.9617 * pow(T - 65.4049, 3) + 1.8709 * pow(T - 65.4049, 2) + 0.0727 * (T - 65.4049) + 0.0031;
    }
    else if (p <= P2)
    {
        scalar cpP1 = 1706.9617 * pow(T - 65.4049, 3) + 1.8709 * pow(T - 65.4049, 2) + 0.0727 * (T - 65.4049) + 0.0031;
        scalar cpP2 = 1878.2529 * pow(T - 93.7469, 3) + 15.5442 * pow(T - 93.7469, 2) + 0.6391 * (T - 93.7469) + 0.0276;
        scalar cpWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        cpValue = cpP1 + cpWeight * (cpP2 - cpP1);
    }
    else
    {
        cpValue = 1878.2529 * pow(T - 93.7469, 3) + 15.5442 * pow(T - 93.7469, 2) + 0.6391 * (T - 93.7469) + 0.0276;;
    }

    return cpValue;
}


template<class Thermo>
inline Foam::scalar Foam::H2_lTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return this->Cp(p, T)*mu(p, T)*rPr_;
}


template<class Thermo>
inline Foam::scalar Foam::H2_lTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return mu(p, T)*rPr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::H2_lTransport<Thermo>::operator+=
(
    const H2_lTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        rPr_ = 1.0/(Y1/rPr_ + Y2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::H2_lTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::H2_lTransport<Thermo> Foam::operator+
(
    const H2_lTransport<Thermo>& ct1,
    const H2_lTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return H2_lTransport<Thermo>
        (
            t,
            0,
            ct1.rPr_
        );
    }
    else
    {
        scalar Y1 = ct1.Y()/t.Y();
        scalar Y2 = ct2.Y()/t.Y();

        return H2_lTransport<Thermo>
        (
            t,
            1.0/(Y1/ct1.rPr_ + Y2/ct2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::H2_lTransport<Thermo> Foam::operator*
(
    const scalar s,
    const H2_lTransport<Thermo>& ct
)
{
    return H2_lTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        1.0/ct.rPr_
    );
}


// ************************************************************************* //

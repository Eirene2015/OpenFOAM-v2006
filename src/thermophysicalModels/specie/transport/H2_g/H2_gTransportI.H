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
// #include <cmath> // Include the C++ math library for pow

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::H2_gTransport<Thermo>::H2_gTransport
(
    const Thermo& t,
    const scalar Pr
)
:
    Thermo(t),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::H2_gTransport<Thermo>::H2_gTransport
(
    const word& name,
    const H2_gTransport& ct
)
:
    Thermo(name, ct),
    rPr_(ct.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::H2_gTransport<Thermo>>
Foam::H2_gTransport<Thermo>::clone() const
{
    return autoPtr<H2_gTransport<Thermo>>::New(*this);
}


template<class Thermo>
inline Foam::autoPtr<Foam::H2_gTransport<Thermo>>
Foam::H2_gTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<H2_gTransport<Thermo>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::H2_gTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{

    //- Fitted Curve For P = 80000 Pa: mu*1e6 = 5.3001 * (T - 75.5049)^3 + 0.0686 * (T - 75.5049)^2 
    //- Fitted Curve For P = 1000000 Pa: mu*1e6 = 7.8401 * (T - 103.8469)^3 + 0.0503 * (T - 103.8469)^2 + 0.0004 * (T - 103.8469)
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar muValue = 0.0;

    if (p <= P1)
    {
        muValue = (5.3001 * Foam::pow(T - 75.5049, 3) + 0.0686 * Foam::pow(T - 75.5049, 2))*1e-6;
    }
    else if (p <= P2)
    {
        scalar muP1 = ( 5.3001 * Foam::pow(T - 75.5049, 3) + 0.0686 * Foam::pow(T - 75.5049, 2) ) * 1e-6;
        scalar muP2 = ( 7.8401 * Foam::pow(T - 103.8469,3) + 0.0503 * Foam::pow(T - 103.8469, 2) + 0.0004 * (T - 103.8469) ) * 1e-6;
        scalar muWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        muValue = muP1 + muWeight * (muP2 - muP1);
    }
    else
    {
        muValue = ( 7.8401 * Foam::pow(T - 103.8469,3) + 0.0503 * Foam::pow(T - 103.8469, 2) + 0.0004 * (T - 103.8469) ) * 1e-6;
    }

    return muValue;
}


template<class Thermo>
inline Foam::scalar Foam::H2_gTransport<Thermo>::cp
(
    const scalar p,
    const scalar T
) const
{
    //- TODO: i have to check to see what is the p value in openfoam. (is it pressure or prgh?)


    //- Fitted Curve P = 800000 Pa: 1110.1966 * pow(T - 75.5049, 3) - 3.8299 * pow(T - 75.5049, 2) + 0.1524 * pow(T - 75.5049) - 0.0048
    //- Fitted Curve P = 150000 Pa: 1151.8260 * Foam::pow(T - 80.9446, 3) - 6.0014 * Foam::pow(T - 80.9446, 2) + 0.2508 * Foam::pow(T - 80.9446) - 0.0086
    //- Fitted Curve P = 200000 Pa: 1179.3649 * Foam::pow(T - 83.7258, 3) - 7.5313 * Foam::pow(T - 83.7258, 2) + 0.3284 * Foam::pow(T - 83.7258) - 0.0120
    //- Fitted Curve P = 400000 Pa: 1285.0742 * Foam::pow(T - 91.3327, 3) - 14.2465 * Foam::pow(T - 91.3327, 2) + 0.7310 * Foam::pow(T - 91.3327) - 0.0314
    //- Fitted Curve P = 500000 Pa: 1338.6252 * Foam::pow(T - 94.0950, 3) - 18.1936 * Foam::pow(T - 94.0950, 2) + 1.0047 * Foam::pow(T - 94.0950) - 0.0459
    //- Fitted Curve P = 1000000 Pa: 1646.8479 * Foam::pow(T - 103.8469, 3) - 48.6356 * Foam::pow(T - 103.8469, 2) + 3.7652 * Foam::pow(T - 103.8469) - 0.2203
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 150000;
    scalar P3 = 200000;
    scalar P4 = 400000;
    scalar P5 = 500000;
    scalar P6 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar cpValue = 0.0;

    if (p <= P1)
    {
        cpValue = 1110.1966 * pow(T - 75.5049, 3) - 3.8299 * pow(T - 75.5049, 2) + 0.1524 * (T - 75.5049) - 0.0048;
    }
    else if (p <= P2)
    {
        scalar cpP1 = 1110.1966 * pow(T - 75.5049, 3) - 3.8299 * pow(T - 75.5049, 2) + 0.1524 * (T - 75.5049) - 0.0048;
        scalar cpP2 = 1151.8260 * Foam::pow(T - 80.9446, 3) - 6.0014 * Foam::pow(T - 80.9446, 2) + 0.2508 * (T - 80.9446) - 0.0086;
        scalar cpWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        cpValue = cpP1 + cpWeight * (cpP2 - cpP1);
    }
    else if (p <= P3)
    {
        scalar cpP2 = 1151.8260 * Foam::pow(T - 80.9446, 3) - 6.0014 * Foam::pow(T - 80.9446, 2) + 0.2508 * (T - 80.9446) - 0.0086;
        scalar cpP3 = 1179.3649 * Foam::pow(T - 83.7258, 3) - 7.5313 * Foam::pow(T - 83.7258, 2) + 0.3284 * (T - 83.7258) - 0.0120;
        scalar cpWeight = (p - P2) / (P3 - P2);  // Linear interpolation weight
        cpValue = cpP2 + cpWeight * (cpP3 - cpP2);
    }
    else if (p <= P4)
    {
        scalar cpP3 = 1179.3649 * Foam::pow(T - 83.7258, 3) - 7.5313 * Foam::pow(T - 83.7258, 2) + 0.3284 * (T - 83.7258) - 0.0120;
        scalar cpP4 = 1285.0742 * Foam::pow(T - 91.3327, 3) - 14.2465 * Foam::pow(T - 91.3327, 2) + 0.7310 * (T - 91.3327) - 0.0314;
        scalar cpWeight = (p - P3) / (P4 - P3);  // Linear interpolation weight
        cpValue = cpP3 + cpWeight * (cpP4 - cpP3);
    }
    else if (p <= P5)
    {
        scalar cpP4 = 1285.0742 * Foam::pow(T - 91.3327, 3) - 14.2465 * Foam::pow(T - 91.3327, 2) + 0.7310 * (T - 91.3327) - 0.0314;
        scalar cpP5 = 1338.6252 * Foam::pow(T - 94.0950, 3) - 18.1936 * Foam::pow(T - 94.0950, 2) + 1.0047 * (T - 94.0950) - 0.0459;
        scalar cpWeight = (p - P4) / (P5 - P4);  // Linear interpolation weight
        cpValue = cpP4 + cpWeight * (cpP5 - cpP4);
    }    
    else if (p <= P6)
    {
        scalar cpP5 = 1338.6252 * Foam::pow(T - 94.0950, 3) - 18.1936 * Foam::pow(T - 94.0950, 2) + 1.0047 * (T - 94.0950) - 0.0459;
        scalar cpP6 = 1646.8479 * Foam::pow(T - 103.8469, 3) - 48.6356 * Foam::pow(T - 103.8469, 2) + 3.7652 * (T - 103.8469) - 0.2203;
        scalar cpWeight = (p - P5) / (P6 - P5);  // Linear interpolation weight
        cpValue = cpP5 + cpWeight * (cpP6 - cpP5);
    }
    else
    {
        cpValue = 1646.8479 * Foam::pow(T - 103.8469, 3) - 48.6356 * Foam::pow(T - 103.8469, 2) + 3.7652 * (T - 103.8469) - 0.2203;
    }

    return cpValue;
}



template<class Thermo>
inline Foam::scalar Foam::H2_gTransport<Thermo>::cv
(
    const scalar p,
    const scalar T
) const
{
    //- TODO: i have to check to see what is the p value in openfoam. (is it pressure or prgh?)

    //- Fitted Curve P = 800000 Pa: 813.4101 * pow(T - 75.5049, 3) - 3.8299 * pow(T - 75.5049, 2) + 0.1524 * (T - 75.5049) - 0.0048
    //- Fitted Curve P = 150000 Pa: 855.0395 * pow(T - 80.9446, 3) - 6.0014 * pow(T - 80.9446, 2) + 0.2508 * (T - 80.9446) - 0.0086
    //- Fitted Curve P = 200000 Pa: 882.5783 * pow(T - 83.7258, 3) - 7.5313 * pow(T - 83.7258, 2) + 0.3284 * (T - 83.7258) - 0.0120
    //- Fitted Curve P = 400000 Pa: 988.2877 * pow(T - 91.3327, 3) - 14.2465 * pow(T - 91.3327, 2) + 0.7310 * (T - 91.3327) - 0.0314
    //- Fitted Curve P = 500000 Pa: 1041.8387 * pow(T - 94.0950, 3) - 18.1936 * pow(T - 94.0950, 2) + 1.0047 * (T - 94.0950) - 0.0459
    //- Fitted Curve P = 1000000 Pa: 1350.0614 * pow(T - 103.8469, 3) - 48.6356 * pow(T - 103.8469, 2) + 3.7652 * (T - 103.8469) - 0.2203
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 150000;
    scalar P3 = 200000;
    scalar P4 = 400000;
    scalar P5 = 500000;
    scalar P6 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar cpValue = 0.0;

    if (p <= P1)
    {
        cpValue = 813.4101 * pow(T - 75.5049, 3) - 3.8299 * pow(T - 75.5049, 2) + 0.1524 * (T - 75.5049) - 0.0048;
    }
    else if (p <= P2)
    {
        scalar cpP1 = 813.4101 * pow(T - 75.5049, 3) - 3.8299 * pow(T - 75.5049, 2) + 0.1524 * (T - 75.5049) - 0.0048;
        scalar cpP2 = 855.0395 * pow(T - 80.9446, 3) - 6.0014 * pow(T - 80.9446, 2) + 0.2508 * (T - 80.9446) - 0.0086;
        scalar cpWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        cpValue = cpP1 + cpWeight * (cpP2 - cpP1);
    }
    else if (p <= P3)
    {
        scalar cpP2 = 855.0395 * pow(T - 80.9446, 3) - 6.0014 * pow(T - 80.9446, 2) + 0.2508 * (T - 80.9446) - 0.0086;
        scalar cpP3 = 882.5783 * pow(T - 83.7258, 3) - 7.5313 * pow(T - 83.7258, 2) + 0.3284 * (T - 83.7258) - 0.0120;
        scalar cpWeight = (p - P2) / (P3 - P2);  // Linear interpolation weight
        cpValue = cpP2 + cpWeight * (cpP3 - cpP2);
    }
    else if (p <= P4)
    {
        scalar cpP3 = 882.5783 * pow(T - 83.7258, 3) - 7.5313 * pow(T - 83.7258, 2) + 0.3284 * (T - 83.7258) - 0.0120;
        scalar cpP4 = 988.2877 * pow(T - 91.3327, 3) - 14.2465 * pow(T - 91.3327, 2) + 0.7310 * (T - 91.3327) - 0.0314;
        scalar cpWeight = (p - P3) / (P4 - P3);  // Linear interpolation weight
        cpValue = cpP3 + cpWeight * (cpP4 - cpP3);
    }
    else if (p <= P5)
    {
        scalar cpP4 = 988.2877 * pow(T - 91.3327, 3) - 14.2465 * pow(T - 91.3327, 2) + 0.7310 * (T - 91.3327) - 0.0314;
        scalar cpP5 = 1041.8387 * pow(T - 94.0950, 3) - 18.1936 * pow(T - 94.0950, 2) + 1.0047 * (T - 94.0950) - 0.0459;
        scalar cpWeight = (p - P4) / (P5 - P4);  // Linear interpolation weight
        cpValue = cpP4 + cpWeight * (cpP5 - cpP4);
    }    
    else if (p <= P6)
    {
        scalar cpP5 = 1041.8387 * pow(T - 94.0950, 3) - 18.1936 * pow(T - 94.0950, 2) + 1.0047 * (T - 94.0950) - 0.0459;
        scalar cpP6 = 1350.0614 * pow(T - 103.8469, 3) - 48.6356 * pow(T - 103.8469, 2) + 3.7652 * (T - 103.8469) - 0.2203;
        scalar cpWeight = (p - P5) / (P6 - P5);  // Linear interpolation weight
        cpValue = cpP5 + cpWeight * (cpP6 - cpP5);
    }
    else
    {
        cpValue = 1350.0614 * pow(T - 103.8469, 3) - 48.6356 * pow(T - 103.8469, 2) + 3.7652 * (T - 103.8469) - 0.2203;
    }

    return cpValue;
}



template<class Thermo>
inline Foam::scalar Foam::H2_gTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return this->Cp(p, T)*mu(p, T)*rPr_;
}


template<class Thermo>
inline Foam::scalar Foam::H2_gTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return mu(p, T)*rPr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::H2_gTransport<Thermo>::operator+=
(
    const H2_gTransport<Thermo>& st
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
inline void Foam::H2_gTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::H2_gTransport<Thermo> Foam::operator+
(
    const H2_gTransport<Thermo>& ct1,
    const H2_gTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return H2_gTransport<Thermo>
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

        return H2_gTransport<Thermo>
        (
            t,
            1.0/(Y1/ct1.rPr_ + Y2/ct2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::H2_gTransport<Thermo> Foam::operator*
(
    const scalar s,
    const H2_gTransport<Thermo>& ct
)
{
    return H2_gTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        1.0/ct.rPr_
    );
}


// ************************************************************************* //

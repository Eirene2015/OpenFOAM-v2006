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

    //- Fitted Curve for P = 80000 Pa: VISCOSITY*1e-6 = 23.4844 * Foam::pow(T - 14.5942, 3) + -2.7867 * Foam::pow(T - 14.5942, 2) + 0.2981 * (T - 14.5942) + -0.0315
    //- Fitted Curve for P = 1000000 Pa: VISCOSITY*1e-6 = 9.0637 * Foam::pow(T - 26.3925, 3) + -0.6848 * Foam::pow(T - 26.3925, 2) + 0.0152 * (T - 26.3925) + -0.0028    
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar muValue = 0.0;

    if (p <= P1)
    {
        muValue = ( 23.4844 * Foam::pow(T - 14.5942, 3) - 2.7867 * Foam::pow(T - 14.5942, 2) + 0.2981 * (T - 14.5942) - 0.0315 ) * 1e-6;
    }
    else if (p <= P2)
    {
        scalar muP1 = ( 23.4844 * Foam::pow(T - 14.5942, 3) - 2.7867 * Foam::pow(T - 14.5942, 2) + 0.2981 * (T - 14.5942) - 0.0315 ) * 1e-6;
        scalar muP2 = ( 9.0637 * Foam::pow(T - 26.3925, 3) - 0.6848 * Foam::pow(T - 26.3925, 2) + 0.0152 * (T - 26.3925) - 0.0028 ) * 1e-6;
        scalar muWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        muValue = muP1 + muWeight * (muP2 - muP1);
    }
    else
    {
        muValue = ( 9.0637 * Foam::pow(T - 26.3925, 3) - 0.6848 * Foam::pow(T - 26.3925, 2) + 0.0152 * (T - 26.3925) - 0.0028 ) * 1e-6;
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

    //- Fitted Curve P = 80000 Pa: 7189.9311 * Foam::pow(T - 14.5942, 3) + 330.0995 * Foam::pow(T - 14.5942, 2) + 32.4909 * (T - 14.5942) - 3.4206;
    //- Fitted Curve P = 1000000 Pa: 13732.7309 * Foam::pow(T - 26.3925, 3) + 1148.9115 * Foam::pow(T - 26.3925, 2) + 127.9205 * (T - 26.3925) + 35.0688;
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar cpValue = 0.0;

    if (p <= P1)
    {
        cpValue = 7189.9311 * Foam::pow(T - 14.5942, 3) + 330.0995 * Foam::pow(T - 14.5942, 2) + 32.4909 * (T - 14.5942) - 3.4206;
    }
    else if (p <= P2)
    {
        scalar cpP1 = 7189.9311 * Foam::pow(T - 14.5942, 3) + 330.0995 * Foam::pow(T - 14.5942, 2) + 32.4909 * (T - 14.5942) - 3.4206;
        scalar cpP2 = 13732.7309 * Foam::pow(T - 26.3925, 3) + 1148.9115 * Foam::pow(T - 26.3925, 2) + 127.9205 * (T - 26.3925) + 35.0688;
        scalar cpWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        cpValue = cpP1 + cpWeight * (cpP2 - cpP1);
    }
    else
    {
        cpValue = 13732.7309 * Foam::pow(T - 26.3925, 3) + 1148.9115 * Foam::pow(T - 26.3925, 2) + 127.9205 * (T - 26.3925) + 35.0688;
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


    //- Fitted Curve P = 80000 Pa: 6893.1445 * Foam::pow(T - 14.5942, 3) + 330.0995 * Foam::pow(T - 14.5942, 2) + 32.4909 * (T - 14.5942) - 3.4206;
    //- Fitted Curve P = 1000000 Pa: 13435.9443 * Foam::pow(T - 26.3925, 3) + 1148.9115 * Foam::pow(T - 26.3925, 2) + 127.9205 * (T - 26.3925) + 35.0688;
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar cvValue = 0.0;

    if (p <= P1)
    {
        cvValue = 6893.1445 * Foam::pow(T - 14.5942, 3) + 330.0995 * Foam::pow(T - 14.5942, 2) + 32.4909 * (T - 14.5942) - 3.4206;
    }
    else if (p <= P2)
    {
        scalar cvP1 = 6893.1445 * Foam::pow(T - 14.5942, 3) + 330.0995 * Foam::pow(T - 14.5942, 2) + 32.4909 * (T - 14.5942) - 3.4206;
        scalar cvP2 = 13435.9443 * Foam::pow(T - 26.3925, 3) + 1148.9115 * Foam::pow(T - 26.3925, 2) + 127.9205 * (T - 26.3925) + 35.0688;
        scalar cvWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        cvValue = cvP1 + cvWeight * (cvP2 - cvP1);
    }
    else
    {
        cvValue = 13435.9443 * Foam::pow(T - 26.3925, 3) + 1148.9115 * Foam::pow(T - 26.3925, 2) + 127.9205 * (T - 26.3925) + 35.0688;
    }

    return cvValue;
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

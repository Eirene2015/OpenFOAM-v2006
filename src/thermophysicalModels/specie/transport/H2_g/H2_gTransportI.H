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

    //- Fitted Curve for P = 80000 Pa: VISCOSITY*1e-6 = 0.9660 * Foam::pow(T - 19.6942, 3) + 0.0631 * Foam::pow(T - 19.6942, 2) + -0.0014 * (T - 19.6942) + 0.0001
    //- Fitted Curve for P = 1000000 Pa: VISCOSITY*1e-6 = 1.9235 * Foam::pow(T - 31.4925, 3) + -0.0300 * Foam::pow(T - 31.4925, 2) + 0.0249 * (T - 31.4925) + -0.0036
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar muValue = 0.0;

    if (p <= P1)
    {
        muValue = (0.9660 * Foam::pow(T - 19.6942, 3) + 0.0631 * Foam::pow(T - 19.6942, 2) - 0.0014 * (T - 19.6942) + 0.0001)*1e-6;
    }
    else if (p <= P2)
    {
        scalar muP1 = ( 0.9660 * Foam::pow(T - 19.6942, 3) + 0.0631 * Foam::pow(T - 19.6942, 2) - 0.0014 * (T - 19.6942) + 0.0001 )*1e-6;
        scalar muP2 = ( 1.9235 * Foam::pow(T - 31.4925, 3) - 0.0300 * Foam::pow(T - 31.4925, 2) + 0.0249 * (T - 31.4925) - 0.0036 ) * 1e-6;
        scalar muWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        muValue = muP1 + muWeight * (muP2 - muP1);
    }
    else
    {
        muValue = ( 1.9235 * Foam::pow(T - 31.4925, 3) - 0.0300 * Foam::pow(T - 31.4925, 2) + 0.0249 * (T - 31.4925) - 0.0036 ) * 1e-6;
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

    //- Fitted Curve P = 80000 Pa: 11713.8604 * Foam::pow(T - 19.6942, 3) - 300.7934 * Foam::pow(T - 19.6942, 2) + 50.4396 * (T - 19.6942) - 4.6363;
    //- Fitted Curve P = 150000 Pa: 12651.7286 * Foam::pow(T - 21.8744, 3) - 498.0551 * Foam::pow(T - 21.8744, 2) + 85.2330 * (T - 21.8744) - 8.1271;
    //- Fitted Curve P = 200000 Pa: 13339.7833 * Foam::pow(T - 23.0104, 3) - 661.8079 * Foam::pow(T - 23.0104, 2) + 116.5643 * (T - 23.0104) - 11.3971;
    //- Fitted Curve P = 400000 Pa: 16542.8938 * Foam::pow(T - 26.1764, 3) - 1671.4210 * Foam::pow(T - 26.1764, 2) + 354.4450 * (T - 26.1764) - 39.3229;
    //- Fitted Curve P = 500000 Pa: 18600.3126 * Foam::pow(T - 27.3428, 3) - 2547.0617 * Foam::pow(T - 27.3428, 2) + 605.1176 * (T - 27.3428) - 71.9698;
    //- Fitted Curve P = 1000000 Pa: 45684.1755 * Foam::pow(T - 31.4925, 3) - 28006.4291 * Foam::pow(T - 31.4925, 2) + 12054.7135 * (T - 31.4925) - 1955.8585;
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
        cpValue = 11713.8604 * Foam::pow(T - 19.6942, 3) - 300.7934 * Foam::pow(T - 19.6942, 2) + 50.4396 * (T - 19.6942) - 4.6363;
    }
    else if (p <= P2)
    {
        scalar cpP1 = 11713.8604 * Foam::pow(T - 19.6942, 3) - 300.7934 * Foam::pow(T - 19.6942, 2) + 50.4396 * (T - 19.6942) - 4.6363;
        scalar cpP2 = 12651.7286 * Foam::pow(T - 21.8744, 3) - 498.0551 * Foam::pow(T - 21.8744, 2) + 85.2330 * (T - 21.8744) - 8.1271;
        scalar cpWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        cpValue = cpP1 + cpWeight * (cpP2 - cpP1);
    }
    else if (p <= P3)
    {
        scalar cpP2 = 12651.7286 * Foam::pow(T - 21.8744, 3) - 498.0551 * Foam::pow(T - 21.8744, 2) + 85.2330 * (T - 21.8744) - 8.1271;
        scalar cpP3 = 13339.7833 * Foam::pow(T - 23.0104, 3) - 661.8079 * Foam::pow(T - 23.0104, 2) + 116.5643 * (T - 23.0104) - 11.3971;
        scalar cpWeight = (p - P2) / (P3 - P2);  // Linear interpolation weight
        cpValue = cpP2 + cpWeight * (cpP3 - cpP2);
    }
    else if (p <= P4)
    {
        scalar cpP3 = 13339.7833 * Foam::pow(T - 23.0104, 3) - 661.8079 * Foam::pow(T - 23.0104, 2) + 116.5643 * (T - 23.0104) - 11.3971;
        scalar cpP4 = 16542.8938 * Foam::pow(T - 26.1764, 3) - 1671.4210 * Foam::pow(T - 26.1764, 2) + 354.4450 * (T - 26.1764) - 39.3229;
        scalar cpWeight = (p - P3) / (P4 - P3);  // Linear interpolation weight
        cpValue = cpP3 + cpWeight * (cpP4 - cpP3);
    }
    else if (p <= P5)
    {
        scalar cpP4 = 16542.8938 * Foam::pow(T - 26.1764, 3) - 1671.4210 * Foam::pow(T - 26.1764, 2) + 354.4450 * (T - 26.1764) - 39.3229;
        scalar cpP5 = 18600.3126 * Foam::pow(T - 27.3428, 3) - 2547.0617 * Foam::pow(T - 27.3428, 2) + 605.1176 * (T - 27.3428) - 71.9698;
        scalar cpWeight = (p - P4) / (P5 - P4);  // Linear interpolation weight
        cpValue = cpP4 + cpWeight * (cpP5 - cpP4);
    }    
    else if (p <= P6)
    {
        scalar cpP5 = 18600.3126 * Foam::pow(T - 27.3428, 3) - 2547.0617 * Foam::pow(T - 27.3428, 2) + 605.1176 * (T - 27.3428) - 71.9698;
        scalar cpP6 = 45684.1755 * Foam::pow(T - 31.4925, 3) - 28006.4291 * Foam::pow(T - 31.4925, 2) + 12054.7135 * (T - 31.4925) - 1955.8585;
        scalar cpWeight = (p - P5) / (P6 - P5);  // Linear interpolation weight
        cpValue = cpP5 + cpWeight * (cpP6 - cpP5);
    }
    else
    {
        cpValue = 45684.1755 * Foam::pow(T - 31.4925, 3) - 28006.4291 * Foam::pow(T - 31.4925, 2) + 12054.7135 * (T - 31.4925) - 1955.8585;
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

    //- Fitted Curve P = 80000 Pa: 11417.0739 * Foam::pow(T - 19.6942, 3) - 300.7934 * Foam::pow(T - 19.6942, 2) + 50.4396 * (T - 19.6942) - 4.6363;
    //- Fitted Curve P = 150000 Pa: 12354.9421 * Foam::pow(T - 21.8744, 3) - 498.0551 * Foam::pow(T - 21.8744, 2) + 85.2330 * (T - 21.8744) - 8.1271;
    //- Fitted Curve P = 200000 Pa: 13042.9967 * Foam::pow(T - 23.0104, 3) - 661.8079 * Foam::pow(T - 23.0104, 2) + 116.5643 * (T - 23.0104) - 11.3971;
    //- Fitted Curve P = 400000 Pa: 16246.1073 * Foam::pow(T - 26.1764, 3) - 1671.4210 * Foam::pow(T - 26.1764, 2) + 354.4450 * (T - 26.1764) - 39.3229;
    //- Fitted Curve P = 500000 Pa: 18303.5261 * Foam::pow(T - 27.3428, 3) - 2547.0617 * Foam::pow(T - 27.3428, 2) + 605.1176 * (T - 27.3428) - 71.9698;
    //- Fitted Curve P = 1000000 Pa: 45387.3889 * Foam::pow(T - 31.4925, 3) - 28006.4291 * Foam::pow(T - 31.4925, 2) + 12054.7135 * (T - 31.4925) - 1955.8585;
    // Define the pressure values and corresponding viscosity equations

    scalar P1 = 80000;
    scalar P2 = 150000;
    scalar P3 = 200000;
    scalar P4 = 400000;
    scalar P5 = 500000;
    scalar P6 = 1000000;

    // Calculate mu based on pressure (p) and temperature (T) using linear interpolation
    scalar cvValue = 0.0;

    if (p <= P1)
    {
        cvValue = 11417.0739 * Foam::pow(T - 19.6942, 3) - 300.7934 * Foam::pow(T - 19.6942, 2) + 50.4396 * (T - 19.6942) - 4.6363;
    }
    else if (p <= P2)
    {
        scalar cvP1 = 11417.0739 * Foam::pow(T - 19.6942, 3) - 300.7934 * Foam::pow(T - 19.6942, 2) + 50.4396 * (T - 19.6942) - 4.6363;
        scalar cvP2 = 12354.9421 * Foam::pow(T - 21.8744, 3) - 498.0551 * Foam::pow(T - 21.8744, 2) + 85.2330 * (T - 21.8744) - 8.1271;
        scalar cvWeight = (p - P1) / (P2 - P1);  // Linear interpolation weight
        cvValue = cvP1 + cvWeight * (cvP2 - cvP1);
    }
    else if (p <= P3)
    {
        scalar cvP2 = 12354.9421 * Foam::pow(T - 21.8744, 3) - 498.0551 * Foam::pow(T - 21.8744, 2) + 85.2330 * (T - 21.8744) - 8.1271;
        scalar cvP3 = 13042.9967 * Foam::pow(T - 23.0104, 3) - 661.8079 * Foam::pow(T - 23.0104, 2) + 116.5643 * (T - 23.0104) - 11.3971;
        scalar cvWeight = (p - P2) / (P3 - P2);  // Linear interpolation weight
        cvValue = cvP2 + cvWeight * (cvP3 - cvP2);
    }
    else if (p <= P4)
    {
        scalar cvP3 = 13042.9967 * Foam::pow(T - 23.0104, 3) - 661.8079 * Foam::pow(T - 23.0104, 2) + 116.5643 * (T - 23.0104) - 11.3971;
        scalar cvP4 = 16246.1073 * Foam::pow(T - 26.1764, 3) - 1671.4210 * Foam::pow(T - 26.1764, 2) + 354.4450 * (T - 26.1764) - 39.3229;
        scalar cvWeight = (p - P3) / (P4 - P3);  // Linear interpolation weight
        cvValue = cvP3 + cvWeight * (cvP4 - cvP3);
    }
    else if (p <= P5)
    {
        scalar cvP4 = 16246.1073 * Foam::pow(T - 26.1764, 3) - 1671.4210 * Foam::pow(T - 26.1764, 2) + 354.4450 * (T - 26.1764) - 39.3229;
        scalar cvP5 = 18303.5261 * Foam::pow(T - 27.3428, 3) - 2547.0617 * Foam::pow(T - 27.3428, 2) + 605.1176 * (T - 27.3428) - 71.9698;
        scalar cvWeight = (p - P4) / (P5 - P4);  // Linear interpolation weight
        cvValue = cvP4 + cvWeight * (cvP5 - cvP4);
    }    
    else if (p <= P6)
    {
        scalar cvP5 = 18303.5261 * Foam::pow(T - 27.3428, 3) - 2547.0617 * Foam::pow(T - 27.3428, 2) + 605.1176 * (T - 27.3428) - 71.9698;
        scalar cvP6 = 45387.3889 * Foam::pow(T - 31.4925, 3) - 28006.4291 * Foam::pow(T - 31.4925, 2) + 12054.7135 * (T - 31.4925) - 1955.8585;
        scalar cvWeight = (p - P5) / (P6 - P5);  // Linear interpolation weight
        cvValue = cvP5 + cvWeight * (cvP6 - cvP5);
    }
    else
    {
        cvValue = 45387.3889 * Foam::pow(T - 31.4925, 3) - 28006.4291 * Foam::pow(T - 31.4925, 2) + 12054.7135 * (T - 31.4925) - 1955.8585;
    }

    return cvValue;
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

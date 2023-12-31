/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "dictionary.H"
#include "primitiveEntry.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::dictionary::reportDefault
(
    const word& keyword,
    const T& deflt,
    const bool added
) const
{
    InfoErr
        << "Dictionary: " << relativeName(true).c_str()
        << " Entry: " << keyword;

    if (added)
    {
        InfoErr
            << " Added";
    }
    InfoErr
        << " Default: " << deflt << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Compare>
Foam::wordList Foam::dictionary::sortedToc(const Compare& comp) const
{
    return hashedEntries_.sortedToc(comp);
}


template<class T>
Foam::entry* Foam::dictionary::add(const keyType& k, const T& v, bool overwrite)
{
    return add(new primitiveEntry(k, v), overwrite);
}


template<class T>
Foam::entry* Foam::dictionary::set(const keyType& k, const T& v)
{
    return set(new primitiveEntry(k, v));
}


template<class T>
T Foam::dictionary::get
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    T val;
    readEntry<T>(keyword, val, matchOpt);
    return val;
}


template<class T, class Predicate>
T Foam::dictionary::getCheck
(
    const word& keyword,
    const Predicate& pred,
    enum keyType::option matchOpt
) const
{
    T val;
    readCheck<T, Predicate>(keyword, val, pred, matchOpt);
    return val;
}


template<class T>
T Foam::dictionary::getCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    T val;
    readCompat<T>(keyword, compat, val, matchOpt);
    return val;
}


template<class T>
T Foam::dictionary::getOrDefault
(
    const word& keyword,
    const T& deflt,
    enum keyType::option matchOpt
) const
{
    const const_searcher finder(csearch(keyword, matchOpt));

    if (finder.good())
    {
        T val;

        ITstream& is = finder.ptr()->stream();
        is >> val;

        checkITstream(is, keyword);

        return val;
    }
    else if (writeOptionalEntries)
    {
        if (writeOptionalEntries > 1)
        {
            FatalIOErrorInFunction(*this)
                << "No optional entry: " << keyword
                << " Default: " << deflt << nl
                << exit(FatalIOError);
        }
        else
        {
            reportDefault(keyword, deflt);
        }
    }

    return deflt;
}


template<class T>
T Foam::dictionary::getOrAdd
(
    const word& keyword,
    const T& deflt,
    enum keyType::option matchOpt
)
{
    const const_searcher finder(csearch(keyword, matchOpt));

    if (finder.good())
    {
        T val;

        ITstream& is = finder.ptr()->stream();
        is >> val;

        checkITstream(is, keyword);

        return val;
    }
    else if (writeOptionalEntries)
    {
        if (writeOptionalEntries > 1)
        {
            FatalIOErrorInFunction(*this)
                << "No optional entry: " << keyword
                << " Default: " << deflt << nl
                << exit(FatalIOError);
        }
        else
        {
            reportDefault(keyword, deflt, true);
        }
    }

    add(new primitiveEntry(keyword, deflt));
    return deflt;
}


template<class T, class Predicate>
T Foam::dictionary::getCheckOrDefault
(
    const word& keyword,
    const T& deflt,
    const Predicate& pred,
    enum keyType::option matchOpt
) const
{
    if (!pred(deflt))
    {
        // Could be as FULLDEBUG instead?
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' with invalid default in dictionary "
            << name()
            << exit(FatalIOError);
    }

    const const_searcher finder(csearch(keyword, matchOpt));

    if (finder.good())
    {
        T val;

        ITstream& is = finder.ptr()->stream();
        is >> val;

        checkITstream(is, keyword);

        if (!pred(val))
        {
            raiseBadInput(is, keyword);
        }

        return val;
    }
    else if (writeOptionalEntries)
    {
        if (writeOptionalEntries > 1)
        {
            FatalIOErrorInFunction(*this)
                << "No optional entry: " << keyword
                << " Default: " << deflt << nl
                << exit(FatalIOError);
        }
        else
        {
            reportDefault(keyword, deflt);
        }
    }

    return deflt;
}


template<class T, class Predicate>
T Foam::dictionary::getCheckOrAdd
(
    const word& keyword,
    const T& deflt,
    const Predicate& pred,
    enum keyType::option matchOpt
)
{
    if (!pred(deflt))
    {
        // Could be as FULLDEBUG instead?
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' with invalid default in dictionary "
            << name()
            << exit(FatalIOError);
    }

    const const_searcher finder(csearch(keyword, matchOpt));

    if (finder.good())
    {
        T val;

        ITstream& is = finder.ptr()->stream();
        is >> val;

        checkITstream(is, keyword);

        if (!pred(val))
        {
            raiseBadInput(is, keyword);
        }

        return val;
    }
    else if (writeOptionalEntries)
    {
        if (writeOptionalEntries > 1)
        {
            FatalIOErrorInFunction(*this)
                << "No optional entry: " << keyword
                << " Default: " << deflt << nl
                << exit(FatalIOError);
        }
        else
        {
            reportDefault(keyword, deflt, true);
        }
    }

    add(new primitiveEntry(keyword, deflt));
    return deflt;
}


template<class T>
bool Foam::dictionary::readEntry
(
    const word& keyword,
    T& val,
    enum keyType::option matchOpt,
    bool mandatory
) const
{
    const const_searcher finder(csearch(keyword, matchOpt));

    if (finder.good())
    {
        ITstream& is = finder.ptr()->stream();
        is >> val;

        checkITstream(is, keyword);

        return true;
    }
    else if (mandatory)
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << name() << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class T, class Predicate>
bool Foam::dictionary::readCheck
(
    const word& keyword,
    T& val,
    const Predicate& pred,
    enum keyType::option matchOpt,
    bool mandatory
) const
{
    const const_searcher finder(csearch(keyword, matchOpt));

    if (finder.good())
    {
        ITstream& is = finder.ptr()->stream();
        is >> val;

        checkITstream(is, keyword);

        if (!pred(val))
        {
            raiseBadInput(is, keyword);
        }

        return true;
    }
    else if (mandatory)
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << name() << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class T>
bool Foam::dictionary::readCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    T& val,
    enum keyType::option matchOpt,
    bool mandatory
) const
{
    const const_searcher finder(csearchCompat(keyword, compat, matchOpt));

    if (finder.good())
    {
        ITstream& is = finder.ptr()->stream();
        is >> val;

        checkITstream(is, keyword);

        return true;
    }
    else if (mandatory)
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << name() << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class T>
bool Foam::dictionary::readIfPresent
(
    const word& keyword,
    T& val,
    enum keyType::option matchOpt
) const
{
    // Read is non-mandatory
    return readEntry<T>(keyword, val, matchOpt, false);
}


template<class T, class Predicate>
bool Foam::dictionary::readCheckIfPresent
(
    const word& keyword,
    T& val,
    const Predicate& pred,
    enum keyType::option matchOpt
) const
{
    // Read is non-mandatory
    return readCheck<T, Predicate>(keyword, val, pred, matchOpt, false);
}


template<class T>
T Foam::dictionary::getOrDefaultCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    const T& deflt,
    enum keyType::option matchOpt
) const
{
    const const_searcher finder(csearchCompat(keyword, compat, matchOpt));

    if (finder.good())
    {
        T val;

        ITstream& is = finder.ptr()->stream();
        is >> val;

        checkITstream(is, keyword);

        return val;
    }
    else if (writeOptionalEntries)
    {
        if (writeOptionalEntries > 1)
        {
            FatalIOErrorInFunction(*this)
                << "No optional entry: " << keyword
                << " Default: " << deflt << nl
                << exit(FatalIOError);
        }
        else
        {
            reportDefault(keyword, deflt);
        }
    }

    return deflt;
}


template<class T>
bool Foam::dictionary::readIfPresentCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    T& val,
    enum keyType::option matchOpt
) const
{
    // Read is non-mandatory
    return readCompat<T>(keyword, compat, val, matchOpt, false);
}


// ************************************************************************* //

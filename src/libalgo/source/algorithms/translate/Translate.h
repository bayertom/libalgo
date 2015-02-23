// Description: Arithmetic parser using postfix notation

// Copyright (c) 2010 - 2011
// Tomas Bayer
// Charles University in Prague, Faculty of Science
// bayertom@natur.cuni.cz

// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#ifndef Translate_H
#define Translate_H

#include <ostream>
#include <iostream>

//Arithmetic parser converting equation from infix to postfix notation
class Translate
{
        public:
                template <typename T>
                static T translateEq ( const char * equation, const T x, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename T>
                static T translateEq ( const char * equation, const T x, const T y, const bool print_exception = true, std::ostream * output = &std::cout );

                template <typename T>
                static T translateEq ( const char * equation, const T lat, const T lon, const T R, const T a, const T b,  const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exception = true, std::ostream * output = &std::cout );

        private:
                static void infixToPostfix ( const char * infix, char * postfix );

                template <typename T>
                static T parseEquation ( const char * equation, const T x, const T y, const T lat, const T lon, const T R, const T a, const T b,  const T lat0, const T lat1, const T lat2, const T lon0 );

                static void findSequence ( const char * equation, char * operator_text );

};

#include "Translate.hpp"

#endif

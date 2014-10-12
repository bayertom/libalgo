// Description: Ambiguous value error

// Copyright (c) 2010 - 2013
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


#ifndef ErrorMathAmbiguousValue_H
#define ErrorMathAmbiguousValue_H


#include "ErrorMath.h"


//Math error: ambiguous value
template <typename T>
class ErrorMathAmbiguousValue : public ErrorMath <T>
{
        protected:
                T arg;

        public:
                ErrorMathAmbiguousValue <T> ( const char * exception_text_, const char * function_text_, const T arg_ ) :
                        ErrorMath <T> ( exception_text_, function_text_ ), arg ( arg_ ) { }

        public:
                virtual ~ErrorMathAmbiguousValue <T> () throw() {};

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Error::printException ( output );
                        *output << arg << '\n';
                }

                virtual T getArg() const {return arg;}
                virtual short getErrorCode() const { return 10;}
};

#endif

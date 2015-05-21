// Description: Math underflow error class

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


#ifndef ErrorMathUnderflow_H
#define ErrorMathUnderflow_H


#include "ErrorMath.h"


//Math error: number underflow
template <typename T>
class ErrorMathUnderflow : public ErrorMath <T>
{
        protected:
                T arg;

        public:
                ErrorMathUnderflow ( const char * exception_text_, const char * function_text_ , const T &arg_ ) :
                        ErrorMath <T> ( exception_text_, function_text_ ), arg ( arg_ ) {};

        public:
                virtual ~ErrorMathUnderflow() throw() {};

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        ErrorMath <T>::printException ( output );
                        *output << arg << '\n';
                }

                virtual T getArg() const {return arg;}
                virtual short getErrorCode() const { return 18;}
};

#endif

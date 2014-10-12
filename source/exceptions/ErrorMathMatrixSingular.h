// Description: Singular matrix error class

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


#ifndef ErrorMathMatrixSingular_H
#define ErrorMathMatrixSingular_H


#include "ErrorMathMatrix.h"


//Matrix error: singular matrix
template <typename TMatrix>
class ErrorMathMatrixSingular : public ErrorMathMatrix <TMatrix>
{
        public:
                ErrorMathMatrixSingular ( const char * exception_text_, const char * function_text_, const TMatrix &M_ ) :
                        ErrorMathMatrix <TMatrix> ( exception_text_, function_text_, M_ ) {}

        public:
                virtual ~ErrorMathMatrixSingular() throw() {};

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        ErrorMathMatrix<TMatrix>::printException ( output );
                }

                virtual TMatrix getArg() const { return ErrorMathMatrix <TMatrix>::getArg(); }
                virtual short getErrorCode() const { return 14;}
};

#endif

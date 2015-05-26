// Description: Matrix error class, other matrix error classes derived from this class

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


#ifndef ErrorMathMatrix_H
#define ErrorMathMatrix_H


#include "libalgo/src/structures/matrix/Matrix.h"

#include "ErrorMath.h"


//Math matrix error
template <typename TMatrix>
class ErrorMathMatrix : public ErrorMath <TMatrix>
{
        protected:
                TMatrix M;

        public:

                ErrorMathMatrix ( const char * exception_text_, const char * function_text_,  const TMatrix &M_ ) :
                        ErrorMath <TMatrix> ( exception_text_, function_text_ ), M ( M_ ) {}

        public:
                virtual ~ErrorMathMatrix() throw() {};

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        ErrorMath <TMatrix>::printException ( output );
                        *output << "Matrix A, rows count: " << M.rows() << ", cols count: " << M.cols() << '\n';
                        M.print ( output );
                }

                virtual TMatrix getArg( ) const {return M;}
                virtual short getErrorCode() const { return 12;}

};

#endif

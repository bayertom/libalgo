// Description: Downhill simplex optimization method

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


#ifndef SimplexMethod_H
#define SimplexMethod_H

#include "libalgo/src/algorithms/matrixoperations/MatrixOperations.h"

//Downhill simplex optimization method
class SimplexMethod
{

        public:

                template <typename T, typename Function>
                static T NelderMead ( Function function, const Matrix <T> &XMIN, const Matrix <T> &XMAX, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, unsigned int &iterations, const T max_error, const unsigned int max_iterations, std::ostream * output = &std::cout );

        private:

                template <typename T>
                static Matrix <T> createRandSimplex ( const Matrix <T> &XMIN, const Matrix <T> &XMAX );


                template <typename T, typename Function>
                static void shrink ( Function function, Matrix <T> &W, Matrix <T> &X, const Matrix<T> &XMIN, const Matrix <T> &XMAX, Matrix <T> &Y,  Matrix <T> &V, const T Sigma );

                template <typename T>
                static void reflection ( const Matrix <T> &XMIN, const Matrix <T> &XMAX, const unsigned int dim, Matrix <T> &X );

};

#include "SimplexMethod.hpp"

#endif

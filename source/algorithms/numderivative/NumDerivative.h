// Description: Compute numeric derivative using Stirling method, function is defined using functors

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


#ifndef NumDerivative_H
#define NumDerivative_H


#include "libalgo/source/structures/matrix/Matrix.h"

/*
//Set derived equation
typedef enum
{
	FunctionF1 = 0,			//Drivative of the first coordinate function
	FunctionF2,			//Derivative of the second coordinate function
} TDerivativeFunction;
*/

//Set type of the derivative
typedef enum
{
	FirstDerivative = 0,
	SecondDerivative
}TDerivativeType;


//Set variable for the derivarive
typedef enum
{
        VariableX1 = 0,		//Derivative according to the first variable
        VariableX2,
        VariableX3,
        VariableX4,
        VariableX5,
        VariableX6,
        VariableX7,
        VariableX8,			//Derivative according to the last variable
} TDerivativeVariable;



//Compute numerical derivative using Stirling method
class NumDerivative
{
        public:

                template <typename T, typename Function>
		static T getDerivative(Function function, const Matrix <T> &args, const TDerivativeType deriv_type, const TDerivativeVariable deriv_var, const T deriv_step, const bool print_exceptions = false);


        private:

                template <typename T, typename Function>
		static void computeFunctionValues(Function function, const Matrix <T> &args, const TDerivativeVariable deriv_var, const T deriv_step, T * values, T * fvalues, const bool print_exceptions = false);


                template <typename T>
		static T computeStirlingFormula(T * values, T * fvalues, const TDerivativeType deriv_type, const T deriv_step);

};

#include "NumDerivative.hpp"

#endif

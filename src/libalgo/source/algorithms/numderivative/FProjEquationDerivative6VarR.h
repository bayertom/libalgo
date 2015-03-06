// Description: Functor, compute numeric derivative of the projection equation using the Stirling formula
// Derivatives are: R, latp, lonp, lat0, lon0 (used for NLPS solution)

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



#ifndef FProjEquationDerivative6VarR_H
#define FProjEquationDerivative6VarR_H


#include "libalgo/source/structures/matrix/Matrix.h"

#include "libalgo/source/algorithms/carttransformation/CartTransformation.h"


//Functor, performs partial derivative of the projection equation using the Stirling formula
// Derivatives are: R, latp, lonp, lat0, lon0, c
template <typename T>
class FProjEquationDerivative6VarR
{

private:
	//Map projection parameters
	const char * equation;
	const T lat;
	const T lon;
	const T a;
	const T b;
	const T lat1;
	const T lat2;
	const TTransformedLongtitudeDirection trans_lon_dir;


public:

	FProjEquationDerivative6VarR(const char * equation_, const T lat_, const T lon_, const T a_, const T b_, const T lat1_, const T lat2_, const TTransformedLongtitudeDirection trans_lon_dir_) :
		equation(equation_), lat(lat_), lon(lon_), a(a_), b(b_), lat1(lat1_), lat2(lat2_), trans_lon_dir(trans_lon_dir_) {}

	T operator () (const Matrix <T> &arg)
	{
		//Reduce lon
		const T lon_red = CartTransformation::redLon0(lon, arg(0, 4) * 180 / M_PI);

		//Convert ( lat, lon ) -> ( lat, lon)_trans
		const T lat_trans = CartTransformation::latToLatTrans(lat, lon_red, arg(0, 1) * 180 / M_PI, arg(0, 2) * 180 / M_PI);
		const T lon_trans = CartTransformation::lonToLonTrans(lat, lon_red, lat_trans, arg(0, 1) * 180 / M_PI, arg(0, 2) * 180 / M_PI, trans_lon_dir);

		//Compute partial derivative of the map projection equation
		T res = ArithmeticParser::parseEq(equation, lat_trans, lon_trans, arg(0, 0), a, b, arg(0, 5), arg(0, 3) * 180 / M_PI, lat1, lat2, false);
		return res;
	}
};

#endif

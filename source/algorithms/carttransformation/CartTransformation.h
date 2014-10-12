// Description: transformation of the geographic coordinates (lat, lon) to (lat, lon)_trans. Conversion (lat, lon) to (x, y).

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


#ifndef CartTransformation_H
#define CartTransformation_H

#include <math.h>

#include "libalgo/source/const/Const.h"

//Forward declaration
template <typename T>
class Point3DCartesian;

template <typename T>
class Point3DGeographic;

template <typename T>
class Projection;

//Set a direction of the transformed longitude: positive value in east (normal) or west (reversed) direction
//from meridian passing through the cartographic pole
typedef enum
{
        NormalDirection = 1,	//Azimuthal, conic projections: positive values in east direction from meridian passing the pole
        ReversedDirection,	//Azimuthal, conic projections: positive values in west direction from meridian passing the pole
        EquatorDirection,	//Cylindrical projections: positive values in north direction from cartographic equator
} TTransformedLongtitudeDirection;


//Basic cartographic transformations
class CartTransformation
{
        public:

                template <typename T>
                static T latToLatTrans ( const Point3DGeographic <T> *p, const Point3DGeographic <T> *pole );

                template <typename T>
                static T lonToLonTrans ( const Point3DGeographic <T> *p, const T lat_trans, const Point3DGeographic <T> *pole, const TTransformedLongtitudeDirection lon_direction = NormalDirection );

                //template <typename T>
                //static T latTransToLat ( const Point3DGeographic <T> *p, const Point3DGeographic <T> *pole );

                //template <typename T>
                //static T lonTransToLon ( const Point3DGeographic <T> *p, const T lat, const Point3DGeographic <T> *pole );

                template <typename T>
                static T latLonToX ( const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exception = true );

                template <typename T>
                static T latLonToY ( const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exception = true );

        public:
                template <typename T>
                static T redLon0 ( const T lon, const T lon0 ) { return ( lon - lon0 < MIN_LON ? 360.0 + ( lon - lon0 ) : ( lon - lon0 > MAX_LON ? ( lon - lon0 ) - 360 : lon - lon0 ) );}

                template <typename T>
                static T latToLatTrans ( const T lat, const T lon, const T latp, const T lonp );

                template <typename T>
                static T lonToLonTrans ( const T lat, const T lon, const T lat_trans, const T latp, const T lonp, const TTransformedLongtitudeDirection lon_direction = NormalDirection );

                //template <typename T>
                //static T latTransToLat ( const T lat_trans, const T lon_trans, const T latp, const T lonp );

                //template <typename T>
                //static T lonTransToLon ( const T lat_trans, const T lon_trans, const T lat, const T latp, const T lonp );

                template <typename T>
                static T latLonToX ( const char * equation_x,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T c, const T lat0, const T lat1, const T lat2, const bool print_exception = true );

                template <typename T>
                static T latLonToY ( const char * equation_x,  const T lat, const T lon, const T R, const T a, const T b, const T dy, const T c, const T lat0, const T lat1, const T lat2, const bool print_exception = true );

                template <typename T>
                static void wgsToJTSK ( const Point3DGeographic <T> *p1, Point3DCartesian <T> * p2 );
};

#include "CartTransformation.hpp"

#endif

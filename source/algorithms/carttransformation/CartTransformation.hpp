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


#ifndef CartTransformation_HPP
#define CartTransformation_HPP

#include <cmath>

#include "libalgo/source/const/Const.h"

#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/point/Point3DCartesianProjected.h"
#include "libalgo/source/structures/projection/Projection.h"

#include "libalgo/source/algorithms/arithmeticparser/ArithmeticParser.h"

#include "libalgo/source/exceptions/ErrorMathInvalidArgument.h"
#include "libalgo/source/exceptions/ErrorMathOverflow.h"
#include "libalgo/source/exceptions/ErrorParse.h"
#include "libalgo/source/exceptions/ErrorBadData.h"


template <typename T>
T CartTransformation::latToLatTrans ( const Point3DGeographic <T> *p, const Point3DGeographic <T> *pole )
{
        //Transform latitude  ( lat, lon ) -> ( lat_tans ) using a cartographic pole (latp, lonp)
        return latToLatTrans ( p->getLat(), p->getLon(), pole->getLat(), pole->getLon() );
}


template <typename T>
T CartTransformation::lonToLonTrans ( const Point3DGeographic <T> *p, const T lat_trans, const Point3DGeographic <T> *pole, const TTransformedLongtitudeDirection lon_direction )
{
        //Transform longitude  ( lat, lon, lat_trans ) -> ( lon_tans ) using a cartographic pole (latp, lonp)
        return lonToLonTrans ( p->getLat(), p->getLon(), lat_trans, pole->getLat(), pole->getLon(), lon_direction );
}


template <typename T>
T CartTransformation::latLonToX ( const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exception )
{
        //Compute x coordinate of the point P[lat, lon] in specific projection
        return latLonToX ( proj->getXEquat(),  p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDx(), proj->getC(),
                           proj->getLat0(), proj->getLat1(), proj->getLat2(), print_exception ) ;
}


template <typename T>
T CartTransformation::latLonToY ( const Point3DGeographic <T> *p, const Projection <T> *proj, const bool print_exception )
{
        //Compute y coordinate of the point P[lat, lon] in specific projection
        return latLonToY ( proj->getYEquat(),  p->getLat(), p->getLon(), proj->getR(), proj->getA(), proj->getB(), proj->getDy(), proj->getC(),
                           proj->getLat0(), proj->getLat1(), proj->getLat2(), print_exception ) ;
}


template <typename T>
T CartTransformation::latToLatTrans ( const T lat, const T lon, const T latp, const T lonp )
{
        //Transform latitude  ( lat, lon ) -> ( lat_tans ) using a cartographic pole (latp, lonp)

        //Throw exception: bad lat
        if ( fabs ( lat ) > MAX_LAT )
        {
                throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: ", "can not convert lat to lat_trans, lat > +- Pi/2", lat );
        }

        //Throw exception: bad lon
        if ( fabs ( lon ) > MAX_LON )
        {
                throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: ", "can not convert lon to lon_trans, lon > +- Pi", lon );
        }

        //Projection in normal aspect
        if ( fabs ( MAX_LAT - latp ) < ANGLE_ROUND_ERROR )
        {
                return lat;
        }

        //Projection in oblique aspect
        else
        {
                //Same coordinates as the cartographic pole, singular point
                if ( ( fabs ( lon - lonp ) < ANGLE_ROUND_ERROR ) && ( fabs ( lat - latp ) < ANGLE_ROUND_ERROR ) )
                {
                        return MAX_LAT;
                }

                //Compute latitude
                T lat_trans_asin = sin ( lat * M_PI / 180.0 ) * sin ( latp * M_PI / 180.0 ) + cos ( lat * M_PI / 180.0 ) * cos ( latp * M_PI / 180.0 ) * cos ( ( lonp - lon ) * M_PI / 180.0 );

                //Throw exception
                if ( ( lat_trans_asin > 1.0 + ARGUMENT_ROUND_ERROR )  || ( lat_trans_asin < -1.0 - ARGUMENT_ROUND_ERROR ) )
                {
                        throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: ", "can not convert lat to lat_trans, asin(arg), arg = ", lat_trans_asin );
                }

                //Correct latitude
                if ( lat_trans_asin > 1.0 )
                {
                        return MAX_LAT;
                }

                //Correct latitude
                if ( lat_trans_asin < -1.0 )
                {
                        return MIN_LAT;
                }

                T test = asin ( lat_trans_asin ) * 180.0 / 3.14;

                //Compute transformed latitude
                return  asin ( lat_trans_asin ) * 180.0 / M_PI;
        }
}


template <typename T>
T CartTransformation::lonToLonTrans ( const T lat, const T lon, const T lat_trans, const T latp, const T lonp, const TTransformedLongtitudeDirection lon_direction )
{
        //Transform longitude  ( lat, lon, lat_trans ) -> ( lon_tans ) using a cartographic pole (latp, lonp)
        //Normal direction (positive value in east direction) or reversed  direction (positive value in west direction, suitable for JTSK - Czech)

        //Throw exception: bad lat
        if ( fabs ( lat ) > MAX_LAT )
        {
                throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: ", " can not convert lon to lon_trans, lat > +- Pi/2", lat );
        }

        //Throw exception: bad lon
        if ( fabs ( lon ) > MAX_LON )
        {
                throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: ", " can not convert lon to lon_trans, lon > +- Pi", lat );
        }

        //Projection in normal position
        if ( fabs ( MAX_LAT - latp ) < ANGLE_ROUND_ERROR )
        {
                return lon;
        }

        //Singular point
        if ( ( lonp >= 0 ) && ( lon == lonp + MIN_LON ) )
        {
                throw ErrorBadData ( "ErrorBadData: ", " can not convert lon to lon_trans, singular point (lon = lonp - 180)." );
        }

        //Singular point
        if ( ( lonp < 0 ) && ( lon == lonp + MAX_LON ) )
        {
                throw ErrorBadData ( "ErrorBadData: ", " can not convert lon to lon_trans, singular point (lon = lonp + 180)." );
        }

        //Projection in oblique aspect
        else
        {
                //Same coordinates as the cartographic pole: singular point
                if ( ( fabs ( lon - lonp ) < ANGLE_ROUND_ERROR ) && ( fabs ( lat - latp ) < ANGLE_ROUND_ERROR ) )
                {
                        throw ErrorBadData ( "ErrorBadData: ", " can not convert lon to lon_trans, singular point (lat=latp, lon=lonp)." );
                }

                //First computation of the longitude
                T lon_trans_asin = sin ( ( lonp - lon ) * M_PI / 180.0 ) * cos ( lat * M_PI / 180.0 ) / cos ( lat_trans * M_PI / 180.0 );

                //Opposite pole, singular point: unable to compute lon_trans
                if ( lon_trans_asin > 1.0 + ARGUMENT_ROUND_ERROR  || lon_trans_asin < -1.0 - ARGUMENT_ROUND_ERROR )
                {
                        throw ErrorBadData ( "ErrorBadData: ", " can not convert lon to lon_trans, singular point (lat=-latp, lon=lonp)." );
                }

                //Correct longitude
                if ( lon_trans_asin > 1 )
                {
                        lon_trans_asin = 1;
                }

                //Correct longitude
                else if ( lon_trans_asin < -1 )
                {
                        lon_trans_asin = -1;
                }

                //Compute lon_temp
                const T lon_temp = fabs ( asin ( lon_trans_asin ) * 180 / M_PI );

                //Second calculation of the longitude
                const T lon_trans_acos = ( -cos ( latp * M_PI / 180 ) * sin ( lat * M_PI / 180 ) + sin ( latp * M_PI / 180 ) * cos ( lat * M_PI / 180 ) * cos ( ( lonp - lon ) * M_PI / 180 ) ) / cos ( lat_trans * M_PI / 180 );

                //(0, Pi/2)
                if ( ( lon_trans_asin >= 0.0 ) && ( lon_trans_acos >= 0.0 ) )
                {
                        if ( lon_direction == NormalDirection ) return -lon_temp;
                        else if ( lon_direction == ReversedDirection ) return lon_temp;
                        else return lon_temp - 90.0;

                        //return ( lon_direction == NormalDirection ? -lon_temp : lon_temp );
                }

                //(Pi/2, Pi)
                if ( ( lon_trans_asin >= 0.0 ) && ( lon_trans_acos <= 0.0 ) )
                {
                        if ( lon_direction == NormalDirection ) return lon_temp - 180.0;
                        else if ( lon_direction == ReversedDirection ) return 180.0 - lon_temp;
                        else return 90.0 - lon_temp;

                        //return ( lon_direction == NormalDirection ? lon_temp - 180.0 : 180.0 - lon_temp );
                }

                //(Pi, 3/2 * Pi)
                if ( ( lon_trans_asin <= 0.0 ) && ( lon_trans_acos <= 0.0 ) )
                {
                        if ( lon_direction == NormalDirection ) return 180.0 - lon_temp;
                        else if ( lon_direction == ReversedDirection ) return lon_temp - 180.0;
                        else return lon_temp + 90.0;

                        //return ( lon_direction == NormalDirection ? 180.0 - lon_temp : lon_temp - 180.0 );
                }

                //(3/2 * Pi, 2 * Pi)
                if ( lon_direction == NormalDirection ) return lon_temp;
                else if ( lon_direction == ReversedDirection ) return - lon_temp;
                else return -lon_temp - 90.0;

                //return ( lon_direction == NormalDirection ? lon_temp : - lon_temp );
        }
}

/*
template <typename T>
T CartTransformation::latTransToLat ( const T lat_trans, const T lon_trans, const T latp, const T lonp )
{
	//Transform latitude  ( lat_trans, lon_trans ) -> ( lat ) using a cartographic pole (lapt, lonp)
	//Throw exception: bad lat
        if ( fabs ( lat_trans ) > MAX_LAT )
        {
                throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: ", "can not convert lat to lat_trans, lat > +- Pi/2", lat );
        }

        //Throw exception: bad lon
        if ( fabs ( lon_trans ) > MAX_LON )
        {
                throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: ", "can not convert lon to lon_trans, lon > +- Pi", lon );
        }

        //Projection in normal position
        if ( fabs ( MAX_LAT - latp ) < ANGLE_ROUND_ERROR )
        {
                return lat_trans;
        }

        //Projection in oblique position
        else
        {

	}

	return 0.0;
}


template <typename T>
T CartTransformation::lonTransToLon ( const T lat_trans, const T lon_trans, const T lat, const T latp, const T lonp )
{
	        //Transform longitude  ( lat_trans, lon_trans, lat ) -> ( lon ) using a cartographic pole (latp, lonp)

	 //Throw exception: bad lat
        if ( fabs ( lat_trans ) > MAX_LAT )
        {
                throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: ", " can not convert lon to lon_trans, lat > +- Pi/2", lat );
        }

        //Throw exception: bad lon
        if ( fabs ( lon_trans ) > MAX_LON )
        {
                throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: ", " can not convert lon to lon_trans, lon > +- Pi", lat );
        }

        //Projection in normal position
        if ( fabs ( MAX_LAT - latp ) < ANGLE_ROUND_ERROR )
        {
                return lon_trans;
        }

	/*
        //Singular point
        if ( ( lonp >= 0 ) && ( lon_trans == lonp + MIN_LON ) )
        {
                throw ErrorBadData ( "ErrorBadData: ", " can not convert lon to lon_trans, singular point (lon = lonp - 180)." );
        }

        //Singular point
        if ( ( lonp < 0 ) && ( lon == lonp + MAX_LON ) )
        {
                throw ErrorBadData ( "ErrorBadData: ", " can not convert lon to lon_trans, singular point (lon = lonp + 180)." );
        }

        //Projection in oblique position
        else
	{
	}

	return 0.0;

}
*/


template <typename T>
T CartTransformation:: latLonToX ( const char * equation_x,  const T lat, const T lon, const T R, const T a, const T b, const T dx, const T c,  const T lat0, const T lat1, const T lat2, const bool print_exception )
{
        //Compute x coordinate of the point P[lat, lon] in specific projection
        const T x_total = ArithmeticParser::parseEq ( equation_x,  lat, lon, R, a, b, c, lat0, lat1, lat2, print_exception ) + dx;

        //Throw exception
        if ( ( x_total > MAX_POINT_COORDINATE ) || ( x_total < -MAX_POINT_COORDINATE ) )
        {
                throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not compute x coordinate, x coordinate > MAX_POINT_COORDINATE, ", " MAX_POINT_COORDINATE  = 1.0e +09", x_total );
        }

        //Result
        return  x_total;
}


template <typename T>
T CartTransformation::latLonToY ( const char * equation_y,  const T lat, const T lon, const T R, const T a, const T b, const T dy, const T c, const T lat0, const T lat1, const T lat2, const bool print_exception )
{
        //Compute y coordinate of the point P[lat, lon]  in specific projection
        const T y_total = ArithmeticParser::parseEq ( equation_y,  lat, lon, R, a, b, c, lat0, lat1, lat2, print_exception ) + dy;

        //Throw exception
        if ( ( y_total > MAX_POINT_COORDINATE ) || ( y_total < -MAX_POINT_COORDINATE ) )
        {
                throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not compute y coordinate, y coordinate > MAX_POINT_COORDINATE, ", " MAX_POINT_COORDINATE  = 1.0e +09", y_total );
        }

        //Result
        return  y_total;
}


template <typename T>
void CartTransformation::wgsToJTSK ( const Point3DGeographic <T> *p1, Point3DCartesian <T> * p2 )
{
        //Deg => rad
        const T PI =  4.0 * atan ( 1.0 ), Ro = 57.29577951;

        //WGS-84
        const T A_WGS = 6378137.0000, B_WGS = 6356752.3142;
        const T E2_WGS = ( A_WGS * A_WGS - B_WGS * B_WGS ) / ( A_WGS * A_WGS );

        //Bessel
        const T A_Bes = 6377397.1550, B_Bes = 6356078.9633;
        const T E2_Bes = ( A_Bes * A_Bes - B_Bes * B_Bes ) / ( A_Bes * A_Bes ), E_Bes = sqrt ( E2_Bes );

        //Scale, Translation, Rotation
        const T m =  -3.5623e-6;
        const T X0 = -570.8285, Y0 = -85.6769, Z0 = -462.8420;
        const T OMX = 4.9984 / 3600 / Ro, OMY = 1.5867 / 3600 / Ro, OMZ = 5.2611 / 3600 / Ro;

        //JTSK
        const T  FI0 = 49.5 / Ro;
        const T U0 = ( 49.0 + 27.0 / 60 + 35.84625 / 3600 ) / Ro;
        const T ALFA = 1.000597498372;
        const T LA_FERRO = ( 17.0 + 40.0 / 60 ) / Ro;
        const T K = 0.9965924869;
        const T R = 6380703.6105;
        const T UK = ( 59.0 + 42.0 / 60 + 42.6969 / 3600 ) / Ro;
        const T VK = ( 42.0 + 31.0 / 60 + 31.41725 / 3600 ) / Ro;
        const T Ro0 = 1298039.0046;
        const T S0 = 78.5 / Ro;

        //Point parameters
        T W = sqrt ( 1 - E2_WGS * sin ( p1->getLat() / Ro ) * sin ( p1->getLat() / Ro ) );
        T M = A_WGS * ( 1 - E2_WGS ) / ( W * W * W );
        T N = A_WGS / W;

        //Transformation (B,L,H)WGS => (X,Y,Z)WGS
        T X_WGS = ( N + p1->getH() ) * cos ( p1->getLat() / Ro ) * cos ( p1->getLon() / Ro );
        T Y_WGS = ( N + p1->getH() ) * cos ( p1->getLat() / Ro ) * sin ( p1->getLon() / Ro );
        T Z_WGS = ( N * ( 1 - E2_WGS ) + p1->getH() ) * sin ( p1->getLat() / Ro );

        //Transformation (X,Y,Z)WGS =>(X,Y,Z)Bes
        T X_Bes = X0 + ( m + 1 ) * ( X_WGS + Y_WGS * OMZ - Z_WGS * OMY );
        T Y_Bes = Y0 + ( m + 1 ) * ( -X_WGS * OMZ + Y_WGS + Z_WGS * OMX );
        T Z_Bes = Z0 + ( m + 1 ) * ( X_WGS * OMY - Y_WGS * OMX + Z_WGS );

        //Transformation (X,Y,Z)Bes => (BLH)Bess
        T rad = sqrt ( X_Bes * X_Bes + Y_Bes * Y_Bes );
        T la = 2 * atan ( Y_Bes / ( rad + X_Bes ) );
        T p = atan ( ( A_Bes * Z_Bes ) / ( B_Bes * rad ) );
        T t = ( Z_Bes + E2_Bes * A_Bes * A_Bes / B_Bes * pow ( sin ( p ), 3 ) ) / ( rad - E2_Bes * A_Bes * pow ( cos ( p ), 3 ) );
        T fi = atan ( t );
        T H = sqrt ( 1 + t * t ) * ( rad - A_Bes / sqrt ( 1 + ( 1 - E2_Bes ) * t * t ) );

        //Transformation (fi,la)Bes => (u,v)sphere  (Gauss conformal projection)
        la = la + LA_FERRO;
        T ro = 1 / K * pow ( tan ( fi / 2 + PI / 4 ) * pow ( ( 1 - E_Bes * sin ( fi ) ) / ( 1 + E_Bes * sin ( fi ) ), E_Bes / 2 ), ALFA );
        T u = 2 * atan ( ro ) - PI / 2;
        T v = ALFA * la;

        //Transformation (u,v)sphere => (s,d)sphere
        T s = asin ( sin ( UK ) * sin ( u ) + cos ( UK ) * cos ( u ) * cos ( VK - v ) );
        T d = asin ( sin ( VK - v ) * cos ( u ) / cos ( s ) );

        //Transformation (s,d)sphere => (Ro,Eps)plane (Lambert conformal projection)
        T n = sin ( S0 );
        T Ro_JTSK = Ro0 * pow ( ( tan ( ( S0 / 2 + PI / 4 ) ) / ( tan ( s / 2 + PI / 4 ) ) ), n );
        T eps_JTSK = n * d;

        //(Ro, eps) => (x,y)
        T X_JTSK = Ro_JTSK * cos ( eps_JTSK );
        T Y_JTSK = Ro_JTSK * sin ( eps_JTSK );

        //Set computed parameters to the projected point
        if ( p1->getPointLabel() != NULL )
        {
                p2->setPointLabel ( p1->getPointLabel() );
        }

        p2->setX ( X_JTSK );
        p2->setY ( Y_JTSK );
        p2->setZ ( H );
}


#endif

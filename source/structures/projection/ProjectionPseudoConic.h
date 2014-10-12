// Description: Pseudoconic projection, derived from Projection

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

#ifndef ProjectionPseudoConic_H
#define ProjectionPseudoConic_H


#include "Projection.h"


//Pseudoconic projection
template <typename T>
class ProjectionPseudoConic : virtual public Projection <T>
{
        private:
                T lat0;

        public:
                ProjectionPseudoConic() : Projection <T> (), lat0 ( 45 ) {}
                ProjectionPseudoConic ( const T R_, const T lat0_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_,  const char * projection_name_ ) :
                        Projection <T> ( R_, lon0_, dx_, dy_, c_, x_equat_, y_equat_, projection_name_ ), lat0 ( lat0_ ) {}
                virtual ~ProjectionPseudoConic() {};

        public:
                virtual Point3DGeographic <T> getCartPole() const {return Point3DGeographic <T> ( MAX_LAT, 0.0 );}
                virtual T getLat0() const {return lat0;}
                virtual T getLat1() const {return lat0;}
                virtual T getLat2() const {return lat0;}
                virtual T getA() const {return this->R;}
                virtual T getB() const {return this->R;}

                virtual TMinMax <T> getLatPInterval () const {return TMinMax <T> ( MAX_LAT, MAX_LAT );}
                virtual TMinMax <T> getLonPInterval () const {return TMinMax <T> ( 0.0, 0.0 );}
                virtual TMinMax <T> getLat0Interval () const {return TMinMax <T> ( MIN_LAT0, MAX_LAT0 );}
                virtual TMinMax <T> getLatPIntervalH ( const TMinMax <T> &lat ) const {return getLatPInterval();}
                virtual TMinMax <T> getLonPIntervalH ( const TMinMax <T> &lon ) const {return getLonPInterval();}
                virtual TTransformedLongtitudeDirection getLonDir () const { return ( TTransformedLongtitudeDirection ) 4;}

                virtual void setCartPole ( const Point3DGeographic <T> & pole )  {}
                virtual void setLat0 ( const T lat0_ ) {lat0 = lat0_;}
                virtual void setLat1 ( const T lat1 ) {}
                virtual void setLat2 ( const T lat2 ) {}
                virtual void setA ( const T a ) {}
                virtual void setB ( const T b ) {}
                virtual void setLonDir ( const TTransformedLongtitudeDirection lon_dir_ ) {}

                virtual void getShortCut ( char * shortcut ) const { strcpy ( shortcut, "PsConi" ); }
                virtual ProjectionPseudoConic <T> *clone() const {return new ProjectionPseudoConic <T> ( *this );}
                virtual void print ( std::ostream * file = &std::cout ) const {}
};

#endif

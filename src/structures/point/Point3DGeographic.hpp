// Description: 3D geographic point

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


#ifndef Point3DGeographic_HPP
#define Point3DGeographic_HPP

#include <string.h>

#include "libalgo/src/const/Const.h"


//Set point ID  to 0, initialize static variable
template <typename T>
unsigned int Point3DGeographic <T>::points_geo_id_counter = 0;


template <typename T>
Point3DGeographic <T> ::Point3DGeographic ( const char * point_label_, const T lat_, const T lon_, const T H_ ) : point_id ( points_geo_id_counter ++ ), lat ( lat_ ), lon ( lon_ ), H ( H_ )
{
        if ( point_label_ != NULL )
        {
                point_label = new char [ strlen ( point_label_ ) + 1 ];
                strcpy ( point_label, point_label_ );
        }

        else
        {
                point_label = NULL;
        }
}


template <typename T>
Point3DGeographic <T> ::Point3DGeographic ( const Point3DGeographic <T> *p ) : point_id ( p->point_id ), lat ( p->lat ), lon ( p->lon ), H ( p->H )
{
        if ( p->point_label != NULL )
        {
                point_label = new char [ strlen ( p->point_label ) + 1 ];
                strcpy ( point_label, p->point_label );
        }

        else
        {
                point_label = NULL;
        }
}


template <typename T>
Point3DGeographic <T> ::Point3DGeographic ( const Point3DGeographic <T> &p ) : point_id ( p.point_id ), lat ( p.lat ), lon ( p.lon ), H ( p.H )
{
        if ( p.point_label != NULL )
        {
                point_label = new char [ strlen ( p.point_label ) + 1 ];
                strcpy ( point_label, p.point_label );
        }

        else
        {
                point_label = NULL;
        }
}


template <typename T>
Point3DGeographic <T> & Point3DGeographic <T> ::operator = ( const Point3DGeographic <T> &p )
{
        if ( this != &p )
        {
                point_id = p.point_id;

                if ( point_label != NULL )
                {
                        delete [] point_label;
                        point_label = NULL;
                }

                if ( p.point_label != NULL )
                {
                        point_label = new char [ strlen ( p.point_label ) + 1 ];
                        strcpy ( point_label, p.point_label );
                }

                lat = p.lat;
                lon = p.lon;
                H = p.H;
        }

        return *this;
}


template <typename T>
bool Point3DGeographic <T> ::operator == ( const Point3DGeographic <T> &p ) const
{
        return ( lat - p.lat ) * ( lat - p.lat ) + ( lon - p.lon ) * ( lon - p.lon ) < MIN_POSITION_DIFF * MIN_POSITION_DIFF;
}


template <typename T>
Point3DGeographic<T>:: ~Point3DGeographic()
{
        if ( point_label != NULL )
        {
                delete point_label;
                point_label = NULL;
        }
}


template <typename T>
void Point3DGeographic <T>::print ( std::ostream * output ) const
{
        //Print point
        *output << std::fixed << std::setprecision ( 3 );
        *output << point_id << "   "  << lat << "   " << lon  << "   " << H << '\n';
}


template <typename T>
void Point3DGeographic <T>::setPointLabel ( const char * point_label_ )
{
        if ( point_label != NULL )
        {
                delete [] point_label;
        }

        if ( point_label_ != NULL )
        {
                point_label = new char [ strlen ( point_label_ ) + 1 ];
                strcpy ( point_label, point_label_ );
        }
}

#endif

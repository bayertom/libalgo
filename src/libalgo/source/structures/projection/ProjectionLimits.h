// Description: Limits in cartographic pole coordinates (lat_pole, lon_pole) and lat of the undistorted parallel

// Copyright (c) 2010 - 2011
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

#ifndef ProjectionLimits_H
#define ProjectionLimits_H


#include "Projection.h"

#include "Limits.h"


//Abstract class: limits of lat/lon
template <typename T>
class ProjectionLimits: virtual public Projection <T>, public Limits <T>
{

        public:
                ProjectionLimits() : Projection <T>(), Limits <T>() {}
                ProjectionLimits ( const T latp_min_, const T latp_max_, const T lonp_min_, const T lonp_max_, const T lat0_min_, const T lat0_max_ ) :
                        Projection <T> (), Limits <T> ( latp_min_, latp_max_, lonp_min_, lonp_max_, lat0_min_, lat0_max_ ) {}
                ProjectionLimits ( const T R_, const T lon0_, const T dx_, const T dy_, const char * x_equat_, const char * y_equat_, const char * projection_name_, const T latp_min_, const T latp_max_, const T lonp_min_, const T lonp_max_, const T lat0_south_min_, const T lat0_south_max_, const T lat0_north_min_, const T lat0_north_max_ ) :
                        Projection <T> ( R_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_ ), Limits <T> ( latp_min_, latp_max_, lonp_min_, lonp_max_, lat0_south_min_, lat0_south_max_, lat0_north_min_, lat0_north_max_ ) {}

        public:
                virtual ~ProjectionLimits() = 0;

};

#include "ProjectionLimits.hpp"

#endif


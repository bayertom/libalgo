// Description: Conic projection with constraints
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

#ifndef ProjectionConicLimits_H
#define ProjectionConicLimits_H


#include "ProjectionConic.h"
#include "ProjectionLimits.h"


//Conic projection: limits of lat/lon
template <typename T>
class ProjectionConicLimits : public ProjectionConic <T>, public ProjectionLimits <T>
{
        public:

                ProjectionConicLimits() :  Projection <T> (), ProjectionConic <T> (), ProjectionLimits <T> ( -90, 90, -180, 180, 10, 80 ) {};
                ProjectionConicLimits ( const T R_, const T lat0_, const T lat1_, const T lat2_, const T latp_, const T lonp_, const T lon0_, const T dx_, const T dy_, const char * x_equat_, const char * y_equat_, const char * projection_name_ ) :
                        Projection <T> ( R_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_ ), ProjectionConic <T> ( R_, lat0_, lat1_, lat2_, latp_, lonp_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_ ),
                        ProjectionLimits <T> ( R_, latp_, lonp_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_, -90, 90, -180, 180, 10, 80 ) {}
                virtual ~ProjectionConicLimits() {}

        public:
                virtual ProjectionConicLimits <T> *clone() const {return new ProjectionConicLimits <T> ( *this );}
};

#endif

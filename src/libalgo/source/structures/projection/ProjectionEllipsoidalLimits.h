// Description: Ellipsoidal projection with constraints

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

#ifndef ProjectionEllipsoidalLimits_H
#define ProjectionEllipsoidalLimits_H


#include "ProjectionEllipsoidal.h"
#include "ProjectionLimits.h"


//Ellipsoidal projection: limits of lat/lon
template <typename T>
class ProjectionEllipsoidalLimits : public ProjectionEllipsoidal <T>, public ProjectionLimits <T>
{
        public:

                ProjectionEllipsoidalLimits () : Projection <T> (), ProjectionEllipsoidal <T> (), ProjectionLimits <T> ( 90, 90, 0, 0, 0, 80 ) {}
                ProjectionEllipsoidalLimits ( const T R_, const T a_, const T b_, const T lat0_, const T lon0_, const T dx_, const T dy_, const char * x_equat_, const char * y_equat_,  const char * projection_name_ ) :
                        Projection <T> ( R_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_ ), ProjectionEllipsoidal <T> ( R_, a_, b_, lat0_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_ ),
                        ProjectionLimits <T> ( R_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_, 90, 90, 0, 0, 0, 80 ) {}
                virtual ~ProjectionEllipsoidalLimits() {}

        public:
                virtual ProjectionEllipsoidalLimits <T> *clone() const {return new ProjectionEllipsoidalLimits <T> ( *this );}

};


#endif


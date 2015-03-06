// Description: Miscellaneous projection with constaints

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

#ifndef ProjectionMiscellaneousLimits_H
#define ProjectionMiscellaneousLimits_H


#include "ProjectionMiscellaneous.h"
#include "ProjectionLimits.h"


//Miscellaneous projection: limits of lat/lon
template <typename T>
class ProjectionMiscellaneousLimits : public ProjectionMiscellaneous <T>, public ProjectionLimits <T>
{
        public:
                ProjectionMiscellaneousLimits() : Projection <T> (), ProjectionMiscellaneous <T> (), ProjectionLimits <T> ( 90, 90, 0, 0, 0, 0 ) {}
                ProjectionMiscellaneousLimits ( const T R_, const T latp_, const T lonp_, const T lon0_, const T dx_, const T dy_, const char * x_equat_, const char * y_equat_,  const char * projection_name_ ) :
                        Projection <T> ( R_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_ ), ProjectionMiscellaneous <T> ( R_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_ ),
                        ProjectionLimits <T> ( R_, lon0_, dx_, dy_, x_equat_, y_equat_, projection_name_, 90, 90, 0, 0 ) {}
                virtual ~ProjectionMiscellaneousLimits() {}

        public:
                virtual ProjectionMiscellaneousLimits <T> *clone() const {return new ProjectionMiscellaneousLimits <T> ( *this );}

};

#endif
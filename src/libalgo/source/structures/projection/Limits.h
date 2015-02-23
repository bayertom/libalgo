// Description: Abstract class  storing limits of the latitude / longitude for cartographic pole and undistorted parallel,
// set for each cartographic projection

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


#ifndef Limits_H
#define Limits_H

//Abstract class  storing limits of the latitude / longitude for cartographic pole and undistorted parallel
template <typename T>
class Limits
{
        protected:
                T latp_min;         //Min latitude of the cartographic pole
                T latp_max;         //Max latitude of the cartographic pole
                T lonp_min;         //Min longitude of the cartographic pole
                T lonp_max;         //Max longitude of the cartographic pole
                T lat0_min;	    //Min latitude of the undistorted parallel
                T lat0_max;         //Max latitude of the undistorted parallel (southern hemisphere)


        public:
                Limits() : latp_min ( -90 ), latp_max ( 90 ), lonp_min ( -180 ), lonp_max ( 180 ), lat0_min ( 0 ), lat0_max ( 90 ) {}
                Limits ( const T latp_min_, const T latp_max_, const T lonp_min_, const T lonp_max_, const T lat0_min_, const T lat0_max_ ) :
                        latp_min ( latp_min_ ), latp_max ( latp_max_ ), lonp_min ( lonp_min_ ), lonp_max ( lonp_max_ ), lat0_min ( lat0_min_ ), lat0_max ( lat0_max_ ) {}
                virtual ~Limits() = 0;

        public:
                T getLatPoleMin() const {return latp_min;}
                T getLatPoleMax() const {return latp_max;}
                T getLonPoleMin() const {return lonp_min;}
                T getLonPoleMax() const {return lonp_max;}
                T getLat0Min() const {return lat0_min;}
                T getLat0Max() const {return lat0_max;}

                void setPoleLimts ( const T latp_min_, const T latp_max_, const T lonp_min_, const T lonp_max_ )
                { latp_min = latp_min_; latp_max = latp_max_; lonp_min = lonp_min_; lonp_max = lonp_max_ ;}
                void setLat0Limits ( const T lat0_min_, const T lat0_max_ ) {lat0_min = lat0_min_; lat0_max = lat0_max_;}

        public:
                virtual Limits <T> *clone() const = 0;
};


template <typename T>
Limits<T>::~Limits () {}


#endif

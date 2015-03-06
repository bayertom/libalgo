// Description: Functor, compute matrix V of squares of residuals for cartometric analysis
// Method: Simplex Method, minimize CND

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


#ifndef FAnalyzeProjV2SCND_H
#define FAnalyzeProjV2SCND_H

template <typename T>
class Sample;

#include "libalgo/source/structures/projection/Sample.h"

//#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


//Forward declarations
template <typename T>
class Projection;


template <typename T>
struct TAnalysisParameters;


template <typename T>
class FAnalyzeProjV2SCND
{
        private:

                Container <Node3DCartesian <T> *> &nl_test;
                Container <Point3DGeographic <T> *> &pl_reference;
                typename TMeridiansList <T> ::Type &meridians;
                typename TParallelsList <T> ::Type &parallels;
                const Container <Face <T> *> &faces_test;
                Projection <T> *proj;
                const TAnalysisParameters <T> &analysis_parameters;
                const TProjectionAspect aspect;
                Sample <T> &sample_res;
                unsigned int & created_samples;
                std::ostream * output;

        public:

                FAnalyzeProjV2SCND ( Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
                                  const Container <Face <T> *> &faces_test_, Projection <T> *proj_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, std::ostream * output_ )
                        : nl_test ( nl_test_ ), pl_reference ( pl_reference_ ), meridians ( meridians_ ), parallels ( parallels_ ), faces_test ( faces_test_ ),  proj ( proj_ ),  analysis_parameters ( analysis_parameters_ ), aspect ( aspect_ ), sample_res ( sample_res_ ),
                          created_samples ( created_samples_ ), output ( output_ ) {}

                void operator () ( Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis )
                {

                        //Compute squares of residuals V = v' * W *v
                        const unsigned int m = nl_test.size();
                        const unsigned int m1 = X.rows();
                        const unsigned int m2 = V.rows();

                        //Get lat0 min and lat0 max
                        const T lat0_min = proj->getLat0Interval().min_val;
                        const T lat0_max = proj->getLat0Interval().max_val;

                        //Process all simplex vertices
                        for ( unsigned int i = 0; i < m1; i++ )
                        {

                                //Normal aspect: lat0, lon0
                                if ( aspect == NormalAspect )
                                {
                                        //Correct R, lat0, lon0
                                        if ( X ( i, 0 ) < 0.0 ) X ( i, 0 ) = fabs ( X ( i, 0 ) );

                                        //Subtract period
                                        if ( fabs ( X ( i, 3 ) ) > MAX_LAT ) X ( i, 3 ) = fmod ( X ( i, 3 ), 90 );

                                        if ( fabs ( X ( i, 4 ) ) > MAX_LON ) X ( i, 4 ) = fmod ( X ( i, 4 ), 180 );

                                        //Set to interval
                                        if ( X ( i, 3 ) < lat0_min ) X ( i, 3 ) = lat0_min;

                                        if ( X ( i, 3 ) > lat0_max ) X ( i, 3 ) = lat0_max;

                                        //Set c inside the interval
                                        if ( fabs ( X ( i, 5 ) ) < 0 ) X ( i, 5 ) = - X ( i, 5 );

                                        if ( fabs ( X ( i, 5 ) ) > MAX_C ) X ( i, 5 ) = 2 * MAX_C - X ( i, 5 );

                                }

                                //Transverse aspect: lonp, lat0
                                else  if ( aspect == TransverseAspect )

                                {
                                        //Correct R, lonp, lat0
                                        if ( X ( i, 0 ) < 0.0 ) X ( i, 0 ) = fabs ( X ( i, 0 ) );

                                        //Subtract period
                                        if ( fabs ( X ( i, 2 ) ) > MAX_LON ) X ( i, 2 ) = fmod ( X ( i, 2 ), 180 );

                                        if ( fabs ( X ( i, 3 ) ) > MAX_LAT ) X ( i, 3 ) = fmod ( X ( i, 3 ), 90 );

                                        //Set to interval
                                        if ( X ( i, 3 ) < lat0_min ) X ( i, 3 ) = lat0_min;

                                        if ( X ( i, 3 ) > lat0_max ) X ( i, 3 ) = lat0_max;

                                        //Set c inside the interval
                                        if ( fabs ( X ( i, 5 ) ) < 0 ) X ( i, 5 ) = - X ( i, 5 );

                                        if ( fabs ( X ( i, 5 ) ) > MAX_C ) X ( i, 5 ) = 2 * MAX_C - X ( i, 5 );

                                }

                                //Oblique aspect: latp, lonp, lat0
                                else if ( aspect == ObliqueAspect )
                                {
                                        //Correct R, latp, lonp, lat0
                                        if ( X ( i, 0 ) < 0.0 ) X ( i, 0 ) = fabs ( X ( i, 0 ) );

                                        //Subtract period
                                        if ( fabs ( X ( i, 1 ) ) > MAX_LAT )  X ( i, 1 ) = fmod ( X ( i, 1 ), 90 );

                                        if ( fabs ( X ( i, 2 ) ) > MAX_LON )  X ( i, 2 ) = fmod ( X ( i, 2 ), 180 );

                                        if ( fabs ( X ( i, 3 ) ) > MAX_LAT )  X ( i, 3 ) = fmod ( X ( i, 3 ), 90 );

                                        //Set lat0 inside the interval
                                        if ( X ( i, 3 ) < lat0_min || X ( i, 3 ) > lat0_max ) X ( i, 3 ) = 0.5 * ( lat0_min + lat0_max );

                                        //Set lonp to zero, if latp = 90
                                        if ( fabs ( X ( i, 1 ) - MAX_LAT ) < 3.0 )
                                        {
                                                X ( i, 2 ) = 0.0;
                                        }

                                        //Set c inside the interval
                                        if ( fabs ( X ( i, 5 ) ) < 0 ) X ( i, 5 ) = - X ( i, 5 );

                                        if ( fabs ( X ( i, 5 ) ) > MAX_C ) X ( i, 5 ) = 2 * MAX_C - X ( i, 5 );

                                }


                                //Set properties to the projection: ommit estimated radius, additional constants dx, dy
                                // They will be estimated again using the transformation
                                Point3DGeographic <T> cart_pole ( X ( i, 1 ), X ( i, 2 ) );
                                proj->setCartPole ( cart_pole );
                                proj->setLat0 ( X ( i, 3 ) );
                                proj->setLon0 ( X ( i, 4 ) );
                                proj->setDx ( 0.0 );
                                proj->setDy ( 0.0 );
                                proj->setC ( X ( i, 5 ) );
				proj->setR(X(i, 0));

                                //Compute analysis for one sample
                                if ( compute_analysis )
                                {
                                        try
                                        {
                                                //Compute analysis
                                                try
                                                {
                                                        CartAnalysis::computeAnalysisForOneSample ( nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, sample_res, false, created_samples, output );
                                                }

                                                //Throw exception
                                                catch ( Error & error )
                                                {
                                                        if ( analysis_parameters.print_exceptions )
                                                        {
                                                                //Print error and info about projection properties
                                                                error.printException ( output );
                                                                *output << "proj = " << proj->getProjectionName() << "  latp = " << proj->getCartPole().getLat() << "  lonp = " << proj->getCartPole().getLon() << "  lat0 = " << proj->getLat0() << '\n';
                                                        }
                                                }

                                                //Get index list of the sample
                                                TIndexList non_singular_points_indices = sample_res.getNonSingularPointsIndices();
                                                TIndexList k_best_points_indices = sample_res.getKBestPointsIndices();

                                                //Change weights in W matrix: weights of singular points or outliers are 0, otherwise they are 1
                                                unsigned int index_k_best_points = 0, n_k_best = k_best_points_indices.size(), n_points = pl_reference.size();
                                                int index_point = -1;

                                                //Set initial index of a point
                                                if ( ( m > 0 ) && ( n_k_best > 0 ) )
                                                        index_point = non_singular_points_indices [ k_best_points_indices [index_k_best_points++] ];

                                                //Process all points
                                                for ( int j = 0; ( j < n_points ) && ( n_k_best > 0 ); j++ )
                                                {
                                                        //Set weight of point to 1 (it is not an outlier nor singular)
                                                        if ( j == index_point )
                                                        {
                                                                W ( index_point, index_point ) = 1.0; W ( index_point + m, index_point + m ) = 1.0;

                                                                if ( index_k_best_points < n_k_best ) index_point = non_singular_points_indices [ k_best_points_indices [index_k_best_points++] ];
                                                        }

                                                        //Set weight of point to zero (it is an outlier or singular)
                                                        else
                                                        {
                                                                W ( j, j ) = 0.0; W ( j + m, j + m ) = 0.0;
                                                        }
                                                }
                                        }

                                        //Throw error
                                        catch ( Error & error )
                                        {
                                                if ( analysis_parameters.print_exceptions ) error.printException();
                                        }
                                }

                                //Compute coordinate differences (residuals): items of V matrix
                                Container <Node3DCartesianProjected <T> *> nl_projected_temp;

                                for ( unsigned int j = 0; j < m; j++ )
                                {
                                        //Get type of the direction
                                        TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

                                        //Reduce lon
                                        const T lon_red = CartTransformation::redLon0 ( pl_reference [j]->getLon(), X ( i, 4 ) );

                                        //Convert geographic point to oblique aspect
                                        T lat_trans = 0.0, lon_trans = 0.0, x = 0, y = 0;

                                        try
                                        {
                                                lat_trans = CartTransformation::latToLatTrans ( pl_reference [j]->getLat(), lon_red, X ( i, 1 ),  X ( i, 2 ) );
                                                lon_trans = CartTransformation::lonToLonTrans ( pl_reference [j]->getLat(), lon_red, lat_trans, X ( i, 1 ),  X ( i, 2 ), trans_lon_dir );

                                                //Compute x, y coordinates
                                                x =  ArithmeticParser::parseEq ( proj->getXEquat(), lat_trans, lon_trans, X ( i, 0 ), proj->getA(), proj->getB(), X ( i, 5 ), X ( i, 3 ), proj->getLat1(), proj->getLat2(), false );
                                                y =  ArithmeticParser::parseEq ( proj->getYEquat(), lat_trans, lon_trans, X ( i, 0 ), proj->getA(), proj->getB(), X ( i, 5 ), X ( i, 3 ), proj->getLat1(), proj->getLat2(), false );
                                        }

                                        //Bad conversion: singular point
                                        catch ( Error &error )
                                        {
                                                //Disable point from analysis: set weight to zero
                                                W ( j, j ) = 0; W ( j + m, j + m ) = 0;
                                        }

                                        //Create new cartographic point
                                        Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T> ( x, y );

                                        //Add point to the list
                                        nl_projected_temp.push_back ( n_projected );
                                }

                                //Computer centers of mass
                                unsigned int n_points = 0;
                                T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;

                                for ( unsigned int j = 0; j < m; j++ )
                                {
                                        //Use only non singular points
                                        if ( W ( j, j ) != 0.0 )
                                        {
                                                x_mass_test += nl_test[j]->getX();
                                                y_mass_test += nl_test[j]->getY();

                                                x_mass_reference += nl_projected_temp[j]->getX();
                                                y_mass_reference += nl_projected_temp[j]->getY();

                                                n_points++;
                                        }
                                }

                                x_mass_test = x_mass_test / n_points;
                                y_mass_test = y_mass_test / n_points;
                                x_mass_reference = x_mass_reference / n_points;
                                y_mass_reference = y_mass_reference / n_points;

                                //Compute coordinate differences (residuals): estimated - input
                                Matrix <T> v ( 2 * m, 1 );

                                for ( unsigned int j = 0; j < m; j++ )
                                {
                                        //Use only non singular points
                                        if ( W ( j, j ) != 0.0 )
                                        {
                                                const T vx = ( ( nl_projected_temp [j]->getX() - x_mass_reference ) - ( nl_test [j]->getX() - x_mass_test ) );
                                                const T vy = ( ( nl_projected_temp [j]->getY() - y_mass_reference ) - ( nl_test [j]->getY() - y_mass_test ) );

                                                //Input is the best vector
                                                if ( ( m1 == 1 ) && ( m2 > 1 ) )
                                                {
                                                        V ( j, 0 ) = vx;
                                                        V ( j + m, 0 ) = vy;
                                                }

                                                //Input is a simplex
                                                else
                                                {
                                                        v ( j, 0 ) = vx;
                                                        v ( j + m, 0 ) = vy;
                                                }
                                        }
                                }

                                //Compute squares of residuals: only if input is a simplex (nega)
                                if ( ( ( m1 > 1 ) || ( m2 == 1 ) ) )
                                {
                                        V ( i, 0 ) = MatrixOperations::norm ( MatrixOperations::trans ( v ) * W * v );
                                }

                                //Set DX, DY, compute it from the first vertex of the simplex
                                if ( i == 0 )
                                {
                                        sample_res.setDx ( x_mass_test - x_mass_reference );
                                        sample_res.setDy ( y_mass_test - y_mass_reference );
                                }
                        }
                }


};


#endif

// Description: Functor, compute matrix V of residuals for cartometric analysis, the cross nearest distance
// Method: NLSP, M6 (without rotation)

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


#ifndef FAnalyzeProjV2CND_H
#define FAnalyzeProjV2CND_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


template <typename T>
class Projection;


template <typename T>
class FAnalyzeProjV2CND
{
private:

		Container <Node3DCartesian <T> *> &nl_test;
		Container <Point3DGeographic <T> *> &pl_reference;
		typename TMeridiansList <T> ::Type &meridians;
		typename TParallelsList <T> ::Type &parallels;
		const Container <Face <T> *> &faces_test;
		Projection <T> *proj;
		T R_def;
		const TAnalysisParameters <T> &analysis_parameters;
		const TProjectionAspect aspect;
		Sample <T> &sample_res;
		unsigned int & created_samples;
		std::ostream * output;

public:

		FAnalyzeProjV2CND(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
				const Container <Face <T> *> &faces_test_, Projection <T> *proj_, const T R_def_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, std::ostream * output_)
				: nl_test(nl_test_), pl_reference(pl_reference_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), R_def(R_def_), analysis_parameters(analysis_parameters_), aspect(aspect_), sample_res(sample_res_),
				created_samples(created_samples_), output(output_) {}

		void operator () (Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis = true)
		{

				//Compute parameters of the V matrix: residuals
				const unsigned int m = nl_test.size();

				//Get lat0 min and lat0 max
				const T lat0_min = proj->getLat0Interval().min_val;
				const T lat0_max = proj->getLat0Interval().max_val;

				//Normal aspect: lat0, lon0
				if (aspect == NormalAspect)
				{
						//Correct R, lat0, lon0
						if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

						//Subtract period
						if (fabs(X(3, 0)) > MAX_LAT) X(3, 0) = fmod(X(3, 0), 90);

						if (fabs(X(4, 0)) > MAX_LON)
								X(4, 0) = fmod(X(4, 0), 180);

						//Set to interval
						if (X(3, 0) < lat0_min) X(3, 0) = lat0_min;

						if (X(3, 0) > lat0_max) X(3, 0) = lat0_max;
				}

				//Transverse aspect: lonp, lat0
				else  if (aspect == TransverseAspect)

				{
						//Correct R, lonp, lat0
						if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

						//if (X(1, 0) != 0.0) X(1, 0) = 0.0;

						//Subtract period
						if (fabs(X(2, 0)) > MAX_LON) X(2, 0) = fmod(X(2, 0), 180);

						if (fabs(X(3, 0)) > MAX_LAT) X(3, 0) = fmod(X(3, 0), 90);

						//Set to interval
						if (X(3, 0) < lat0_min) X(3, 0) = lat0_min;

						if (X(3, 0) > lat0_max) X(3, 0) = lat0_max;
				}

				//Oblique aspect: latp, lonp, lat0
				else if (aspect == ObliqueAspect)
				{
						//Correct R, latp, lonp, lat0
						if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

						//Subtract period
						if (fabs(X(1, 0)) > MAX_LAT)  X(1, 0) = fmod(X(1, 0), 90);

						if (fabs(X(2, 0)) > MAX_LON)  X(2, 0) = fmod(X(2, 0), 180);

						if (fabs(X(3, 0)) > MAX_LAT)  X(3, 0) = fmod(X(3, 0), 90);

						//Set lat0 inside the interval
						if (X(3, 0) < lat0_min || X(3, 0) > lat0_max) X(3, 0) = 0.5 * (lat0_min + lat0_max);

						//Set lonp to zero, if latp = 90
						if (fabs(X(1, 0) - MAX_LAT) < 1.0)  X(2, 0) = 0.0;
				}

				//Set properties to the projection: ommit estimated radius, additional constants dx, dy
				// They will be estimated again using the transformation
				Point3DGeographic <T> cart_pole(X(1, 0), X(2, 0));
				proj->setR(R_def);
				proj->setCartPole(cart_pole);
				proj->setLat0(X(3, 0));
				proj->setLon0(X(4, 0));
				proj->setDx(0.0);
				proj->setDy(0.0);
				proj->setC(0.0);
				proj->setR(X(0, 0));

				//Compute analysis for one sample
				if (compute_analysis)
				{
						try
						{
								//Compute analysis
								try
								{
										CartAnalysis::computeAnalysisForOneSample(nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, sample_res, false, created_samples, output);
								}

								//Throw exception
								catch (Error & error)
								{
										if (analysis_parameters.print_exceptions)
										{
												//Print error and info about projection properties
												error.printException(output);
												*output << "proj = " << proj->getProjectionName() << "  latp = " << proj->getCartPole().getLat() << "  lonp = " << proj->getCartPole().getLon() << "  lat0 = " << proj->getLat0() << "  c = " << proj->getC() << '\n';
										}
								}

								//Get index list of the sample
								TIndexList non_singular_points_indices = sample_res.getNonSingularPointsIndices();
								TIndexList k_best_points_indices = sample_res.getKBestPointsIndices();

								//Change weights in W matrix: weights of singular points or outliers are 0, otherwise they are 1
								unsigned int index_k_best_points = 0, n_k_best = k_best_points_indices.size(), n_points = pl_reference.size();
								int index_point = (n_k_best > 0 ? non_singular_points_indices[k_best_points_indices[index_k_best_points++]] : -1);

								for (int i = 0; (i < n_points) && (n_k_best > 0); i++)
								{
										//Set weight of point to 1 (it is not an outlier nor singular)
										if (i == index_point)
										{
												W(index_point, index_point) = 1.0; W(index_point + n_points, index_point + n_points) = 1.0;

												if (index_k_best_points < n_k_best) index_point = non_singular_points_indices[k_best_points_indices[index_k_best_points++]];
										}

										//Set weight of point to zero (it is an outlier or singular)
										else
										{
												W(i, i) = 0.0; W(i + n_points, i + n_points) = 0.0;
										}
								}
						}

						//Throw error
						catch (Error & error)
						{
								if (analysis_parameters.print_exceptions) error.printException();
						}
				}

				//Set Earth radius
				//proj->setR(X(0, 0));

				//Compute coordinate differences (residuals): items of V matrix
				Container <Node3DCartesianProjected <T> *> nl_projected_temp;

				for (unsigned int i = 0; i < m; i++)
				{
						//Get type of the direction
						TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

						//Reduce lon
						const T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), X(4, 0));

						T lat_trans = 0.0, lon_trans = 0.0, x = 0, y = 0;

						try
						{
								//Convert geographic point to oblique position: use a normal direction of converted longitude
								lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), lon_red, X(1, 0), X(2, 0));
								lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), lon_red, lat_trans, X(1, 0), X(2, 0), trans_lon_dir);

								//Compute x, y coordinates
								x = ArithmeticParser::parseEq(proj->getXEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(5, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);
								y = ArithmeticParser::parseEq(proj->getYEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(5, 0), X(3, 0), proj->getLat1(), proj->getLat2(), false);
						}

						catch (Error &error)
						{
								//Disable point from analysis: set weight to zero
								W(i, i) = 0; W(i + m, i + m) = 0;
						}

						//Create new cartographic point
						Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T>(x, y);

						//Add point to the list
						nl_projected_temp.push_back(n_projected);
				}

				//Computer centers of mass for both systems P, P'
				unsigned int n_points = 0;
				T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;

				for (unsigned int i = 0; i < m; i++)
				{
						//Use only non singular points
						if (W(i, i) != 0.0)
						{
								x_mass_test += nl_test[i]->getX();
								y_mass_test += nl_test[i]->getY();

								x_mass_reference += nl_projected_temp[i]->getX();
								y_mass_reference += nl_projected_temp[i]->getY();

								n_points++;
						}
				}

				//Centers of mass
				x_mass_test = x_mass_test / n_points;
				y_mass_test = y_mass_test / n_points;
				x_mass_reference = x_mass_reference / n_points;
				y_mass_reference = y_mass_reference / n_points;
;
				Container <Node3DCartesian <T> *> nl_test_red;
				Container <Node3DCartesianProjected <T> *> nl_projected_red;
				for (unsigned int i = 0; i < m; i++)
				{
						Node3DCartesian <T> *point_test = new Node3DCartesianProjected <T>(nl_test[i]->getX() - x_mass_test, nl_test[i]->getY() - y_mass_test);
						Node3DCartesianProjected <T> *point_proj = new Node3DCartesianProjected <T>(nl_projected_temp[i]->getX() - x_mass_reference, nl_projected_temp[i]->getY() - y_mass_reference);

						nl_test_red.push_back(point_test);
						nl_projected_red.push_back(point_proj);
				}

				/*
				typename Point1::Type dist = 0;
				const unsigned int n1 = list1.size();
				const unsigned int n2 = list2.size();

				//Throw exception
				if (n1 != n2 || n1 == 0 || n2 == 0)
				{
						throw ErrorBadData("ErrorBadData: can not compute cross nearest diastance,", " both lists have different size.");
				}
				*/
				//Create KD trees for neighbours searching
				KDTree2D <Node3DCartesian <T> > kd_tree1;
				KDTree2D <Node3DCartesianProjected <T> > kd_tree2;
				kd_tree1.createKDTree2D(const_cast < Container <Node3DCartesian <T> *> &> (nl_test_red));
				kd_tree2.createKDTree2D(const_cast < Container <Node3DCartesianProjected <T> * > &> (nl_projected_red));


				//Compute coordinate differences (residuals): estimated - input
				for (unsigned int i = 0; i < m; i++)
				{
						//Use only non singular points
						if (W(i, i) != 0.0)
						{
								Node3DCartesian <T> * point1 = nl_test_red[i];
								Node3DCartesianProjected <T>  * point2 = nl_projected_red[i];
								
								//Find nearest point to point from first dataset in KD tree of the second dataset
								Node3DCartesian <T>  * nearest1 = kd_tree1.findNN(point2);

								//Find nearest point to point from second dataset in KD tree of the first dataset
								Node3DCartesianProjected <T> * nearest2 = kd_tree2.findNN(point1);

								T dx = (point1->getX() - nearest2->getX() - nearest1->getX() + point2->getX() );
								T dy = (point1->getY() - nearest2->getY() - nearest1->getY() + point2->getY() );
								
								V(i, 0) = dx;
								V(i + m, 0) = dy;
						}
				}

				//Set DX, DY
				sample_res.setDx(x_mass_test - x_mass_reference);
				sample_res.setDy(y_mass_test - y_mass_reference);

				//V.print();
		}

};


#endif

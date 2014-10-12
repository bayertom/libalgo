// Description: Functor, compute matrix V of residuals for cartometric analysis
// Method: Simplex method, M7 (5 determined parameters, without rotation, without radius)

// Copyright (c) 2010 - 2014
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


#ifndef FAnalyzeProjV4S_H
#define FAnalyzeProjV4S_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


template <typename T>
class Projection;


template <typename T>
class FAnalyzeProjV4S
{
private:

	Container <Node3DCartesian <T> *> &nl_test;
	Container <Point3DGeographic <T> *> &pl_reference;
	typename TMeridiansList <T> ::Type &meridians;
	typename TParallelsList <T> ::Type &parallels;
	const Container <Face <T> *> &faces_test;
	Projection <T> *proj;
	T &R_def;
	T &q1;
	T &q2;
	const TAnalysisParameters <T> &analysis_parameters;
	const TProjectionAspect aspect;
	Sample <T> &sample_res;
	unsigned int & created_samples;
	std::ostream * output;

	unsigned int &iter;

public:

	FAnalyzeProjV4S(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
		const Container <Face <T> *> &faces_test_, Projection <T> *proj_, T &R_def_, T &q1_, T &q2_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, unsigned int &iter_, std::ostream * output_)
		: nl_test(nl_test_), pl_reference(pl_reference_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), R_def(R_def_), q1(q1_), q2(q2_), analysis_parameters(analysis_parameters_), aspect(aspect_), sample_res(sample_res_),
		created_samples(created_samples_), iter(iter_), output(output_) {}


	void operator () (Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, Matrix <T> &W, const bool compute_analysis = true)
	{

		//Compute parameters of the V matrix: residuals
		const unsigned int m = nl_test.size();
		const unsigned int m1 = X.rows();
		const unsigned int m2 = V.rows();

		//Get lat0 min and lat0 max
		const T lat0_min = proj->getLat0Interval().min_val;
		const T lat0_max = proj->getLat0Interval().max_val;

		//Process all simplex vertices
		for (unsigned int i = 0; i < m1; i++)
		{

			//Normal aspect: lat0, lon0
			if (aspect == NormalAspect)
			{
				//Subtract period
				if (fabs(X(i, 2)) > MAX_LAT) X(i, 2) = fmod(X(i, 2), 90);

				if (fabs(X(i, 3)) > MAX_LON)
					X(i, 3) = fmod(X(i, 3), 180);

				//Set to interval
				if (X(i, 2) < lat0_min) X(i, 2) = lat0_min;

				if (X(i, 2) > lat0_max) X(i, 2) = lat0_max;
			}

			//Transverse aspect: lonp, lat0
			else  if (aspect == TransverseAspect)
			{
				//Subtract period
				if (fabs(X(i, 1)) > MAX_LON) X(i, 1) = fmod(X(i, 1), 180);

				if (fabs(X(i, 2)) > MAX_LAT) X(i, 2) = fmod(X(i, 2), 90);

				//Set to interval
				if (X(i, 2) < lat0_min) X(i, 2) = lat0_min;

				if (X(i, 2) > lat0_max) X(i, 2) = lat0_max;


				//Set lon0
				if (fabs(X(i, 3)) > MAX_LON)  X(i, 3) = fmod(X(i, 1), 180);
			}

			//Oblique aspect: latp, lonp, lat0
			else if (aspect == ObliqueAspect)
			{
				//Subtract period
				if (fabs(X(i, 0)) > MAX_LAT)  X(i, 0) = fmod(X(i, 0), 90);

				if (fabs(X(i, 1)) > MAX_LON)  X(i, 1) = fmod(X(i, 1), 180);

				if (fabs(X(i, 2)) > MAX_LAT)  X(i, 2) = fmod(X(i, 2), 90);

				//Set lat0 inside the interval
				if (X(i, 2) < lat0_min || X(i, 2) > lat0_max) X(i, 2) = 0.5 * (lat0_min + lat0_max);

				//Set lonp to zero, if latp = 90
				if (fabs(X(i, 0) - MAX_LAT) < 1.0)  X(i, 1) = 0.0;

				//Set lonp to zero, if latp = 90

				if (fabs(X(i, 0) - MAX_LAT) < 1)
				{
					X(i, 0) = 90.0;
					X(i, 1) = 0.0;
				}

				else if (fabs(X(i, 0)) < 1)
				{
					//X(0, 0) = 0;
					//X(1, 0) = 90;
				}

				//Set lon0
				if (fabs(X(i, 3)) > MAX_LON)  X(i, 3) = fmod(X(i, 1), 180);
				X(i, 3) = 0;
			}

			//Set properties to the projection: ommit estimated radius, additional constants dx, dy
			// They will be estimated again using the transformation
			Point3DGeographic <T> cart_pole(X(i, 0), X(i, 1));
			proj->setR(R_def);
			proj->setCartPole(cart_pole);
			proj->setLat0(X(i, 2));
			proj->setLon0(X(i, 3));
			proj->setDx(0.0);
			proj->setDy(0.0);
			proj->setC(X(i, 4));


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

					for (int j = 0; (j < n_points) && (n_k_best > 0); j++)
					{
						//Set weight of point to 1 (it is not an outlier nor singular)
						if (j == index_point)
						{
							W(index_point, index_point) = 1.0; W(index_point + n_points, index_point + n_points) = 1.0;

							if (index_k_best_points < n_k_best) index_point = non_singular_points_indices[k_best_points_indices[index_k_best_points++]];
						}

						//Set weight of point to zero (it is an outlier or singular)
						else
						{
							W(j, j) = 0.0; W(j + n_points, j + n_points) = 0.0;
						}
					}
				}

				//Throw error
				catch (Error & error)
				{
					if (analysis_parameters.print_exceptions) error.printException();
				}
			}

			//Compute coordinate differences (residuals): items of V matrix
			Container <Node3DCartesianProjected <T> *> nl_projected_temp;

			for (unsigned int j = 0; j < m; j++)
			{
				//Get type of the direction
				TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

				//Reduce lon
				const T lon_red = CartTransformation::redLon0(pl_reference[j]->getLon(), X(i, 3));

				T lat_trans = 0.0, lon_trans = 0.0, x = 0.0, y = 0.0;

				try
				{
					//Convert geographic point to oblique position: use a normal direction of converted longitude
					lat_trans = CartTransformation::latToLatTrans(pl_reference[j]->getLat(), lon_red, X(i, 0), X(i, 1));
					lon_trans = CartTransformation::lonToLonTrans(pl_reference[j]->getLat(), lon_red, lat_trans, X(i, 0), X(i, 1), trans_lon_dir);

					//Compute x, y coordinates
					x = ArithmeticParser::parseEq(proj->getXEquat(), lat_trans, lon_trans, R_def, proj->getA(), proj->getB(), X(i, 4), X(i, 2), proj->getLat1(), proj->getLat2(), false);
					y = ArithmeticParser::parseEq(proj->getYEquat(), lat_trans, lon_trans, R_def, proj->getA(), proj->getB(), X(i, 4), X(i, 2), proj->getLat1(), proj->getLat2(), false);
				}

				catch (Error &error)
				{
					//Disable point from analysis: set weight to zero
					W(j, j) = 0; W(j + m, j + m) = 0;
				}

				//Create new cartographic point
				Node3DCartesianProjected <T> *n_projected = new Node3DCartesianProjected <T>(x, y);

				//Add point to the list
				nl_projected_temp.push_back(n_projected);
			}

			//Apply transformation
			TTransformationKeyHelmert2D <T> key_helmert;
			HelmertTransformation2D::getTransformKey(nl_test, nl_projected_temp, key_helmert);

			//std::cout << key_helmert.c1 << "   " << key_helmert.c2 << "   " << atan2(key_helmert.c2, key_helmert.c1) * 57.3;


			//Computer centers of mass for both systems P, P'
			unsigned int n_points = 0;
			T x_mass_test = 0.0, y_mass_test = 0.0, x_mass_reference = 0.0, y_mass_reference = 0.0;

			for (unsigned int j = 0; j < m; j++)
			{
				//Use only non singular points
				if (W(j, j) != 0.0)
				{
					x_mass_test += nl_test[j]->getX();
					y_mass_test += nl_test[j]->getY();

					x_mass_reference += nl_projected_temp[j]->getX();
					y_mass_reference += nl_projected_temp[j]->getY();

					n_points++;
				}
			}

			//Centers of mass
			x_mass_test = x_mass_test / n_points;
			y_mass_test = y_mass_test / n_points;
			x_mass_reference = x_mass_reference / n_points;
			y_mass_reference = y_mass_reference / n_points;

			//Compute scale using the least squares adjustment: h = inv (A'WA)A'WL, weighted Helmert transformation
			T sum_xy_1 = 0, sum_xy_2 = 0, sum_xx_yy = 0;
			for (unsigned int j = 0; j < n_points; j++)
			{
				sum_xy_1 = sum_xy_1 + (nl_test[j]->getX() - x_mass_test) * W(j, j) * (nl_projected_temp[j]->getX() - x_mass_reference) +
					(nl_test[j]->getY() - y_mass_test) * W(j, j) * (nl_projected_temp[j]->getY() - y_mass_reference);
				sum_xy_2 = sum_xy_2 + (nl_test[j]->getY() - y_mass_test) * W(j, j) * (nl_projected_temp[j]->getX() - x_mass_reference) -
					(nl_test[j]->getX() - x_mass_test) * W(j, j) * (nl_projected_temp[j]->getY() - y_mass_reference);
				sum_xx_yy = sum_xx_yy + (nl_projected_temp[j]->getX() - x_mass_reference) * W(j, j) * (nl_projected_temp[j]->getX() - x_mass_reference) +
					(nl_projected_temp[j]->getY() - y_mass_reference) * W(j, j) * (nl_projected_temp[j]->getY() - y_mass_reference);
			}

			//Transformation ratios
			q1 = sum_xy_1 / sum_xx_yy;
			q2 = sum_xy_2 / sum_xx_yy;

			//Rotation
			const T alpha = atan2(q2, q1) * 180.0 / M_PI;
			//std::cout << "rot = " << alpha;
			//std::cout << "   scale = " << sqrt(q1*q1 + q2*q2);

			//Compute coordinate differences (residuals): estimated - input
			Matrix <T> v(2 * m, 1);
			for (unsigned int j = 0; j < m; j++)
			{
				//Use only non singular points
				if (W(j, j) != 0.0)
				{
					const T vx = (q1 * (nl_projected_temp[j]->getX() - x_mass_reference) - q2 * (nl_projected_temp[j]->getY() - y_mass_reference) - (nl_test[j]->getX() - x_mass_test));
					const T vy = (q2 * (nl_projected_temp[j]->getX() - x_mass_reference) + q1 * (nl_projected_temp[j]->getY() - y_mass_reference) - (nl_test[j]->getY() - y_mass_test));

					///Input is the best vector
					if ((m1 == 1) && (m2 > 1))
					{
						V(j, 0) = vx;
						V(j + m, 0) = vy;
					}

					//Input is a simplex
					else
					{
						v(j, 0) = vx;
						v(j + m, 0) = vy;
					}
				}
			}

			//Compute squares of residuals: only if input is a simplex (nega)
			if (((m1 > 1) || (m2 == 1)))
			{
				V(i, 0) = MatrixOperations::norm(MatrixOperations::trans(v) * W * v);
			}

			//Set DX, DY, compute it from the first vertex of the simplex
			if (i == 0)
			{
				//Compute DX, DY
				const T dx = x_mass_test - x_mass_reference * q1 + y_mass_reference * q2;
				const T dy = y_mass_test - x_mass_reference * q2 - y_mass_reference * q1;

				sample_res.setDx(dx);
				sample_res.setDy(dy);

				//Set rotation
				sample_res.setRotation(alpha);

				//Perform scaling
				R_def *= sqrt(q1 * q1 + q2 * q2);
				q1 = q1 / R_def;
				q2 = q2 / R_def;
			}

			iter++;
		}

		
	}


};


#endif

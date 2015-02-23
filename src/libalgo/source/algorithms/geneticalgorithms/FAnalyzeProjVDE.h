// Description: Functor, compute matrix V of squares of residuals for cartometric analysis
// Method: Differential evolution

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


#ifndef FAnalyzeProjVDE_H
#define FAnalyzeProjVDE_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"


template <typename T>
class Projection;


template <typename T>
class FAnalyzeProjVDE
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

	unsigned int &iter;

public:

	FAnalyzeProjVDE(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
		const Container <Face <T> *> &faces_test_, Projection <T> *proj_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, unsigned int &iter_, std::ostream * output_)
		: nl_test(nl_test_), pl_reference(pl_reference_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), analysis_parameters(analysis_parameters_), aspect(aspect_), sample_res(sample_res_),
		created_samples(created_samples_), iter(iter_), output(output_) {}

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
			if (fabs(X(0, 3)) > MAX_LAT) X(0, 3) = fmod(X(0, 3), 90);

			if (fabs(X(0, 4)) > MAX_LON) X(0, 4) = fmod(X(0, 4), 180);

			//Set to interval
			if (X(0, 3) < lat0_min) X(0, 3) = lat0_min;

			if (X(0, 3) > lat0_max) X(0, 3) = lat0_max;
		}

		//Transverse aspect: lonp, lat0
		else  if (aspect == TransverseAspect)

		{
			//Correct R, lonp, lat0
			if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

			//Subtract period
			if (fabs(X(0, 2)) > MAX_LON) X(0, 2) = fmod(X(0, 2), 180);

			if (fabs(X(0, 3)) > MAX_LAT) X(0, 3) = fmod(X(0, 3), 90);

			//Set to interval
			if (X(0, 3) < lat0_min || X(0, 3) > lat0_max) X(0, 3) = 0.5 * (lat0_min + lat0_max);
		}

		//Oblique aspect: latp, lonp, lat0
		else if (aspect == ObliqueAspect)
		{
			//Correct R, latp, lonp, lat0
			if (X(0, 0) < 0.0) X(0, 0) = fabs(X(0, 0));

			//Subtract period
			if (fabs(X(0, 1)) > MAX_LAT)  X(0, 1) = fmod(X(0, 1), 90);

			if (fabs(X(0, 2)) > MAX_LON)  X(0, 2) = fmod(X(0, 2), 180);

			if (fabs(X(0, 3)) > MAX_LAT)  X(0, 3) = fmod(X(0, 3), 90);

			//Set lat0 inside the interval
			if (X(0, 3) < lat0_min || X(0, 3) > lat0_max) X(0, 3) = 0.5 * (lat0_min + lat0_max);

			//Set lonp to zero, if latp = 90
			if (fabs(X(0, 1) - MAX_LAT) < 1.0)
			{
				X(0, 1) = 90.0;
				X(0, 2) = 0.0;
			}

			//Set lon0
			X(0, 4) = 0;
		}

		//Set properties to the projection: ommit estimated radius, additional constants dx, dy
		// They will be estimated again using the transformation
		Point3DGeographic <T> cart_pole(X(0, 1), X(0, 2));
		proj->setR(X(0, 0));
		proj->setCartPole(cart_pole);
		proj->setLat0(X(0, 3));
		proj->setLon0(X(0, 4));
		proj->setDx(X(0, 5));
		proj->setDy(X(0, 6));
		proj->setC(X(0, 7));

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
						*output << "proj = " << proj->getProjectionName() << "  latp = " << proj->getCartPole().getLat() << "  lonp = " << proj->getCartPole().getLon() << "  lat0 = " << proj->getLat0() << '\n';
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

		//Compute new coordinates
		for (unsigned int i = 0; i < m; i++)
		{
			T x_new, y_new;

			try
			{
				//Get type of the direction
				TTransformedLongtitudeDirection trans_lon_dir = proj->getLonDir();

				//Reduce lon
				const T lon_red = CartTransformation::redLon0(pl_reference[i]->getLon(), X(0, 4));

				//Convert geographic point to oblique position: use a normal direction of converted longitude
				const T lat_trans = CartTransformation::latToLatTrans(pl_reference[i]->getLat(), lon_red, X(0, 1), X(0, 2));
				const T lon_trans = CartTransformation::lonToLonTrans(pl_reference[i]->getLat(), lon_red, lat_trans, X(0, 1), X(0, 2), trans_lon_dir);

				//Compute new coordinates: add shifts
				Y(i, 0) = ArithmeticParser::parseEq(proj->getXEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(0, 7), X(0, 3), proj->getLat1(), proj->getLat2(), false) + X(0, 5);
				Y(i + m, 0) = ArithmeticParser::parseEq(proj->getYEquat(), lat_trans, lon_trans, X(0, 0), proj->getA(), proj->getB(), X(0, 7), X(0, 3), proj->getLat1(), proj->getLat2(), false) + X(0, 6);
			}

			//Throw exception: bad conversion, a singular point
			catch (Error & error)
			{
				x_new = 1.0;
				y_new = 0.0;
			}

			//Compute coordinate differences (residuals): estimated - input
			V(i, 0) = (Y(i, 0) - nl_test[i]->getX());
			V(i + m, 0) = (Y(i + m, 0) - nl_test[i]->getY());
		}

		iter++;
	}

};


#endif
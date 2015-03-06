// Description: Functor, compute matrix V of residuals for cartometric analysis, lon solved using the bisection
// Method: NLSP, M7 (5 determined parameters, without rotation, without radius)

// Copyright (c) 2010 - 2015
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


#ifndef FAnalyzeProjV4B_H
#define FAnalyzeProjV4B_H


#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/cartanalysis/CartAnalysis.h"
#include "libalgo/source/algorithms/outliers/Outliers.h"
#include "libalgo/source/algorithms/bisection/Bisection.h"
#include "libalgo/source/algorithms/bisection/FV2Bisection.h"

template <typename T>
class Projection;


template <typename T>
class FAnalyzeProjV4B
{
	private:

		Container <Node3DCartesian <T> *> &nl_test;
		Container <Point3DGeographic <T> *> &pl_reference;
		typename TMeridiansList <T> ::Type &meridians;
		typename TParallelsList <T> ::Type &parallels;
		const Container <Face <T> *> &faces_test;
		Projection <T> *proj;
		T &R;
		T &q1;
		T &q2;
		const TAnalysisParameters <T> &analysis_parameters;
		const TProjectionAspect aspect;
		Sample <T> &sample_res;
		unsigned int & created_samples;
		unsigned int &res_evaluation;
		TMEstimatorsWeightFunction me_function;
		T k;
		Matrix <unsigned int> &I;
		std::ostream * output;

	public:

		FAnalyzeProjV4B(Container <Node3DCartesian <T> *> &nl_test_, Container <Point3DGeographic <T> *> &pl_reference_, typename TMeridiansList <T> ::Type &meridians_, typename TParallelsList <T> ::Type &parallels_,
		const Container <Face <T> *> &faces_test_, Projection <T> *proj_, T &R_est_, T &q1_, T &q2_, const TAnalysisParameters <T> & analysis_parameters_, const TProjectionAspect aspect_, Sample <T> &sample_res_, unsigned int & created_samples_, unsigned int &res_evaluation_, const TMEstimatorsWeightFunction &me_function_, const T k_, Matrix <unsigned int> &I_, std::ostream * output_)
		: nl_test(nl_test_), pl_reference(pl_reference_), meridians(meridians_), parallels(parallels_), faces_test(faces_test_), proj(proj_), R(R_est_), q1(q1_), q2(q2_), analysis_parameters(analysis_parameters_), aspect(aspect_), sample_res(sample_res_), 
		created_samples(created_samples_), res_evaluation(res_evaluation_), me_function(me_function_), k(k_), I(I_), output(output_) {}

		T operator () (Matrix <T> &X)
		{
			//Compute residuals, call simple wrapper
			return evaluateResiduals(X, nl_test, pl_reference, meridians, parallels, faces_test, proj, R, q1, q2, analysis_parameters, aspect, sample_res, created_samples, res_evaluation, me_function, k, I, output);
		}
			
};

/*
template <typename T>
T evaluateResiduals(const Matrix <T> &X, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels,
	const Container <Face <T> *> &faces_test, Projection <T> *proj, T &R, T &q1, T &q2, const TAnalysisParameters <T> & analysis_parameters, const TProjectionAspect aspect, Sample <T> &sample_res, unsigned int & created_samples, 
	unsigned int &res_evaluation, const TMEstimatorsWeightFunction &me_function, const T k, Matrix <unsigned int> &I, std::ostream * output)
{
	//Simple wrapper calling method from FAnalyzeProj4
	const unsigned int m = nl_test.size(), n = X.rows();
	Matrix <T> W(2 * m, 2 * m, 0.0, 1), V(2 * m, 1), Y(2 * m, 1);

	//Call function from FAnalyzeProjV4
	evaluateResiduals(X, Y, V, W, nl_test, pl_reference, meridians, parallels, faces_test, proj, R, q1, q2, analysis_parameters, aspect, sample_res, created_samples, res_evaluation, me_function, k, I, output);

	//Compute objective function value
	const T fval = norm(trans(V) * V);

	return fval;
}
*/

#endif

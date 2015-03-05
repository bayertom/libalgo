// Description: Functor, compute matrix V of residuals for cartometric analysis
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


#ifndef FAnalyzeProjV2B_HPP
#define FAnalyzeProjV2B_HPP
/*

template <typename T>
T FAnalyzeProjV2B<T>:: operator () (Matrix <T> &X)
{
	//Compute residuals, call simple wrapper
	return evaluateResiduals(X, nl_test, pl_reference, meridians, parallels, faces_test,proj, analysis_parameters, aspect, sample_res, created_samples, res_evaluation, me_function, k, I, output);
}


template <typename T>
inline T evaluateResiduals(const Matrix <T> &X, Container <Node3DCartesian <T> *> &nl_test, Container <Point3DGeographic <T> *> &pl_reference, typename TMeridiansList <T> ::Type &meridians, typename TParallelsList <T> ::Type &parallels,
	const Container <Face <T> *> &faces_test, Projection <T> *proj, const TAnalysisParameters <T> & analysis_parameters, const TProjectionAspect aspect, Sample <T> &sample_res, unsigned int & created_samples, unsigned short &res_evaluation, const TMEstimatorsWeightFunction &me_function, 
	const T k, Matrix <unsigned int> &I, std::ostream * output)
{
	//Simple wrapper calling method from FAnalyzeProjV2
	const unsigned int m = nl_test.size(), n = X.rows();
	Matrix <T> W(2 * m, 2 * m, 0.0, 1), V(2 * m, 1), Y(2 * m, 1);

	//Call function from FAnalyzeProjV2
	evaluateResiduals(X, Y, V, W, nl_test, pl_reference, meridians, parallels, faces_test, proj, analysis_parameters, aspect, sample_res, created_samples, res_evaluation, me_function, k, I, output);

	//Compute objective function value
	const T fval = norm(trans(V) * V);

	return fval;
}

*/
#endif

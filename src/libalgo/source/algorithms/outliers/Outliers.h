// Description: Outliers detection using the least squares and M-estimators

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


#ifndef Outliers_H
#define Outliers_H

#include "libalgo/source/structures/list/Container.h"

#include "libalgo/source/algorithms/transformation/Transformation2D.h"

//Forward declarations

//Several equations related to 2D transformation
class Outliers
{
	public:
		template <typename Point1, typename Point2, typename TKey>
                static void findOutliersLTS ( const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source, Container <Point1 *> &global_points_dest,
                                Container <Point2 *> &local_points_dest, TKey & min_key,  typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs, const typename Point1::Type perc_ratio = 0.8 );

                template <typename Point1, typename Point2, typename TKey>
                static void findOutliersIRLS ( const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source,
                                Container <Point1 *> &global_points_dest, Container <Point2 *> &local_points_dest,  TKey & min_key, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs );

		template <typename Point1, typename Point2, typename FunctionJ>
                static void findOutliersME ( const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source,
                                Container <Point1 *> &global_points_dest, Container <Point2 *> &local_points_dest, Matrix <typename Point1::Type> &X, FunctionJ function_j, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs, const unsigned int MAX_ITER );


	private:
		template <typename T>
                static void createKBestPairsOfPoints ( const TAccuracyCharacteristics <T> &deviations, typename TDevIndexPairs<T>::Type & point_pairs, const float perc_ratio );



};

#include "Outliers.hpp"

#endif
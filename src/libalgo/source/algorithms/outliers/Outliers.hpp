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


#ifndef Outliers_HPP
#define Outliers_HPP

template <typename Point1, typename Point2, typename TKey>
void Outliers::findOutliersLTS ( const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source, Container <Point1 *> &global_points_dest,
                Container <Point2 *> &local_points_dest, TKey & min_key, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs, const typename Point1::Type perc_ratio )
{
        //Removing outliers using the iteratively LTS
        //Outliers are dectected using weighted transformation and Danish method
        
        const unsigned int n_global_points_source = global_points_source.size();
        unsigned int n_iterations = 10 * n_global_points_source;
        bool best_key_longer = false;

        //Throw exception: not enough points
        if ( n_global_points_source < 3 )
        {
                throw ErrorBadData ( "ErrorBadData: not enough global points to be optimized. ", "Can not find optimal transformation key." );
        }

        //Initialize min error
        typename Point1::Type min_error = MAX_FLOAT;

        //Initialize random number generator
        srand ( ( unsigned ) time ( 0 ) );

        //Perform n iterations to find the best key
        unsigned int i = 0, best_key_size = 2, pairs_actual_size = 0;
        TAccuracyCharacteristics <typename Point1::Type> accuracy;

        do
        {
                //Get 2 randomly generated different indices
                int index2 = 0;
                int index1 = rand() % n_global_points_source ;

                do
                {
                        index2 = rand() % n_global_points_source;
                }
                while ( index1 == index2 );

                //Create point pairs
                typename TDevIndexPairs <typename Point1::Type>::Type pairs, pairs_actual;
                pairs_actual.push_back ( std::make_pair ( 0.0, index1 ) );
                pairs_actual.push_back ( std::make_pair ( 1.0, index2 ) );

                //Rearrange points: use only fist k-best points
                rearrangePoints ( global_points_source, local_points_source, global_points_dest, local_points_dest, pairs_actual );

                //Create keys and errors
                TKey key, key_actual;
                typename Point1::Type error = MAX_FLOAT, error_actual = 0.1 * error;

                //Performs iterations for this key until it convergates
                for ( unsigned int j = 0 ; j < 2, error_actual < error; j++ )
                {
                        //Assign old values
                        error = error_actual;
                        key = key_actual;
                        pairs = pairs_actual;

                        //Compute transformation for the key of 2 points
                        Container <Point1 *> transformed_points;
                        getTransformKey ( global_points_dest, local_points_dest, key_actual );
                        transform ( global_points_source, local_points_source, &transformed_points, key_actual );

                        //Compute deviation
                        TAccuracyCharacteristics <typename Point1::Type> accuracy;
                        accuracy = getAccuracyCharacteristics ( global_points_source, local_points_source, &transformed_points, key_actual );

                        //Find new k-best pairs and sort
                        createKBestPairsOfPoints ( accuracy, pairs_actual, perc_ratio );

                        //Rearrange points: use only fist k-best points
                        rearrangePoints ( global_points_source, local_points_source, global_points_dest, local_points_dest, pairs_actual );

                        //Compute transformation only k_best local points using k-best key
                        Container <Point1 *> transformed_points_dest;
                        transform ( global_points_dest, local_points_dest, transformed_points_dest, key_actual );

                        //Get accuracy characteristic for k-best global points and k-best transformed points using k-best key
                        TAccuracyCharacteristics <typename Point1::Type> accuracy_actual2 = getAccuracyCharacteristics ( global_points_dest, local_points_dest, &transformed_points_dest, key_actual );
                        error_actual = accuracy_actual2.std_dev;
                }

                //We found a better key: remember deviation (error), key and pairs
                if ( error_actual < min_error )
                {
                        min_error = error_actual;
                        min_key = key;
                        min_pairs = pairs;
                }

        }
        while ( ++i < n_iterations && perc_ratio != 1.0 );

        //Sort pairs points indices
        std::sort ( min_pairs.begin(), min_pairs.end(), sortPointPairsByIndices <typename Point1::Type> () );

        //Rearrange points: use only first k-best points
        rearrangePoints ( global_points_source, local_points_source, global_points_dest, local_points_dest, min_pairs );
}


template <typename Point1, typename Point2, typename TKey>
void Outliers::findOutliersIRLS ( const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source,
                Container <Point1 *> &global_points_dest, Container <Point2 *> &local_points_dest, TKey & min_key, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs )
{
        //Removing outliers using the iteratively IRLS
        //Outliers are dectected using weighted transformation and Danish method
        //Method published by Wieser and Brunner, 2000
        const unsigned int n_global_points_source = global_points_source.size();
        typename Point1::Type MIN_DEV_DIFFERENCE = 0.0001, std_dev_old = 0;

        //Accuracy parameters
        TAccuracyCharacteristics <typename Point1::Type> accuracy;

        //Throw exception: not enough points
        if ( n_global_points_source < 3 )
        {
                throw ErrorBadData ( "ErrorBadData: not enough global points to be optimized. ", "Can not find optimal transformation key." );
        }

        //Initialize weights
        typename TWeights <typename Point1::Type> ::Type weights ( n_global_points_source, 1 );

        //First iteration: compute as unweighted triangulation
        Container <Point1 *> transformed_points_init;
        Transformation2D::getTransformKey ( global_points_source, local_points_source, weights, min_key );
        Transformation2D::transform ( global_points_source, local_points_source, transformed_points_init, min_key );
        accuracy = Transformation2D::getAccuracyCharacteristics ( global_points_source, local_points_source, transformed_points_init, min_key, weights );

        //Compute initial q_ee
        typename TItemsList <typename Point1::Type>::Type q_ee ( n_global_points_source, 1 ), q_ee_inv ( n_global_points_source, 1 );

        for ( unsigned int i = 0; i < n_global_points_source; i++ )
        {
                q_ee[i] = ( 1.0 / weights[i] - accuracy.q_xx[i] );
        }

        //Initialize accuracy
        accuracy.std_dev = MAX_FLOAT;

        //Perform iterations until no change of standard deviation
        do
        {
                //Remember old deviation
                std_dev_old = accuracy.std_dev;

                //Invert Q_ee matrix
                for ( unsigned int i = 0; i < n_global_points_source; i++ )
                        q_ee_inv[i] =  1.0 / q_ee[i];

                //Compute weighted transformation
                Container <Point1 *> transformed_points;
                Transformation2D::getTransformKey ( global_points_source, local_points_source, q_ee_inv, min_key );
                Transformation2D::transform ( global_points_source, local_points_source, transformed_points, min_key );

                //Compute accuracy characteristics for weighted transformation
                accuracy = Transformation2D::getAccuracyCharacteristics ( global_points_source, local_points_source, transformed_points, min_key, q_ee_inv );

                //Remove outliers using Danish method
                for ( unsigned int i = 0; i < n_global_points_source; i++ )
                {
                        //Compute new Q_ee
                        q_ee[i] = ( 1.0 / weights[i] - accuracy.q_xx[i] );

                        //Danish method: decrease weight for outliers
                        if ( ( accuracy.res[i].res_xy > 2.0 * accuracy.std_dev * sqrt ( q_ee[i] ) ) /*&& (accuracy.res[i].res_xy > 200 )*/ )
                        {
                                weights[i] *=  exp ( - ( accuracy.res[i].res_xy / ( 2.0 * accuracy.std_dev * sqrt ( q_ee[i] ) ) ) );
                        }
                }

        }
        while ( fabs ( std_dev_old - accuracy.std_dev ) > MIN_DEV_DIFFERENCE );

        //Create best pairs and add them to the key: removing outliers having decreased weights
        for ( unsigned int i = 0; i < n_global_points_source; i++ )
        {
                if ( weights[i] == 1.0 )
                {
                        min_pairs.push_back ( std::make_pair ( accuracy.res[i].res_xy, i ) );
                }
        }

        //Rearrange points: add best points to both dest lists
        Transformation2D::rearrangePoints ( global_points_source, local_points_source, global_points_dest, local_points_dest, min_pairs );

        //Compute transformation key using non-weighted key
        Transformation2D::getTransformKey ( global_points_dest, local_points_dest, min_key );
}


template <typename Point1, typename Point2, typename FunctionJ>
void Outliers::findOutliersME ( const Container <Point1 *> &global_points_source, const Container <Point2 *> &local_points_source,
                                Container <Point1 *> &global_points_dest, Container <Point2 *> &local_points_dest, Matrix <typename Point1::Type> &X, FunctionJ function_j, typename TDevIndexPairs <typename Point1::Type>::Type & min_pairs, const unsigned int MAX_ITER )
{
        //Removing outliers using the M-estimated
        //Outliers are dectected using Huber function
        /*
	const unsigned int n = global_points_source.size();

	//Create matrix Y
	Matrix <typename Point1::Type> Y( 2 * n, 1), A( 2 * n + 1), A0( 2 * n + 1), E ( 2 * n + 1 ), W (2 * n, 2 * n , 0, 1);

	for ( unsigned int i = 0; i < n; i++  )
	{
		Y(i, 0) = global_points_source[i].getX();
		Y(i + n, 0) = global_points_source[i].getY();
		A0(i, 0) = local_points_source[i].getX();
		A0(i + n, 0) = local_points_source[i].getY();
	}

	//Initialize solution
	Matrix <typename Point1::Type> XRES = X;

	//Compute Jacobian matrix
	Matrix <typename Point1::Type> J (2 * n, 6);
	function_j ( X, J );

	//Initialize residuals
	E = ( Y - A0 );

	//Perform iterations
	for ( unsigned int i = 0; i < 10; i++ )
	{
		const Matrix <typename Point1::Type> Var = trans ( E ) * E;
		const typename Point1::Type s_n = sqrt ( Var(0,0) / n );

		//Compute new weights of observations
		for ( unsigned int j = 0; j < 2 * n; j++ )
		{
			//Outside the interval
			typename Point1::Type psi;
			if ( fabs( E( i, 0 ) ) > 1.5 * sigma )
			{
				psi = 1.5 * ( ( E ( i, 0) > 0 ) - ( E ( i, 0 ) < 0 ) ) / s_n;
			}

			//Inside interval
			else
			{
				psi =  E ( i, 0) / s_n;
			}

			//Compute weight
			W(j,j) = psi / ( E ( i, 0) / s_n );
		}

		//Solve the itertively re-weighted least squares
		Matrix <typename Point1::Type> DX = inv ( trans(J) * W * J) * trans(J) * W * YY;

		XRES = XRES + DX;

		//Compute the residuals
		E = ( Y - A0 ) - J * DX;
	}
	
        */
}


template <typename T>
void Outliers::createKBestPairsOfPoints ( const TAccuracyCharacteristics <T> &deviations, typename TDevIndexPairs <T>::Type & point_pairs, const float perc_ratio )
{
        //Create pairs of all points, sort by standard deviation and erase n-k points with greatest outliers
        const unsigned int n = deviations.res.size();

        //Clear old pairs
        point_pairs.clear();

        //Create new pairs
        for ( unsigned int j = 0; j < n ; j++ )
        {
                point_pairs.push_back ( std::make_pair ( deviations.res[j].res_xy, j ) );
        }

        //Sort pairs by standard deviation
        std::sort ( point_pairs.begin(), point_pairs.end(), sortPointPairsByResiduals <T> () );

        //Let only first k-items
        point_pairs.erase ( point_pairs.begin() + n * perc_ratio, point_pairs.end() );
}


#endif

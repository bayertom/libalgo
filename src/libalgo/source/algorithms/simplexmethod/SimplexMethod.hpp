// Description: Downhill simplex optimalization method

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


#ifndef SimplexMethod_HPP
#define SimplexMethod_HPP

#include "FAnalyzeProjV2S.h"

//Set namespace
using namespace MatrixOperations;


template <typename T, typename Function>
T SimplexMethod::NelderMead ( Function function, const Matrix <T> &XMIN, const Matrix <T> &XMAX, Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, unsigned int &iterations, const T max_error, const unsigned int max_iterations, std::ostream * output )
{
        //Compute Nelder-Mead method for a function
        //Algorithm based on Lagarias, Reeds, Wright, 1998, SIAM
        const T Rho = 1.0, Xi = 2.0, Gama = 0.5, Sigma = 0.5;

        //Create random simplex
        Matrix <T> XX = createRandSimplex ( XMIN, XMAX );

        //Get dimensions
        unsigned int m = XX.rows(), n = XX.cols(), m1 = W.rows();

        //Create matrices for the simplex operations
        Matrix <unsigned int> IX ( m, 1 );
        Matrix <T> VV ( m, 1 ), VR ( 1, 1 ), VS ( 1, 1 ), VE ( 1, 1 ), VCO ( 1, 1 ), VCI ( 1, 1 ), YR ( m1, 1 ), YS ( m1, 1 ),
		YE(m1, 1), YCO(m1, 1), YCI(m1, 1), WR(m1, m1), WS(m1, m1), WE(m1, m1), WCO(m1, m1), WCI(m1, m1);
	Matrix <T> YBEST(m1, 1), VBEST(1, 1), WBEST(m1, m1);

        //Set iterations to 0
        iterations = 0;

        //Compute initial function values
        function ( XX, Y, VV, W, true );

        //Sort residuals in ascending order
        sort ( VV, IX );

        //Change order of XX rows: row i to row IX (i)
        Matrix <T> XS ( m, n );

        for ( unsigned int i = 0; i < m; i++ )
                XS.row ( XX.row ( IX ( i, 0 ) ), i );

        //Assign sorted matrix to XX
        XX = XS;

        //Perform Nelder-Mead algorithm
        do
        {
                //Compute centroid
                Matrix <T> XC = sumCols ( XX ( 0, n - 1, 0, n - 1 ) ) * 1.0 / n;

                //Compute reflection point and residuals
                Matrix <T> XR = XC * ( 1.0 + Rho ) - XX ( n, n, 0, n - 1 ) * Rho;

                //Reflection into the search space
                reflection ( XMIN, XMAX, n, XR );

                //Compute residuals
                function ( XR, YR, VR, WR, true );
                const T  fr = VR ( 0, 0 );

                //A reflection point acceptable
                if ( ( VV ( 0, 0 ) <= fr ) && ( fr < VV ( n - 1, 0 ) ) )
                {
                        XX.submat ( XR, n, 0 );
                        VV ( n, 0 ) = fr;
                }

                //Expansion of the simplex
                else if ( fr < VV ( 0, 0 ) )
                {
                        //Compute expanded point
                        Matrix <T> XE = XC * ( 1.0 + Rho * Xi ) - XX ( n, n, 0, n - 1 ) * Rho * Xi;

                        //Reflection into the search space
                        reflection ( XMIN, XMAX, n, XE );

                        //Compute residuals
                        function ( XE, YE, VE, WE, true );
                        const T  fe = VE ( 0, 0 );

                        //An expanded point is acceptable
                        if ( fe < fr )
                        {
                                XX.submat ( XE, n, 0 );
                                VV ( n, 0 ) = fe;
                        }

                        //An expanded point is not acceptable (use a reflected one)
                        else
                        {
                                XX.submat ( XR, n, 0 );
                                VV ( n, 0 ) = fr;
                        }
                }

                //Outside contraction of the simplex
                else if ( ( VV ( n - 1, 0 ) <= fr ) && ( fr < VV ( n , 0 ) ) )
                {
                        //Compute outside contracted point
                        Matrix <T> XCO = XC * ( 1.0 + Rho * Gama ) - XX ( n, n, 0, n - 1 ) * Rho * Gama;

                        //Reflection into the search space
                        reflection ( XMIN, XMAX, n, XCO );

                        //Compute residuals
                        function ( XCO, YCO, VCO, WCO, true );
                        const T  fco = VCO ( 0, 0 );

                        //An outside contracted point is acceptable
                        if ( fco < VV ( n , 0 ) )
                        {
                                XX.submat ( XCO, n, 0 );
                                VV ( n, 0 ) = fco;
                        }

                        //An outside contracted point is not acceptable: shrink a simplex
                        else
                        {
                                shrink ( function, W, XX, XMIN, XMAX, Y, VV, Sigma );

                                //Increment iterations
                                iterations ++;
                        }
                }

                //Inside contraction of the simplex
                else
                {
                        //Compute outside contracted point
                        Matrix <T> XCI = XC * ( 1.0 - Gama ) + XX ( n, n, 0, n - 1 ) * Gama;

                        //Reflection into the search space
                        reflection ( XMIN, XMAX, n, XCI );

                        //Compute residuals
                        function ( XCI, YCI, VCI, WCI, true );
                        const T  fci = VCI ( 0, 0 );

                        //An inside contracted point is acceptable
                        if ( fci < VV ( n , 0 ) )
                        {
                                XX.submat ( XCI, n, 0 );
                                VV ( n, 0 ) = fci;
                        }

                        //An inside contracted point is not acceptable: shrink a simplex
                        else
                        {
                                shrink ( function, W, XX, XMIN, XMAX, Y, VV,  Sigma );

                                //Increment iterations
                                iterations ++;
                        }
                }

                //Perform perturbation
                //perturbation ( A, B, m, XX );

                //Sort residuals in ascending order
                sort ( VV, IX );

                //Change order of XX rows: row i to row IX (i)
                for ( unsigned int i = 0; i < m; i++ )
                        XS.row ( XX.row ( IX ( i, 0 ) ), i );

                //Assign sorted matrix to XX
                XX = XS;

                //Increment iterations
                iterations ++;

                //std::cout << iterations << " ";
                //XX.print();

		//VV.print();
                //std::cout << iterations << " ";

		Matrix <T> XBEST = XX(0, 0, 0, n - 1);
		//function(XBEST, YBEST, VBEST, W, false);  //OK, EQDC
		function(XX, YR, VV, W, false);  //OK, MERC

                if ( iterations % 100 == 0 )
                {
                        std::cout.flush();
                        std::cout << ".";
			//XX.print();
                }
        }
        while ( ( iterations < max_iterations ) && ( fabs ( VV ( 0, 0 ) - VV ( n, 0 ) ) > max_error ) );

        //std::cout << "iter=" << iterations;

        //Get minimum
        X = XX ( 0, 0, 0, n - 1 );

        //Compute residuals for the found solution
        function ( X, Y, V, W, true );

        //X.print();
        //V.print();

        //Print residuals
        //*output << " X:"; X.print ( output );
        //*output << "Residuals: "; V.print ( output );
        //*output << "Iterations: " << iterations << '\n';

        //Return squares of residuals
        return (T) norm ( trans ( V ) * W * V ) ;
}


template <typename T>
Matrix <T> SimplexMethod::createRandSimplex ( const Matrix <T> &XMIN, const Matrix <T> &XMAX )
{
        //Create random simplex
        const unsigned int dim = XMIN.cols();

        //Create matrix for simplex
        Matrix <T> XX ( dim + 1, dim );

        //Initialize random number generator
        srand ( ( unsigned ) time ( 0 ) );

	//Compute difference
	const Matrix <T> DX = XMAX - XMIN;

        //Create random simplex
        for ( unsigned int i = 0; i < dim + 1; i++ )
        {
                for ( unsigned int j = 0; j < dim; j++ )
                {
                        const T rand_num = XMIN ( 0, j ) + DX ( 0, j ) * rand() / ( RAND_MAX + 1.0 );
			if (rand_num < -1000)
				double ccc = -6;
                        XX ( i, j ) = rand_num;
                }
        }

        return XX;
}



template <typename T, typename Function>
void SimplexMethod::shrink ( Function function, Matrix <T> &W, Matrix <T> &XX, const Matrix <T> &XMIN, const Matrix <T> &XMAX, Matrix <T> &Y,  Matrix <T> &VV, const T Sigma )
{
        //Shrink a simplex
        const unsigned int m = XX.rows(), n = XX.cols(), m1 = W.rows();

        //Create matrices
        Matrix <T> V1 ( 1, 1 ), VSH ( 1, 1 ), Y1 ( m1, 1 ), YSH ( m1, 1 ), W1 ( m1, m1 ), WSH ( m1, m1 );

        //Get first point of the simplex (best)
        Matrix <T> X1 = XX ( 0, 0, 0, n - 1 );

        //Compute residuals
        function ( X1, Y1, V1, W1, true );

        //Actualize VV matrix
        VV ( 0, 0 ) = V1 ( 0, 0 );

        //Shrink remaining simplex points
        for ( unsigned int i = 1; i < n + 1 ; i++ )
        {
                //Compute shrink point
                Matrix <T> XSH = XX ( 0, 0, 0, n - 1 ) + ( XX ( i, i, 0, n - 1 ) - XX ( 0, 0, 0, n - 1 ) ) * Sigma;

                reflection ( XMIN, XMAX, n, XSH );

                //Compute residuals
                function ( XSH, YSH, VSH, WSH, true );

                //Set submatrix
                XX.submat ( XSH, i, 0 );
                VV ( i, 0 ) = VSH ( 0, 0 );
        }
}


template <typename T>
void SimplexMethod::reflection ( const Matrix <T> &XMIN, const Matrix <T> &XMAX, const unsigned int dim, Matrix <T> &X )
{
        //Set elements of the simplex into the search space represented by the n-dimensional cuboid
        for ( unsigned int j = 0; j < dim; j++ )
        {
		while ((X(0, j) < XMIN(0, j)) || (X(0, j) > XMAX(0, j)))
                {
			//XMIN == XMAX
			if (XMAX(0, j) - XMIN(0, j) < MAX_FLOAT_OPER_ERROR )
			{
				X(0, j) = XMIN(0, j);
				break;
			}

			//Left form the lower bound
			else if (X(0, j) > XMAX(0, j))
			{
				X(0, j) = 2 * XMAX(0, j) - X(0, j);
			}

			//Right to the upper bound
			else if (X(0, j) < XMIN(0, j))
			{
				X(0, j) = 2 * XMIN(0, j) - X(0, j);
			}
                }
        }
	
}

#endif

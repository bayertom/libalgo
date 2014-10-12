// Description: Non Linear Least Squares algorithms

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


#ifndef NonLinearLeastSquares_HPP
#define NonLinearLeastSquares_HPP

#include <cmath>

#include "libalgo/source/exceptions/ErrorBadData.h"


//Set namespace
using namespace MatrixOperations;


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::GN ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, const T max_error, const unsigned short max_iterations )
{
        //Solving non Linear Least Squares using the Gaussian-Newton method
        //Faster than backtracking method
        unsigned short iterations = 0;
        T cost_old = MAX_FLOAT;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), V ( m, 1 ), dX ( n, 1 ), dX2 ( n, 1 );

        //Assign matrix
        Matrix <T> X_Old = X;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Perform iterations
        while ( iterations < max_iterations )
        {
                //Compute new J matrix: linearize solution using the Taylor approximation
                function_j ( X, J );

                //Stop computation
                if ( sum2 ( trans ( J ) * W * V ) < max_error )
                        break;

                //Compute Minimum Weighteed Least Squares using qr decomposition
                //Jacobian J = [ d_R, d_latp, d_lonp, d_lat0, d_lon0, d_dx, d_dy]
                //dX = mlsqr ( J, W, V ) * ( -1.0 );
                dX = pinv1 ( trans ( J ) * W * J ) * trans ( J ) * W * V * ( -1.0 );

                //Compute new X
                X = X + dX;

                //Compute new V matrix
                function_v ( X, Y, V, W );

                X.print();

                //Increment iterations
                iterations++;
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::GND ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y,  Matrix <T> &V, unsigned short &iterations, const T max_error, const unsigned short max_iterations, std::ostream * output )
{
        //Solving Non linear Least Squares using the damped Gaussian-Newton method
        T cost_old = MAX_FLOAT;

        //Set iterations to 0
        iterations = 0;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), V2 ( m, 1 ), dX ( n, 1 ), J_new ( m, n );

        //Assign matrix
        Matrix <T> Y2 = Y;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute matrices
        function_j ( X, J );
        Matrix <T> G = trans ( J ) * W * V;
        Matrix <T> F = trans ( V ) * W * V;

        //Perform iterations
        while ( iterations < max_iterations )
        {
                //Compute new J matrix: linearize solution using the Taylor approximation
                function_j ( X, J );

                //Compute the direction 
                dX = mlsqr ( J, W, V ) * ( -1.0 );

                //Compute new trial X
                Matrix <T> X2 = X + dX;

                //Compute new trial V matrix (residuals)
                function_v ( X2, Y2, V2, W );

                //Apply a bisection to find a damping factor t
                const T alpha = 0.1 ;
                float t = 1.0;

                while ( ( sum2 ( V2 ) > sum2 ( V ) + sum ( trans ( V ) * J *  dX * t * alpha * 2.0 ) )  && ( t > 1.0e-20 ) )
                {
                        //Step t bisection
                        t /= 2;

                        //Compute new X2
                        X2 = X + dX * t;

                        //Compute new V matrix: do not change parameters in one iteration step
                        function_v ( X2, Y2, V2, W, false );
                }

                //Compute new X using damping factor t
                X = X + dX * t;

                //Compute V matrix
                function_v ( X, Y, V, W );

                //Compute new J matrix
                function_j ( X, J );

                //Compute new residuals and gradient
                Matrix <T> F_new = trans ( V ) * W * V;
                Matrix <T> G_new = trans ( J ) * W * V;

                //Terminal condition
                if ( ( norm ( G ) < max_error ) || ( fabs ( F_new ( 0, 0 ) - F ( 0, 0 ) ) < 1.0e-10 * std::max ( 1.0 , F ( 0, 0 ) ) ) ||
                                ( F ( 0, 0 ) <  max_error ) )
                        break;

                //Assign old values
                F = F_new;
                G = G_new;

                //Increment iterations
                iterations++;


                //std::cout << iterations << " ";
                //dX.print();
                //X.print();
                //F.print();
                //V.print();
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        std::cout << "iter:" << iterations << '\n';
        X.print();
        X.print ( output );
        //T rr = norm ( trans (V) * V );
        //std::cout << "res =" << rr << '\n';

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::LM ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, const T max_error, const unsigned short max_iterations )
{
        //Solving non Linear Least Squares with using the Levenberq-Marquard algorithm
        unsigned short iterations = 0;
        T cost_old = MAX_FLOAT;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), V ( m, 1 ), V2 ( m, 1 ), dX ( n, 1 ), E ( n, n, 0.0, 1.0 );

        //Initialize matrices
        Matrix <T> X_new = X, Y_new = Y, V_new = V, W_new = W;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute initial J matrix
        function_j ( X, J );

        Matrix <T> B = trans ( J ) * W * J;
        Matrix <T> G = trans ( J ) * W * V;

        //Get diagonal matrix D
        Matrix <T> D = diag ( diag ( B ) );
        //Matrix <T> D = diag( diag ( trans ( J ) * J ) );

        //Compute initial solution
        const T Ro_min = 0.25, Ro_max = 0.75;
        T lambda = 1.0, lambdac = 0.75;

        //Compute s
        T s = norm ( trans ( V ) * W * V );

        //Perform iterations
        while ( ( iterations < max_iterations ) )
        {
                //Compute Minimum Weighteed Least Squares using SVD
                dX = pinvs ( B + D * lambda ) * G ;

                //Compute new trial X
                X_new = X - dX;

                //Compute new Y matrix and residual matrix V
                function_v ( X_new, Y_new, V_new, W_new );

                //Stop computation
                if ( sum2 ( trans ( J ) * W * V ) < max_error )
                        break;

                //Compute predicted solution
                const T s_new = norm ( trans ( V_new ) * W * V_new );
                const Matrix <T> PS = trans ( dX ) * ( G * 2.0 - B * dX );
                const T p_sol = PS ( 0, 0 );

                //Compute Ro
                const T Ro = fabs ( ( s - s_new ) / p_sol );

                //New cost is outside interval, greater than Ro
                if ( Ro > Ro_max )
                {
                        //Decrease parameter
                        lambda /= 2;

                        if ( lambda < lambdac ) lambda = 0.0;
                }

                //New cost is outside interval, lower than Ro
                else if ( Ro < Ro_min )
                {

                        //Compute muliplication parameter
                        const Matrix <T> dVV = trans ( dX )  * G ;
                        T nu = 2.0 + ( s_new - s ) / dVV ( 0, 0 );

                        //And correct if outside interval
                        if ( nu < 2 ) nu = 2.0;
                        else if ( nu > 10.0 ) nu = 10.0;

                        //Compute new lambda
                        if ( lambda == 0.0 )
                        {
                                Matrix <T> test = pinvs ( B );
                                lambdac = 1.0 / max ( abs ( diag ( pinvs ( B ) ) ) );
                                lambda  = lambdac;
                                nu /= 2;
                        }

                        //And multiply lambda
                        lambda *= nu;
                }

                //We found a better solution, remember it
                if ( s_new < s )
                {
                        s = s_new;
                        X = X_new;
                        V = V_new;
                        W = W_new;

                        //Compute new J matrix: linearize solution using the Taylor approximation
                        function_j ( X, J );

                        B = trans ( J ) * W * J;
                        G = trans ( J ) * W * V;
                }

                //Increment iterations
                iterations++;

                X.print();

                std::cout << iterations << " ";
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        //X.print();

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::LM2 ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, const T max_error, const unsigned short max_iterations )
{
        //Solving non Linear Least Squares with using the Levenberq-Marquard algorithm
        unsigned short iterations = 0;
        T cost_old = MAX_FLOAT, lambda = 0.01, CR2, lambda_UP_fac = 11.0, lambda_DN_fac = 9.0;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), V ( m, 1 ), dX ( n, 1 ), E ( n, n, 0.0, 1.0 );

        //Initialize matrices
        Matrix <T> X2 = X, Y2 = Y, V2 = V, W2 = W;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute initial J matrix
        function_j ( X, J );

        //Compute new criterium
        CR2 = norm ( trans ( V ) * W * V );

        //Compute B and BB matrices
        Matrix <T> B = trans ( J ) * W * J;
        Matrix <T> G = trans ( J ) * W * V;

        //Assign old values
        T CR2_old = CR2;

        //Perform iterations
        while ( ( iterations < max_iterations ) )
        {
                Matrix <T> D = diag ( diag ( B ) );

                //Compute Minimum Weighteed Least Squares using qr decomposition
                dX = pinv1 ( B + D * lambda ) * G * ( -1.0 );

                //Compute new trial X
                X2 = X + dX;

                //Compute new Y matrix and residual matrix V
                function_v ( X2, Y2, V2, W2 );

                //Stop computation
                if ( sum2 ( trans ( J ) * W * V2 ) < max_error )
                        break;

                //Compute criteria (predicted solution)
                const T CR2_try = norm ( trans ( V2 ) * W * V2 );

                //Compute Ro
                const T Rho1 = CR2 - CR2_try;
                const Matrix <T> Rho2 = trans ( dX ) * ( dX * lambda + G ) * 2.0;
                const T Rho = fabs ( ( CR2 - CR2_try ) / Rho2 ( 0, 0 ) );

                //We found a better solution
                if ( Rho > 1.0e-2 )
                {
                        CR2 = CR2_try;
                        X = X2;

                        //Compute Jacobi matrix
                        function_j ( X, J );

                        //Compute new Y matrix and residual matrix V
                        function_v ( X, Y, V, W );

                        //Compute new B and G matrices
                        B = trans ( J ) * W * J;
                        G = trans ( J ) * W * V;

                        //Update lambda
                        lambda = min ( lambda * lambda_UP_fac, 1.e7 );
                }

                else
                {
                        //Actualize criterium
                        CR2 = CR2_old;

                        //Update lambda
                        lambda = max ( lambda / lambda_DN_fac, 1.e-7 );
                }

                X.print();

                //Increment iterations
                iterations++;

                std::cout << iterations << " ";
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        //X.print();

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::LM3 ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, const T max_error, const unsigned short max_iterations )
{
        //Solving non Linear Least Squares with using the Levenberq-Marquard algorithm
        unsigned short iterations = 0;
        T cost_old = MAX_FLOAT, lambda = 0.01, CR2, lambda_UP_fac = 11.0, lambda_DN_fac = 9.0;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), V ( m, 1 ), dX ( n, 1 ), E ( n, n, 0.0, 1.0 );

        //Initialize matrices
        Matrix <T> X2 = X, Y2 = Y, V2 = V, W2 = W, XN = X;

        //Compute initial J matrix
        function_j ( X, J );

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute residuals
        T f = norm ( trans ( V ) * W * V );

        //Set parameters
        T nu0 = 0.001, nfun = 1, ngrad = 1, numh = 0, numf = 1, numg = 1;

        //Compute gradient
        Matrix <T> G =  trans ( J ) * W * V ;

        //Compute parameter
        T nu = norm ( G ), nup = nu;

        //Perform iterations
        while ( ( iterations < max_iterations ) )
        {
                //Compute Minimum Weighteed Least Squares using qr decomposition
                dX = pinv1 ( trans ( J ) * W * J + E * nu ) * G;

                //Compute new trial X
                X2 = X - dX;

                //Compute trust region
                int itr = 0;
                lmTrust ( function_v, X, X2, f, J, G, nu, nu0, XN, nup, itr );

                if ( itr  > 20 )
                {

                }

                //Assign results
                X = XN;
                nu = nup;
                numf = numf + itr;

                //At least one iteration performed
                if ( itr > 1 )
                {
                        //Compute initial J matrix
                        function_j ( X, J );

                        //Compute initial V matrix (residuals)
                        function_v ( X, Y, V, W );

                        //Compute residuals
                        f = norm ( trans ( V ) * W * V );

                        //Compute gradiend
                        G = trans ( J ) * W * V ;

                        numf++; numg++;
                }

                //Stop computation
                //if ( sum2 ( trans ( J ) * W * V2 ) < max_error )
                //        break;

                X.print();

                //Increment iterations
                iterations++;

                std::cout << iterations << " ";
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        //X.print();

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}


template <typename T, typename FunctionV>
void NonLinearLeastSquares::lmTrust ( FunctionV function_v,  const Matrix <T> X, Matrix <T> X2, const T &f, const Matrix <T> &J, const Matrix <T> &G, T nu, const T nu0, Matrix <T> &XN, T &nup, int &itr )
{
        int numft = 0, numgt = 0;
        T mu0 = 0.1,  mulow = 0.25,  muhigh = 0.75,  mup = 2.0,  mdown = 0.5;

        //Initialize solution
        Matrix <T> Z = X;

        //Set iterations to 0
        itr = 0;

        while ( ( norm ( Z - X ) < 0.001 ) && ( itr < 20 ) )
        {
                const unsigned short m = J.rows(), n = J.cols();
                Matrix <T> Y2 ( m, 1 ), V2 ( m, 1 ), W ( m, m ), W2 ( m, m ), E ( n, n, 0.0, 1.0 );

                //Compute V matrix (residuals)
                function_v ( X2, Y2, V2, W2 );

                //Compute residuals
                const T f2 = norm ( trans ( V2 ) * W2 * V2 );

                Matrix <T> dXX = X2 - X;
                T ared = f - f2;
                Matrix <T> GDX = trans ( G ) * dXX;
                T pred = - 0.5 * GDX ( 0, 0 );

                //Compute ratio
                T rat = fabs ( ared / pred );

                //Case 1
                if ( rat < mu0 )
                {
                        nu = std::max ( nu * mup, nu0 );

                        //Compute new dX
                        Matrix <T> dX = pinv1 ( trans ( J ) * W * J + E * nu ) * G ;

                        X2 = X - dX;
                }

                else if ( rat < mulow )
                {
                        Z = X2;
                        nu = std::max ( nu * mup, nu0 );
                }

                else
                {
                        Z = X2;

                        if ( rat > muhigh )
                        {
                                nu = mdown * nu;

                                if ( nu < nu0 ) nu = 0.0;
                        }
                }

                itr++;
        }

        //Actualize values
        nup = nu;
        XN = Z;
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::LM4 ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, unsigned short &iterations, const T max_error, const unsigned short max_iterations, std::ostream * output )
{
        //Solving non Linear Least Squares with using the Levenberq-Marquard algorithm
        //Initial mu value from  Kelley's "Iterative Methods for Optimization", SIAM 1999.
        T cost_old = MAX_FLOAT;

        //Set iterations to 0
        iterations = 0;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), dX ( n, 1 ), E ( n, n, 0.0, 1.0 );

        //Initialize matrices
        Matrix <T> Y_new = Y, V_new = V, W_new = W, J_new = J;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute initial J matrix
        function_j ( X, J );

        //Compute matrices
        Matrix <T> H = trans ( J ) * W * J;
        Matrix <T> G = trans ( J ) * W * V;

        //Compute norm
        Matrix <T> F = trans ( V ) * W * V;
        T f_n = norm ( F );

        //Initialize Levenberg parameters
        const T lambda_min = 0.001, mu0 = 0.1, mulow = 0.25, muhigh = 0.75, omup = 2, omdown = 0.5;
        T lambda =  norm ( trans ( G ) * G ) / norm ( trans ( G ) * H  * G );
        //lambda = 0.0001 * max ( diag ( H ) );

        //Perform iterations
        while ( iterations < max_iterations )
        {
                //Compute Minimum Weighteed Least Squares using qr decomposition
                dX = pinv1 ( H + E * lambda ) * G * ( -1.0 );

                //Compute new trial X
                Matrix <T> X_new = X + dX;

                //Compute new V matrix and residual matrix V
                function_v ( X_new, Y_new, V_new, W_new );

                //Compute new J matrix
                function_j ( X_new, J_new );

                //Compute predicted solution parameters
                const Matrix <T> G_new = trans ( J_new ) * W_new * V_new;
                const Matrix <T> H_new = trans ( J_new ) * W_new * J_new;
                const Matrix <T> F_new = trans ( V_new ) * W_new * V_new ;
                const T f_n_new = norm ( F_new );

                //Compute reduction parameters from actual and predicted solution
                const T RA = f_n - f_n_new;
                Matrix <T> RP =  trans ( G ) * dX * ( -1.0 ) - trans ( dX ) * H_new * dX * 0.5;
                const T Ro = RA / RP ( 0, 0 );

                //Update trust region radius: decrease radius
                if ( Ro < mulow )
                {
                        lambda = std::max ( omup * lambda, lambda_min ) ;
                }

                //Update trust region radius: increase radius
                else if ( Ro > muhigh )
                {

                        lambda = omdown * lambda;
                        //lambda = std::min( omup * lambda, 100.0 ) ;
                }

                //Lambda is too small, do not use
                //if ( lambda < lambda_min )
                //        lambda = 0;
                //J new solution acceptable: Ro is not too small
                if ( Ro > mu0 )
                {
                        //Stop computation
                        if ( ( norm ( G ) < max_error ) || ( fabs ( F_new ( 0, 0 ) - F ( 0, 0 ) ) < 1.0e-10 * std::max ( 1.0 , F ( 0, 0 ) ) ) ||
                                        ( F ( 0, 0 ) <  max_error ) )
                                break;

                        //Assign matrices
                        F = F_new;
                        X = X_new;
                        f_n = f_n_new;
                        G = G_new;
                        H = H_new;
                }

                //Increment iterations
                iterations++;
                //std::cout << lambda;
                //std::cout << iterations << " ";

                //X.print();
        }

        //Compute final values in V
        function_v ( X, Y, V, W );
        //V.print();
        X.print();
        X.print ( output );

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::LM5 ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, unsigned short &iterations, const T max_error, const unsigned short max_iterations, std::ostream * output )
{
        //Solving non Linear Least Squares with using the Levenberq-Marquard algorithm
        //Initial mu value from  Kelley's "Iterative Methods for Optimization", SIAM 1999.
        //Algorithm by Bruun Nielsen
        T cost_old = MAX_FLOAT;

        //Set iterations to 0
        iterations = 0;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), dX ( n, 1 ), E ( n, n, 0.0, 1.0 );

        //Initialize matrices
        Matrix <T> Y_new = Y, V_new = V, W_new = W, J_new = J;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute initial J matrix
        function_j ( X, J );

        //Compute matrices
        Matrix <T> H = trans ( J ) * W * J;
        Matrix <T> G = trans ( J ) * W * V;
        Matrix <T> F = trans ( V ) * W * V * 0.5;

        //Compute norm
        T ng = norm ( G );

        //Initialize Levenberg parameters
        T nu = 2.0;
        T mu = std::max ( norm ( trans ( G ) * G ) / norm ( trans ( G ) * H  * G ), 0.001 );
        mu = 0.0001 * max ( diag ( H ) );

        //Perform iterations
        while ( iterations < max_iterations )
        {
                //Compute Minimum Weighteed Least Squares using qr decomposition
                dX = pinv1 ( H + E * mu ) * G * ( -1.0 );

                //*output << "iter = " << iterations;
                //std::cout << "iter = " << iterations;
                //dX.print();
                //dX.print(output);

                //Compute new trial X
                Matrix <T> X_new = X + dX;

                //X_new.print();
                //std::cout << mu << ' ';
                //Compute dL
                Matrix <T> dLM = ( trans ( dX ) * ( dX * mu - G ) ) * 0.5;
                const T dL = dLM ( 0, 0 );

                //Compute new V matrix and residual matrix V
                function_v ( X_new, Y_new, V_new, W_new );

                //Compute new J matrix
                function_j ( X_new, J_new );

                //Compute new residuals and difference
                const Matrix <T> F_new =  trans ( V_new ) * W_new * V_new * 0.5;
                const T dF = F ( 0, 0 ) - F_new ( 0, 0 );

                //Update solution: new solution is better
                if ( ( dL > 0 ) && ( dF > 0 ) )
                {
                        //Stop computation
                        //std::cout << "g-norm " << norm(G) << '\n';
                        if ( ( norm ( G ) < max_error ) || ( fabs ( F_new ( 0, 0 ) - F ( 0, 0 ) ) < 1.0e-10 * std::max ( 1.0 , F ( 0, 0 ) ) ) ||
                                        ( F ( 0, 0 ) <  max_error ) )
                                break;

                        //Update solution
                        X = X_new;
                        F = F_new;
                        J = J_new;
                        V = V_new;
                        W = W_new;

                        //Actualize matrices
                        H = trans ( J ) * W * J;
                        G = trans ( J ) * W * V;
                        ng = norm ( G );

                        //Actualize mu: dammping criterium
                        mu = mu * std::max ( 1.0 / 3, 1.0 - pow ( ( 2.0 * dF / dL - 1.0 ), 3 ) );
                        nu = 2.0;
                        //X.print(output);
                        //X.print();
                }

                //Old solution is better, decrese mu, nu
                else
                {
                        mu = mu * nu;
                        nu = 2 * nu;
                }

                //Increment iterations
                iterations++;

                //std::cout << lambda;
                //std::cout << iterations << " ";
                //std::cout << "res= \n";
                // X.print();
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        //T rr = log ( norm ( trans (V) * V ));
        //std::cout << "res =" << rr << '\n';

        //std::cout << "iter:" << iterations << '\n';
        X.print();
        //std::cout << "residuals";
        //V.print();
        X.print ( output );

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::LM6 ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y, Matrix <T> &V, unsigned short &iterations, const T max_error, const unsigned short max_iterations, std::ostream * output )
{
        //Solving non Linear Least Squares with using the Levenberq-Marquard algorithm
        //Initial mu value from  Kelley's "Iterative Methods for Optimization", SIAM 1999.
        //Algorithm by Bruun Nielsen
        T cost_old = MAX_FLOAT;

        //Set iterations to 0
        iterations = 0;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), dX ( n, 1 ), E ( n, n, 0.0, 1.0 );

        //Initialize matrices
        Matrix <T> Y_new = Y, V_new = V, W_new = W, J_new = J;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute initial J matrix
        function_j ( X, J );

        //Compute matrices
        Matrix <T> H = trans ( J ) * W * J;
        Matrix <T> G = trans ( J ) * W * V;
        Matrix <T> F = trans ( V ) * W * V * 0.5;

        //Compute norm
        T ng = norm ( G );

        //Initialize Levenberg parameters
        T nu = 2.0, max_mu = 1.0e6;
        T mu = std::max ( norm ( trans ( G ) * G ) / norm ( trans ( G ) * H  * G ), 0.001 );
        mu = 0.0001 * max ( diag ( H ) );

        //Perform iterations
        while ( iterations < max_iterations )
        {
                //Compute Minimum Weighteed Least Squares using qr decomposition
                dX = pinv1 ( H + E * mu ) * G * ( -1.0 );

                //*output << "iter = " << iterations;
                //std::cout << "iter = " << iterations;
                //dX.print();
                //dX.print(output);

                //Compute new trial X
                Matrix <T> X_new = X + dX;

                //X_new.print();
                //std::cout << mu << ' ';
                //Compute dL
                Matrix <T> dLM = ( trans ( dX ) * ( dX * mu - G ) ) * 0.5;
                const T dL = dLM ( 0, 0 );

                //Compute new V matrix and residual matrix V
                function_v ( X_new, Y_new, V_new, W_new );

                //Compute new J matrix
                function_j ( X_new, J_new );

                //Compute new residuals and difference
                const Matrix <T> F_new =  trans ( V_new ) * W_new * V_new;
                const T dF = F ( 0, 0 ) - F_new ( 0, 0 );
                const T Ro = dF / dL;

                //New solution is better, decrese mu, nu
                if ( Ro > 0.75 )
                {
                        //Actualize mu: dammping criterium
                        mu = mu * std::max ( 1.0 / 3, 1.0 - pow ( ( 2.0 * dF / dL - 1.0 ), 3 ) );
                        nu = 2.0;
                }

                //Old solution is better, increase mu, nu
                else if ( Ro < 0.25 )
                {
                        mu = std::min ( mu * nu, max_mu );
                        nu = 2 * nu;
                }

                //Update solution
                if ( Ro > 0.1 )
                {
                        //Stop computation
                        if ( ( norm ( G ) < max_error ) || ( fabs ( F_new ( 0, 0 ) - F ( 0, 0 ) ) < 1.0e-10 * std::max ( 1.0 , F ( 0, 0 ) ) ) ||
                                        ( F ( 0, 0 ) <  max_error ) )
                                break;

                        //Update solution
                        X = X_new;
                        F = F_new;
                        J = J_new;
                        V = V_new;
                        W = W_new;

                        //Actualize matrices
                        H = trans ( J ) * W * J;
                        G = trans ( J ) * W * V;
                        ng = norm ( G );

                        //X.print();
                }

                X_new.print();
                //Increment iterations
                iterations++;

                //std::cout << lambda;
                //std::cout << iterations << " ";
                //std::cout << "res= \n";
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        //T rr = log ( norm ( trans (V) * V ));
        //std::cout << "res =" << rr << '\n';

        //std::cout << "iter:" << iterations << '\n';
        X.print();
        //std::cout << "residuals";
        //V.print();
        X.print ( output );

        //Return squares of residuals
        return norm ( trans ( V ) * W * V );
}


template <typename T, typename FunctionJ, typename FunctionV, typename FunctionC>
T NonLinearLeastSquares::BFGSH ( FunctionJ function_j, FunctionV function_v, FunctionC function_c,  Matrix <T> &W, Matrix <T> &X, Matrix <T> &Y,  Matrix <T> &V, unsigned short &iterations, const T alpha, const T nu, const T max_error, const unsigned short max_iterations, const T max_diff, std::ostream * output )
{
        //Solving non Linear Least Squares with using the hybrid BFGS algorithm
        //Algorithm by L Luksan
	//Default values: nu = 0.0001, alpha = 0.0001;
        T cost_old = MAX_FLOAT;

        //Create matrices
        const unsigned short m = W.rows(), n = X.rows();
        Matrix <T> J ( m, n ), V2 ( m, 1 ),  dX ( n, 1 ), E ( n, n, 0.0, 1.0 );

        //Assign matrix
        Matrix <T> Y2 = Y;

        //Compute initial V matrix (residuals)
        function_v ( X, Y, V, W );

        //Compute initial J matrix
        function_j ( X, J );

        //Compute matrices
        Matrix <T> B = trans ( J ) * W * J;
        Matrix <T> G = trans ( J ) * W * V;
        Matrix <T> F = trans ( V ) * W * V;

        //Initialize BFGS parameter
        Matrix <T> B_new = B, G_new = G, F_new = F;

        //Set iterations to 0
        iterations = 0;
	//X.print();
	//J.print();
	//V.print();

        //Perform iterations
        while ( iterations < max_iterations )
        {
                //Compute new dX
		//dX = pinvs ( B ) * G * ( -1.0 );
		dX = mlsqr ( J, W, V ) * ( -1.0 );

		//dX.print();
                //Compute new trial X
                Matrix <T> X2 = X + dX;

                //Compute new trial V matrix (residuals)
                function_v ( X2, Y2, V2, W );

                //Apply back-step using a bisection
                const T t_min = 1.0e-10; 
		T t = 1.0;

                while ( ( sum2 ( V2 ) > sum2 ( V ) + sum ( trans ( V ) * J *  dX ) * t * alpha * 2.0 ) && ( t > t_min ) )
                {
                        //Step t bisection
                        t /= 2;

                        //Compute new X2
                        X2 = X + dX * t;

                        //Compute new V matrix: do not change parameters in one iteration step
                        function_v ( X2, Y2, V2, W, false );
                }

		
		//Compute new X using back-step method
                X = X + dX * t;

                //Compute new V matrix and residual matrix V
                function_v ( X, Y, V, W );

                //Compute new J matrix
                function_j ( X, J );

                //Compute new residuals and gradient
                F_new = trans ( V ) * V;
                G_new = trans ( J ) * W * V;		
		
		/*std::cout << norm(G) << '\n';
		std::cout << fabs(F_new(0, 0) - F(0, 0)) << '\n';
		std::cout << max_diff * std::max(1.0, F(0, 0)) << '\n';
		std::cout << F(0, 0) << '\n';
		std::cout << F_new(0, 0) << '\n';*/
		
		//Terminal condition
                if ( ( norm ( G ) < max_error ) || ( fabs ( F_new ( 0, 0 ) - F ( 0, 0 ) ) < 1.0 * max_diff * std::max ( 1.0 , F ( 0, 0 ) ) ) ||
                                ( F ( 0, 0 ) <  max_error ) )
                        break;
				
                //Compute selection criterium
                const T dF = ( ( F ( 0, 0 ) - F_new ( 0, 0 ) ) / F ( 0, 0 ) );

                //Compute Hessian matrix in a common way
                if ( dF > nu )
                {
                        B_new = trans ( J ) * W * J;
                }

                //Use BFGS method
                else
                {
                        //Compute solution difference
                        Matrix <T> d = dX * t;

                        //Compute gradient difference
                        Matrix <T> y = G_new - G;

                        //Matrix <T> B2 = trans ( J ) * W * J;

                        //Compute new Hessian matrix using BFGS
                        const Matrix <T> D1 =  trans ( y ) * d;
                        const Matrix <T> D2 =  trans ( d ) * B * d;

                        Matrix <T> dB = y * trans ( y ) / D1 ( 0, 0 ) - B * d * trans ( B * d ) / D2 ( 0, 0 );
                        B_new = B + dB;
                }

                //Assign values
                B = B_new;
                G = G_new;
                F = F_new;

		//X.print(output);
		//*output << iterations << '\t' << F_new(0, 0) << '\n';
		//X.print();

                //Increment iterations
                iterations++;
        }

        //Compute final values in V
        function_v ( X, Y, V, W );

        //std::cout << "iter:" << iterations << '\n';
        //X.print(output);
        //X.print(output);
        //V.print();
        //Matrix <T> RES = ( trans ( V ) * V );
        //RES.print();
        //std::cout << "sum " << (sum(W));
        //Return squares of residuals
        return norm ( trans ( V ) * W * V ) / sum ( W );
}


#endif

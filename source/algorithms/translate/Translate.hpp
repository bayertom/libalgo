// Description: Arithmetic parser using postfix notation

// Copyright (c) 2010 - 2011
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


#ifndef Translate_HPP
#define Translate_HPP

#include <stack>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <ctype.h>


#include "libalgo/source/const/Const.h"

#include "libalgo/source/exceptions/ErrorMathInvalidArgument.h"
#include "libalgo/source/exceptions/ErrorMathOverflow.h"
#include "libalgo/source/exceptions/ErrorMathZeroDevision.h"
#include "libalgo/source/exceptions/ErrorMathRange.h"
#include "libalgo/source/exceptions/ErrorMathMatrixSquare.h"
#include "libalgo/source/exceptions/ErrorParse.h"
#include "libalgo/source/exceptions/Error.h"


template <typename T>
T Translate::translateEq ( const char * equation, const T x, const bool print_exception, std::ostream * output )
{
        //Translation equation
        T res = 0;
        char postfix[4096];

        try
        {
                //Translate infix notation to postfix notation
                infixToPostfix ( equation, postfix );

                //Parse equation
                res = parseEquation ( postfix, x, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );

                return res;
        }

        //Throw exception
        catch ( Error & error )
        {
                if ( print_exception )
                {
                        error.printException ( output );
                }

                throw;
        }
}


template <typename T>
T Translate::translateEq ( const char * equation, const T x, const T y, const bool print_exception, std::ostream * output )
{
        //Translation equation
        T res = 0;
        char postfix[4096];

        try
        {
                //Translate infix notation to postfix notation
                infixToPostfix ( equation, postfix );

                //Parse equation
                res = parseEquation ( postfix, x, y, 0, 0, 0, 0, 0, 0, 0, 0, 0 );

                return res;
        }

        //Throw exception
        catch ( Error & error )
        {
                //Print exeption
                if ( print_exception )
                {
                        error.printException ( output );
                }

                //Throw exception
                throw;
        }
}


template <typename T>
T Translate::translateEq ( const char * equation,  const T lat, const T lon, const T R, const T a, const T b, const T lat0, const T lat1, const T lat2, const T lon0, const bool print_exception, std::ostream * output )
{
        //Translation cartographic equation
        T res = 0;
        char postfix[4096];

        try
        {
                //Translate infix notation to postfix notation
                infixToPostfix ( equation, postfix );

                //Parse equation
                res = parseEquation ( postfix, T ( 0 ), T ( 0 ), lat, lon, R, a, b,  lat0, lat1, lat2, lon0 );

                return res;
        }

        //Throw exception
        catch ( Error & error )
        {
                //Print exeption
                if ( print_exception )
                {
                        error.printException();
                }

                //Throw exception
                throw;
        }
}


template <typename T>
T Translate::parseEquation ( const char * equation, const T x, const T y, const T lat, const T lon, const T R, const T a, const T b,  const T lat0, const T lat1, const T lat2, const T lon0 )

{
        //Parse postfix notation
        T result = 0;
        std::stack <T> operands;

        //Process postfix notation
        while ( *equation != '\0' )
        {
                //NUMBER
                if ( isdigit ( *equation ) )
                {
                        char number_text[32];

                        //Find number and convert to T
                        findSequence ( equation, number_text );
                        T number = atof ( number_text );

                        //Move pointer
                        equation += strlen ( number_text );

                        //Add into stack
                        operands.push ( number );
                }

                //Function or variable
                else if ( isalpha ( *equation ) )
                {
                        char function_text[32];
                        findSequence ( equation, function_text );

                        //Move pointer
                        equation += strlen ( function_text );

                        //SIN(x)
                        if ( strcmp ( function_text, "sin" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing: ", "sin(x)" );

                                //Pop one operand
                                T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "sin(x), x > MAX.", op );

                                //Result
                                result = sin ( op * M_PI / 180 );

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //COS(x)
                        else if ( strcmp ( function_text, "cos" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing: ", "cos(x)." );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow  <T> ( "ErrorMathOverflow: can not parse equation: ", "cos(x), x > MAX." );

                                //Result
                                result = cos ( op * M_PI / 180 );

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //TAN(x)
                        else if ( strcmp ( function_text, "tg" ) == 0 || strcmp ( function_text, "tan" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing: ", "tan(x)" );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "tan(x), x > MAX.", op );
                                else if ( fabs ( cos ( op * M_PI / 180 ) ) < MIN_FLOAT )
                                {
                                        //throw ErrorMathMatrixSquare ( "ERangError: can not parse equation ", "tan(x), x = ", 1, 1);
                                        throw ErrorMathRange <T> ( "ErrorMathRange: can not parse equation ", "tan(x), x = ", op );
                                }

                                //Result
                                result = tan ( op * M_PI / 180 );

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //COT(x)
                        else if ( strcmp ( function_text, "cotg" ) == 0 || strcmp ( function_text, "cot" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing: ", "cotg(x)" );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exceptions
                                if ( op > MAX_FLOAT )
                                {
                                        throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "cotg(x), x > MAX." );
                                }

                                else if ( fabs ( sin ( op * M_PI / 180 ) ) < MIN_FLOAT )
                                {
                                        //throw ErrorMathMatrixSquare ( "ERangError: can not parse equation ", "tan(x), x = ", 1, 1);
                                        throw ErrorMathRange <T> ( "ErrorMathRange: can not parse equation ", "cotg(x), x = ", op );
                                }

                                //Result
                                result = 1 / tan ( op * M_PI / 180 );

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //ASIN(x)
                        else if ( strcmp ( function_text, "asin" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing: ", "asin(x)" );

                                //Pop one operand
                                T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "asin(x), x > MAX.", op );

                                //Throw exception: fabs (abn_acos) - 1 > MIN_FLOAT
                                if ( ( op > 1 + ARGUMENT_ROUND_ERROR ) || ( op < - 1 - ARGUMENT_ROUND_ERROR ) )
                                {
                                        throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: can not parse equation ", "asin(x), x = ", op );
                                }

                                //Correct round errors
                                else if ( op > 1 )
                                {
                                        op = 1;
                                }

                                else if ( op < - 1 )
                                {
                                        op = -1;
                                }

                                //Result
                                result = asin ( op ) * 180 / M_PI;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //ACOS(x)
                        else if ( strcmp ( function_text, "acos" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing: ", "acos(x)" );

                                //Pop one operand
                                T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "acos(x), x > MAX." );

                                //Throw exception: fabs (abn_acos) - 1 > MIN_FLOAT
                                if ( ( op > 1 + ARGUMENT_ROUND_ERROR ) || ( op < - 1 - ARGUMENT_ROUND_ERROR ) )
                                {
                                        throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: can not parse equation", "acos(x), x = ", op );
                                }

                                //Correct round errors
                                else if ( op > 1 )
                                {
                                        op = 1;
                                }

                                else if ( op < - 1 )
                                {
                                        op = -1;
                                }

                                //Result
                                result = acos ( op ) * 180 / M_PI;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //ATAN(x)
                        else if ( strcmp ( function_text, "atan" ) == 0 || strcmp ( function_text, "atg" )  == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing: ", "atan(x)." );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "atan(x), x > MAX." );

                                //Result
                                result = atan ( op ) * 180 / M_PI;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //LN(x)
                        else if ( strcmp ( function_text, "ln" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, ", "ln(x)" );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "ln(x), x > MAX.", op );
                                else if ( op <= MIN_FLOAT ) throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: can not parse equation ", "ln(x), x = ", op );

                                //Result
                                result = log ( op );

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //LOG(x)
                        else if ( strcmp ( function_text, "log" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: ", function_text );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "log(x), x > MAX.", op );
                                else if ( op <= MIN_FLOAT ) throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: can not parse equation ", "log(x), x = ", op );

                                //Result
                                result = log10 ( op );

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //EXP(x)
                        else if ( strcmp ( function_text, "exp" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation. argument missing: ", function_text );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > 3 * MAX_FLOAT_EXPONENT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "exp^x, x > MAX.", op );

                                //Result
                                result = exp ( op );

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //SQR(x)
                        else if ( strcmp ( function_text, "sqr" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing: ", function_text );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > sqrt ( MAX_FLOAT ) ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "sqr(x), x > MAX.", op );

                                //Result
                                result = op * op;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //SQRT(x)
                        else if ( strcmp ( function_text, "sqrt" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing: ", function_text );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "sqr(x), x > MAX.", op );

                                if ( op < 0 ) throw ErrorMathInvalidArgument <T> ( "ErrorMathInvalidArgument: can not parse equation ", "sqr(x), x = ", op );

                                //Result
                                result = sqrt ( op );

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //ABS(x)
                        else if ( strcmp ( function_text, "abs" ) == 0 )
                        {
                                //Empty stack ?
                                if ( operands.empty() ) throw ErrorParse ( "ErrorParse: can not parse equation, argument missing:  ", function_text );

                                //Pop one operand
                                const T op = operands.top();

                                //Remove operand from the stack
                                operands.pop();

                                //Exception
                                if ( op > MAX_FLOAT ) throw ErrorMathOverflow  <T> ( "ErrorMathOverflow: can not parse equation ", "abs(x), x > MAX.", op );

                                //Result
                                result = fabs ( op );

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //VARIABLE x
                        else if ( function_text[0] == 'x' )
                        {
                                //Get result
                                result = x;

                                //Add result to the stack
                                operands.push ( result );
                        }


                        //VARIABLE y
                        else if ( function_text[0] == 'y' )
                        {
                                //Get result
                                result = y;

                                //Add result to the stack
                                operands.push ( result );
                        }


                        //VARIABLE RO
                        else if ( strcmp ( function_text, "RO" )  == 0 || strcmp ( function_text, "Ro" ) == 0 )
                        {
                                //Get result
                                result = 180 / M_PI;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //VARIABLE PI
                        else if ( strcmp ( function_text, "Pi" ) == 0 || strcmp ( function_text, "PI" ) == 0 )
                        {
                                //Get result
                                result = M_PI;

                                //Add result to the stack
                                operands.push ( result );
                        }


                        //VARIABLE R
                        else if ( function_text[0] == 'R' )
                        {
                                //Get result
                                result = R;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //VARIABLE lon
                        else if ( strcmp ( function_text, "lon" ) == 0 || strcmp ( function_text, "lam" )  == 0 || strcmp ( function_text, "v" ) == 0 || strcmp ( function_text, "V" ) == 0 )
                        {
                                //Get result
                                result = lon;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //VARIABLE lat
                        else if ( strcmp ( function_text, "lat" ) == 0 || strcmp ( function_text, "phi" ) == 0 || strcmp ( function_text, "u" ) == 0 || strcmp ( function_text, "U" ) == 0 )
                        {
                                //Get result
                                result = lat;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //VARIABLE lat0
                        else if ( strcmp ( function_text, "lat0" ) == 0 || strcmp ( function_text, "u0" ) == 0 || strcmp ( function_text, "U0" ) == 0 )
                        {
                                //Get result
                                result = lat0;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //VARIABLE lat1
                        else if ( strcmp ( function_text, "lat1" ) == 0 )
                        {
                                //Get result
                                result = lat1;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //VARIABLE lat2
                        else if ( strcmp ( function_text, "lat2" ) == 0 )
                        {
                                //Get result
                                result = lat2;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //VARIABLE lat
                        else if ( strcmp ( function_text, "lon0" ) == 0 )
                        {
                                //Get result
                                result = lon0;

                                //Add result to the stack
                                operands.push ( result );
                        }

                        //Unknown variable
                        else
                        {
                                //Throw exception
                                throw ErrorParse ( "EParsError: ", function_text );
                        }
                }

                //OPERATORS
                else if ( *equation == '^' || *equation == '*' || *equation == '/' || *equation == '+' || *equation == '-' )
                {
                        //Empty stack ?
                        if ( operands.empty() ) throw ErrorParse ( "ErrorParse: ", "Invalid first argument for operation +, -, *, /." );

                        //Pop the first operand
                        const  T op1 = operands.top();

                        //Remove operand from the stack
                        operands.pop();

                        //Empty stack ?
                        if ( operands.empty() ) throw ErrorParse ( "ErrorParse: ", "Invalid second argument for operation +, -, *, /." );

                        //Pop the second operand
                        const T op2 = operands.top();

                        //Remove operand from the stack
                        operands.pop();

                        //Compute result (Power)
                        if ( *equation == '^' )
                        {
                                //Convergate to zero
                                if ( op1 < -MAX_FLOAT_EXPONENT )  result = 0;

                                //Overflows
                                if ( op1 > MAX_FLOAT_EXPONENT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "x^y, number > MAX", op1 );

                                if ( op2 > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "x^y, exponent > MAX.", op2 );

                                result = pow ( op2, op1 );
                        }

                        //Compute result (Multiply)
                        else if ( *equation == '*' )
                        {
                                //Exception
                                if ( op1 > MAX_FLOAT || op2 > MAX_FLOAT )
                                {
                                        throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "x * y, number > MAX", op1 );
                                }

                                result = op2 * op1;
                        }

                        //Compute result (Divide)
                        else if ( *equation == '/' )
                        {
                                //Exception
                                if ( fabs ( op1 ) < MIN_FLOAT )
                                {
                                        throw  ErrorMathZeroDevision <T> ( "ErrorMathDivisonByZero: can not parse equation ", "x / y, y = 0.", op1 );
                                }

                                result = op2 / op1;
                        }

                        //Compute result (Add)
                        else if ( *equation == '+' )
                        {
                                //Exception
                                if ( fabs ( op1 ) > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "x + y, result > MAX", op1 );

                                if ( fabs ( op2 ) > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "x + y, result > MAX", op2 );


                                result = op2 + op1;
                        }

                        //Compute result (Subtract)
                        else if ( *equation == '-' )
                        {
                                //Exception
                                if ( fabs ( op1 ) > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "x - y, number > MAX", op1 );

                                if ( fabs ( op2 ) > MAX_FLOAT ) throw ErrorMathOverflow <T> ( "ErrorMathOverflow: can not parse equation ", "x - y, number > MAX", op2 );

                                result = op2 - op1;
                        }

                        //Add result to the stack
                        operands.push ( result );

                        equation ++;
                }

                //Space, TAB
                else if ( *equation == ' ' || *equation == '\t' )
                {
                        equation ++;
                }

                //Illegal character
                else
                {
                        //Error
                        throw ErrorParse ( "ErrorParse: ", "Illegal character in equation, parsing stopped." );
                }
        }

        //Get result from the stack
        if ( !operands.empty() )
        {
                //Get the result
                result = operands.top();

                //Remove result from the stack
                operands.pop();

                //Is this result correct?
                if ( operands.empty() )
                {
                        //Result is correct
                        return result;
                }

                //Something in the stack
                else
                {
                        //Error in equation
                        char oper[64];
                        sprintf ( oper, "%20.4f", operands.top() );

                        throw ErrorParse ( "ErrorParse: can not parse equation, bad argument: ", oper );
                }
        }

        //Empty stack before result calculation
        else
        {
                throw ErrorParse ( "ErrorParse: can not parse equation, ", " no equation." );
        }
}

#endif

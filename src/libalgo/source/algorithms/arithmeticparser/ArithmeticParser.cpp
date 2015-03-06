// Description: Arithmetic parser using postfix notation based on modified Shunting-yard algorithm

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


#include "ArithmeticParser.h"


//List of variables
const char * vars[] =
{
        "x",
        "y",
        "R",
        "a",
        "b",
        "c",
        "u",
        "v",
        "lat",
        "lon",
        "phi",
        "lam",
        "u0",
        "lat0",
        "phi0",
        "u1",
        "lat1",
        "phi1",
        "u2",
        "lat2",
        "phi2",
        "v0",
        "lon0",
        "lam0",
	"theta"
};


//List of constants
const char * consts [] =
{
        "pi",
        "Pi",
        "PI",
        "e",
        "E",
        "Ro",
        "RO"
};


//List of functions
const char * functs [] =
{
        "sin",
        "cos",
        "tan",
        "tg",
        "cot",
        "cotg",
        "asin",
        "acos",
        "atan",
        "ln",
        "log",
        "exp",
        "sqr",
        "sqrt",
        "abs",
        "sign"
};


void ArithmeticParser::init ( TVarConsFunctMap & vars_list, TVarConsFunctMap & consts_list, TVarConsFunctMap & functs_list )
{
        //Add variables, constants, functions to the map to enable the search

        //Variables
        vars_list[vars[v_x]] = v_x;
        vars_list[vars[v_y]] = v_y;
        vars_list[vars[v_R]] = v_R;
        vars_list[vars[v_a]] = v_a;
        vars_list[vars[v_b]] = v_b;
        vars_list[vars[v_c]] = v_c;
        vars_list[vars[v_u]] = v_u;
        vars_list[vars[v_v]] = v_v;
        vars_list[vars[v_lat]] = v_lat;
        vars_list[vars[v_lon]] = v_lon;
        vars_list[vars[v_phi]] = v_phi;
        vars_list[vars[v_lam]] = v_lam;
        vars_list[vars[v_u0]] = v_u0;
        vars_list[vars[v_lat0]] = v_lat0;
        vars_list[vars[v_phi0]] = v_phi0;
        vars_list[vars[v_u1]] = v_u1;
        vars_list[vars[v_lat1]] = v_lat1;
        vars_list[vars[v_phi1]] = v_phi1;
        vars_list[vars[v_u2]] = v_u2;
        vars_list[vars[v_lat2]] = v_lat2;
        vars_list[vars[v_phi2]] = v_phi2;
        vars_list[vars[v_v0]] = v_v0;
        vars_list[vars[v_lon0]] = v_lon0;
        vars_list[vars[v_lam0]] = v_lam0;
	vars_list[vars[v_theta]] = v_theta;

        //Constants
        consts_list[consts[c_pi]] = c_pi;
        consts_list[consts[c_Pi]] = c_Pi;
        consts_list[consts[c_PI]] = c_PI;
        consts_list[consts[c_e]] = c_e;
        consts_list[consts[c_pi]] = c_E;
        consts_list[consts[c_Ro]] = c_Ro;
        consts_list[consts[c_RO]] = c_RO;

        //Functions
        functs_list[functs[f_sin]] = f_sin;
        functs_list[functs[f_cos]] = f_cos;
        functs_list[functs[f_tan]] = f_tan;
        functs_list[functs[f_tg]] = f_tg;
        functs_list[functs[f_cot]] = f_cot;
        functs_list[functs[f_cotg]] = f_cotg;
        functs_list[functs[f_asin]] = f_asin;
        functs_list[functs[f_acos]] = f_acos;
        functs_list[functs[f_atan]] = f_atan;
        functs_list[functs[f_ln]] = f_ln;
        functs_list[functs[f_log]] = f_log;
        functs_list[functs[f_exp]] = f_exp;
        functs_list[functs[f_sqr]] = f_sqr;
        functs_list[functs[f_sqrt]] = f_sqrt;
        functs_list[functs[f_abs]] = f_abs;
        functs_list[functs[f_sign]] = f_sign;
}


void ArithmeticParser::infixToPostfix ( const char * infix, char * postfix, TPlusMinusOperatorTypes & plus_minus_types_postfix )
{
        //Convert equation from infix notation to postfix notation
        TPlusMinusOperatorType plus_minus_oper_next = UnaryOperator;
        std::stack < std::string> operators;

        //List of variables, constants, functions: initialize
        TVarConsFunctMap vars_list, consts_list, functs_list;
        init ( vars_list, consts_list, functs_list );

        //Stack of binary or unary +- operators
        std::stack <TPlusMinusOperatorType> plus_minus_types_temp;

        //int index = 0;

        //Process all characters of the infix equation
        while ( *infix != '\0' )
        {
                //std::cout << index++ << '\n';

                //Get text
                char text[32];

                //Get Sequence
                findSequence ( &infix, text );

                //Space, tab or last character
                if ( ( *text == '\t' ) || ( *text == ' ' ) || ( *text == '\0' ) )
                {
                        continue;
                }

                //Number
                else if ( isdigit ( ( unsigned char ) *text ) )
                {
                        //Add space
                        *postfix++ = ' ';
                        *postfix = '\0';

                        //Add to the postfix notation
                        strcat ( postfix, text );

                        //Move pointer about length of the text
                        postfix += strlen ( text );

                        //Set next operator +- as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //Variable or const
                else if ( vars_list.find ( std::string ( text ) ) != vars_list.end() || consts_list.find ( std::string ( text ) ) != consts_list.end() )
                {
                        //Add space
                        *postfix++ = ' ';
                        *postfix = '\0';

                        //Add to the postfix notation
                        strcat ( postfix, text );

                        //Move pointer about length of the text
                        postfix += strlen ( text );

                        //Set next operator +- as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //Left bracket
                else if ( *text == '(' )
                {
                        //Add to the stack
                        operators.push ( std::string ( text ) );

                        //Set next operator +- as unary
                        plus_minus_oper_next = UnaryOperator;
                }

                //Right bracket
                else if ( *text == ')' )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                char c[32];

                                //Run until operator with lowest priority found
                                do
                                {
                                        //Error, stop parsing
                                        if ( operators.empty() )
                                        {
                                                //Miss "(" bracket;
                                                throw ErrorParse ( "ErrorParse: can not parse equation,", "missing ( ." );
                                        }

                                        //New operator on the top of the stack
                                        strcpy ( c, operators.top().c_str() );

                                        //Type of the plus minus operator on the top of the stack
                                        TPlusMinusOperatorType plus_minus;
                                        plus_minus = ( !plus_minus_types_temp.empty() ? plus_minus_types_temp.top() : UnaryOperator );

                                        //Add new operator to the postfix notation
                                        if ( *c != '(' ) //Do not add bracket
                                        {
                                                //Add space before function string
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add + - type operator to the output and remove from the stack
                                                if ( ( *c == '+' ) || ( *c == '-' ) )
                                                {
                                                        //Add + - operator type to the output
                                                        plus_minus_types_postfix.push_back ( plus_minus );

                                                        //Remove + - operator type on the top of the stack
                                                        plus_minus_types_temp.pop();
                                                }

                                                //Add new operator to the postfix notation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );
                                        }

                                        //Remove operator on the top of the stack
                                        operators.pop();

                                }
                                while ( *c != '(' ) ;
                        }

                        //Error
                        else
                        {
                                //Miss first "(" bracket;
                                throw ErrorParse ( "ErrorParse: can not parse equation, ", " missing ( ." );
                        }

                        //Set next +- operator as binary
                        plus_minus_oper_next = BinaryOperator;
                }


                //***************************************************************************

                //Functions (Priority Level 4)
                else if ( functs_list.find ( std::string ( text ) ) != functs_list.end() )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                //Run until operator with lowest priority found
                                for ( ;; )
                                {
                                        //Operator on the top  of the stack
                                        const char * c = operators.top().c_str();

                                        //Found operator with lower priority, stop
                                        if ( ( *c == '^' ) || ( *c == '*' ) || ( *c == '/' ) || ( *c == '+' ) || ( *c == '-' ) || ( *c == '(' ) )
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add new operator to the postfix notation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }

                        //Add text to the stack
                        operators.push ( text );

                        //Set next +- operator as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //***************************************************************************

                //Power (Priority level 3)
                else if ( *text == '^' )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                //Run until operator with lowest priority found
                                for ( ;; )
                                {
                                        //Operator on the top  of the stack
                                        const char * c = operators.top().c_str();

                                        //Found operator with lower priority, stop
                                        if ( ( *c == '*' ) || ( *c == '/' ) || ( *c == '+' ) || ( *c == '-' ) || ( *c == '(' ) )
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add new operator to the postfix notation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }

                        //No operator with lowest priority found or empty stack
                        operators.push ( std::string ( text ) );

                        //Set next operator + - as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //***************************************************************************

                //Multiple, divide (Priority level 2)
                else if ( ( *text == '*' ) || ( *text == '/' ) )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                //Run until operator with lowest priority found
                                for ( ;; )
                                {
                                        //Operator on the top  of the stack
                                        const char * c = operators.top().c_str();

                                        //Found operator with lower priority, stop
                                        if ( ( *c == '+' ) || ( *c == '-' ) || ( *c == '(' ) )
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add new operator to the postfix notation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }

                        //No operator with lowest priority found or empty stack
                        operators.push ( std::string ( text ) );

                        //Set next +- operator as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //***************************************************************************

                //Add, subtract (Priority level 1)
                else if ( ( *text == '+' ) || ( *text == '-' ) )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                //Run until operator with lowest priority found
                                for ( ;; )
                                {
                                        //Operator on the top  of the stack
                                        const char * c = operators.top().c_str();

                                        //Type of the plus minus operator on the top of the stack
                                        TPlusMinusOperatorType plus_minus;
                                        plus_minus = ( !plus_minus_types_temp.empty() ? plus_minus_types_temp.top() : UnaryOperator );

                                        //Found operator with lower priority
                                        if ( *c == '(' )
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add new operator to the postfix notationnotation
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Add plus minus type to the output and remove from the stack
                                                if ( ( *c == '+' ) || ( *c == '-' ) )
                                                {
                                                        //Remove plus minus type of the operator from the top of the stack
                                                        plus_minus_types_temp.pop();

                                                        //Add plus minus the type to the output
                                                        plus_minus_types_postfix.push_back ( plus_minus );
                                                }

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }

                        //No operator with lowest priority found or stack empty
                        operators.push ( std::string ( text ) );

                        //Set the next + - operator as unary or binary
                        plus_minus_types_temp.push ( plus_minus_oper_next );

                        //Set next +- operator as binary
                        plus_minus_oper_next = BinaryOperator;
                }

                //Throw exception
                else
                {
                        //Unknown variable or function
                        throw ErrorParse ( "ErrorParse: can not parse equation, unknown function or variable. ", text );
                }
        }

        //Add rest of the stack to the postfix
        while ( ! operators.empty() )
        {
                //Get operator on the top
                char c[32];

                //New operator on the top of the stack
                strcpy ( c, operators.top().c_str() );

                //Type of the plus minus operator on the top of the stack
                TPlusMinusOperatorType plus_minus;
                plus_minus = ( !plus_minus_types_temp.empty() ? plus_minus_types_temp.top() : UnaryOperator );

                //Add new operator to the postfix notation
                if ( *c != '(' )
                {
                        //Add space before function string
                        *postfix++ = ' ';
                        *postfix = '\0';

                        //Add new operator to the postfix notation
                        strcat ( postfix, c );
                        postfix += strlen ( c );

                        //Add plus minus type to the output and remove from the stack
                        if ( ( *c == '+' ) || ( *c == '-' ) )
                        {
                                //Remove plus minus type of the operator from the top of the stack
                                plus_minus_types_postfix.push_back ( plus_minus );

                                //Remove type of the +- operator on the top of the stack
                                plus_minus_types_temp.pop();
                        }
                }

                //Error, stop
                else
                {
                        //Miss ")" bracket
                        throw ErrorParse ( "ErrorParse: can not parse equation, ", "misssing ) ." );
                        break;
                }

                //Remove operator on the top of the stack
                operators.pop();

        }

        //End postfix notation
        *postfix = '\0';
}


void ArithmeticParser::findSequence ( const char ** equation, char * operator_text )
{
        //Create first possible valid sequence of characters in the postfix notation
        //	0 = space .................. A
        //	1 = text ................... B
        //	2 = decimal separator ...... A
        //	3 = arithmetic operator .... A
        do
        {
                //Actual char and increment equation char
                const int c = * ( * equation ) ++;

                //Compare actual char end next char
                //Alphanumeric and non-alphanumeric chars or vice versa
                //Detect changes: [0->1], [1->0], [1->2], [1->3], [2->1], [3->1]
                if ( ( bool ) isalnum ( ( unsigned char ) ** equation ) != ( bool ) isalnum ( ( unsigned char ) c ) )
                {
                        //Detect changes: [0->1], [1->2], [1->3], [2->1], [3->1]
                        if ( !isspace ( c ) )
                        {
                                *operator_text++ = c;

                                //Detect changes: [0->1], [1->3], [3->1], stops for decimal separator
                                if ( ( c != '.' ) && ( c != ',' ) && ( **equation != '.' ) && ( **equation != ',' ) )
                                {
                                        break;
                                }
                        }
                }

                //Both alphanumeric chars or both non-alphanumeric chars
                //Detect changes: [0->0], [0->2x], [0->3], [1->1], [2->0], [2->2x], [2->3x], [3->0], [3->2x], [3->3x]
                else
                {
                        // Detect changes [0->3], [2->3x] [3,3x]
                        if ( ( c == '+' ) || ( c == '-' ) || ( c == '*' ) || ( c == '/' ) || ( c == '^' ) || ( c == '(' ) || ( c == ')' ) )
                        {
                                *operator_text++ = c;
                                break;
                        }

                        //Detect changes [0->2x], [1->1], [2->2x], [3->0], [3->2x]
                        else if ( !isspace ( c ) )
                        {
                                *operator_text++ = c;
                        }
                }

        }
        while ( **equation != '\0' );

        //Correctly end char with \0
        *operator_text = '\0';
}




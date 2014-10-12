#include "Translate.h"

void Translate::infixToPostfix ( const char * infix, char * postfix )
{
        //Convert equation from infix notation to postfix notation
        std::stack < std::string > operators;

        //Process all characters of the infix equation
        while ( *infix != '\0' )
        {
                //Get operand
                char operand[32];

                //Get Sequence
                findSequence ( infix, operand );

                //Move pointer
                infix += strlen ( operand );

                //Number
                if ( isdigit ( ( unsigned char ) *operand ) )
                {
                        //Add space
                        *postfix++ = ' ';
                        *postfix = '\0';

                        //Add to the postfix notation
                        strcat ( postfix, operand );

                        //Move pointer about length of the operand
                        postfix += strlen ( operand );
                }

                //Variables x, y, R, lon, lat, lat0, lat1, lat2, lon0, phi, lam
                else if ( *operand == 'x' || *operand == 'y' || *operand == 'R' || *operand == 'u' || *operand == 'v' || *operand == 'U' || *operand == 'V'
                                || strcmp ( operand, "lon" ) == 0 || strcmp ( operand, "lat" ) == 0 || strcmp ( operand, "phi" ) == 0 || strcmp ( operand, "la" ) == 0
                                || strcmp ( operand, "lam" ) == 0 || strcmp ( operand, "lat0" ) == 0 || strcmp ( operand, "u0" ) == 0 || strcmp ( operand, "U0" ) == 0
                                || strcmp ( operand, "lat1" ) == 0 || strcmp ( operand, "lat2" ) == 0 || strcmp ( operand, "lon0" ) == 0 || strcmp ( operand, "Pi" ) == 0 || strcmp ( operand, "PI" ) == 0
                                || strcmp ( operand, "Ro" ) == 0 || strcmp ( operand, "RO" ) == 0 )
                {
                        //Add space
                        *postfix++ = ' ';
                        *postfix = '\0';

                        //Add to the postfix notation
                        strcat ( postfix, operand );

                        //Move pointer about length of the operand
                        postfix += strlen ( operand );
                }

                //Left bracket
                else if ( *operand == '(' )
                {
                        //Add to the stack
                        operators.push ( ( std::string ( operand ) ) );
                }

                //Right bracket
                else if ( *operand == ')' )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                char c[32];

                                //Run until operator with lowest priority found
                                do
                                {
                                        //Error, stop
                                        if ( operators.empty() )
                                        {
                                                //Miss "(" bracket;
                                                throw ErrorParse ( "ErrorParse: can not parse equation,", "missing ( ." );
                                        }

                                        //New operator on the top of the stack
                                        strcpy ( c, operators.top().c_str() );

                                        //Add operator to the postfix
                                        if ( *c != '(' ) //Do not add bracket
                                        {
                                                //Add space before function string
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add new operator
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

                }


                //***************************************************************************

                //Functions (Priority Level 4)
                else if ( strcmp ( operand , "sin" ) == 0 || strcmp ( operand , "cos" ) == 0 || strcmp ( operand , "tan" ) == 0 || strcmp ( operand , "tg" ) == 0
                                || strcmp ( operand , "cotg" ) == 0 || strcmp ( operand , "cot" ) == 0 || strcmp ( operand , "ln" ) == 0 || strcmp ( operand , "log" ) == 0
                                || strcmp ( operand , "exp" ) == 0 || strcmp ( operand , "sqrt" ) == 0 || strcmp ( operand , "abs" ) == 0 )
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
                                        if ( *c == '^' || *c == '*' || *c == '/' || *c == '+' || *c == '-' || *c == '(' )
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add operator to the postfix
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }

                        //Add operand to the stack
                        operators.push ( operand );
                }

                //***************************************************************************

                //Power (Priority level 3)
                else if ( *operand == '^' )
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
                                        if ( *c == '*' || *c == '/' || *c == '+' || *c == '-' || *c == '(' )
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add operator to the postfix
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
                        operators.push ( std::string ( operand ) );
                }

                //***************************************************************************

                //Multiple, divide (Priority level 2)
                else if ( *operand == '*' || *operand == '/' )
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
                                        if ( *c == '+' || *c == '-' || *c == '(' )
                                        {
                                                break;
                                        }

                                        //Found operator with higher or same priority
                                        else
                                        {
                                                //Add space
                                                *postfix++ = ' ';
                                                *postfix = '\0';

                                                //Add operator to the postfix
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
                        operators.push ( std::string ( operand ) );
                }

                //***************************************************************************

                //Add, subtract (Priority level 1)

                else if ( *operand == '+' || *operand == '-' )
                {
                        //If the stack not empty
                        if ( !operators.empty() )
                        {
                                //Run until operator with lowest priority found
                                for ( ;; )
                                {
                                        //Operator on the top  of the stack
                                        const char * c = operators.top().c_str();

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

                                                //Add operator to the postfix
                                                strcat ( postfix, c );
                                                postfix += strlen ( c );

                                                //Remove operator from the top of the stack
                                                operators.pop();
                                        }

                                        //If stack empty, stop
                                        if ( operators.empty() ) break;
                                }
                        }

                        //No operator with lowest priority found or stack empty
                        operators.push ( std::string ( operand ) );
                }
        }

        //Add rest of the stack to the postfix
        while ( ! operators.empty() )
        {
                //Get operator on the top
                char c[32];

                //New operator on the top of the stack
                strcpy ( c, operators.top().c_str() );

                //Add operator to the postfix
                if ( *c != '(' )
                {
                        //Add space before function string
                        *postfix++ = ' ';
                        *postfix = '\0';

                        //Add new operator
                        strcat ( postfix, c );
                        postfix += strlen ( c );
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


void Translate::findSequence ( const char * equation, char * operator_text )
{
        //Create valid sequence of characters (due to the postfix notation)

        do //Create first possible valid sequence
        {
                //Actual char
                const int c = *equation;

                //Increment pointer
                equation++;

                //Type of the alphanumeric character has been changed (text <-> operator)
                if ( ( bool ) isalnum ( ( unsigned char ) *equation ) != ( bool ) isalnum ( ( unsigned char ) c ) )
                {
                        // space, TAB <-> text and space <-> number represent an illegal changes
                        if ( ( c != ' ' ) && ( c != '\t' ) )
                        {
                                *operator_text++ = c;

                                //For decimal number do not stop number parsing
                                if ( ( c != '.' ) && ( c != ',' ) && ( *equation != '.' ) && ( *equation != ',' ) )
                                {
                                        break;
                                }
                        }
                }

                //Type of the character has not been changed or operator <-> operator
                else
                {
                        // Operator <-> operator
                        if ( c == '+' || c == '-' || c == '*' || c == '/' || c == '^' || c == '('  || c == ')' )
                        {
                                *operator_text++ = c;
                                break;
                        }

                        //Not a space nor TAB, alphanumeric character (text <-> text, number <-> number), has not been changed
                        else if ( ( c != ' ' ) && ( c != '\t' ) )
                        {
                                *operator_text++ = c;
                        }
                }

        }
        while ( *equation != '\0' );

        //End char
        *operator_text = '\0';
}


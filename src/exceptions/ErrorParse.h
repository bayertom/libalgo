// Description: Parse equation error class

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


#ifndef ErrorParse_H
#define ErrorParse_H


#include "Error.h"


//Parser error: infix to postfix notation
class ErrorParse : public Error
{
        protected:
                char * function_text;

        public:

                ErrorParse ( const char * exception_text_, const char * function_text_ ) : Error ( exception_text_ )
                {
                        if ( function_text_ != NULL )
                        {
                                function_text = new char [ strlen ( function_text_ ) + 1 ];
                                strcpy ( function_text, function_text_ );
                        }

                        else
                        {
                                function_text = NULL;
                        }
                }


                ErrorParse ( const ErrorParse & e ) : Error ( e.exception_text )
                {
                        if ( e.function_text != NULL )
                        {
                                function_text = new char [ strlen ( e.function_text ) + 1 ];
                                strcpy ( function_text, e.function_text );
                        }

                        else
                        {
                                function_text = NULL;
                        }
                }


                virtual ~ErrorParse() throw()
                {
                        if ( function_text != NULL )
                        {
                                delete [] function_text;
                                function_text = NULL;
                        }
                }

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Error::printException ( output );
                        *output << function_text << '\n';
                }

                virtual short getErrorCode() const { return 20;}
};

#endif

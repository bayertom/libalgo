// Description: Exception class, other error classes derived from this class

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


#ifndef Error_H
#define Error_H


#include <exception>
#include <ostream>
#include <iostream>
#include <string.h>


//Error, inherited from std::exception
class Error : public std::exception
{
        protected:
                char * exception_text;

        public:
                Error() : exception_text ( NULL ) {}

                Error ( const char * exception_text_ )
                {
                        if ( exception_text_ != NULL )
                        {
                                exception_text = new char [ strlen ( exception_text_ ) + 1 ];
                                strcpy ( exception_text, exception_text_ );
                        }

                        else
                        {
                                exception_text = NULL;
                        }
                }


                Error ( const Error & e )
                {
                        if ( e.exception_text != NULL )
                        {
                                exception_text = new char [ strlen ( e.exception_text ) + 1 ];
                                strcpy ( exception_text, e.exception_text );
                        }

                        else
                        {
                                exception_text = NULL;
                        }
                }


                virtual ~Error () throw()
                {
                        if ( exception_text != NULL )
                        {
                                delete [] exception_text;
                                exception_text = NULL;
                        }
                }

        public:
                const char * const getExceptionText() const {return exception_text;}

        public:
                virtual void printException ( std::ostream * output = &std::cout ) const { *output << exception_text << '\n'; }
                virtual short getErrorCode() const { return 1;}
};


#endif


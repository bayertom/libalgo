// Description: Math error class, other math error classes derived from this class

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


#ifndef ErrorMath_H
#define ErrorMath_H


#include "Error.h"


//Math error
template <typename T>
class ErrorMath : public Error
{
        protected:
                char * math_text;

        public:

                ErrorMath ( const char * exception_text_, const char * math_text_ ) : Error ( exception_text_ )
                {
                        if ( math_text_ != NULL )
                        {
                                math_text = new char [ strlen ( math_text_ ) + 1 ];
                                strcpy ( math_text, math_text_ );
                        }

                        else
                        {
                                math_text = NULL;
                        }
                }


                ErrorMath ( const ErrorMath & e ) : Error ( e.exception_text )
                {
                        if ( e.math_text != NULL )
                        {
                                math_text = new char [ strlen ( e.math_text ) + 1 ];
                                strcpy ( math_text, e.math_text );
                        }

                        else
                        {
                                math_text = NULL;
                        }
                }


                virtual ~ErrorMath() throw()
                {
                        if ( math_text != NULL )
                        {
                                delete [] math_text;
                                math_text = NULL;
                        }
                }

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Error::printException ( output );
                        *output << math_text << '\n';
                }

                virtual T getArg() const = 0;

                virtual short getErrorCode() const { return 9;}
};


#endif

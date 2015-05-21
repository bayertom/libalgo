// Description: Error bad casting class

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


#ifndef ErrorBadCast_H
#define ErrorBadCast_H


#include "Error.h"

//Bad cast error
class ErrorBadCast : public Error
{
        protected:
                char * type_text;

        public:
                ErrorBadCast ( const char * exception_text_, const char * type_text_ ) : Error ( exception_text_ )
                {
                        if ( type_text_ != NULL )
                        {
                                type_text = new char [ strlen ( type_text_ ) + 1 ];
                                strcpy ( type_text, type_text_ );
                        }

                        else
                        {
                                type_text = NULL;
                        }
                }


                ErrorBadCast ( const ErrorBadCast & e ) : Error ( e.exception_text )
                {
                        if ( e.type_text != NULL )
                        {
                                type_text = new char [ strlen ( e.type_text ) + 1 ];
                                strcpy ( type_text, e.type_text );
                        }

                        else
                        {
                                type_text = NULL;
                        }
                }


                virtual ~ErrorBadCast() throw()
                {
                        if ( type_text != NULL )
                        {
                                delete [] type_text;
                                type_text = NULL;
                        }
                }

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Error::printException ( output );
                        *output << type_text << '\n';
                }

                virtual short getErrorCode() const { return 3;}
};

#endif

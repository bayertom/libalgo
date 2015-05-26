// Description: Error bad data class

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


#ifndef ErrorBadData_H
#define ErrorBadData_H


#include "Error.h"


//Bad data error
class ErrorBadData : public Error
{
        protected:
                char * data_text;

        public:

                ErrorBadData ( const char * exception_text_, const char * data_text_ ) : Error ( exception_text_ )
                {
                        if ( data_text_ != NULL )
                        {
                                data_text = new char [ strlen ( data_text_ ) + 1 ];
                                strcpy ( data_text, data_text_ );
                        }

                        else
                        {
                                data_text = NULL;
                        }
                }


                ErrorBadData ( const ErrorBadData & e ) : Error ( e.exception_text )
                {
                        if ( e.data_text != NULL )
                        {
                                data_text = new char [ strlen ( e.data_text ) + 1 ];
                                strcpy ( data_text, e.data_text );
                        }

                        else
                        {
                                data_text = NULL;
                        }
                }


                virtual ~ErrorBadData() throw()
                {
                        if ( data_text != NULL )
                        {
                                delete [] data_text;
                                data_text = NULL;
                        }
                }

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Error::printException ( output );
                        *output << data_text << '\n';
                }

                virtual short getErrorCode() const  { return 4;}
};


#endif

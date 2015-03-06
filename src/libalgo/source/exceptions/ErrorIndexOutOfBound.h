// Description: Error index ouf of bound class

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


#ifndef ErrorIndexOutOfBound_H
#define ErrorIndexOutOfBound_H


#include "Error.h"

//Index out of bound error
class ErrorIndexOutOfBound : public Error
{
        protected:

                char * index_text;

        public:

                ErrorIndexOutOfBound ( const char * exception_text_, const char * index_text_ ) : Error ( exception_text_ )
                {
                        if ( index_text_ != NULL )
                        {
                                index_text = new char [ strlen ( index_text_ ) + 1 ];
                                strcpy ( index_text, index_text_ );
                        }

                        else
                        {
                                index_text = NULL;
                        }
                }


                ErrorIndexOutOfBound ( const ErrorIndexOutOfBound & e ) : Error ( e.exception_text )
                {
                        if ( e.index_text != NULL )
                        {
                                index_text = new char [ strlen ( e.index_text ) + 1 ];
                                strcpy ( index_text, e.index_text );
                        }

                        else
                        {
                                index_text = NULL;
                        }
                }


                virtual ~ErrorIndexOutOfBound() throw()
                {
                        if ( index_text != NULL )
                        {
                                delete [] index_text;
                                index_text = NULL;
                        }
                }

        public:


                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Error::printException ( output );
                        *output << index_text << '\n';
                }

                virtual short getErrorCode() const { return 8;}
};



#endif

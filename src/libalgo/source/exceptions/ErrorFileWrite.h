// Description: Error file write class

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


#ifndef ErrorFileWrite_H
#define ErrorFileWrite_H


#include "Error.h"

//File write error
class ErrorFileWrite : public Error
{
        protected:
                char * file_text;

        public:

                ErrorFileWrite ( const char * exception_text_, const char * file_text_ ) : Error ( exception_text_ )
                {
                        if ( file_text_ != NULL )
                        {
                                file_text = new char [ strlen ( file_text_ ) + 1 ];
                                strcpy ( file_text, file_text_ );
                        }

                        else
                        {
                                file_text = NULL;
                        }
                }


                ErrorFileWrite ( const ErrorFileWrite & e ) : Error ( e.exception_text )
                {
                        if ( e.file_text != NULL )
                        {
                                file_text = new char [ strlen ( e.file_text ) + 1 ];
                                strcpy ( file_text, e.file_text );
                        }

                        else
                        {
                                file_text = NULL;
                        }
                }


                virtual ~ErrorFileWrite() throw()
                {
                        if ( file_text != NULL )
                        {
                                delete [] file_text;
                                file_text = NULL;
                        }
                }

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Error::printException ( output );
                        *output << file_text << '\n';
                }

                virtual short getErrorCode() const { return 7;}
};

#endif

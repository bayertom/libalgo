// Description: Error bad Output class

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


#ifndef ErrorBadOutput_H
#define ErrorBadOutput_H


#include "Error.h"


//Bad Output error
class ErrorBadOutput : public Error
{
        protected:
                char * output_text;

        public:

                ErrorBadOutput ( const char * exception_text_, const char * output_text_ ) : Error ( exception_text_ )
                {
                        if ( output_text_ != NULL )
                        {
                                output_text = new char [ strlen ( output_text_ ) + 1 ];
                                strcpy ( output_text, output_text_ );
                        }

                        else
                        {
                                output_text = NULL;
                        }
                }


                ErrorBadOutput ( const ErrorBadOutput & e ) : Error ( e.exception_text )
                {
                        if ( e.output_text != NULL )
                        {
                                output_text = new char [ strlen ( e.output_text ) + 1 ];
                                strcpy ( output_text, e.output_text );
                        }

                        else
                        {
                                output_text = NULL;
                        }
                }


                virtual ~ErrorBadOutput() throw()
                {
                        if ( output_text != NULL )
                        {
                                delete [] output_text;
                                output_text = NULL;
                        }
                }

        public:

                virtual void printException ( std::ostream * output = &std::cout ) const
                {
                        Error::printException();
                        *output << output_text << '\n';
                }

                virtual short getErrorCode() const  { return 5;}
};


#endif

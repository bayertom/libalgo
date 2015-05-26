// Description: Load text file, items will be denominated into cols

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


#ifndef File_H
#define File_H

#include <vector>
#include <list>
#include <string>

//New user types
typedef std::vector < std::string >  TFileLines;
typedef std::vector < std::vector < std::string > > TFileWords;

//Input operation with file
class File
{
        public:
                static TFileLines loadFileToLines ( const char * file_name );
                static TFileWords loadFileToWords ( const char * file_name );

        private:

                static void commaToDot ( char ** text );
};

#endif

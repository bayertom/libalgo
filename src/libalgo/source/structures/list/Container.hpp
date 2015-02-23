// Description: Class representing a container of items (container of points, edges, projections)

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


#ifndef Container_HPP
#define Container_HPP

#include <algorithm>

#include "libalgo/source/structures/point/Point3DGeographic.h"
#include "libalgo/source/structures/point/Node3DCartesian.h"

#include "libalgo/source/structures/line/Edge.h"

#include "libalgo/source/structures/projection/ProjectionAzimuthal.h"
#include "libalgo/source/structures/projection/ProjectionConic.h"
#include "libalgo/source/structures/projection/ProjectionCylindrical.h"
#include "libalgo/source/structures/projection/ProjectionEllipsoidal.h"
#include "libalgo/source/structures/projection/ProjectionMiscellaneous.h"
#include "libalgo/source/structures/projection/ProjectionPolyConic.h"
#include "libalgo/source/structures/projection/ProjectionPseudoAzimuthal.h"
#include "libalgo/source/structures/projection/ProjectionPseudoConic.h"
#include "libalgo/source/structures/projection/ProjectionPseudoCylindrical.h"

#include "IndexLists.h"

#include "libalgo/source/algorithms/round/Round.h"

#include "libalgo/source/io/File.h"

#include "libalgo/source/comparators/sortPointsByID.h"

#include "libalgo/source/exceptions/Error.h"
#include "libalgo/source/exceptions/ErrorBadData.h"


template <typename Point, const TDestructable destructable>
template <typename CompSort, typename CompEqual>
void Container <Point, destructable>::removeDuplicateElements ( typename TItemsList <Point>::Type ::iterator it_begin, typename TItemsList <Point>::Type ::iterator it_end, CompSort comp_sort, CompEqual comp_equal )
{
        //Remove duplicate points (Generic container)
        std::sort ( it_begin, it_end, comp_sort );
        typename TItemsList <Point>::Type ::iterator i_new_end = std::unique ( it_begin, it_end, comp_equal );

        //Erase duplicate items
        this->items.erase ( i_new_end, this->items.end() );

        //Sort points by ID
        std::sort ( it_begin, this->items.end(), sortPointsByID <Point> () );

        //Renumber all sorted points
        GenericContainer2 <Point, destructable>::updateIDOfItems();
}


template <typename Point, const TDestructable destructable>
template <TDimension dim>
void Container <Point, destructable> ::load ( const char * file,  const bool print_exception, std::ostream * output )
{
        //Perform loading of points (Generic container)
        loadPoints ( file, Dimension <dim>(),  print_exception, output );
}


template <typename Point, const TDestructable destructable>
void Container <Point, destructable>::toIndexList ( TIndexList & il )
{
        //Export indices of points to index list
        GenericContainer2 <Point, destructable>::updateIDOfItems();

        //Find min point id
        unsigned int point_id_min = ( *std::min_element ( this->items.begin(), this->items.end(),
                                      sortPointsByID <Point> ( ) ) ).getPointID();

        //Subtract point_id_min from each point id
        for ( unsigned int i = 0; i < this->items.size(); i++ )
        {
                il.push_back ( this->items[i].getPointID() - point_id_min );
        }
}


template <typename Point, const TDestructable destructable>
template <TDimension dim>
void Container <Point, destructable> ::loadPoints ( const char * file, Dimension <dim>, const bool print_exception, std::ostream * output )
{
        //Load points from file: load points of dimension 3 (Generic container)
        try
        {
                //Load file
                const TFileWords file_content = File::loadFileToWords ( file );

                //Process file
                for ( unsigned int i = 0; i < file_content.size(); i++ )
                {
                        // 3D Geographic point in DD MM SS mode: <point_label> <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        if ( file_content[i].size() == 8 )
                        {
                                this->items.push_back ( Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ) + atof ( file_content[i][2].c_str() ) / 60.0  + atof ( file_content[i][3].c_str() ) / 3600.0,
                                                                atof ( file_content[i][4].c_str() ) + atof ( file_content[i][5].c_str() ) / 60.0  + atof ( file_content[i][6].c_str() ) / 3600.0, atof ( file_content[i][7].c_str() ) ) );
                        }

                        //3D Geographic point in DD MM SS mode: <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        else if ( file_content[i].size() == 7 )
                        {
                                this->items.push_back ( Point ( atof ( file_content[i][0].c_str() ) + atof ( file_content[i][1].c_str() ) / 60.0  + atof ( file_content[i][2].c_str() ) / 3600.0,
                                                                atof ( file_content[i][3].c_str() ) + atof ( file_content[i][4].c_str() ) / 60.0  + atof ( file_content[i][5].c_str() ) / 3600.0, atof ( file_content[i][6].c_str() ) ) );
                        }

                        //Common 3D point with label
                        else if ( file_content[i].size() == 4 )
                        {
                                this->items.push_back ( Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ), atof ( file_content[i][2].c_str() ), atof ( file_content[i][3].c_str() ) ) );
                        }

                        //Common 3D point without label
                        else if ( file_content[i].size() == 3 )
                        {
                                this->items.push_back ( Point ( atof ( file_content[i][0].c_str() ), atof ( file_content[i][1].c_str() ), atof ( file_content[i][2].c_str() ) ) );
                        }

                        //Throw exception
                        else throw ErrorBadData ( "ErrorBadData: unknown data format. ", file );
                }
        }

        //Some error during points processing has appeared
        catch ( Error & error )
        {
                //Print error and do not add projection to the list
                if ( print_exception )
                {
                        error.printException();
                }

                //Clear nodes
                this->items.clear();

                //Throw exception
                throw;
        }
}


template <typename Point, const TDestructable destructable>
void Container <Point, destructable> ::loadPoints ( const char * file, Dimension <Dim2D>,  const bool print_exception, std::ostream * output )
{
        //Load points from file: overloaded function for points of dimension 2 (Generic container)
        try
        {
                //Load file
                const TFileWords file_content = File::loadFileToWords ( file );

                //Process file
                for ( unsigned int i = 0; i < file_content.size(); i++ )
                {
                        // 2D Geographic point in DD MM SS mode: <point_label> <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        if ( file_content[i].size() == 7 )
                        {
                                this->items.push_back ( Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ) + atof ( file_content[i][2].c_str() ) / 60.0  + atof ( file_content[i][3].c_str() ) / 3600.0,
                                                                atof ( file_content[i][4].c_str() ) + atof ( file_content[i][5].c_str() ) / 60.0  + atof ( file_content[i][6].c_str() ) / 3600.0 ) );
                        }

                        //2D Geographic point in DD MM SS mode: <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        else if ( file_content[i].size() == 6 )
                        {
                                this->items.push_back ( Point ( atof ( file_content[i][0].c_str() ) + atof ( file_content[i][1].c_str() ) / 60.0  + atof ( file_content[i][2].c_str() ) / 3600.0,
                                                                atof ( file_content[i][3].c_str() ) + atof ( file_content[i][4].c_str() ) / 60.0  + atof ( file_content[i][5].c_str() ) / 3600.0 ) );
                        }

                        //Common 2D point with label
                        else if ( file_content[i].size() == 3 )
                        {
                                this->items.push_back ( Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ), atof ( file_content[i][2].c_str() ) ) );
                        }

                        //Common 2D point without label
                        else if ( file_content[i].size() == 2 )
                        {
                                this->items.push_back ( Point ( atof ( file_content[i][0].c_str() ), atof ( file_content[i][1].c_str() ) ) );
                        }

                        //Throw exception
                        else throw ErrorBadData ( "ErrorBadData: unknown data format. ", file );
                }
        }

        //Some error during points processing has appeared
        catch ( Error & error )
        {
                //Print error and do not add projection to the list
                if ( print_exception )
                {
                        error.printException();
                }

                //Clear nodes
                this->items.clear();

                //Throw exception
                throw;
        }
}


//*********************************************************************************************************************
//Partial specialization for Point * (List of points)
//*********************************************************************************************************************


template <typename Point, const TDestructable destructable>
Container <Point*, destructable > & Container <Point*, destructable >:: operator = ( const Container <Point*, destructable > &source )
{
        //Operator = (Partial specialization for Point *)
        GenericContainer2 <Point*, destructable >::operator = ( source );

        return *this;
}


template <typename Point, const TDestructable destructable>
template <const TDestructable destructable2>
Container <Point*, destructable > & Container <Point*, destructable >:: operator = ( const Container <Point*, destructable2 > &source )
{
        //Operator = (Partial specialization for Point *)
        GenericContainer2 <Point*, destructable>::operator = ( source );

        return *this;
}


template <typename Point, const TDestructable destructable>
template <TDimension dim>
void Container <Point *, destructable> ::load ( const char * file,  const bool print_exception, std::ostream * output )
{
        //Perform loading of points (Partial specialization for Point *)
        loadPoints ( file, Dimension <dim>(),  print_exception, output );
}


template <typename Point, const TDestructable destructable>
template <TDimension dim>
void Container <Point *, destructable> ::loadPoints ( const char * file, Dimension <dim>,  const bool print_exception, std::ostream * output )
{
        //Load points from file: load points of dimension 3 (Partial specialization for Point *)
        try
        {
                //Load file
                const TFileWords file_content = File::loadFileToWords ( file );

                //Process file
                for ( unsigned int i = 0; i < file_content.size(); i++ )
                {
                        // 3D Geographic point in DD MM SS mode: <point_label> <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        if ( file_content[i].size() == 8 )
                        {
                                this->items.push_back ( new Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ) + atof ( file_content[i][2].c_str() ) / 60.0  + atof ( file_content[i][3].c_str() ) / 3600.0,
                                                                    atof ( file_content[i][4].c_str() ) + atof ( file_content[i][5].c_str() ) / 60.0  + atof ( file_content[i][6].c_str() ) / 3600.0, atof ( file_content[i][7].c_str() ) ) );
                        }

                        //3D Geographic point in DD MM SS mode: <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        else if ( file_content[i].size() == 7 )
                        {
                                this->items.push_back ( new Point ( atof ( file_content[i][0].c_str() ) + atof ( file_content[i][1].c_str() ) / 60.0  + atof ( file_content[i][2].c_str() ) / 3600.0,
                                                                    atof ( file_content[i][3].c_str() ) + atof ( file_content[i][4].c_str() ) / 60.0  + atof ( file_content[i][5].c_str() ) / 3600.0, atof ( file_content[i][6].c_str() ) ) );
                        }

                        //Common 3D point with label
                        else if ( file_content[i].size() == 4 )
                        {
                                this->items.push_back ( new Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ), atof ( file_content[i][2].c_str() ), atof ( file_content[i][3].c_str() ) ) );
                        }

                        //Common 3D point without label
                        else if ( file_content[i].size() == 3 )
                        {
                                this->items.push_back ( new Point ( atof ( file_content[i][0].c_str() ), atof ( file_content[i][1].c_str() ), atof ( file_content[i][2].c_str() ) ) );
                        }

                        //Throw exception
                        else throw ErrorBadData ( "ErrorBadData: unknown data format. ", file );
                }
        }

        //Some error during points processing has appeared
        catch ( Error & error )
        {
                //Print error and do not add projection to the list
                if ( print_exception )
                {
                        error.printException();
                }

                //Clear nodes
                this->items.clear();

                //Throw exception
                throw;
        }
}


template <typename Point, const TDestructable destructable>
void Container <Point *, destructable>::toIndexList ( TIndexList & il )
{
        //Export points to index list
        GenericContainer2 <Point*, destructable>::updateIDOfItems();

        //Get min point ID
        unsigned int point_id_min = ( *std::min_element ( this->items.begin(), this->items.end(),
                                      sortPointsByID <Point> () ) )->getPointID();

        //Subtract point_id_min from each point id
        for ( unsigned int i = 0; i < this->items.size(); i++ )
        {
                il.push_back ( this->items[i]->getPointID() - point_id_min );
        }
}


template <typename Point, const TDestructable destructable>
void Container <Point *, destructable> ::loadPoints ( const char * file, Dimension <Dim2D>,  const bool print_exception, std::ostream * output )
{
        //Load points from file: overloaded function for points of dimension 2 (Partial specialization for Point *)
        try
        {
                //Load file
                const TFileWords file_content = File::loadFileToWords ( file );

                //Process file
                for ( unsigned int i = 0; i < file_content.size(); i++ )
                {
                        // 2D Geographic point in DD MM SS mode: <point_label> <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        if ( file_content[i].size() == 7 )
                        {
                                this->items.push_back ( new Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ) + atof ( file_content[i][2].c_str() ) / 60.0  + atof ( file_content[i][3].c_str() ) / 3600.0,
                                                                    atof ( file_content[i][4].c_str() ) + atof ( file_content[i][5].c_str() ) / 60.0  + atof ( file_content[i][6].c_str() ) / 3600.0 ) );
                        }

                        //2D Geographic point in DD MM SS mode: <DD_lat> <MM_lat> <SS_lat> <DD_lon> <MM_lon> <SS_lon> <H>
                        else if ( file_content[i].size() == 6 )
                        {
                                this->items.push_back ( new Point ( atof ( file_content[i][0].c_str() ) + atof ( file_content[i][1].c_str() ) / 60.0  + atof ( file_content[i][2].c_str() ) / 3600.0,
                                                                    atof ( file_content[i][3].c_str() ) + atof ( file_content[i][4].c_str() ) / 60.0  + atof ( file_content[i][5].c_str() ) / 3600.0 ) );
                        }

                        //Common 2D point with label
                        else if ( file_content[i].size() == 3 )
                        {
                                this->items.push_back ( new Point ( file_content[i][0].c_str(), atof ( file_content[i][1].c_str() ), atof ( file_content[i][2].c_str() ) ) );
                        }

                        //Common 2D point without label
                        else if ( file_content[i].size() == 2 )
                        {
                                this->items.push_back ( new Point ( atof ( file_content[i][0].c_str() ), atof ( file_content[i][1].c_str() ) ) );
                        }

                        //Throw exception
                        else throw ErrorBadData ( "ErrorBadData: unknown data format. ", file );
                }
        }

        //Some error during points processing has appeared
        catch ( Error & error )
        {
                //Print error and do not add projection to the list
                if ( print_exception )
                {
                        error.printException();
                }

                //Clear nodes
                this->items.clear();

                //Throw exception
                throw;
        }
}


template <typename Point, const TDestructable destructable>
template <typename CompSort, typename CompEqual>
void Container <Point *, destructable>::removeDuplicateElements ( typename TItemsList <Point *>::Type ::iterator it_begin, typename TItemsList <Point*>::Type ::iterator it_end, CompSort comp_sort, CompEqual comp_equal )
{
        //Remove duplicate points (Partial specialization for Item *)
        if ( ( it_begin != it_end ) && ( it_begin != this->items.end() ) )
        {
                //Sort items
                std::sort ( it_begin, it_end, comp_sort );

                //Remove duplicate items including a correct destruction
                unsigned int index_last_unique = it_begin - this->items.begin(), index_start = index_last_unique;;
                unsigned int end = this->items.size() - ( this->items.end() - it_end );

                for ( unsigned int i = index_last_unique ; i < end ; i++ )
                {
                        //Compare adjacent items
                        if ( ( *this->items[i] != * this->items[index_last_unique] ) )
                        {
                                //Swap duplicate item with actual item (duplicate items will be moved rights to index_last_unique)
                                Point * temp = this->items[++index_last_unique];
                                this->items[index_last_unique] = this->items[i];
                                this->items[i] = temp;
                        }
                }

                //Call destructor for items rights to index_last_unique (all duplicate elements)
                for ( unsigned int i = index_last_unique + 1; i < end ; i++ )
                {
                        if ( this->items[i] != NULL ) delete this->items[i];
                }

                //Resize items to unique index
                this->items.resize ( index_last_unique + 1 );

                //Sort points by ID
                std::sort ( this->items.begin() += index_start, this->items.begin() += std::min ( index_last_unique + 1, end ), sortPointsByID <Point> () );

                //Renumber all points
                GenericContainer2 <Point*, destructable>::updateIDOfItems();
        }
}


//*********************************************************************************************************************
//Partial specialization for Edge (List of Edges)
//*********************************************************************************************************************

template <typename Point, const TDestructable destructable>
void Container <Edge <Point> , destructable>:: load ( const char * file, Container <Point, destructable> *pl, const bool print_exception )
{
        //Load edges from file (Partial specialization)
        try
        {
                //Load file
                const TFileWords file_content = File::loadFileToWords ( file );

                //Process file
                for ( unsigned int i = 0; i < file_content.size(); i++ )
                {
                        //Add node to list (insert only lines with correct number of items)
                        if ( file_content[i].size() == 2 )
                        {
                                //Get indices
                                const unsigned int index1 = atoi ( file_content[i][0].c_str() );
                                const unsigned int index2 = atoi ( file_content[i][1].c_str() );

                                //Add edge to the list
                                this->items.push_back ( Edge <Point> ( ( *pl ) [index1 - 1], ( *pl ) [index2 - 1] ) );
                        }
                }
        }

        //Some error during processing Projection has appeared
        catch ( Error & error )
        {
                //Print error and do not add projection to the list
                if ( print_exception )
                {
                        error.printException();
                }

                //Clear nodes
                this->clear();

                //Throw exception
                throw;
        }
}


//*********************************************************************************************************************
//Partial specialization for Face (List of Faces)
//*********************************************************************************************************************

template <typename T, const TDestructable destructable>
Container <Face<T> *, destructable > & Container <Face<T>*, destructable >:: operator = ( const Container <Face<T>*, destructable > &source )
{
        //Operator = (Partial specialization for Point *)
        GenericContainer2 <Face<T> *, destructable >::operator = ( source );

        return *this;
}



//*********************************************************************************************************************
//Partial specialization for VoronoiCell (List of Voronoi Cells)
//*********************************************************************************************************************

template <typename T, const TDestructable destructable>
Container <VoronoiCell<T> *, destructable > & Container <VoronoiCell<T>*, destructable >:: operator = ( const Container <VoronoiCell<T>*, destructable > &source )
{
        //Operator = (Partial specialization for Point *)
        GenericContainer2 <VoronoiCell<T> *, destructable >::operator = ( source );

        return *this;
}



//*********************************************************************************************************************
//Partial specialization for Projection (List of Projections)
//*********************************************************************************************************************

template <typename T, const TDestructable destructable>
void Container <Projection <T> *, destructable> ::load ( const char * file, const bool print_exception )
{

        //Load all projections from the specified file
        unsigned int non_empty_lines = 0;
        char terminator[BUFF];
        char line_number[BUFF];
        Point3DGeographic <T> pole;

        //Indicators
        bool set_projection_type = false;
        bool set_name = false;
        bool set_x_equation = false;
        bool set_y_equation = false;
        bool set_lat_pole = false;
        bool set_lon_pole = false;
        bool set_lat0 = false;
        bool set_lat1 = false;
        bool set_lat2 = false;
        bool set_lon0 = false;
        bool set_lon_dir = false;
        bool set_dx = false;
        bool set_dy = false;
        bool set_r = false;
        bool set_a = false;
        bool set_b = false;

        terminator[0] = '\0';

        //Load all items from file
        TFileLines file_content = File::loadFileToLines ( file );

        //Process all cartographic projection
        Projection <T> * proj = NULL;

        unsigned int i = 0;

        for ( ; i < file_content.size(); i++ )
        {
                try
                {
                        //Get actual line of the file
                        const char * line_text = file_content[i].c_str();

                        //Process line with terminator
                        if ( strcmp ( line_text, "<projection>" ) == 0 || strcmp ( line_text, "<Projection>" )  == 0 || strcmp ( line_text, "<name>" ) == 0 || strcmp ( line_text, "<Name>" ) == 0 ||
                                        strcmp ( line_text, "<x>" ) == 0 || strcmp ( line_text, "<X>" ) == 0 || strcmp ( line_text, "<y>" ) == 0 || strcmp ( line_text, "<Y>" ) == 0 ||
                                        strcmp ( line_text, "<lat_pole>" ) == 0 || strcmp ( line_text, "<Lat_pole>" ) == 0 || strcmp ( line_text, "<lon_pole>" ) == 0 ||
                                        strcmp ( line_text, "<Lon_pole>" ) == 0 || strcmp ( line_text, "<lat0>" ) == 0 || strcmp ( line_text, "<Lat0>" ) == 0 || strcmp ( line_text, "<lat1>" ) == 0 ||
                                        strcmp ( line_text, "<Lat1>" ) == 0 || strcmp ( line_text, "<lat2>" ) == 0 || strcmp ( line_text, "<Lat2>" ) == 0 || strcmp ( line_text, "<lon0>" ) == 0 ||
                                        strcmp ( line_text, "<Lon0>" )  == 0 || strcmp ( line_text, "<lon_dir>" ) == 0 || strcmp ( line_text, "<Lon_dir>" )  == 0 || strcmp ( line_text, "<dx>" ) == 0 ||
                                        strcmp ( line_text, "<Dx>" ) == 0 || strcmp ( line_text, "<dy>" )  == 0 || strcmp ( line_text, "<Dy>" ) == 0 || strcmp ( line_text, "<r>" ) == 0 ||
                                        strcmp ( line_text, "<R>" ) == 0 || strcmp ( line_text, "<a>" ) == 0 || strcmp ( line_text, "<A>" ) == 0 || strcmp ( line_text, "<b>" ) == 0 || strcmp ( line_text, "<B>" ) == 0 )
                        {

                                //Next non empty line after terminator
                                if ( non_empty_lines == 1 || terminator[0] == '\0' )
                                {
                                        //Copy buffer without white spaces to the char
                                        strcpy ( terminator, line_text );

                                        //Reset non-empty lines count
                                        non_empty_lines = 0;
                                }

                                //More than one terminators
                                else
                                {
                                        //Get line number
                                        sprintf ( line_number, "%d", i + 1 );

                                        //Throw exception
                                        throw ErrorBadData ( "ErrorBadData: error in projection file, following two terminators, line ", line_number );
                                }
                        }

                        //Process line without terminator
                        else if ( non_empty_lines == 0 )
                        {
                                //Test projection type
                                if ( strcmp ( terminator, "<projection>" ) * strcmp ( terminator, "<Projection>" ) == 0 )
                                {
                                        //Definition of new projection starts, but definition of old projection still has not been completed
                                        if ( set_projection_type )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, projection type has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_projection_type = true;

                                                //Create new azimuthal projection and set properties
                                                if ( strcmp ( file_content[i].c_str(), "azimuthal" ) == 0 ||  strcmp ( file_content[i].c_str(), "Azimuthal" ) == 0 )
                                                {
                                                        proj = new ProjectionAzimuthal <T> ();
                                                }

                                                //Create new conic projection and set properties
                                                else if ( strcmp ( file_content[i].c_str(), "conic" ) == 0 || strcmp ( file_content[i].c_str(), "Conic" ) == 0 )
                                                {
                                                        proj = new ProjectionConic <T> ();
                                                }

                                                //Set projection type: cylindrical
                                                else if ( strcmp ( file_content[i].c_str(), "cylindrical" ) == 0 || strcmp ( file_content[i].c_str(), "Cylindrical" ) == 0 )
                                                {
                                                        proj = new ProjectionCylindrical <T> ();
                                                }

                                                //Set projection type: ellipsoidal
                                                else if ( strcmp ( file_content[i].c_str(), "ellipsoidal" ) == 0 || strcmp ( file_content[i].c_str(), "Ellipsoidal" ) == 0 )
                                                {
                                                        proj = new ProjectionEllipsoidal <T> ();
                                                }

                                                //Set projection type: miscellaneous
                                                else if ( strcmp ( file_content[i].c_str(), "miscellaneous" ) == 0 || strcmp ( file_content[i].c_str(), "Miscellaneous" ) == 0 )
                                                {
                                                        proj = new ProjectionMiscellaneous <T> ();
                                                }

                                                //Set projection type: polyconic
                                                else if ( strcmp ( file_content[i].c_str(), "polyconic" ) == 0 || strcmp ( file_content[i].c_str(), "Polyconic" ) == 0 )
                                                {
                                                        proj = new ProjectionPolyconic <T> ();
                                                }

                                                //Set projection type: pseudoazimuthal
                                                else if ( strcmp ( file_content[i].c_str(), "pseudoazimuthal" ) == 0 || strcmp ( file_content[i].c_str(), "Pseudoazimuthal" ) == 0 )
                                                {
                                                        proj = new ProjectionPseudoAzimuthal <T> ();
                                                }

                                                //Set projection type: pseudoconic
                                                else if ( strcmp ( file_content[i].c_str(), "pseudoconic" ) == 0 || strcmp ( file_content[i].c_str(), "Pseudoconic" ) == 0 )
                                                {
                                                        proj = new ProjectionPseudoConic <T> ();
                                                }

                                                //Set projection type: pseudocylindrical
                                                else if ( strcmp ( file_content[i].c_str(), "pseudocylindrical" ) == 0 || strcmp ( file_content[i].c_str(), "Pseudocylindrical" ) == 0 )
                                                {
                                                        proj = new ProjectionPseudoCylindrical <T> ();
                                                }

                                                //Unknown projection, throw exception
                                                else
                                                {
                                                        //Get line number
                                                        sprintf ( line_number, "%d", i + 1 );

                                                        //Throw exception
                                                        throw ErrorBadData ( "ErrorBadData: error in projection file, bad projection type, line ", line_number );
                                                }
                                        }
                                }

                                //Name of the projection
                                else if ( strcmp ( terminator, "<name>" ) * strcmp ( terminator, "<Name>" ) == 0 )
                                {
                                        //Projection name has already been defined
                                        if ( set_name )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, projection name has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_name = true;
                                                proj->setProjectionName ( line_text );
                                        }
                                }

                                //Test X equation
                                else if ( strcmp ( terminator, "<x>" ) * strcmp ( terminator, "<X>" ) == 0 )
                                {
                                        //X equation has already been defined
                                        if ( set_x_equation )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, x equation has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_x_equation = true;
                                                proj->setXEquat ( line_text );
                                        }
                                }

                                //Test Y equation
                                else if ( strcmp ( terminator, "<y>" ) * strcmp ( terminator, "<Y>" ) == 0 )
                                {
                                        //Y equation has already been defined
                                        if ( set_y_equation )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, y equation has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_y_equation = true;
                                                proj->setYEquat ( line_text );
                                        }
                                }

                                //Test Lat pole
                                else if ( strcmp ( terminator, "<lat_pole>" ) * strcmp ( terminator, "<Lat_pole>" ) == 0 )
                                {
                                        //Lat pole has already been defined
                                        if ( set_lat_pole )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, latitude of the pole has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_lat_pole = true;
                                                pole.setLat ( atof ( line_text ) );

                                                proj->setCartPole ( pole );
                                        }
                                }

                                //Test Lon pole
                                else if ( strcmp ( terminator, "<lon_pole>" ) * strcmp ( terminator, "<Lon_pole>" ) == 0 )
                                {
                                        //Lon pole has already been defined
                                        if ( set_lon_pole )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, longitude of the pole has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_lon_pole = true;
                                                pole.setLon ( atof ( line_text ) );

                                                proj->setCartPole ( pole );
                                        }
                                }

                                //Test Lat0
                                else if ( strcmp ( terminator, "<lat0>" ) * strcmp ( terminator, "<Lat0>" ) == 0 )
                                {
                                        //Lat0 has already been defined
                                        if ( set_lat0 )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, lat0 has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_lat0 = true;
                                                proj->setLat0 ( atof ( line_text ) );
                                        }
                                }

                                //Test Lat1
                                else if ( strcmp ( terminator, "<lat1>" ) * strcmp ( terminator, "<Lat1>" ) == 0 )
                                {
                                        //Lat1 has already been defined
                                        if ( set_lat1 )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, lat1 has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_lat1 = true;
                                                proj->setLat1 ( atof ( line_text ) );
                                        }
                                }

                                //Test Lat2
                                else if ( strcmp ( terminator, "<lat2>" ) * strcmp ( terminator, "<Lat2>" ) == 0 )
                                {
                                        //Lat2 has already been defined
                                        if ( set_lat2 )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, lat2 has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_lat2 = true;
                                                proj->setLat2 ( atof ( line_text ) );
                                        }
                                }

                                //Test lon0
                                else if ( strcmp ( terminator, "<lon0>" ) * strcmp ( terminator, "<Lon0>" ) == 0 )
                                {
                                        //Lon0 has already been defined
                                        if ( set_lon0 )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );;

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, lon0 has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_lon0 = true;
                                                proj->setLon0 ( atof ( line_text ) );
                                        }
                                }

                                //Test lon_dir
                                else if ( strcmp ( terminator, "<lon_dir>" ) * strcmp ( terminator, "<Lon_dir>" ) == 0 )
                                {
                                        //Lon_dir has already been defined
                                        if ( set_lon_dir )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );;

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, londir has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_lon_dir = true;
                                                proj->setLonDir ( ( TTransformedLongtitudeDirection ) atoi ( line_text ) );
                                        }
                                }

                                //Test dx
                                else if ( strcmp ( terminator, "<dx>" ) * strcmp ( terminator, "<Dx>" ) == 0 )
                                {
                                        //dx has already been defined
                                        if ( set_dx )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, dx has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_dx = true;
                                                proj->setDx ( atof ( line_text ) );
                                        }
                                }

                                //Test dy
                                else if ( strcmp ( terminator, "<dy>" ) * strcmp ( terminator, "<Dy>" ) == 0 )
                                {
                                        //dy has already been defined
                                        if ( set_dy )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, dy has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_dy = true;
                                                proj->setDy ( atof ( line_text ) );
                                        }
                                }

                                //Test R
                                else if ( strcmp ( terminator, "<r>" ) * strcmp ( terminator, "<R>" ) == 0 )
                                {
                                        //R has already been defined
                                        if ( set_r )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, R has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_r = true;
                                                proj->setR ( atof ( line_text ) );
                                        }
                                }

                                //Test a
                                else if ( strcmp ( terminator, "<a>" ) * strcmp ( terminator, "<A>" ) == 0 )
                                {
                                        //a has already been defined
                                        if ( set_a )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, a has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_a = true;
                                                proj->setA ( atof ( line_text ) );
                                        }
                                }

                                //Test b
                                else if ( strcmp ( terminator, "<b>" ) * strcmp ( terminator, "<B>" ) == 0 )
                                {
                                        //b has already been defined
                                        if ( set_b )
                                        {
                                                //Get line number
                                                sprintf ( line_number, "%d", i + 1 );

                                                //Throw exception
                                                throw ErrorBadData ( "ErrorBadData: error in projection file, b has already been defined, line ", line_number );
                                        }

                                        else
                                        {
                                                set_b = true;
                                                proj->setB ( atof ( line_text ) );
                                        }
                                }

                                //Bad terminator (none of above)
                                else
                                {
                                        //Get line number
                                        sprintf ( line_number, "%d", i + 1 );

                                        //Throw exception
                                        throw ErrorBadData ( "ErrorBadData: error in projection file, bad terminator <>, line ", line_number );
                                }

                                //Increment non empty lines indicator
                                non_empty_lines ++;
                        }

                        //More than 1 non empty line with attributes
                        else
                        {
                                //Get line number
                                sprintf ( line_number, "%d", i + 1 );

                                //Throw exception
                                throw ErrorBadData ( "ErrorBadData: error in projection file, following two attribute lines, line ", line_number );
                        }


                        //All items were set?
                        if ( set_projection_type && set_name && set_x_equation && set_y_equation && set_lat_pole && set_lon_pole && set_lat0 &&
                                        set_lat1 && set_lat2 && set_lon0 && set_lon_dir && set_dx && set_dy && set_r && set_a && set_b )
                        {

                                //Add projection to the list
                                this->items.push_back ( proj );

                                //Reset all indicators
                                set_projection_type = false;
                                set_name = false;
                                set_x_equation = false;
                                set_y_equation = false;
                                set_lat_pole = false;
                                set_lon_pole = false;
                                set_lat0 = false;
                                set_lat1 = false;
                                set_lat2 = false;
                                set_lon0 = false;
                                set_lon_dir = false;
                                set_dx = false;
                                set_dy = false;
                                set_r = false;
                                set_a = false;
                                set_b = false;

                                //Set projection to null
                                proj = NULL;
                        }
                }

                //Delete uncomplete projection
                catch ( Error & error )
                {
                        error.printException();

                        if ( proj != NULL )
                        {
                                delete proj;
                        }

                        throw;
                }
        }
}


#endif

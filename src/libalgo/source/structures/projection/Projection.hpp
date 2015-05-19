// Description: General cartographic projection with pure wirtual methods
// Foe evaluations, its equations in the postfix notation are stored

// Copyright (c) 2010 - 2015
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


#ifndef Projection_HPP
#define Projection_HPP

#include <string.h>
#include <stdio.h>

template <typename T>
class ProjectionAzimuthal;

template <typename T>
class ProjectionConic;

template <typename T>
class ProjectionCylindrical;

template <typename T>
class ProjectionEllipsoidal;

template <typename T>
class ProjectionPolyconic;

template <typename T>
class ProjectionPseudoAzimuthal;

template <typename T>
class ProjectionPseudoConic;

template <typename T>
class ProjectionPseudoCylindrical;

template <typename T>
class ProjectionMiscellaneous;


template <typename T>
Projection <T> ::Projection(const T R_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_, const char * x_equat_postfix_, const char * y_equat_postfix_)
	: R(R_), lon0(lon0_), dx(dx_), dy(dy_), c(c_)
{
	projection_family = NULL;
	projection_name = NULL;

        if ( x_equat_ != NULL )
        {
                x_equat = new char [ strlen ( x_equat_ ) + 1 ];
                strcpy ( x_equat, x_equat_ );
        }

	else
	{
		x_equat = NULL;
	}

	if (y_equat_ != NULL)
	{
		y_equat = new char[strlen(y_equat_) + 1];
		strcpy(y_equat, y_equat_);
	}

        else
        {
               
                y_equat = NULL;
        }

	if (x_equat_postfix_ != NULL)
	{
		x_equat_postfix = new char[strlen(x_equat_postfix_) + 1];
		strcpy(x_equat_postfix, x_equat_postfix_);
	}

	else
	{
		x_equat_postfix = NULL;
	}

	if (y_equat_postfix_ != NULL)
	{
		y_equat_postfix = new char[strlen(y_equat_postfix_) + 1];
		strcpy(y_equat_postfix, y_equat_postfix_);
	}

	else
	{

		y_equat_postfix = NULL;
	}
}


template <typename T>
Projection <T> ::Projection(const T R_, const T lon0_, const T dx_, const T dy_, const T c_, const char * x_equat_, const char * y_equat_, const char * x_equat_postfix_, const char * y_equat_postfix_, const char * projection_family_, const char * projection_name_)
	: R(R_), lon0(lon0_), dx(dx_), dy(dy_), c(c_)
{
	if (x_equat_ != NULL)
	{
		x_equat = new char [ strlen ( x_equat_ ) + 1 ];
		strcpy(x_equat, x_equat_);
	}

	else
	{
		x_equat = NULL;
	}

	if (y_equat_ != NULL)
	{
		y_equat = new char[strlen(y_equat_) + 1];
		strcpy(y_equat, y_equat_);
	}

	else
	{
		y_equat = NULL;
	}

	if (x_equat_postfix_ != NULL)
	{
		x_equat_postfix = new char[strlen(x_equat_postfix_) + 1];
		strcpy(x_equat_postfix, x_equat_postfix_);
	}

	else
	{
		x_equat_postfix = NULL;
	}

	if (y_equat_postfix_ != NULL)
	{
		y_equat_postfix = new char[strlen(y_equat_postfix_) + 1];
		strcpy(y_equat_postfix, y_equat_postfix_);
	}

	else
	{
		y_equat_postfix = NULL;
	}

	if (projection_family_ != NULL)
	{
		projection_family = new char[strlen(projection_family_) + 1];
		strcpy(projection_family, projection_family_);
	}

	else
	{
		projection_family = NULL;
	}

	if (projection_name_ != NULL)
	{
		projection_name = new char[strlen(projection_name_) + 1];
		strcpy(projection_name, projection_name_);
	}

	else
	{
		projection_name = NULL;
	}
}


template <typename T>
Projection <T> ::Projection(const Projection <T> &proj) : R(proj.R), lon0(proj.lon0), dx(proj.dx), dy(proj.dy), c(proj.c)
{

	if (proj.x_equat != NULL)
	{
		x_equat = new char[strlen(proj.x_equat) + 1];
		strcpy(x_equat, proj.x_equat);
	}

	else
	{
		x_equat = NULL;
	}

	if (proj.y_equat != NULL)
	{
		y_equat = new char[strlen(proj.y_equat) + 1];
		strcpy(y_equat, proj.y_equat);
	}

	else
	{
		y_equat = NULL;
	}

	if (proj.x_equat_postfix != NULL)
	{
		x_equat_postfix = new char[strlen(proj.x_equat_postfix) + 1];
		strcpy(x_equat_postfix, proj.x_equat_postfix);
	}

	else
	{
		x_equat_postfix = NULL;
	}

	if (proj.y_equat_postfix != NULL)
	{
		y_equat_postfix = new char[strlen(proj.y_equat_postfix) + 1];
		strcpy(y_equat_postfix, proj.y_equat_postfix);
	}

	else
	{
		y_equat_postfix = NULL;
	}

	if (proj.projection_family != NULL)
	{
		projection_family = new char[strlen(proj.projection_family) + 1];
		strcpy(projection_family, proj.projection_family);
	}

	else
	{
		projection_family = NULL;
	}

	if (proj.projection_name != NULL)
	{
		projection_name = new char[strlen(proj.projection_name) + 1];
		strcpy(projection_name, proj.projection_name);
	}

	else
	{
		projection_name = NULL;
	}

}


template <typename T>
Projection <T> ::Projection(const Projection <T> *proj) : R(proj->R), lon0(proj->lon0), dx(proj->dx), dy(proj->dy), c(proj->c)
{

	if (proj->x_equat != NULL)
	{
		x_equat = new char[strlen(proj->x_equat) + 1];
		strcpy(x_equat, proj->x_equat);
	}

	else
	{
		x_equat = NULL;
	}

	if (proj->y_equat != NULL)
	{
		y_equat = new char[strlen(proj->y_equat) + 1];
		strcpy(y_equat, proj->y_equat);
	}

	else
	{
		y_equat = NULL;
	}

	if (proj->x_equat_postfix != NULL)
	{
		x_equat_postfix = new char[strlen(proj->x_equat_postfix) + 1];
		strcpy(x_equat_postfix, proj->x_equat_postfix);
	}

	else
	{
		x_equat_postfix = NULL;
	}

	if (proj->y_equat_postfix != NULL)
	{
		y_equat_postfix = new char[strlen(proj->y_equat_postfix) + 1];
		strcpy(y_equat_postfix, proj->y_equat_postfix);
	}

	else
	{
		y_equat_postfix = NULL;
	}

	if (proj->projection_family != NULL)
	{
		projection_family = new char[strlen(proj->projection_family) + 1];
		strcpy(projection_family, proj->projection_family);
	}

	else
	{
		projection_family = NULL;
	}

	if (proj->projection_name != NULL)
	{
		projection_name = new char[strlen(proj->projection_name) + 1];
		strcpy(projection_name, proj->projection_name);
	}

	else
	{
		projection_name = NULL;
	}
}


template <typename T>
Projection <T>:: ~Projection()
{
        if ( x_equat != NULL )
        {
                delete [] x_equat;
                x_equat = NULL;
        }

        if ( y_equat != NULL )
        {
                delete [] y_equat;
                y_equat = NULL;
        }

	if (x_equat_postfix != NULL)
	{
		delete[] x_equat_postfix;
		x_equat_postfix = NULL;
	}

	if (y_equat_postfix != NULL)
	{
		delete[] y_equat_postfix;
		y_equat_postfix = NULL;
	}

	if (projection_family != NULL)
	{
		delete[] projection_family;
		projection_family = NULL;
	}

        if ( projection_name != NULL )
        {
                delete [] projection_name;
                projection_name = NULL;
        }
}


template <typename T>
void Projection <T> ::setXEquat ( const char * x_equat_ )
{
        if ( x_equat_ != NULL )
        {
                x_equat = new char [ strlen ( x_equat_ ) + 1 ];
                strcpy ( x_equat, x_equat_ );
        }

        else
        {
                x_equat = NULL;
        }
}


template <typename T>
void Projection <T> :: setYEquat ( const char * y_equat_ )
{
        if ( y_equat_ != NULL )
        {
                y_equat = new char [ strlen ( y_equat_ ) + 1 ];
                strcpy ( y_equat, y_equat_ );
        }

        else
        {
                y_equat = NULL;
        }
}


template <typename T>
void Projection <T> ::setXEquatPostfix(const char * x_equat_postfix_)
{
	if (x_equat_postfix_ != NULL)
	{
		x_equat_postfix_ = new char[strlen(x_equat_postfix_) + 1];
		strcpy(x_equat_postfix, x_equat_postfix_);
	}

	else
	{
		x_equat_postfix = NULL;
	}
}


template <typename T>
void Projection <T> ::setYEquatPostfix(const char * y_equat_postfix_)
{
	if (y_equat_postfix_ != NULL)
	{
		y_equat_postfix = new char[strlen(y_equat_postfix_) + 1];
		strcpy(y_equat_postfix, y_equat_postfix_);
	}

	else
	{
		y_equat_postfix = NULL;
	}
}


template <typename T>
void Projection <T> ::setProjectionFamily(const char * projection_family_)
{
	if (projection_family_ != NULL)
	{
		projection_family = new char[strlen(projection_family_) + 1];
		strcpy(projection_family, projection_family_);
	}

	else
	{
		projection_family = NULL;
	}
}


template <typename T>
void Projection <T> ::setProjectionName ( const char * projection_name_ )
{
        if ( projection_name_ != NULL )
        {
                projection_name = new char [ strlen ( projection_name_ ) + 1 ];
                strcpy ( projection_name, projection_name_ );
        }

        else
        {
                projection_name = NULL;
        }
}


template <typename T>
void Projection <T> ::XEquatToPostfix()
{
	//Convert the infix notation to the postfix notation
	char x_equat_postfix_[MAX_TEXT_LENGTH];

	try
	{
		if (x_equat_postfix != NULL)
			delete x_equat_postfix;

		//Parse X equation
		ArithmeticParser::infixToPostfix(x_equat, x_equat_postfix_);
		x_equat_postfix = new char[strlen(x_equat_postfix_) + 1];
		strcpy(x_equat_postfix, x_equat_postfix_);
	}

	//Throw exception
	catch (ErrorMath <T> & error)
	{
		//Throw exception
		throw;
	}

	//Throw exception
	catch (Error & error)
	{
		//Throw exception
		throw error;
	}
}


template <typename T>
void Projection <T> ::YEquatToPostfix()
{
	//Convert the infix notation to the postfix notation
	char y_equat_postfix_[MAX_TEXT_LENGTH];

	try
	{
		//Convert the infix notation to the postfix notation
		if (y_equat_postfix != NULL)
			delete y_equat_postfix;

		//Parse Y equation
		ArithmeticParser::infixToPostfix(y_equat, y_equat_postfix_);
		y_equat_postfix = new char[strlen(y_equat_postfix_) + 1];
		strcpy(y_equat_postfix, y_equat_postfix_);
	}

	//Throw exception
	catch (ErrorMath <T> & error)
	{
		//Throw exception
		throw;
	}

	//Throw exception
	catch (Error & error)
	{
		//Throw exception
		throw error;
	}
}


#endif

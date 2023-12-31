/**
 * \file    Matrix.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    December, 2019
 *
 * \brief   Declares the Matrix class (3 by 3 only)
 * \details A 3 by 3 matrix is represented as a vertical stack of 3 row
 *          vectors.
 */

#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED 1

#include <iostream>
#include <iomanip>
#include <cmath> 
#include <list>
#include <sstream>

class Vector;

class Matrix
{
	public:
		/**
		 * \brief Default constructor; makes a 0 matrix.
		 */
		Matrix();

		/**
		 * \brief Constructor
		 * \param rowOne Row vectors of the matrix.
		 */
		Matrix(Vector rowOne, Vector rowTwo, Vector rowThree);

		/**
		 * \brief Copy constructor
		 * \param right The matrix to copy from
		 */
//		Matrix(Matrix& right);

		/**
		 * \brief Makes the current matrix into a diagonal matrix. 
		 */
		void diagonal(double top, double mid, double bot);

		/**
		 * \brief Data getter; gets one row vector at a time
		 * \param index Starting from 0, index of the row to get
		 */
		Vector getRow(int index) const;

		/**
		 * \brief Matrix - vector dot operation; returns a vector
		 */
		Vector dot(Vector& right);

		/**
		 * \brief Matrix  arithmetics
		 */
		Matrix operator+(const Matrix& right);
		Matrix operator-(const Matrix& right);
		Matrix operator/(const double denom);
		Matrix operator*(const double mult);

	private:
        // these have to be pointers because Vector is forward declared
		Vector* rowOne_;
		Vector* rowTwo_;
		Vector* rowThree_;


};

#endif

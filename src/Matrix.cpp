/**
 * \file    Matrix.cc
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    December, 2019
 *
 * \brief   Implements the mathematical Matrix class. Provides vector tensor
 *          products and other arithmetic operations.
 * 
 */

#include "Matrix.hpp"
#include "Vector.hpp"

Matrix::Matrix()
{
    rowOne_ = new Vector();
    rowTwo_ = new Vector();
    rowThree_ = new Vector();
}

Matrix::Matrix(Vector rowOne, Vector rowTwo, Vector rowThree)
{
	rowOne_ = new Vector(rowOne);
	rowTwo_ = new Vector(rowTwo);
	rowThree_ = new Vector(rowThree);
}

/*
Matrix::Matrix(Matrix& right)
		: rowOne_(right.getRow(0)), rowTwo_(right.getRow(1)), 
		rowThree(right.getRow(2))
{
	// nothing to do here
}
*/

void Matrix::diagonal(double top, double mid, double bot)
{
	rowOne_ = new Vector(top, 0,   0);
	rowTwo_ = new Vector(0,   mid, 0);
	rowThree_ = new Vector(0,   0,   bot);
	return;
}

Vector Matrix::getRow(int index) const
{
    Vector rtn;
	switch(index) {
		case 0:
//            std::cerr << "case 0";
            rtn =  *rowOne_;
            break;
		case 1:
//            std::cerr << "case 1";
            rtn =  *rowTwo_;
            break;
		case 2:
//            std::cerr << "case 2";
            rtn =  *rowThree_;
            break;
        default:
            std::cerr << "invalid row index at Matrix::getRow" << std::endl;
	}
    return rtn;
}

Vector Matrix::dot(Vector& right)
{
	double x = rowOne_->dot(right);
	double y = rowTwo_->dot(right);
	double z = rowThree_->dot(right);

	Vector result(x, y, z);
	return result;
}


Matrix Matrix::operator+(const Matrix& right)
{
    Vector rowOne = getRow(0) + right.getRow(0);
    Vector rowTwo = getRow(1) + right.getRow(1);
    Vector rowThree = getRow(2) + right.getRow(2);
    Matrix result(rowOne, rowTwo, rowThree);
    return result;
}

Matrix Matrix::operator-(const Matrix& right)
{
    Vector rowOne = getRow(0) - right.getRow(0);
    Vector rowTwo = getRow(1) - right.getRow(1);
    Vector rowThree = getRow(2) - right.getRow(2);
    Matrix result(rowOne, rowTwo, rowThree);
    return result;
}

Matrix Matrix::operator/(const double denom)
{
	Vector rowOne = getRow(0) / denom;
    Vector rowTwo = getRow(1) / denom;
    Vector rowThree = getRow(2) / denom;
    Matrix result(rowOne, rowTwo, rowThree);

	return result;
}

Matrix Matrix::operator*(const double mult)
{
	Vector rowOne = getRow(0) * mult;
    Vector rowTwo = getRow(1) * mult;
    Vector rowThree = getRow(2) * mult;
    Matrix result(rowOne, rowTwo, rowThree);

	return result;
}

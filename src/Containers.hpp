/**
 * \file    Containers.h
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    June, 2020
 *
 * \brief   Declares the VecDoub and MatDoub classes
 * 
 */

#ifndef CONTAINER_H_INCLUDED
#define CONTAINER_H_INCLUDED 1

#include <exception>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <vector>

typedef std::vector<double> VecDoub;

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
//    os << "[";
    for (int i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << std::endl;
    return os;
}

/**
 * \brief Data container class for 2D arrays
 * \details Row major. Matrices are stored as a set of row vectors.
 */
class MatDoub
{
	public:
		/**
		 * \brief Constructs a nw x nh array, filled with value a
		 * \param nw Number of horizontal elements (# of cols)
		 * \param nh Number of vertical elements (# of rows)
		 * \param a  Default to 0. Fill array with constant values.
		 */
		MatDoub(int nw, int nh, double a = 0);

        /**
         *\brief Constructs a single row matrix from a vector
         */
        MatDoub(VecDoub& vec);
    
		/**
		 * \brief Element access
		 * \param w  Index along row
		 * \param h Index along column
		 * \return An editable reference to matrix element
		 */
		double& at(int w, int h);

		/**
		 * \brief Get an editable reference to a certain row
		 */
		VecDoub& getRow(int h);
    
        /**
         *\brief Get a copy of the given column
         *\note This is not editable. The matrix is row major so we'll have to return a copy
         */
        VecDoub getCol(int w);

		/**
		 * \brief get number of elements in a row
		 */
		int getW();

		/**
		 * \brief get number of elements in a col
		 */
		int getH();
    
        /**
         *\brief Resizes the matrix to the new specified size
         *\param nw Number of horizontal elements
         *\param nh Number of vertical elements
         */
        void resize(int nw, int nh);
    
        /**
         *\brief Arbitrary matrix vector multiplication
         *\param right Vector to be dotted with. Length of right must match nw_
         */
        VecDoub dot(VecDoub& right);

		/**
         * \brief Define how a vector is printed
         */
        friend std::ostream& operator<<(std::ostream &os, const MatDoub& v);

	private:
		// Default constructor disabled.
		MatDoub();

		std::vector<std::vector<double>> data_;
		int nw_;
		int nh_;
};

inline std::ostream& operator<<(std::ostream &os, const MatDoub& m) {
    // printing one row at a time
    for (int ih = 0; ih < m.nh_; ++ih){
        for (int iw = 0; iw < m.nw_; ++iw){
            os << m.data_.at(ih).at(iw);
            if (iw != m.nw_ - 1){
                os << ",";
            } else {
                os << std::endl;
            }
        }
    }
    return os;
}

/**
 *\brief Housekeeping. Allows for exiting without memory leaks
 */
struct CtnrException : public std::exception {
    int code_; ///< Exit code
    CtnrException(int code): code_(code) {} // trivial constructor
    
    const char * what () const throw () {
        std::cerr << "CtnrException encountered: " << std::endl;
        switch (code_){
            case 0: std::cerr << "dimension mismatch in MatDoub.dot" << std::endl;
                exit(0); //normal exit
        }
        return "Uncaught exceptions";
   }
};



#endif // CONTAINER_H_INCLUDED

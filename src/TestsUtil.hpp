/**
 * \file    TestsUtil.hpp
 * \author  Xin Zhang, Princeton Plasma Physics Laboratory
 * \date    October, 2020
 *
 * \brief   Declares the TestsUtil class, a set of testing utilities
 */

#ifndef testsutil_H_INCLUDED
#define testsutil_H_INCLUDED 1


class TestsUtil
{
public:
    /**
    *\brief Check if two double are almost equal
    */
	inline bool doubEqual(double d1, double d2)
	{
		return fabs(d1 - d2) < PRECISION;
	}

    /**
     *\brief Check if two VecDoub are almost equal
     */
    inline bool vecEqual(VecDoub& v1, VecDoub& v2)
    {
    	auto it2 = v2.begin();
	    if (v1.size() != v2.size()){
	        return 0;
	    }
	    bool rtn(1);
	    for (auto it = v1.begin(); it != v1.end(); ++it) {
	        bool equal = ( fabs(*it - *it2) < PRECISION );
	        rtn = (rtn && equal);
	        ++it2;
	    }
	    return rtn;
    }

    /**
     * \brief Check if two MatDoub are almost equal
     */
    inline bool matEqual(MatDoub& m1, MatDoub& m2)
    {
    	bool rtn(1);
	    if (m1.getW() != m2.getW() || m1.getH() != m2.getH()){
	        return 0;
	    }
	    for (int i = 0; i < m1.getH(); ++i){
	        VecDoub v1 = m1.getRow(i);
	        VecDoub v2 = m2.getRow(i);
	        bool equal = vecEqual(v1, v2);
	        rtn = (rtn && equal);
	    }
	    return rtn;
    }

    
    inline VecDoub linspace(double start, double end, int num)
    {
        VecDoub out(num);
        double d = (end - start) / (num - 1);
        for(int i = 0; i<num; i++){
            out.at(i) = start + i * d;
        }
        return out;
    }
    
    inline VecDoub randArray(double min, double max, int num)
    {
        VecDoub out(num);
        std::default_random_engine generator(int(time(NULL)));
        // std::uniform_real_distribution<double> distribution(min,max);
        double mean = (min + max) / 2;
        double sigma = (mean - min)/2;
        std::normal_distribution<double> distribution(mean, sigma);
        for(int i = 0; i<num; i++){
            double val(min - 2);
            while (val < min || val > max){
                val = distribution(generator);
            }
            out.at(i) = val;
        }
        return out;
    }

    inline VecDoub meshgrid(std::vector<double>& rr, std::vector<double>& zz)
    {
        int nr = rr.size();
        int nz = zz.size();

        VecDoub out(nr * nz * 2);
        for(int i = 0; i < nr; i++){
            for(int j = 0; j < nz; j++){
                int index = (nz * i + j) * 2;
                out.at(index)    = rr.at(i);
                out.at(index +1) = zz.at(j);
            }
        }
        return out;
    }
    
    inline void meshgridMat(std::vector<double>& rr, std::vector<double>& zz, MatDoub& rout, MatDoub& zout)
    {
        int nr = rr.size();
        int nz = zz.size();

        MatDoub out(nr, nz);
        for(int ir = 0; ir < nr; ir++){
            for(int jz = 0; jz < nz; jz++){
                rout.at(ir, jz)    = rr.at(ir);
                zout.at(ir, jz)    = zz.at(jz);
            }
        }
        return;
    }
    
    inline void result(bool passed, std::string& what)
    {
        std::cerr << "Running test for: " << what << std::endl;

        if (passed){
            std::cerr << " ---- passed." << std::endl;
        } else {
            std::cerr << " ---- failed." << std::endl;
        }
    }
};

#endif

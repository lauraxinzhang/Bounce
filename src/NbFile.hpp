#include <iostream>
#include <netcdf>
#include <vector>

//using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

struct NbFile
{
    /**
     * \brief Data member. The NUBEAM output file to be read
     */
    NcFile * dataFile_;
    
    /**
     *\brief Constructor. Open the NcFile for read-only.
     */
    NbFile(std::string dirIn, std::string shotID, int cpuNo, int cpuTot);
    
    /**
     * \brief Data member. A list of whether each "particle" needs to be updated
     * TODO: ???? are we updating NUBEAM particle lists?
     */
//    std::vector<bool> update_;
    
    /**
     *\brief Read a single variable from the File
     *\param varname Name of the variable to be read
     */
    std::vector<double> readNBVar(std::string varname);
    
    /**
     *\brief Read a single variable from the file
     *\param varname Nameof variable to be read
     *\note Same as the previous function. Overloaded for different data type
     */
    std::vector<double> readNBVar(const char* varname);
    
    /**
     * \brief (TODO) update variable in the file
     */
//    void updateNBVar(std::string varname);
    
    /**
     * \brief Overloaded for different data type
     */
//    void updateNBVar(const char* varname);
    
};

inline NbFile::NbFile(std::string dirIn, std::string shotID, int cpuNo, int cpuTot)
{
    std::string cpu = "cpu" + std::to_string(cpuNo) + "_" + std::to_string(cpuTot);
    std::string fname = shotID + "_nbi_ptcl_state_" + cpu + ".cdf";
    
    try{
        std::cerr << "Opening netCDF file: " << fname << std::endl;
        dataFile_ = new NcFile(dirIn + fname, NcFile::read);
        //dataFile_ = &dataFile;
        std::cerr << "netCDF file opened." << std::endl;
    } catch (NcException& e) {
        std::cerr << "Failed to open file. Error thrown by NC: " << std::endl;
        std::cerr << e.what() << std::endl;
        exit(0);
    }
}

inline std::vector<double> NbFile::readNBVar(std::string varname)
{
    try {
        NcVar var = dataFile_->getVar(varname);
        if ( var.isNull() ) {
                std::cout << "var is null " << std::endl;
                exit(0);
        }
        // if data does exist
        std::cout << "got variable w/ dimension" << var.getDimCount() << std::endl;
        std::vector<NcDim> dims = var.getDims();
        size_t size0 = dims.at(0).getSize(); // this is always 2 for NUBEAM
        size_t size1 = dims.at(1).getSize(); // size of buffer
        std::cout << size0 << ',' << size1 << std::endl;
        std::vector<size_t> start {0, 0};
        std::vector<size_t> count {1, size1};
        double array[size1]; // data containers must be pre-allocated
        var.getVar(start, count, array); // get data value, store into npts
        std::vector<double> rtn(array, array + size1);
        std::cout << "got array of size: " << rtn.size() << std::endl;
        return rtn;
    } catch (NcException& e) {
        std::cerr << "Failed to get variable: " << varname << std::endl;
        std::cerr << e.what() << std::endl;
        exit(0);
    }
}

inline std::vector<double> NbFile::readNBVar(const char* varname)
{
    std::string name(varname);
    return readNBVar(name);
}

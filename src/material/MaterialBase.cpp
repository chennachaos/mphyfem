#include "MaterialBase.h"
#include <boost/algorithm/string.hpp>
#include <unordered_map>
#include "util.h"


using namespace std;



MaterialBase::MaterialBase()
{
  matData.resize(100);
  MIXED_ELEMENT = false;

  ResiMagnfield.resize(3);
  ApplMagnfield.resize(3);

  rho0 = 0.0;

  tis = "Galpha";
  spectralRadius = 0.0;
  td.resize(100);
  td.setZero();
}





MaterialBase::~MaterialBase()
{
}




void MaterialBase::printData()
{

    return;
}





int MaterialBase::readInput(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    int ii = 0;

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);
        cout << line << endl;

        if(line[0] == '}') break;


        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            if(stringlist[0] == "id")
            {
                id = stoi(stringlist[1]);
            }
            else if(stringlist[0] == "density")
            {
                rho0 = stod(stringlist[1]);
            }
            else if(stringlist[0] == "data")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

                matData.resize(stringlist.size());

                ii = 0;
                for(auto& str: stringlist)
                {
                    boost::trim(str);
                    matData[ii++] = stod(str);
                }

                printVector(matData);
            }
            else if(stringlist[0] == "data_deviatoric")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

                data_deviatoric.resize(stringlist.size());

                ii = 0;
                for(auto& str: stringlist)
                {
                    boost::trim(str);
                    data_deviatoric[ii++] = stod(str);
                }

                printVector(data_deviatoric);
            }
            else if(stringlist[0] == "data_volumetric")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
                for(auto& str: stringlist)  boost::trim(str);

                Utype = stoi(stringlist[0]);

                Kinv  = stod(stringlist[1]);
            }
            else if(stringlist[0] == "data_viscoelastic")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

                data_viscoelastic.resize(stringlist.size());

                ii = 0;
                for(auto& str: stringlist)
                {
                    boost::trim(str);
                    data_viscoelastic[ii++] = stod(str);
                }

                printVector(data_viscoelastic);
            }
            else if(stringlist[0] == "data_elecfield")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

                data_elecfield.resize(stringlist.size());

                ii = 0;
                for(auto& str: stringlist)
                {
                    boost::trim(str);
                    data_elecfield[ii++] = stod(str);
                }

                printVector(data_elecfield);
            }
            else if(stringlist[0] == "data_magnfield")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

                data_magnfield.resize(stringlist.size());

                ii = 0;
                for(auto& str: stringlist)
                {
                    boost::trim(str);
                    data_magnfield[ii++] = stod(str);
                }

                printVector(data_magnfield);
            }
            else if(stringlist[0] == "resi_magnfield")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

                if( stringlist.size() < 3)
                {
                    throw runtime_error("Insufficient data for residual magnetic field in MaterialBase::readInput ...");
                }

                ResiMagnfield.resize(3);

                ii = 0;
                for(auto& str: stringlist)
                {
                    boost::trim(str);
                    ResiMagnfield[ii++] = stod(str);
                }

                printVector(ResiMagnfield);
            }
            else if(stringlist[0] == "appl_magnfield")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

                if( stringlist.size() < 3)
                {
                    throw runtime_error("Insufficient data for applied magnetic field in MaterialBase::readInput ...");
                }

                ApplMagnfield.resize(3);

                ii = 0;
                for(auto& str: stringlist)
                {
                    boost::trim(str);
                    ApplMagnfield[ii++] = stod(str);
                }

                printVector(ApplMagnfield);
            }
            else if(stringlist[0] == "timescheme")
            {
                tis = stringlist[1];
            }
            else
            {
                throw runtime_error("MaterialBase::readInput ... invalid data type");
            }
        }
    }

    return 0;
}



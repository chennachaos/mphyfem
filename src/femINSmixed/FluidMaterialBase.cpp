#include "FluidMaterialBase.h"
#include <boost/algorithm/string.hpp>
#include <unordered_map>
#include "util.h"


using namespace std;



FluidMaterialBase::FluidMaterialBase()
{
  matData.resize(100);

  rho = 0.0;
  mu  = 0.1;
}





FluidMaterialBase::~FluidMaterialBase()
{
}




void FluidMaterialBase::printData()
{

    return;
}





int FluidMaterialBase::readInput(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    unordered_map<string,int>   map_keys = {
                        {"rho",                1},
                        {"mu",                 2},
                        {"data",               3},
                        };

    //getline(infile,line);    boost::trim(line);

    int ii = 0;

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        if(line[0] == '}') break;

        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            cout << line << endl;

            boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            unordered_map<string,int>::const_iterator got = map_keys.find(stringlist[0]);

            if( got == map_keys.end() )
            {
                throw runtime_error("Key not found in FluidMaterialBase::readInput ...");
                //return -1;
            }

            cout << " got->first  = " << got->first  << endl;
            cout << " got->second = " << got->second << endl;

            switch(got->second)
            {
                case 1:                                     // density

                    rho = stod(stringlist[1]);

                break;

                case 2:                                     // dynamic viscosity

                    mu = stod(stringlist[1]);

                break;

                default:

                break;

            }
        }
    }

    return 0;
}



#include "ElementTypeData.h"
#include <boost/algorithm/string.hpp>
#include <unordered_map>
#include <fstream>
#include "utilitiesmaterial.h"
#include "util.h"
#include "myIDmaps.h"

ElementTypeData::ElementTypeData()
{
    modeltype = "3D";
    thickness = 1.0;
    sss  = 2;
}



ElementTypeData::~ElementTypeData()
{
}


int ElementTypeData::readData(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    cout << line << endl;

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}')
        {
            getline(infile,line);
            break;
        }


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            if(stringlist[0] == "id")
            {
                id = stoi(stringlist[1]);
            }
            else if(stringlist[0] == "type")
            {
                elname = stringlist[1];
                
                elmTypeNameNum = getElementID_Standard(elname);
            }
            else if(stringlist[0] == "model")
            {
                modeltype = stringlist[1];

                if(modeltype == "pstress")
                  sss = 1;
                else if (modeltype == "pstrain")
                  sss = 2;
                else if (modeltype == "axsy")
                  sss = 3;
                else
                {
                  throw runtime_error("Model type not available in ElementTypeData::readData");
                }
            }
            else if(stringlist[0] == "thickness")
            {
                thickness = stod(stringlist[1]);
            }
            else if(stringlist[0] == "data")
            {
                boost::algorithm::split(stringlist, stringlist[1], boost::is_any_of(" "), boost::token_compress_on);
                for(auto& str: stringlist)  boost::trim(str);

                for(int i=0; i<stringlist.size(); i++)
                  elemData.push_back(stod(stringlist[i]));
                
                printVector(elemData);
            }
            else
            {
                throw runtime_error("Option not available in ElementTypeData::readData");
            }
        }
    }

    return 0;
}



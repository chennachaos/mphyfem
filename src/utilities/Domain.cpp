
#include "Domain.h"
#include <iostream>
#include "util.h"
#include <boost/algorithm/string.hpp>

using namespace std;


Domain::Domain()
{
}



Domain::~Domain()
{

}



Domain::Domain(string& label_, int tag_, int ndim_)
{
    label = label_;
    tag   = tag_;
    ndim  = ndim_;
}



int Domain::readData(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        //cout << line << endl;

        if(line[0] == '}')
        {
            getline(infile,line);
            break;
        }

        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            if(stringlist[0] == "Material")
            {
              matId = stoi(stringlist[1]) - 1;
            }
            else if(stringlist[0] == "Element")
            {
              elemId = stoi(stringlist[1]) - 1;
            }
            else
            {
                throw runtime_error("Option not available in Domain::readData");
            }
        }
    }

    return 0;
}



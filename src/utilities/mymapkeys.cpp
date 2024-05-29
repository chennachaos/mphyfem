
#include "mymapkeys.h"
#include <unordered_map>
#include <stdexcept>

using namespace std;


int  getDOFfromString(string& matkey)
{
    unordered_map<string,int>   map_keys_dofs = {
                        {"Xvelocity",       0},
                        {"Yvelocity",       1},
                        {"Zvelocity",       2},
                        //
                        {"phi",             0},
                        {"eta",             1},
                        //
                        {"VX",              0},
                        {"VY",              1},
                        {"VZ",              2},
                        //
                        {"Xdisplacement",   0},
                        {"Ydisplacement",   1},
                        {"Zdisplacement",   2},
                        //
                        {"UX",              0},
                        {"UY",              1},
                        {"UZ",              2},
                        };

    unordered_map<string,int>::const_iterator got = map_keys_dofs.find(matkey);

    if( got == map_keys_dofs.end() )
    {
        throw runtime_error("Key not found in getDOFfromString ...");
    }

    return got->second;
}



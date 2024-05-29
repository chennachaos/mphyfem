#ifndef incl_ElementTypeData_h
#define incl_ElementTypeData_h

#include <vector>
#include <iostream>

using namespace std;



class ElementTypeData
{
  public:

    //member variables
    int             id, elmTypeNameNum, sss;
    string          elname, modeltype;
    double          thickness, bodyForce[3];
    vector<double>  elemData;

    //member functions

    ElementTypeData();

    virtual ~ElementTypeData();

    int  getId()
    { return id; }

    string  getElementName()
    { return elname; }

    int  getElemTypeNameNum()
    { return  elmTypeNameNum; }

    int  getModeltypeNum()
    { return sss; }

    string  getModeltype()
    { return modeltype; }

    double  getThickness()
    { return thickness; }

    vector<double>  getData()
    { return  elemData; }

    int readData(ifstream& infile, string& line);

};

#endif


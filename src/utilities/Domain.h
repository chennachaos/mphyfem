
#ifndef incl_Domain_h
#define incl_Domain_h


#include <vector>
#include <string>

using namespace std;


class Domain
{ 
  public:
    string  label;

    int     ndim, tag, nElem, matId, elemId, ndof;

    Domain();

    Domain(string& label_, int tag_, int ndim_);

    ~Domain();

    void set(double dt_, double dtMin_, double dtMax_);

    int readData(ifstream& infile, string& line);


    void setTag(int  idd)
    {
      tag = idd;    return;
    }

    int  getTag()
    { return  tag;}

    void setLabel(string&  namee)
    {
      label = namee;    return;
    }

    string getLabel()
    {
      return label;
    }

    void setDimension(int  dim)
    {
      ndim = dim;    return;
    }

    int getDimension()
    {
      return ndim;
    }

    int getDomainTag()
    {
      return tag;
    }

    void setElementId(int  idd)
    {
      elemId = idd;    return;
    }

    int getElementType()
    {
      return elemId;
    }

    void setMatlId(int  idd)
    {
      matId = idd;    return;
    }

    int getMatlId()
    {
      return matId;
    }

};

#endif


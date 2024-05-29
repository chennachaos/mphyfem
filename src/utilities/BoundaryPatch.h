#ifndef incl_BoundaryPatch_h
#define incl_BoundaryPatch_h

#include <vector>
#include <iostream>
#include "myMathFunction.h"
#include "util.h"
#include <fstream>
#include <memory>


using namespace std;



class BoundaryPatch
{
  public:

    //member variables
    int             tag, ndim, nElem, nNode;
    string          label;
    bool            outputflag;
    myPoint         centroid;

    vector<int>              nodeNums;         // list of node numbers on the boundary patch
    vector<vector<int> >     elemConn;         // element-node connectivity for the boundary elements
    vector<vector<int> >     forAssyVecAll;    // element-dof connectivity for assembly of force due to traction

    vector<string>           BCType;           // Dirichlet or Neumann or Robin
    vector<string>           valuetype;        // specified or wall or traction
    vector<int>              dof_specified;    // specified DOF
    vector<string>           expressions;      // expressions for specified boundary conditions
    vector<int>              timeFunctionNums; // time function numbers

    ofstream  forcedata;

    //typedef  tuple<int, int, double>  data_dof;
    //vector<data_dof>  DirichletBCs;

    //member functions

    BoundaryPatch();

    BoundaryPatch(string& label_, int tag_, int ndim_);

    ~BoundaryPatch();

    void setTag(int tag_)
    { tag=tag_;  return; }

    int  getTag()
    { return tag; }

    void setNdim(int ndim_)
    { ndim = ndim_; return;}

    void  setLabel(string& labl)
    { label = labl; return; }

    string  getLabel()
    { return label; }

    int  getNumberOfElements()
    { return nElem; }

    int  getNumberOfNodes()
    { return nNode; }

    bool  getOutputFlag()
    { return outputflag; }

    void  setOutputFlag();

    void addElement(vector<int> & elnodeNums);

    int readData(ifstream& infile, string& line);

    void printData();

    void processData();
    
    void updateNodeNumbers(vector<int>& node_map_get_new);

    //void setSpecifiedDOFs(vector<int>& node_map_get_new, vector<vector<bool> >& NodeTypeOld, vector<int>& dofs_specified);

    //void setBoundaryConditions(vector<int>& node_map_get_old, vector<myPoint>& node_coords, VectorXd& solnApplied);

};

#endif


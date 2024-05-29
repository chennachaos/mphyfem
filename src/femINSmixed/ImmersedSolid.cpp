
#include "ImmersedSolid.h"
#include "SolutionData.h"
#include "GeomDataLagrange.h"
#include "ImmersedIntegrationElement.h"
#include <boost/algorithm/string.hpp>
#include "BasisFunctionsLagrange.h"


int ImmersedSolid::count = 0;


ImmersedSolid::ImmersedSolid()
{
  id = count++;

  nElem = 0;
  PENALTY = 1.0e2;
  BC_ENFORCE_TYPE = BC_ENFORCE_TYPE_LAGRANGE;

  isNitsche = false;
  NitscheFact = 1.0;

  uGrid    =  vtkSmartPointer<vtkUnstructuredGrid>::New();

  polyDataVTK =  vtkSmartPointer<vtkPolyData>::New();

  totalForce.resize(6);
  totalForce.setZero();
}




ImmersedSolid::~ImmersedSolid()
{
  for(vector<ImmersedIntegrationElement*>::iterator pObj = ImmIntgElems.begin();
       pObj != ImmIntgElems.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  ImmIntgElems.clear(); // Purge the contents so no one tries to delete them again

  --count;
}




void  ImmersedSolid::readInput(string& fname)
{
    cout << " ImmersedSolid::readInput \n\n " << endl;

    ifstream  infile(fname);

    if(infile.fail())
    {
       cout << " Could not open immersed elements file " << endl;
       exit(-1);
    }

    string line;
    vector<string>  stringlist;

    while(getline(infile,line))
    {
        //if (line == "") continue;

        //std::cout << line << std::endl;

        //cout << " size       = " << line.size() << endl;

        boost::trim(line);
 
        //cout << " size       = " << line.size() << endl;
        //cout << " first word = " << line[0] << endl;

        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            //std::cout << line << std::endl;

            if(line.compare(string("Dimension")) == 0)
            {
                cout << "Dimension" << endl;
                readDimension(infile, line);
            }
            else if(line.compare(string("Nodes")) == 0)
            {
                cout << "Nodes" << endl;
                readNodes(infile, line);
            }
            else if(line.compare(string("Immersed Boundary Elements")) == 0)
            {
                cout << "Immersed Boundary Elements" << endl;
                readImmersedBoundaryElements(infile, line);
            }
            else if(line.compare(string("Prescribed Motion")) == 0)
            {
                cout << "Prescribed Motion" << endl;
                readPrescribedMotion(infile, line);
            }
            else
            {
                cout << "key =  " <<  line << endl;
                throw runtime_error("Key not found in ImmersedSolid::readInput...");
                //return -1;
            }
      }
    }

    infile.close();

    cout << " Immersed solid file has been read successfully \n\n " << endl;

    return;
}



void ImmersedSolid::readDimension(ifstream& infile, string& line)
{
    // read the value of dimension
    infile >> ndim ;

    cout << "ndim = " << ndim << endl;

    assert( (ndim > 0) && (ndim < 4) );

    return;
}



void ImmersedSolid::readNodes(ifstream& infile, string& line)
{
    // read the number of nodes
    infile >> nNode ;

    cout << "nNode = " << nNode << endl;

    nodePosOrig.resize(nNode);
    nodePosCur.resize(nNode);

    int  ind;
    if(ndim == 1)
    {
        for(int ii=0; ii<nNode; ii++)
        {
            infile  >>  ind >>  nodePosOrig[ii][0];

            nodePosCur[ii][0] = nodePosOrig[ii][0];
        }
    }
    else if(ndim == 2)
    {
        for(int ii=0; ii<nNode; ii++)
        {
            infile  >>  ind >>  nodePosOrig[ii][0] >>  nodePosOrig[ii][1];

            nodePosCur[ii][0] = nodePosOrig[ii][0];
            nodePosCur[ii][1] = nodePosOrig[ii][1];
        }
    }
    else
    {
        for(int ii=0; ii<nNode; ii++)
        {
            infile  >>  ind >>  nodePosOrig[ii][0] >>  nodePosOrig[ii][1] >>  nodePosOrig[ii][2];

            nodePosCur[ii][0] = nodePosOrig[ii][0];
            nodePosCur[ii][1] = nodePosOrig[ii][1];
            nodePosCur[ii][2] = nodePosOrig[ii][2];
        }
    }

    specVelocity.resize(nNode);
    for(int ii=0; ii<nNode; ii++)
    {
        specVelocity[ii].setZero();
    }    
    //cout.precision(17);
    //for(int ii=0; ii<nNode; ii++)
      //cout << ii << '\t' << setprecision(12) << nodePosData[ii][0]  << '\t' << setprecision(12) << nodePosData[ii][1]  << '\t' << setprecision(12) << nodePosData[ii][2] << endl;
    //cout << endl;    cout << endl;

    return;
}


void  ImmersedSolid::readImmersedBoundaryElements(ifstream& infile, string& line)
{
    // read the number of elements
    getline(infile,line);    boost::trim(line);

    nElem = stoi(line);

    cout << "nElem = " << nElem << endl;

    vector<string>  stringlist;
    bool activeflag = true;
    vector<int>   nodeNums;

    ImmIntgElems.resize(nElem);

    for(int ee=0; ee<nElem; ee++)
    {
        getline(infile,line);
        boost::trim(line);
        cout << line << endl;
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        //stringlist.erase(std::remove(stringlist.begin(), stringlist.end(), " "), stringlist.end());

        int ind = stringlist.size();
        cout << "ind = " << ind << '\t' << stoi(stringlist[1]) << endl;

        activeflag = (stoi(stringlist[1]) == 1);

        nodeNums.resize(ind-2);
        for(int j=2; j<ind; j++)
            nodeNums[j-2] = stoi(stringlist[j]) - 1;

        //printVector(enodeNums);
        
        ImmIntgElems[ee] = new ImmersedIntegrationElement;

        ImmIntgElems[ee]->nodeNums = nodeNums;

        if(!activeflag)
          ImmIntgElems[ee]->deactivate();

        ImmIntgElems[ee]->initialiseDOFvalues();
    }
    
    npElem = nodeNums.size();

    return;
}






void  ImmersedSolid::readPrescribedMotion(ifstream& infile, string& line)
{
    // read {
    getline(infile,line);    boost::trim(line);

    vector<string>  stringlist;

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;
    }

    return;
}


void ImmersedSolid::computePointAtGP(int elnum, double gpcoord, myPoint& pt)
{
  int  ii, jj;
  vector<double>  N(npElem);

  Lagrange_BasisFuns1D(npElem-1, gpcoord, &N[0]);

  vector<int> nodeNums = ImmIntgElems[elnum]->nodeNums;
  pt.setZero();
  for(ii=0; ii<npElem; ii++)
  {
    for(jj=0; jj<ndim; jj++)
      pt[jj] += (nodePosCur[nodeNums[ii]][jj] * N[ii] );
  }

  return;
}





void ImmersedSolid::reset()
{

  return;
}





void  ImmersedSolid::setDataForOutput(vector<vector<int> >& vectemp)
{
    OutputData.resize(vectemp.size());
    for(int ii=0; ii<vectemp.size(); ii++)
    {
      OutputData[ii] = vectemp[ii];
    }

    return;
}



void  ImmersedSolid::computeCentroid(int index)
{
    centroid.setZero();

    switch(index)
    {
      case 0: // original configuration

        for(int aa=0; aa<nNode; aa++)
          centroid += nodePosOrig[aa];

      break;

      case 1: // configuration at t_{n+1}

        for(int aa=0; aa<nNode; aa++)
          centroid += nodePosCur[aa];

      break;
    }

    centroid /= nNode;

    bool printOut = 0;
    if(printOut)
    {
      printf("\n  Centroid of the immersed body ... %5d \n", id);
      printf(" X-dir \t Y-dir \t Z-dir \n ");
      if(ndim == 2)
        printf(" %12.6f \t  %12.6f \t  %12.6f \n", centroid[0], centroid[1], 0.0);
      else
        printf(" %12.6f \t  %12.6f \t  %12.6f \n", centroid[0], centroid[1], centroid[2]);
    }

    return;
}




void  ImmersedSolid::computeAABB(int index)
{
    bbox.initialize();

    int  ii=0, aa=0;
    double  val=0.0;

    switch(index)
    {
      case 0: // original configuration

        for(aa=0; aa<nNode; aa++)
        {
          for(ii=0; ii<ndim; ii++)
          {
            val = nodePosOrig[aa][ii];

            if( val < bbox.minBB[ii]  )
              bbox.minBB[ii] = val ;

            if( val > bbox.maxBB[ii]  )
              bbox.maxBB[ii] = val ;
          }
        }

      break;

      case 1: // configuration at t_{n+1}

        for(aa=0; aa<nNode; aa++)
        {
          for(ii=0; ii<ndim; ii++)
          {
            val = nodePosCur[aa][ii];

            if( val < bbox.minBB[ii]  )
              bbox.minBB[ii] = val ;

            if( val > bbox.maxBB[ii]  )
              bbox.maxBB[ii] = val ;
          }
        }

      break;
    }

    return;
}



void  ImmersedSolid::computeTotalForce()
{
    totalForce.resize(6);
    totalForce.setZero();
    if(isBoundaryConditionTypeLagrange())
    {
      computeCentroid(1);

      for(int aa=0;aa<ImmIntgElems.size();aa++)
      {
        ImmIntgElems[aa]->integrateForceAndMoment(totalForce, centroid);
      }
    }

    bool printOut = 0;
    if(printOut)
    {
      printf("\n  Sum of the forces on immersed body ... %5d \n", id);
      printf(" X-dir \t Y-dir \t Z-dir \n ");
      printf(" %12.6f \t  %12.6f \t  %12.6f \n\n\n", totalForce[0], totalForce[1], totalForce[2]);
      printf("\n  Sum of the moments on immersed body ... %5d \n", id);
      printf(" X-dir \t Y-dir \t Z-dir \n ");
      printf(" %12.6f \t  %12.6f \t  %12.6f \n\n\n", totalForce[3], totalForce[4], totalForce[5]);
    }

  return;
}



int ImmersedSolid::updatePointPositions(double timeNow)
{
    // rigid body translation

    double  dispX = 0.3*sin(2.0*PI*1.0*timeNow), dispY=0.0;
    centroid[0]=4.0, centroid[1]=0.2;

    for(int ii=0; ii<nNode; ii++)
    {
        nodePosCur[ii][0] = nodePosOrig[ii][0] + dispX;
        nodePosCur[ii][1] = nodePosOrig[ii][1] + dispY;
    }

/*
    // rigid body rotation

    double  xx, yy, theta = (PI/4.0)*sin(2.0*PI*1.0*timeNow);
    double  cst = cos(theta), snt = sin(theta);
    centroid[0]=4.0, centroid[1]=0.2;

    for(int ii=0; ii<nNode; ii++)
    {
        xx = nodePosOrig[ii][0] - centroid[0];
        yy = nodePosOrig[ii][1] - centroid[1];

        nodePosCur[ii][0] = centroid[0] + xx*cst + yy*snt;
        nodePosCur[ii][1] = centroid[1] - xx*snt + yy*cst;
    }
*/

    return 0;
}






int ImmersedSolid::updatePointPositions(int ndof, VectorXd& disp)
{
    int ii, jj, ind;
    for(ii=0; ii<nNode; ii++)
    {
        ind = ii*ndof;
        for(jj=0; jj<ndim; jj++)
          nodePosCur[ii][jj] = nodePosOrig[ii][jj] + disp[ind+jj];
    }

    return 0;
}







int ImmersedSolid::updatePointVelocities(int ndof, VectorXd& velo)
{
    int ii, jj, ind;
    for(ii=0; ii<nNode; ii++)
    {
        ind = ii*ndof;
        for(jj=0; jj<ndim; jj++)
          specVelocity[ii][jj] = velo[ind+jj];
    }

    return 0;
}























/* Node ordering for the quadratic triangular element
 3
 | |
 |   |
 6     5
 |       |
 |         |
 1-----4-----2
*/

/* Node ordering for the quadratic quadrilateral element
 4-----7-----3
 |     |     |
 |     |     |
 8-----9-----6
 |     |     |
 |     |     |
 1-----5-----2
*/



#include "BernsteinElem3DFace.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "BasisFunctionsBernstein.h"


using namespace std;

extern MyTime myTime;
extern List<TimeFunction> timeFunction;



BernsteinElem3DFace::BernsteinElem3DFace()
{
  ndim   = 3;
  ndof   = 3;
}

BernsteinElem3DFace::~BernsteinElem3DFace()
{
}


void BernsteinElem3DFace::prepareElemData()
{
    npElem = nodeNums.size();

    //assert(nodeNums.size() == npElem);

    vals2project.resize(npElem);

    nlbfU = npElem;
    nsize = nlbfU*ndof;

    globalDOFnums.resize(nsize);

    int ii, jj, kk=0, ind;
    for(ii=0; ii<nodeNums.size(); ii++)
    {
      ind = nodeNums[ii]*ndof;
      for(jj=0;jj<ndof;jj++)
        globalDOFnums[kk++] = ind+jj;
    }

    //primvar.setDim(nsize);
    //resi.setDim(nsize);
    //printVector(globalDOFnums);

    // set the element property variables

    //cout <<  elmType <<  '\t' <<  matType <<  endl;

    elmDat = &(SolnData->ElemProp[elmType]->data[0]);
    matDat = &(SolnData->MatlProp[matType]->data[0]);

    nGP   = (int) elmDat[0] ;

    thick = 1.0;


    degree = 2;
    if( (npElem == 3) || (npElem == 4) )
    {
      degree = 1;
    }

    ELEM_SHAPE = ELEM_SHAPE_QUAD_BERNSTEIN;
    if( (npElem == 3) || (npElem == 6) )
    {
      ELEM_SHAPE = ELEM_SHAPE_TRIA_BERNSTEIN;
    }

    return;
}



double BernsteinElem3DFace::computeVolume(bool init)
{
  double  dvol, Jac, param[2];

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

  int   ii, gp;

  double xNode[npElem], yNode[npElem], xx, yy;
  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
  }

  vector<double>  gausspoints1, gausspoints2, gaussweights;

  getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    elemVol=0.0;
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          if(init)
            GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          else
            GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);

          xx = yy= 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            xx += N[ii]*xNode[ii];
            yy += N[ii]*yNode[ii];
          }

          if(axsy)
            dvol *= 2.0*PI*xx;

      elemVol += dvol;
  }//gp

  //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

  return elemVol;
}




int BernsteinElem3DFace::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    return 0;
}



int BernsteinElem3DFace::calcResidual(VectorXd& Flocal)
{
    return 0;
}






int BernsteinElem3DFace::calcLoadVector(VectorXd& Flocal)
{
    int  ii, jj, TI, TIp1, TIp2, gp;
    double  param[3], trac[3], N[npElem];
    double  Jac, dvol, specVal, curvature;

    int  config = (int) elmDat[1];
    double tx = elmDat[2];
    double ty = elmDat[3];
    double tz = elmDat[4];
    double gamma = elmDat[5];

    double normal[3], tangent1[3], tangent2[3];

    double  xNode[npElem], yNode[npElem], zNode[npElem], xx, yy, zz;
    if(config)
    {
      for(ii=0; ii<npElem; ii++)
      {
        xNode[ii] = GeomData->NodePosCur[nodeNums[ii]][0];
        yNode[ii] = GeomData->NodePosCur[nodeNums[ii]][1];
        zNode[ii] = GeomData->NodePosCur[nodeNums[ii]][2];
      }
    }
    else
    {
      for(ii=0; ii<npElem; ii++)
      {
        xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
        yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
        zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
      }
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
      getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);
    else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);
    else
    {
        throw runtime_error("Element type no defined...");
    }


    if(Flocal.rows() != nsize)
      Flocal.resize(nsize);
    Flocal.setZero();


    elemVol=0.0;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        //cout << " gp = " << param[0] << '\t' << param[1] << '\t' << npElem << '\t' << nlbfU << endl;

        if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
          BernsteinBasisFunsFaceQuad(degree, param, xNode, yNode, zNode, N, normal, tangent1, tangent2, curvature, Jac);
        else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
          BernsteinBasisFunsFaceTria(degree, param, xNode, yNode, zNode, N, normal, tangent1, tangent2, curvature, Jac);

        dvol = gaussweights[gp]*Jac;
        elemVol += dvol;

        //cout << " normal = " << normal[0] << '\t' << normal[1] << '\t' << normal[2] << '\t' << Jac << '\t' << dvol << endl;

        xx = yy = zz = 0.0;
        for(ii=0; ii<nlbfU; ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
          zz += N[ii]*zNode[ii];
        }

        //  specified traction values
        trac[0] = (tx*normal[0]+ty*tangent1[0]+tz*tangent2[0])*dvol;
        trac[1] = (tx*normal[1]+ty*tangent1[1]+tz*tangent2[1])*dvol;
        trac[2] = (tx*normal[2]+ty*tangent1[2]+tz*tangent2[2])*dvol;

        // specified surface tension
        //trac[0] += (2.0*curvature*gamma*normal[0])*dvol;
        //trac[1] += (2.0*curvature*gamma*normal[1])*dvol;
        //trac[2] += (2.0*curvature*gamma*normal[2])*dvol;

        //cout << " Funcs  = " << N[0] << '\t' << N[1] << '\t' << N[2] << '\t' << N[3] << '\t' << N[4] << '\t' << N[5] << endl;
        //cout << " trac   = " << trac[0] << '\t' << trac[1] << '\t' << trac[2] << endl;

        for(ii=0; ii<npElem; ii++)
        {
          TI   = ii*ndof;
          TIp1 = TI+1;
          TIp2 = TI+2;

          Flocal(TI)   += N[ii]*trac[0];
          Flocal(TIp1) += N[ii]*trac[1];
          Flocal(TIp2) += N[ii]*trac[2];
        }
    } // for(gp=0; gp<nGP; gp++)
    //cout << " Volume = " << elemVol << endl;

    //printVector(Flocal);

    return 0;
}



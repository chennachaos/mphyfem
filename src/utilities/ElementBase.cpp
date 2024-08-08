
#include "ElementBase.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "BasisFunctionsBernstein.h"
#include "ExactSolutionsElasticity.h"
#include "utilitiesmaterial.h"


extern MyTime           myTime;
extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern bool debug;

using namespace std;



ElementBase::ElementBase()
{
  //if(debug)  cout << "     ElementBase: constructor ...\n\n";

  elenum = elemcount++;

  nlbfU = ndof = nsize = nivGP = nGP = elenum = 0;
  processorId = 0;

}





ElementBase::~ElementBase()
{
  //if(debug)  cout << "     ElementBase: destructor ...\n\n";

}




void ElementBase::prepareElemData()
{
    //printVector(nodeNums);
    globalDOFnums.resize(nsize);

    int ii, jj, kk=0, ind;
    for(ii=0; ii<nodeNums.size(); ii++)
    {
      ind = nodeNums[ii]*ndof;
      for(jj=0;jj<ndof;jj++)
        globalDOFnums[kk++] = ind+jj;
    }

    int  idd = MatlData->getMaterialTypeNameNumber();
    //cout << " idd = " << idd << endl;
    //if( !( (idd == 1) || (idd == 101) || (idd == 102)) && (!finite) )
    //{
        //throw runtime_error("\n\n Hyperelastic materials are not supported with Small-strain option \n\n");
    //}

    initialiseIntVar();

    return;
}




void ElementBase::setnivGP()
{

    return;
}




void ElementBase::initialiseIntVar()
{
    // set nivGP value
    switch(MatlData->getMaterialTypeNameNumber())
    {
        case 101:
        case 103:
        case 104:
        case 105:
        case 106:
        case 2006:
        case 2007:

            nivGP = (int) MatlData->data_viscoelastic[0];

        break;

        default:
            nivGP = 0;
        break;
    }

    //cout << " nGP   = " << nGP << endl;
    //cout << " nivGP = " << nivGP << endl;

    // allocate memory for internal variables
    if(nivGP > 0)
    {
        ivar.initialise(nGP, nivGP);
    }

    return;
}







double  ElementBase::computeGeomOrig(int dir, VectorXd& NN)
{
  double val=0.0;

  for(int ii=0;ii<nlbfU;ii++)
    val += GeomData->NodePosOrig[nodeNums[ii]][dir] * NN[ii];

  return  val;
}


double  ElementBase::computeGeomNew(int dir, VectorXd& NN)
{
  double val=0.0;

  for(int ii=0;ii<nlbfU;ii++)
    val += GeomData->NodePosNew[nodeNums[ii]][dir] * NN[ii];

  return  val;
}



double  ElementBase::computeGeomCur(int dir, VectorXd& NN)
{
  double val=0.0;

  for(int ii=0;ii<nlbfU;ii++)
    val += GeomData->NodePosCur[nodeNums[ii]][dir] * NN[ii];

  return  val;
}


void  ElementBase::computeGeomNew(double* param,  double* geom)
{
    VectorXd  N(nlbfU);

    if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
        BernsteinBasisFunsTria(degree, param[0], param[1], &N(0));
    else if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
        BernsteinBasisFunsQuad(degree, param[0], param[1], &N(0));
    else if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
        BernsteinBasisFunsTetra(degree, param[0], param[1], param[2], &N(0));
    else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
        BernsteinBasisFunsWedge(degree, param[0], param[1], param[2], &N(0));
    else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
        BernsteinBasisFunsHexa(degree, param[0], param[1], param[2], &N(0));

    if(ndim == 2)
    {
      geom[0]=geom[1]=0.0;
      for(int ii=0;ii<nlbfU;ii++)
      {
        geom[0] += GeomData->NodePosNew[nodeNums[ii]][0] * N[ii];
        geom[1] += GeomData->NodePosNew[nodeNums[ii]][1] * N[ii];
      }
    }
    else
    {
      geom[0]=geom[1]=geom[2]=0.0;
      for(int ii=0;ii<nlbfU;ii++)
      {
        geom[0] += GeomData->NodePosNew[nodeNums[ii]][0] * N[ii];
        geom[1] += GeomData->NodePosNew[nodeNums[ii]][1] * N[ii];
        geom[2] += GeomData->NodePosNew[nodeNums[ii]][2] * N[ii];
      }
    }

    return;
}


void  ElementBase::computeDefGrad2D(VectorXd& dN_dx, VectorXd& dN_dy, MatrixXd& F)
{
    double  ux, uy;
    F = Matrix3d::Identity();

    for(int ii=0;ii<nlbfU;ii++)
    {
        int ind = nodeNums[ii]*ndof;

        ux = SolnData->disp[ind];
        uy = SolnData->disp[ind+1];

        F(0,0) += (ux*dN_dx[ii]);
        F(0,1) += (ux*dN_dy[ii]);

        F(1,0) += (uy*dN_dx[ii]);
        F(1,1) += (uy*dN_dy[ii]);
    }

    return;
}




void  ElementBase::computeDefGrad2DPrev(VectorXd& dN_dx, VectorXd& dN_dy, MatrixXd& F)
{
    double  ux, uy;
    F = Matrix3d::Identity();

    for(int ii=0;ii<nlbfU;ii++)
    {
        int ind = nodeNums[ii]*ndof;

        ux = SolnData->dispPrev[ind];
        uy = SolnData->dispPrev[ind+1];

        F(0,0) += (ux*dN_dx[ii]);
        F(0,1) += (ux*dN_dy[ii]);

        F(1,0) += (uy*dN_dx[ii]);
        F(1,1) += (uy*dN_dy[ii]);
    }

    return;
}




void  ElementBase::computeDefGrad2DCur(VectorXd& dN_dx, VectorXd& dN_dy, MatrixXd& F)
{
    double  ux, uy;
    F = Matrix3d::Identity();

    for(int ii=0;ii<nlbfU;ii++)
    {
        int ind = nodeNums[ii]*ndof;

        ux = SolnData->dispCur[ind];
        uy = SolnData->dispCur[ind+1];

        F(0,0) += (ux*dN_dx[ii]);
        F(0,1) += (ux*dN_dy[ii]);

        F(1,0) += (uy*dN_dx[ii]);
        F(1,1) += (uy*dN_dy[ii]);
    }

    return;
}



void  ElementBase::computeDefGrad(VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz, MatrixXd& F)
{
    double  ux, uy, uz;
    F = Matrix3d::Identity();

    for(int ii=0;ii<nlbfU;ii++)
    {
        int ind = nodeNums[ii]*ndof;

        ux = SolnData->disp[ind];
        uy = SolnData->disp[ind+1];
        uz = SolnData->disp[ind+2];

        F(0,0) += (ux*dN_dx[ii]);
        F(0,1) += (ux*dN_dy[ii]);
        F(0,2) += (ux*dN_dz[ii]);

        F(1,0) += (uy*dN_dx[ii]);
        F(1,1) += (uy*dN_dy[ii]);
        F(1,2) += (uy*dN_dz[ii]);

        F(2,0) += (uz*dN_dx[ii]);
        F(2,1) += (uz*dN_dy[ii]);
        F(2,2) += (uz*dN_dz[ii]);
    }

    return;
}




void  ElementBase::computeDefGradPrev(VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz, MatrixXd& F)
{
    double  ux, uy, uz;
    F = Matrix3d::Identity();

    for(int ii=0;ii<nlbfU;ii++)
    {
        int ind = nodeNums[ii]*ndof;

        ux = SolnData->dispPrev[ind];
        uy = SolnData->dispPrev[ind+1];
        uz = SolnData->dispPrev[ind+2];

        F(0,0) += (ux*dN_dx[ii]);
        F(0,1) += (ux*dN_dy[ii]);
        F(0,2) += (ux*dN_dz[ii]);

        F(1,0) += (uy*dN_dx[ii]);
        F(1,1) += (uy*dN_dy[ii]);
        F(1,2) += (uy*dN_dz[ii]);

        F(2,0) += (uz*dN_dx[ii]);
        F(2,1) += (uz*dN_dy[ii]);
        F(2,2) += (uz*dN_dz[ii]);
    }

    return;
}




void  ElementBase::computeDefGradCur(VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz, MatrixXd& F)
{
    double  ux, uy, uz;
    F = Matrix3d::Identity();

    for(int ii=0;ii<nlbfU;ii++)
    {
        int ind = nodeNums[ii]*ndof;

        ux = SolnData->dispCur[ind];
        uy = SolnData->dispCur[ind+1];
        uz = SolnData->dispCur[ind+2];

        F(0,0) += (ux*dN_dx[ii]);
        F(0,1) += (ux*dN_dy[ii]);
        F(0,2) += (ux*dN_dz[ii]);

        F(1,0) += (uy*dN_dx[ii]);
        F(1,1) += (uy*dN_dy[ii]);
        F(1,2) += (uy*dN_dz[ii]);

        F(2,0) += (uz*dN_dx[ii]);
        F(2,1) += (uz*dN_dy[ii]);
        F(2,2) += (uz*dN_dz[ii]);
    }

    return;
}



double  ElementBase::computeDisplacement(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->disp[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}


double  ElementBase::computeDisplacementCur(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->dispCur[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}




double  ElementBase::computeVelocity(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->velo[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}


double  ElementBase::computeVelocityCur(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->veloCur[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}



double  ElementBase::computeAcceleration(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->acce[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}


double  ElementBase::computeAccelerationCur(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->acceCur[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}


double  ElementBase::computeValue(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->pres[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}



double  ElementBase::computeValueIncr(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->presIncr[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}



double  ElementBase::computeValuePrev(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->presPrev[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}




double  ElementBase::computeValueCur(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->dispCur[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}




double  ElementBase::computeValueDot(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->velo[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}



double  ElementBase::computeValueDotCur(int dir, VectorXd& NN)
{
   double val=0.0;
   for(int ii=0;ii<nlbfU;ii++)
     val += SolnData->dispDotCur[nodeNums[ii]*ndof+dir] * NN[ii];

   return  val;
}




double  ElementBase::computeValue2(int dir, VectorXd& NN)
{
    double val=0.0;
    for(int ii=0;ii<nlbfP;ii++)
      //val += SolnData->pres[nodeNums[ii]]*NN[ii];
      val += SolnData->pres[forAssyVecPres[ii]]*NN[ii];

    return  val;
}



double  ElementBase::computeValue2Prev(int dir, VectorXd& NN)
{
    double val=0.0;
    for(int ii=0;ii<nlbfP;ii++)
      //val += SolnData->var2Prev[nodeNums[ii]]*NN[ii];
      val += SolnData->presPrev[forAssyVecPres[ii]]*NN[ii];

    return  val;
}



double  ElementBase::computeValue2Cur(int dir, VectorXd& NN)
{
    double val=0.0;
    for(int ii=0;ii<nlbfP;ii++)
      //val += SolnData->var2Cur[nodeNums[ii]]*NN[ii];
      val += SolnData->presCur[forAssyVecPres[ii]]*NN[ii];

    return  val;
}



double  ElementBase::computeValue3(int dir, VectorXd& NN)
{
    double val=0.0;
    for(int ii=0;ii<nlbfF;ii++)
      val += SolnData->mpot[nodeNums[ii]]*NN[ii];

    return  val;
}



double  ElementBase::computeValue3Prev(int dir, VectorXd& NN)
{
    double val=0.0;
    for(int ii=0;ii<nlbfF;ii++)
      val += SolnData->mpotPrev[nodeNums[ii]]*NN[ii];

    return  val;
}



double  ElementBase::computeValue3Cur(int dir, VectorXd& NN)
{
    double val=0.0;
    for(int ii=0;ii<nlbfF;ii++)
      val += SolnData->mpotCur[nodeNums[ii]]*NN[ii];

    return  val;
}




int  ElementBase::applyDirichletBCs(MatrixXd& Klocal, VectorXd& Flocal)
{
  // add contributions to the rhs vector
  // from nodes with specified displacement BCs

  int aa, bb, ii, jj;
  double fact=0.0;

  for(ii=0; ii<nsize; ii++)
  {
    aa = forAssyVec[ii];
    if(aa == -1) // this DOF has a prescibed value
    {
      fact = SolnData->dispApplied[globalDOFnums[ii]];

      for(jj=0; jj<nsize; jj++)
      {
        bb = forAssyVec[jj];
        if( bb != -1 )
        {
          Flocal(jj) -= Klocal(jj, ii) * fact;
        }
      }
    }
  }

  return 0;
}





int  ElementBase::applyDirichletBCsMixed(int var2_offset, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2)
{
  // add contributions to the rhs vector
  // from nodes with specified displacement BCs

  int aa, bb, ii, jj, sizeP = forAssyVecPres.size();
  double fact=0.0;

  if(var2_offset == 0) // for the explicit scheme
  {
    // applied displacement (velocity) dof for the solid (fluid)
    for(ii=0; ii<nsize; ii++)
    {
      aa = forAssyVec[ii];
      if(aa == -1) // this DOF has a prescibed value
      {
        fact = SolnData->dispApplied[globalDOFnums[ii]];

        for(jj=0; jj<sizeP; jj++)
        {
          //cout << Flocal2(jj) << endl;
          //cout << ii << '\t' << aa << '\t' << globalDOFnums[ii] << '\t' << jj << '\t' << forAssyVecPres[jj] << endl;
          //cout << Flocal2(jj) << '\t' << Kpu(jj, ii) << endl;
          if( forAssyVecPres[jj] != -1 )
          {
            Flocal2(jj) -= Kpu(jj, ii) * fact;
          }
        }
      }
    }

    // applied pressure
    for(ii=0; ii<sizeP; ii++)
    {
      aa = forAssyVecPres[ii];
      if(aa == -1) // this DOF has a prescibed value
      {
        fact = SolnData->presApplied[nodeNums[ii]];

        for(jj=0; jj<nsize; jj++)
        {
          if( forAssyVec[jj] != -1 )
          {
            Flocal1(jj) -= Kup(jj, ii) * fact;
          }
        }

        for(jj=0; jj<sizeP; jj++)
        {
          if( forAssyVecPres[jj] != -1 )
          {
            Flocal2(jj) -= Kpp(jj, ii) * fact;
          }
        }
      }
    }
  }
  else  // for the implicit scheme
  {
    // applied displacement (velocity) dof for the solid (fluid)
    for(ii=0; ii<nsize; ii++)
    {
      aa = forAssyVec[ii];
      if(aa == -1) // this DOF has a prescibed value
      {
        fact = SolnData->dispApplied[globalDOFnums[ii]];

        for(jj=0; jj<nsize; jj++)
        {
          if( forAssyVec[jj] != -1 )
          {
            Flocal1(jj) -= Kuu(jj, ii) * fact;
          }
        }

        for(jj=0; jj<sizeP; jj++)
        {
          if( forAssyVecPres[jj] != -1 )
          {
            Flocal2(jj) -= Kpu(jj, ii) * fact;
          }
        }
      }
    }

    // applied pressure
    for(ii=0; ii<sizeP; ii++)
    {
      aa = forAssyVecPres[ii];
      if(aa == -1) // this DOF has a prescibed value
      {
        fact = SolnData->presApplied[nodeNums[ii]];

        for(jj=0; jj<nsize; jj++)
        {
          if( forAssyVec[jj] != -1 )
          {
            Flocal1(jj) -= Kup(jj, ii) * fact;
          }
        }

        for(jj=0; jj<sizeP; jj++)
        {
          if( forAssyVecPres[jj] != -1 )
          {
            Flocal2(jj) -= Kpp(jj, ii) * fact;
          }
        }
      }
    }
  }

  return 0;
}



int  ElementBase::applyDirichletBCs2field(int schemeType, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, VectorXd& Flocal1, VectorXd& Flocal2, vector<int>& forAssyVar1, vector<int>& forAssyVar2, VectorXd& dispApplied, VectorXd& presApplied)
{
  // add contributions to the rhs vector
  // from nodes with specified displacement BCs

  int aa, bb, ii, jj, size1 = forAssyVar1.size(), size2 = forAssyVar2.size();
  double fact=0.0;

  //printVector(globalDOFnums);

    // applied displacement dof for the solid
    for(ii=0; ii<size1; ii++)
    {
      aa = forAssyVar1[ii];
      if(aa == -1) // this DOF has a prescibed value
      {
        fact = dispApplied[globalDOFnums[ii]];

        for(jj=0; jj<size1; jj++)
        {
          if( forAssyVar1[jj] != -1 )
          {
            Flocal1(jj) -= K11(jj, ii) * fact;
          }
        }

        for(jj=0; jj<size2; jj++)
        {
          if( forAssyVar2[jj] != -1 )
          {
            Flocal2(jj) -= K21(jj, ii) * fact;
          }
        }
      }
    }

    // applied pressure
    for(ii=0; ii<size2; ii++)
    {
      aa = forAssyVar2[ii];
      if(aa == -1) // this DOF has a prescribed value
      {
        fact = presApplied[nodeNums[ii]];

        //cout <<  ii <<  '\t' <<  fact <<  endl;

        for(jj=0; jj<size1; jj++)
        {
          if( forAssyVar1[jj] != -1 )
          {
            Flocal1(jj) -= K12(jj, ii) * fact;
          }
        }

        for(jj=0; jj<size2; jj++)
        {
          if( forAssyVar2[jj] != -1 )
          {
            Flocal2(jj) -= K22(jj, ii) * fact;
          }
        }
      }
    }

    return 0;
}




int  ElementBase::applyDirichletBCs3field(int schemeType, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, MatrixXd& K13, MatrixXd& K31, MatrixXd& K33, VectorXd& Flocal1, VectorXd& Flocal2, VectorXd& Flocal3, vector<int>& forAssyVar1, vector<int>& forAssyVar2, vector<int>& forAssyVar3, VectorXd& dispApplied, VectorXd& presApplied, VectorXd& mpotApplied)
{
  // add contributions to the rhs vector
  // from nodes with specified displacement BCs

  int aa, bb, ii, jj;
  int size1 = forAssyVar1.size(), size2 = forAssyVar2.size(), size3 = forAssyVar3.size();
  double fact=0.0;

  //printVector(globalDOFnums);
  //printVector(nodeNums);
  //printVector(forAssyVecEpot);

    // applied displacement dof for the solid
    for(ii=0; ii<size1; ii++)
    {
      aa = forAssyVar1[ii];
      if(aa == -1) // this DOF has a prescibed value
      {
        fact = dispApplied[globalDOFnums[ii]];

        for(jj=0; jj<size1; jj++)
        {
          if( forAssyVar1[jj] != -1 )
          {
            Flocal1(jj) -= K11(jj, ii) * fact;
          }
        }

        for(jj=0; jj<size2; jj++)
        {
          if( forAssyVar2[jj] != -1 )
          {
            Flocal2(jj) -= K21(jj, ii) * fact;
          }
        }

        for(jj=0; jj<size3; jj++)
        {
          if( forAssyVar3[jj] != -1 )
          {
            Flocal3(jj) -= K31(jj, ii) * fact;
          }
        }
      }
    }

    //printVector(presApplied);

    // applied potential
    for(ii=0; ii<size2; ii++)
    {
      aa = forAssyVar2[ii];
      if(aa == -1) // this DOF has a prescribed value
      {
        fact = presApplied[nodeNums[ii]];

        //cout <<  ii <<  '\t' << nodeNums[ii] <<  '\t' <<  fact <<  endl;

        for(jj=0; jj<size1; jj++)
        {
          if( forAssyVar1[jj] != -1 )
          {
            Flocal1(jj) -= K12(jj, ii) * fact;
          }
        }

        for(jj=0; jj<size2; jj++)
        {
          if( forAssyVar2[jj] != -1 )
          {
            Flocal2(jj) -= K22(jj, ii) * fact;

            //cout << ii << '\t' <<  jj << '\t' << fact << '\t' << K22(jj, ii) <<  '\t' <<  Flocal2(jj) <<  endl;
          }
        }
        //cout << endl;cout << endl;cout << endl;
      }
    }

    // applied pressure
    for(ii=0; ii<size3; ii++)
    {
      aa = forAssyVar3[ii];
      if(aa == -1) // this DOF has a prescribed value
      {
        fact = mpotApplied[nodeNums[ii]];

        //cout <<  ii <<  '\t' <<  fact <<  endl;

        for(jj=0; jj<size1; jj++)
        {
          if( forAssyVar1[jj] != -1 )
          {
            Flocal1(jj) -= K13(jj, ii) * fact;
          }
        }

        for(jj=0; jj<size3; jj++)
        {
          if( forAssyVar3[jj] != -1 )
          {
            Flocal3(jj) -= K33(jj, ii) * fact;
          }
        }
      }
    }

    return 0;
}




int  ElementBase::applyDirichletBCsElecField(MatrixXd& Kff, VectorXd& FlocalF)
{
  // add contributions to the rhs vector
  // from nodes with specified electric potential

  int aa, bb, ii, jj, sizeF = forAssyVecEpot.size();
  double fact=0.0;

  for(ii=0; ii<sizeF; ii++)
  {
    aa = forAssyVecEpot[ii];
    if(aa == -1) // this DOF has a prescribed value
    {
      fact = SolnData->mpotApplied[nodeNums[ii]];

      for(jj=0; jj<sizeF; jj++)
      {
        if( forAssyVecEpot[jj] != -1 )
        {
          FlocalF(jj) -= Kff(jj, ii) * fact;
        }
      }
    }
  }

  return 0;
}



int  ElementBase::calcError(int index, double* val)
{
  if(ndim == 2)
    calcError2D(index, val);
  else
    calcError3D(index, val);

  return 0;
}



int  ElementBase::calcError2D(int index, double* errordata)
{
  int   ii, jj, kk, gp;

  double  param[2], Jac, dvol0, dvol, rad, theta, val, fact;
  double  dispExact[2], dispNum[2], streDev[4], pNum, pExact;
  double  detF, streExact[4], streNum[4], Np[nlbfP], error, numer, denom;

  VectorXd  N(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU), stre(9);
  MatrixXd  F(3,3);

  double xNode[nlbfU], yNode[nlbfU], geom[3];
  for(ii=0;ii<nlbfU;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
  }

  double thickness = ElemTypeData->getThickness();

  vector<double>  matDat = MatlData->matData;

  //ThickCylinder  analy(matDat[0], matDat[1]);
  ThickCylinderComposite  analy(matDat[0], matDat[1]);
  //ElasticityFiniteStrain  analy(matDat[1], matDat[2]);
  //ThickCylinder  analy(sss, matDat[1], matDat[2]);
  double dt   = myTime.dt;
  double tCur = myTime.cur;

  //cout << "nlbfU = " << nlbfU << endl;

  vector<double>  gausspoints1, gausspoints2, gaussweights;
  nGP = getGaussPoints2D(npElem, gausspoints1, gausspoints2, gaussweights);

  errordata[0] = 0.0;
  errordata[1] = 0.0;

  if(index == 0) // L2 norm
  {
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          //cout << param[0] << '\t' << param[1] << endl;

          GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dNu_dx(0), &dNu_dy(0), Jac);

          dvol = gaussweights[gp] * Jac;

          geom[0] = geom[1] = 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            //cout << ii << '\t' << xNode[ii] << '\t' << yNode[ii] << '\t' << N[ii] << endl;
            geom[0] += xNode[ii]*N[ii];
            geom[1] += yNode[ii]*N[ii];
          }

          rad = sqrt(geom[0]*geom[0]+geom[1]*geom[1]);
          theta = atan2(geom[1],geom[0]);

          dispNum[0] = computeDisplacement(0, N);
          dispNum[1] = computeDisplacement(1, N);


          dispExact[0] = analy.dispX(rad, theta);
          dispExact[1] = analy.dispY(rad, theta);

          //cout << rad << '\t' << theta << endl;
          //cout << dispExact[0] << '\t' << dispNum[0] << endl;
          //cout << dispExact[1] << '\t' << dispNum[1] << endl; cout << endl;

          numer = 0.0;
          denom = 0.0;
          for(ii=0; ii<ndim; ii++)
          {
            error = dispExact[ii] - dispNum[ii];
            numer += error*error;

            denom += (dispExact[ii]*dispExact[ii]);
          }

          // absolute error
          fact = numer;

          // relative error
          fact = numer/denom;

          errordata[0] += ( numer * dvol );
          errordata[1] += ( denom * dvol );
    }//gp
  }
  else if(index == 1) // pressure
  {
    for(gp=0;gp<nGP;gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dNu_dx(0), &dNu_dy(0), Jac);

        dvol = gaussweights[gp]*(Jac*thickness);

        geom[0] = geom[1] = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += xNode[ii]*N[ii];
          geom[1] += yNode[ii]*N[ii];
        }
        rad = sqrt(geom[0]*geom[0]+geom[1]*geom[1]);
        theta = atan2(geom[1],geom[0]);

        // evaluate pressure at the quadrature points

        if(nlbfP == 1)
        {
          Np[0] = 1.0;
          pres = presDOF[0];
        }
        else
        {
          if(nlbfP == 3)
            LagrangeBasisFunsTria(nlbfP, param[0], param[1], Np);
          else if(nlbfP == 4)
            LagrangeBasisFunsQuad(nlbfP, param[0], param[1], Np);

          pNum = 0.0;
          for(ii=0; ii<nlbfP; ii++)
            pNum += (Np[ii]*SolnData->pres[nodeNumsPres[ii]]);
        }

        pExact = analy.pressure(rad, theta);

        error = pExact-pres;

        numer = error*error;
        denom = pExact*pExact;

        // absolute error
        fact = numer;

        // relative error
        fact = numer/denom;

        errordata[0] += ( numer * dvol );
        errordata[1] += ( denom * dvol );
    }//gp
  }
  else if(index == 2)
  {
    for(gp=0;gp<nGP;gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dNu_dx(0), &dNu_dy(0), Jac);

        dvol = gaussweights[gp]*(Jac*thickness);

        geom[0] = geom[1] = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += xNode[ii]*N[ii];
          geom[1] += yNode[ii]*N[ii];
        }
        rad = sqrt(geom[0]*geom[0]+geom[1]*geom[1]);
        theta = atan2(geom[1],geom[0]);

        // evaluate pressure at the quadrature points

        if(nlbfP == 1)
        {
          Np[0] = 1.0;
          pres = presDOF[0];
        }
        else
        {
          if(nlbfP == 3)
            LagrangeBasisFunsTria(nlbfP, param[0], param[1], Np);
          else if(nlbfP == 4)
            LagrangeBasisFunsQuad(nlbfP, param[0], param[1], Np);

          pNum = 0.0;
          for(ii=0; ii<nlbfP; ii++)
            pNum += (Np[ii]*SolnData->pres[nodeNumsPres[ii]]);
        }

        computeDefGrad2DCur(dNu_dx, dNu_dy, F);

          // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
        adjust_deformationgradient_2D(sss, finite, F);

        MatlData->computeMechanicalStress(F, pres, stre);

          streNum[0] = stre[0]; // xx
          streNum[1] = stre[4]; // yy
          streNum[2] = stre[1]; // xy
          streNum[3] = stre[8]; // zz

          analy.stresses1(rad, theta, streExact);
          //streExact[0] = analy.stressXX(rad, theta);
          //streExact[1] = analy.stressYY(rad, theta);
          //streExact[2] = analy.stressXY(rad, theta);

          numer = 0.0;
          denom = 0.0;
          for(ii=0; ii<4; ii++)
          {
            //cout << ii << '\t' << streExact[ii] << '\t' << streNum[ii] << endl;
            error = streExact[ii] - streNum[ii];
            numer += error*error;

            denom += (streExact[ii]*streExact[ii]);
          }

          // absolute error
          fact = numer;

          // relative error
          fact = numer/denom;

          errordata[0] += ( numer * dvol );
          errordata[1] += ( denom * dvol );
    }//gp
  }


  return 0;
}



int  ElementBase::calcError3D(int index, double* errordata)
{
  int   err,  isw,  count,  count1, ll, ii, jj, kk, gp;

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

  double  param[3], Jac, dvol0, dvol, rad, theta, phi, val, fact;
  double  dispExact[3], dispNum[3], pNum, pExact;
  double  detF, streExact[6], streNum[6], Np[nlbfP], error, numer, denom;

  double xNode[nlbfU], yNode[nlbfU], zNode[nlbfU], geom[3];
  for(ii=0;ii<nlbfU;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
  }
  //cout << " index = " << index << endl;

  vector<double>  matDat = MatlData->matData;

  //ThickSphere  analy(matDat[0], matDat[1]);
  ThickSphereComposite  analy(matDat[0], matDat[1]);
  //ElasticityFiniteStrain  analy(matDat[1], matDat[2]);
  //ThickCylinder  analy(sss, matDat[1], matDat[2]);
  double dt   = myTime.dt;
  double tCur = myTime.cur;


  vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
  nGP = getGaussPoints3D(npElem, gausspoints1, gausspoints2, gausspoints3, gaussweights);

  errordata[0] = 0.0;
  errordata[1] = 0.0;

  if(index == 0) // L2 norm
  {
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];
          param[2] = gausspoints3[gp];

          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

          dvol = gaussweights[gp] * Jac;

          geom[0] = geom[1] = geom[2] = 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            geom[0] += xNode[ii]*N[ii];
            geom[1] += yNode[ii]*N[ii];
            geom[2] += zNode[ii]*N[ii];
          }

          //rad = sqrt(geom[0]*geom[0]+geom[1]*geom[1]+geom[2]*geom[2]);
          //theta = atan2(geom[1],geom[0]);

          dispNum[0] = computeDisplacement(0, N);
          dispNum[1] = computeDisplacement(1, N);
          dispNum[2] = computeDisplacement(2, N);


          //dispExact[0] = analy.dispX(rad, theta);
          //dispExact[1] = analy.dispY(rad, theta);
          //dispExact[2] = analy.dispZ(rad, theta);
          analy.disp(geom[0], geom[1], geom[2], dispExact);

          //cout << dispExact[0] << '\t' << dispNum[0] << endl;
          //cout << dispExact[1] << '\t' << dispNum[1] << endl;
          //cout << dispExact[2] << '\t' << dispNum[2] << endl;

          numer = 0.0;
          denom = 0.0;
          for(ii=0; ii<ndim; ii++)
          {
            error = dispExact[ii] - dispNum[ii];
            numer += error*error;

            denom += (dispExact[ii]*dispExact[ii]);
          }

          // absolute error
          fact = numer;

          // relative error
          fact = numer/denom;

          errordata[0] += ( numer * dvol );
          errordata[1] += ( denom * dvol );
    }//gp
  }
  else if(index == 1)
  {
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];
          param[2] = gausspoints3[gp];

          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

          dvol = gaussweights[gp] * Jac;

          geom[0] = geom[1] = geom[2] = 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            geom[0] += xNode[ii]*N[ii];
            geom[1] += yNode[ii]*N[ii];
            geom[2] += zNode[ii]*N[ii];
          }

          rad   = sqrt(geom[0]*geom[0] + geom[1]*geom[1] + geom[2]*geom[2]);
          theta = atan(geom[1]/geom[0]);
          phi   = acos(geom[2]/rad);

          // evaluate pressure at the quadrature points

          if(nlbfP == 1)
          {
            Np[0] = 1.0;
            pNum = presDOF[0];
          }
          else
          {
            if(nlbfP == 4)
              LagrangeBasisFunsTetra(nlbfP, param[0], param[1], param[2], Np);
            else if(nlbfP == 8)
              LagrangeBasisFunsHexa(nlbfP, param[0], param[1], param[2], Np);

            pNum = 0.0;
            for(ii=0; ii<nlbfP; ii++)
              pNum += (Np[ii]*SolnData->pres[nodeNumsPres[ii]]);
          }

        pExact = analy.pressure(rad, theta, phi);

        //cout << "pExact =    " << pExact << '\t' << pres << endl;

        error = pExact-pres;

        numer = error*error;
        denom = pExact*pExact;

        // absolute error
        fact = numer;

        // relative error
        fact = numer/denom;

        errordata[0] += ( numer * dvol );
        errordata[1] += ( denom * dvol );
    }//gp
  }
  else if(index == 2)
  {
    MatrixXd  F(3,3);
    VectorXd  stre(9);

    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];
          param[2] = gausspoints3[gp];

          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

          dvol = gaussweights[gp] * Jac;

          geom[0] = geom[1] = geom[2] = 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            geom[0] += xNode[ii]*N[ii];
            geom[1] += yNode[ii]*N[ii];
            geom[2] += zNode[ii]*N[ii];
          }

          rad   = sqrt(geom[0]*geom[0] + geom[1]*geom[1] + geom[2]*geom[2]);
          theta = atan(geom[1]/geom[0]);
          phi   = acos(geom[2]/rad);

          computeDefGradCur(dN_dx, dN_dy, dN_dz, F);

          // evaluate pressure at the quadrature points

          if(nlbfP == 1)
          {
            Np[0] = 1.0;
            pNum = presDOF[0];
          }
          else
          {
            if(nlbfP == 4)
              LagrangeBasisFunsTetra(nlbfP, param[0], param[1], param[2], Np);
            else if(nlbfP == 8)
              LagrangeBasisFunsHexa(nlbfP, param[0], param[1], param[2], Np);

            pNum = 0.0;
            for(ii=0; ii<nlbfP; ii++)
              pNum += (Np[ii]*SolnData->pres[nodeNumsPres[ii]]);
          }


          MatlData->computeMechanicalStress(F, pNum, stre);

          analy.stresses(geom[0], geom[1], geom[2], streExact);

          streNum[0] = stre[0]; //xx
          streNum[1] = stre[4]; //yy
          streNum[2] = stre[8]; //zz
          streNum[3] = stre[1]; //xy
          streNum[4] = stre[5]; //yz
          streNum[5] = stre[2]; //xz

          numer = 0.0;
          denom = 0.0;
          for(ii=0; ii<6; ii++)
          {
            //cout << ii << '\t' << streExact[ii] << '\t' << streNum[ii] << endl;
            error = streExact[ii] - streNum[ii];
            numer += error*error;

            denom += (streExact[ii]*streExact[ii]);
          }

          // absolute error
          fact = numer;

          // relative error
          fact = numer/denom;

          errordata[0] += ( numer * dvol );
          errordata[1] += ( denom * dvol );
    }//gp
  }

  return 0;
}



void ElementBase::diffStiffTest(double ddd, int dig, int dig2, bool gfrmt)
{
/*
    //  function for u-p-J formulation

    int    i, ii, jnd, jj, j, k, ind, bb, nn;
    int  size1 = 2*nlbfU+nlbfP+nlbfP;

    cout << nlbfU << '\t' << nlbfP << '\t' << size1 <<  endl;

    ddd = 1.0e-6;
    double dd[6]  = {-3.*ddd, -2.*ddd, -ddd, ddd, 2.*ddd, 3.*ddd };
    double  r[6*nsize], r2[6*nlbfP], r3[6*nlbfP];

    MatrixXd  Kfull1(size1, size1), Kfull2(size1, size1), Kdiff(size1, size1);
    MatrixXd  Kuu, Kuf, Kfu, Kff;
    VectorXd  Flocal, Flocal2, Flocal3;

    Kfull1.setZero();
    Kfull2.setZero();
    Kdiff.setZero();

  for(ii=0; ii<size1; ii++)
  {
    cout << '\t' ;
    for(jj=0; jj<size1; jj++)
    {
      cout  <<  fixed << Kfull1(ii,jj) << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;

  for(ii=0; ii<size1; ii++)
  {
    cout << '\t' ;
    for(jj=0; jj<size1; jj++)
    {
      cout  <<  fixed << Kfull2(ii,jj) << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;

    cout << " Diff Test for element " << elenum << '\t' << ddd << endl;

//
    // loop over columns of s
    for (jnd=0; jnd<nlbfU; jnd++)                          // nodes
    //for (jnd=0; jnd<1; jnd++)                              // nodes
    {
      for (jj=0; jj<ndof; jj++)  // dof
      {
        j = jnd*ndof + jj;

        nn  = nodeNums[jnd];
        ind = nn*2+jj;

        // loop over perturbations

        for (k=0; k<6; k++)
        {
          cout << jnd << '\t' << jj << '\t' <<  k << '\t' <<  dd[k] << endl;
          // apply pertubation
          SolnData->disp[ind] += dd[k];

          SolnData->dispCur[ind]  =  SolnData->disp[ind];

          GeomData->NodePosNew[nn][jj] = GeomData->NodePosOrig[nn][jj] + SolnData->disp[ind];
          GeomData->NodePosCur[nn][jj] = GeomData->NodePosOrig[nn][jj] + SolnData->dispCur[ind];


          // calculate residual
          //calcStiffnessAndResidual(Kuu, Kuf, Kfu, Kff, Flocal, Flocal2, true);
          //printMatrix(Kpu);
          //printVector(Flocal);
          //printVector(Flocal);
          //printVector(Flocal2);

          // remove pertubation
          SolnData->disp[ind] -= dd[k];

          SolnData->dispCur[ind]  =  SolnData->disp[ind];

          GeomData->NodePosNew[nn][jj] = GeomData->NodePosOrig[nn][jj] + SolnData->disp[ind];
          GeomData->NodePosCur[nn][jj] = GeomData->NodePosOrig[nn][jj] + SolnData->dispCur[ind];


          // loop over rows of s, store residual

          for(i=0; i<nsize; i++)
            r[i*6+k] = Flocal[i];

          r2[k] = -Rt;

          r3[k] = -Rp(0);
        }
        cout <<  " kkkkkkkkkkkkkkk " <<  endl;

        // loop over rows of s
        for (i=0; i<nsize; i++)
        Kfull2(j,i) = (    +       r[i*6+0]
                           -  9. * r[i*6+1]
                           + 45. * r[i*6+2]
                           - 45. * r[i*6+3]
                           +  9. * r[i*6+4]
                           -       r[i*6+5] ) / (60. * ddd);
        cout <<  " kkkkkkkkkkkkkkk " <<  endl;

        Kfull2(j, nsize) = (   +       r2[0]
                                 -  9. * r2[1]
                                 + 45. * r2[2]
                                 - 45. * r2[3]
                                 +  9. * r2[4]
                                 -       r2[5] ) / (60. * ddd);
        cout <<  " kkkkkkkkkkkkkkk " <<  endl;

        Kfull2(j, nsize+nlbfP) = (     +       r3[0]
                                       -  9. * r3[1]
                                       + 45. * r3[2]
                                       - 45. * r3[3]
                                       +  9. * r3[4]
                                       -       r3[5] ) / (60. * ddd);
        cout <<  " kkkkkkkkkkkkkkk " <<  endl;

    }
  }
*/
  cout <<  " mmmmmmmmmmmmm " <<  endl;

/*
  // calculate stiffness
  //calcStiffnessAndResidual(Kuu, Kuf, Kfu, Kff, Flocal, Flocal2, true);
  printMatrix(Kuu);
  printMatrix(Kut);
  printMatrix(Ktu);
  printMatrix(Kup);
  printMatrix(Kpu);
  cout <<  " mmmmmmmmmmmmm " <<  endl;

  // Kuu, Kut, Kup
  for(ii=0; ii<nsize; ii++)
  {
    cout << "ii= " <<  ii << endl;
    for(jj=0;jj<nsize;jj++)
    {
      Kfull1(ii,jj) = Kuu(ii,jj);
    }

    Kfull1(ii,nsize) = Kut(ii,0);

    Kfull1(ii,nsize+nlbfP) = Kup(ii,0);
  }
  cout <<  " ooooooooooooooooo " <<  endl;

  // Ktu, Ktt, Ktp
    for(jj=0;jj<nsize;jj++)
    {
      Kfull1(nsize,jj) = Ktu(0,jj);
    }
    Kfull1(nsize,nsize)   = Ktt;
    Kfull1(nsize,nsize+1) = Kpt;

  // Kpu, Kpt, Kpp
    for(jj=0;jj<nsize;jj++)
    {
      Kfull1(nsize+nlbfP,jj) = Kpu(0,jj);
    }
    Kfull1(nsize+nlbfP,nsize) = Kpt;

  cout.setf(ios::fixed);
  cout.setf(ios::showpoint);
  cout.precision(8);

  cout << "    Analytical Stiffness Matrix   " << endl;
  cout << endl;
  for(ii=0; ii<size1; ii++)
  {
    cout << '\t' ;
    for(jj=0; jj<size1; jj++)
    {
      cout  <<  fixed << Kfull1(ii,jj) << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;


  cout << "    diff Stiffness Matrix   " << endl;
  cout << endl;
  cout << endl;
  for(ii=0; ii<size1; ii++)
  {
    cout << '\t' ;
    for(jj=0; jj<size1; jj++)
    {
      Kdiff(ii,jj) = Kfull1(ii,jj) - Kfull2(ii,jj);
      cout  <<  fixed <<  Kfull2(ii,jj) << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;

  cout << "    difference of Stiffness Matrices   " << endl;
  cout << endl;
  cout << endl;

  for(ii=0; ii<size1; ii++)
  {
    cout << '\t' ;
    for(jj=0; jj<size1; jj++)
    {
      cout  <<  fixed << Kdiff(ii, jj) << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;
*/
  return;
}






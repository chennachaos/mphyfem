
#include "Elem_Solid_2D_Tria10_Mixed1.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"


using namespace std;

extern MyTime myTime;
extern vector<unique_ptr<TimeFunction> > timeFunctions;



Elem_Solid_2D_Tria10_Mixed1::Elem_Solid_2D_Tria10_Mixed1()
{
  ndim = 2;
  ndof = 2;
}



Elem_Solid_2D_Tria10_Mixed1::~Elem_Solid_2D_Tria10_Mixed1()
{
}



void Elem_Solid_2D_Tria10_Mixed1::prepareElemData()
{
  npElem = nodeNums.size();
  nlbfU  = npElem;

  nsize = nlbfU*ndof;


  ELEM_SHAPE = ELEM_SHAPE_TRIA;
  nGP = 7;
  nlbfP = 3;

  nodeNumsPres.resize(nlbfP);
  int ind=nlbfP*elenum;
  for(int ii=0; ii<nlbfP; ii++)
    nodeNumsPres[ii] = ind+ii;

  forAssyVecPres = nodeNumsPres;


  ElementBase::prepareElemData();

  return;
}




int Elem_Solid_2D_Tria10_Mixed1::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
  return 0;
}




int Elem_Solid_2D_Tria10_Mixed1::calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalP)
{
  ElementBase::calcStiffnessAndResidualMixed2D(Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP);

  return 0;
}




void Elem_Solid_2D_Tria10_Mixed1::elementContourplot(int vartype, int varindex, int index)
{
   double outval[50];

   switch(vartype)
   {
       case 1:  // plot total strain
       case 2:  // plot elastic strain
       case 3:  // plot plastic strain

                projectStrain(true, vartype, varindex, index, outval);

              break;

       case 4:  // plot stress

                //intVar2 = intVar1;
                projectStress(true, vartype, varindex, index, outval);

              break;

       case 5:  // plot element internal variables

                projectInternalVariable(true, vartype, varindex, index, outval);

              break;

       default:

              cout  << " Invalid Variable Type to project in 'Elem_Solid_2D_Tria10_Mixed1::projectToNodes'" << endl;
              break;
    }

    assert(vals2project.size() >= 1);

    vals2project[0] = 0.0;
    for(int ii=0; ii<nGP; ii++)
    {
      vals2project[0] += outval[ii];
    }

    vals2project[0] /= nGP;

    return;
}


void Elem_Solid_2D_Tria10_Mixed1::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
{
   double outval[50];

   switch(vartype)
   {
       case 1:  // plot total strain
       case 2:  // plot elastic strain
       case 3:  // plot plastic strain

                projectStrain(extrapolateFlag, vartype, varindex, index, outval);

              break;

       case 4:  // plot stress

                //intVar2 = intVar1;
                projectStress(extrapolateFlag, vartype, varindex, index, outval);

              break;

       case 5:  // plot element internal variables

                projectInternalVariable(extrapolateFlag, vartype, varindex, index, outval);

              break;

       default:

              cout  << "           Invalid Variable Type to project " << endl;
              break;
    }


    assert(vals2project.size() == npElem);

    // project from quadrature points to nodes/control points
    if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
        stressrecovery_extrapolate_Quadrilateral(degree, outval, &vals2project[0]);
    else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
        stressrecovery_extrapolate_Triangle(degree, outval, &vals2project[0]);

    //printVector(vals2project);

    return;
}


void Elem_Solid_2D_Tria10_Mixed1::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    int  err = 0, ii, jj, gp;

    double  detF, dvol, Jac, fact, dvol0, param[2], Np[nlbfP];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), stre(9);

    //double  timeFact = elmDat[10]*timeFunctions[0].getValue();
    MatrixXd  F(3,3), Fn(3,3);
    F.setZero();
    Fn.setZero();
    double dt   = myTime.dt;

    double xNode[9], yNode[9], geom[3];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
        getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);
    else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
        getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);


    for(int gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        geom[0] = geom[1] = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += N[ii]*xNode[ii];
          geom[1] += N[ii]*yNode[ii];
        }

        computeDefGrad2D(dN_dx, dN_dy, F);

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        adjust_deformationgradient_2D(sss, finite, F);


        // basis functions for pressure dof
        if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
            BernsteinBasisFunsQuad(1, param[0], param[1], Np);
        else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
            BernsteinBasisFunsTria(1, param[0], param[1], Np);

        // evaluate pressure at the quadrature points
        pres = 0.0;
        for(ii=0; ii<nlbfP; ii++)
          pres += (Np[ii]*SolnData->pres[nodeNums[ii]]);
          //pres += (Np[ii]*SolnData->var2[forAssyVecPres[ii]]);


        MatlData->computeMechanicalStress(F, pres, stre);

        if(varindex < 9)
           outval[gp] = stre[varindex];
        else if(varindex == 9)
           outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6* stre[1]*stre[1])/2);
        else if(varindex == 10)
          outval[gp] = pres;

        //outval[gp] = detF;

    }//gp

    return;
}



void Elem_Solid_2D_Tria10_Mixed1::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void Elem_Solid_2D_Tria10_Mixed1::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.size() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       //outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}








#include "Elem_Solid_2D_Disp.h"
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



Elem_Solid_2D_Disp::Elem_Solid_2D_Disp()
{
  ndim = 2;
  ndof = 2;
}



Elem_Solid_2D_Disp::~Elem_Solid_2D_Disp()
{
}



void Elem_Solid_2D_Disp::prepareElemData()
{
  npElem = nodeNums.size();
  nlbfU  = npElem;

  nsize = nlbfU*ndof;


  if( (nlbfU == 3) )
  {
    ELEM_SHAPE = ELEM_SHAPE_TRIA;
    nGP = 1;
  }
  else if( (nlbfU == 6) )
  {
    ELEM_SHAPE = ELEM_SHAPE_TRIA;
    nGP = 3;
  }
  else if( (nlbfU == 10) )
  {
    ELEM_SHAPE = ELEM_SHAPE_TRIA;
    nGP = 7;
  }
  else if( (nlbfU == 4) )
  {
    ELEM_SHAPE = ELEM_SHAPE_QUAD;
    nGP = 4;
  }
  else if( (nlbfU == 9) )
  {
    ELEM_SHAPE = ELEM_SHAPE_QUAD;
    nGP = 9;
  }
  else if( (nlbfU == 16) )
  {
    ELEM_SHAPE = ELEM_SHAPE_QUAD;
    nGP = 16;
  }
  else
  {
     throw runtime_error("Elem_Solid_2D_Disp::prepareElemData() ... invalid element type");
  }


  ElementBase::prepareElemData();

  return;
}




int Elem_Solid_2D_Disp::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
  return 0;
}




int Elem_Solid_2D_Disp::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
    int   err, index, ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  detF, fact, dvol, dvol0, Jac, volstrain;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[2], bforce[2]={0.0,0.0}, force[2];
    double  veloCur[2], acceCur[2], sig[2];


    int    sss   = ElemTypeData->getModeltypeNum();
    double thickness = ElemTypeData->getThickness();
    bool   axsy  = (sss == 3);
    bool  finite = MatlData->isFiniteStrain();

    double rho0  = MatlData->getDensity() ;
    double rho   = rho0;

    //bforce[0]   = elmDat[6]*timeFunction[0]->getFactor() ;
    //bforce[1]   = elmDat[7]*timeFunction[0]->getFactor() ;

    double af   = SolnData->td(2);
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;


    double xNode[nlbfU], yNode[nlbfU], xx, yy;
    for(ii=0;ii<nlbfU;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    nGP = getGaussPoints2D(npElem, gausspoints1, gausspoints2, gaussweights);


    MatrixXd  Cmat(9,9), Gc(3,9), F(3,3);
    Cmat.setZero();
    Gc.setZero();
    F.setZero();
    VectorXd  stre(9);
    stre.setZero();


    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(thickness*Jac);

        xx = yy = 0.0;
        for(ii=0;ii<npElem;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
        }

        if(axsy)
          dvol0 *= 2.0*PI*xx;

        dvol  = dvol0;

        computeDefGrad2DCur(dN_dx, dN_dy, F);

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        volstrain = F(0,0) + F(1,1) - 2.0;

        acceCur[0] = computeAccelerationCur(0, N);
        acceCur[1] = computeAccelerationCur(1, N);

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
        adjust_deformationgradient_2D(sss, finite, F);

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(thickness*Jac);
        }

        if(finite)  dvol *= F(2,2);

        MatlData->computeStressAndTangent(true, sss, F, F, pres, stre, Cmat, ivar, gp, myTime.dt);
        //if(err !=0)    return 1;
        //printMatrix(Cmat);


        // Calculate Stiffness and Residual
        //==============================================

        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        //   part 1. -- material part (not necessarily symmetric!!)
        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          Gc(0,0) = bb1 * Cmat(0,0) + bb2 * Cmat(3,0);
          Gc(0,1) = bb1 * Cmat(0,1) + bb2 * Cmat(3,1);
          Gc(0,3) = bb1 * Cmat(0,3) + bb2 * Cmat(3,3);
          Gc(0,4) = bb1 * Cmat(0,4) + bb2 * Cmat(3,4);

          Gc(1,0) = bb1 * Cmat(1,0) + bb2 * Cmat(4,0);
          Gc(1,1) = bb1 * Cmat(1,1) + bb2 * Cmat(4,1);
          Gc(1,3) = bb1 * Cmat(1,3) + bb2 * Cmat(4,3);
          Gc(1,4) = bb1 * Cmat(1,4) + bb2 * Cmat(4,4);

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[1] + bb2*stre[4];

          Flocal(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal(TIp1) += (bb5*force[1] - sig[1]) ;

          for(jj=0; jj<nlbfU; jj++)
          {
              cc1 = dN_dx[jj];
              cc2 = dN_dy[jj];
              cc5 = N[jj];

              TJ   = 2*jj;
              TJp1 = TJ+1;

              // acceleration term
              fact  = bb5*acceFact1*cc5*rho0;

              // material Stiffness contribution
              //fact  += af*(sig[0]*cc1+sig[1]*cc2)*FiniteFact;

              Klocal(TI,   TJ)    += fact ;
              Klocal(TIp1, TJp1)  += fact ;

              Klocal(TI,   TJ)    +=  af*(Gc(0,0) * cc1 + Gc(0,3) * cc2) ;
              Klocal(TI,   TJp1)  +=  af*(Gc(0,1) * cc1 + Gc(0,4) * cc2) ;
              Klocal(TIp1, TJ)    +=  af*(Gc(1,0) * cc1 + Gc(1,3) * cc2) ;
              Klocal(TIp1, TJp1)  +=  af*(Gc(1,1) * cc1 + Gc(1,4) * cc2) ;
          }
        }
    } //gp

    //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

    return 0;
}




void Elem_Solid_2D_Disp::elementContourplot(int vartype, int varindex, int index)
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

              cout  << " Invalid Variable Type to project in 'Elem_Solid_2D_Disp::projectToNodes'" << endl;
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


void Elem_Solid_2D_Disp::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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


void Elem_Solid_2D_Disp::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
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



void Elem_Solid_2D_Disp::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void Elem_Solid_2D_Disp::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}







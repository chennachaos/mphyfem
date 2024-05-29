
#include "Elem_Solid_3D_Disp.h"
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



Elem_Solid_3D_Disp::Elem_Solid_3D_Disp()
{
  ndim = 3;
  ndof = 3;
}



Elem_Solid_3D_Disp::~Elem_Solid_3D_Disp()
{
}



void Elem_Solid_3D_Disp::prepareElemData()
{
  npElem = nodeNums.size();
  nlbfU  = npElem;

  nsize = nlbfU*ndof;


  if( (npElem == 4) )
  {
    ELEM_SHAPE = ELEM_SHAPE_TETRA;
    nGP = 1;
  }
  else if( (npElem == 10) )
  {
    ELEM_SHAPE = ELEM_SHAPE_TETRA;
    nGP = 4;
  }
  else if( (npElem == 11) )
  {
    ELEM_SHAPE = ELEM_SHAPE_TETRA;
    nGP = 11;
  }
  else if( (npElem == 20) )
  {
    ELEM_SHAPE = ELEM_SHAPE_TETRA;
    nGP = 11;
  }
  else if( npElem == 5 )
  {
    ELEM_SHAPE = ELEM_SHAPE_PYRAMID;
    nGP = 1;
  }
  else if( npElem == 13 )
  {
    ELEM_SHAPE = ELEM_SHAPE_PYRAMID;
    nGP = 4;
  }
  else if( (npElem == 6) )
  {
    ELEM_SHAPE = ELEM_SHAPE_WEDGE;
    nGP = 2;
  }
  else if( (npElem == 18) )
  {
    ELEM_SHAPE = ELEM_SHAPE_WEDGE;
    nGP = 9;
  }
  else if( (npElem == 8) )
  {
    ELEM_SHAPE = ELEM_SHAPE_HEXA;
    nGP = 8;
  }
  else if( (npElem == 27) )
  {
    ELEM_SHAPE = ELEM_SHAPE_HEXA;
    nGP = 27;
  }
  else
  {
     throw runtime_error("Elem_Solid_3D_Disp::prepareElemData() ... invalid element type");
  }


  ElementBase::prepareElemData();

  return;
}





int Elem_Solid_3D_Disp::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != nsize) || (Mlocal.cols() != nsize) )
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = 0.0; //elmDat[5]*computeVolume(true)/nlbfU;

      for(int ii=0; ii<nsize; ii++)
      {
        Mlocal(ii,ii) = fact;
      }
    }
    else
    {
      int  nGPt, ii, jj, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;
      double  fact, dvol0, Jac, bb1, cc1, param[3];

      double rho0 = MatlData->getDensity();

      VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

      vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
      if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
      {
        nGPt=10;
        getGaussPointsTetra(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);
      }
      else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
      {
        nGPt=9;
        getGaussPointsWedge(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);
      }
      else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
      {
        nGPt=27;
        getGaussPointsHexa(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);
      }
      else
      {
        throw runtime_error("Element type not defined...");
      }

      for(gp=0; gp<nGPt; gp++)
      {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];
          param[2] = gausspoints3[gp];

          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

          dvol0 = gaussweights[gp]*Jac;

          for(ii=0;ii<nlbfU;ii++)
          {
            bb1 = (dvol0*rho0)*N[ii];

            TI   = 3*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;

            for(jj=0; jj<nlbfU; jj++)
            {
              TJ   = 3*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              fact  = bb1*N[jj];

              Mlocal(TI,   TJ)    += fact ;
              Mlocal(TIp1, TJp1)  += fact ;
              Mlocal(TIp2, TJp2)  += fact ;
            }
          }
      } //gp
      //cout << " elemVol = " << elemVol << endl;
      //printMatrix(Klocal);  printf("\n\n\n");
    }

    return 0;
}





int Elem_Solid_3D_Disp::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
    int   err, index, ii, jj, kk, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double  detF, fact, dvol, dvol0, Jac, volstrain;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[3], bforce[3]={0.0,0.0,0.0}, force[3];
    double  veloCur[3], acceCur[3], sig[3];


    bool  finite = MatlData->isFiniteStrain();

    double rho0  = MatlData->getDensity() ;
    double rho   = rho0;

    //bforce[0]   = elmDat[6]*timeFunction[0]->getFactor() ;
    //bforce[1]   = elmDat[7]*timeFunction[0]->getFactor() ;

    double af   = SolnData->td(2);
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;


    double xNode[npElem], yNode[npElem], zNode[npElem], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    nGP = getGaussPoints3D(npElem, gausspoints1, gausspoints2, gausspoints3, gaussweights);


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
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = gaussweights[gp]*Jac;

        xx = yy = zz = 0.0;
        for(ii=0;ii<npElem;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
          zz += N[ii]*zNode[ii];
        }

        dvol  = dvol0;

        computeDefGradCur(dN_dx, dN_dy, dN_dz, F);

        volstrain = F.trace() - 3.0;
        detF = F.determinant();

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        acceCur[0] = computeAccelerationCur(0, N);
        acceCur[1] = computeAccelerationCur(1, N);
        acceCur[2] = computeAccelerationCur(2, N);

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;
        }

        MatlData->computeStressAndTangent(true, sss, F, F, pres, stre, Cmat, ivar, gp, myTime.dt);
        //if(err !=0)    return 1;
        //printMatrix(Cmat);


        // Calculate Stiffness and Residual
        //==============================================

        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];
        force[2] = rho0*bforce[2];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];
        force[2] -= rho0*acceCur[2];

        //   part 1. -- material part (not necessarily symmetric!!)
        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*dN_dz[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          for(kk=0; kk<9; kk++)
          {
            Gc(0,kk) = bb1 * Cmat(0,kk) + bb2 * Cmat(3,kk) + bb3 * Cmat(6,kk);
            Gc(1,kk) = bb1 * Cmat(1,kk) + bb2 * Cmat(4,kk) + bb3 * Cmat(7,kk);
            Gc(2,kk) = bb1 * Cmat(2,kk) + bb2 * Cmat(5,kk) + bb3 * Cmat(8,kk);
          }

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          sig[0] = bb1*stre[0] + bb2*stre[3] + bb3*stre[6] ;
          sig[1] = bb1*stre[1] + bb2*stre[4] + bb3*stre[7] ;
          sig[2] = bb1*stre[2] + bb2*stre[5] + bb3*stre[8] ;

          Flocal(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal(TIp1) += (bb5*force[1] - sig[1]) ;
          Flocal(TIp2) += (bb5*force[2] - sig[2]) ;

          for(jj=0; jj<nlbfU; jj++)
          {
            cc1 = dN_dx[jj];
            cc2 = dN_dy[jj];
            cc3 = dN_dz[jj];
            cc5 = N[jj];

            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            // acceleration term
            fact  = bb5*acceFact1*cc5*rho0;

            // material Stiffness contribution
            //fact += af*(sig[0]*cc1 + sig[1]*cc2 + sig[2]*cc3)*FiniteFact;

            Klocal(TI,   TJ)    +=  fact;
            Klocal(TIp1, TJp1)  +=  fact;
            Klocal(TIp2, TJp2)  +=  fact;

            Klocal(TI,   TJ)    +=  af*(Gc(0,0) * cc1 + Gc(0,3) * cc2 + Gc(0,6) * cc3) ;
            Klocal(TI,   TJp1)  +=  af*(Gc(0,1) * cc1 + Gc(0,4) * cc2 + Gc(0,7) * cc3) ;
            Klocal(TI,   TJp2)  +=  af*(Gc(0,2) * cc1 + Gc(0,5) * cc2 + Gc(0,8) * cc3) ;

            Klocal(TIp1, TJ)    +=  af*(Gc(1,0) * cc1 + Gc(1,3) * cc2 + Gc(1,6) * cc3) ;
            Klocal(TIp1, TJp1)  +=  af*(Gc(1,1) * cc1 + Gc(1,4) * cc2 + Gc(1,7) * cc3) ;
            Klocal(TIp1, TJp2)  +=  af*(Gc(1,2) * cc1 + Gc(1,5) * cc2 + Gc(1,8) * cc3) ;

            Klocal(TIp2, TJ)    +=  af*(Gc(2,0) * cc1 + Gc(2,3) * cc2 + Gc(2,6) * cc3) ;
            Klocal(TIp2, TJp1)  +=  af*(Gc(2,1) * cc1 + Gc(2,4) * cc2 + Gc(2,7) * cc3) ;
            Klocal(TIp2, TJp2)  +=  af*(Gc(2,2) * cc1 + Gc(2,5) * cc2 + Gc(2,8) * cc3) ;
          }
        }
    } //gp

    //if(elenum == 0)
      //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

    return 0;
}




void Elem_Solid_3D_Disp::elementContourplot(int vartype, int varindex, int index)
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

              cout  << " Invalid Variable Type to project in 'Elem_Solid_3D_Disp::projectToNodes'" << endl;
              break;
    }

    assert(vals2project.size() >= 1);

    vals2project[0] = 0.0;
    for(int ii=0; ii<nGP; ii++)
      vals2project[0] += outval[ii];

    vals2project[0] /= nGP;

    return;
}


void Elem_Solid_3D_Disp::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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

                //intVar = intVarPrev;
                projectStress(extrapolateFlag, vartype, varindex, index, outval);

              break;

       case 5:  // plot element internal variables

                projectInternalVariable(extrapolateFlag, vartype, varindex, index, outval);

              break;

       default:

              cout  << "           Invalid Variable Type to project in Elem_Solid_3D_Disp::projectToNodes " << endl;
              break;
    }

    assert(vals2project.size() == npElem);

    // project from quadrature points to nodes/control points

    if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
      stressrecovery_extrapolate_Tetrahedron(degree, outval, &vals2project[0]);
    else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
      stressrecovery_extrapolate_Wedge(degree, outval, &vals2project[0]);
    else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
      stressrecovery_extrapolate_Hexahedron(degree, outval, &vals2project[0]);
    else
    {
        throw runtime_error("Element type not defined...");
    }

    //printVector(vals2project);

    return;
}




void Elem_Solid_3D_Disp::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    double  detF, fact, dvol, dvol0, Jac, param[3];
 
    MatrixXd  F(3,3), Fn(3,3), Cmat(9,9), stresMat(3,3);
    F.setZero();
    Fn.setZero();
    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU), Np(nlbfP), stre(9);

    int   err,  isw,  count,  count1, ll, ii, jj, gp;

    double dt = myTime.dt;

    double xNode[nlbfU], yNode[nlbfU], zNode[nlbfU], xx, yy, zz;
    for(ii=0;ii<nlbfU;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
      getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);
    else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
      getGaussPointsWedge(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);
    else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
      getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);
    else
    {
        throw runtime_error("Element type not defined...");
    }


    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        computeDefGradPrev(dN_dx, dN_dy, dN_dz, Fn);
        computeDefGrad(dN_dx, dN_dy, dN_dz, F);

        detF = F.determinant();

        // basis functions for pressure dof
        if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
          BernsteinBasisFunsTetra(1, param[0], param[1], param[2], &Np(0));
        else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
          BernsteinBasisFunsWedge(1, param[0], param[1], param[2], &Np(0));
        else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
          BernsteinBasisFunsHexa(1, param[0], param[1], param[2], &Np(0));


        // evaluate pressure at the quadrature points
        pres = computeValue2(0, Np);

        MatlData->computeStressAndTangent(true, sss, Fn, F, pres, stre, Cmat, ivar, gp, dt);
        //MatlData->computeMechanicalStress(F, pres, stre);
        //vector2matrix(stre, stresMat);

        //stresMat = detF*stresMat*F.inverse().transpose();

        //matrix2vector(stresMat, stre);

        if(varindex < 9)
           outval[gp] = stre[varindex];
        else if(varindex == 9)
           outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6.0*(stre[1]*stre[1]+stre[2]*stre[2]+stre[5]*stre[5]))/2.0);
        else if(varindex == 10)
           outval[gp] = pres;
    } //gp

    return;
}




void Elem_Solid_3D_Disp::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void Elem_Solid_3D_Disp::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}






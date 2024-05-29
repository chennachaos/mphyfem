
/* Node ordering for the quadratic edge element
 1-----3-----2
*/


#include "BernsteinElem2DEdge3Node.h"
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



BernsteinElem2DEdge3Node::BernsteinElem2DEdge3Node()
{
  ndim   = 2;
  degree = 2;
  npElem = 3;
  nlbfU  = 3;
  nlbfF  = 0;
  nlbfP  = 0;
  ndof   = 2;
  nsize  = npElem*ndof;
}


BernsteinElem2DEdge3Node::~BernsteinElem2DEdge3Node()
{
}


void BernsteinElem2DEdge3Node::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double BernsteinElem2DEdge3Node::computeVolume(bool init)
{
  double  dvol, Jac, param[2];

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

  int   ii, gp;

  double xNode[6], yNode[6], xx, yy;
  for(ii=0;ii<npElem;ii++)
  {
    //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
    //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

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
            GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          else
            GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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




int BernsteinElem2DEdge3Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1, TJ, TJp1;
    int  kk, ind1, ind2;

    MatrixXd  Mlocal(nsize, nsize);
    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);
    VectorXd  dispC(nsize), velC(nsize), accC(nsize);

    double  F[4], detF, F33, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  stre[4], cc[4][4], bc[2][4], param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];

    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
      //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    for(ii=0;ii<npElem;ii++)
    {
      ind1 = ndof*ii;
      ind2 = nodeNums[ii]*ndof;

      for(kk=0;kk<ndof;kk++)
      {
        dispC(ind1+kk)  =  SolnData->var1Cur[ind2+kk];
        velC(ind1+kk)   =  SolnData->var1DotCur[ind2+kk];
        accC(ind1+kk)   =  SolnData->var1DotDotCur[ind2+kk];
      }
    }

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();
    Mlocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    elemVol=0.0;
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol0 = gaussweights[gp]*(Jac*thick);
          dvol  = dvol0;
          elemVol += dvol;

          xx = yy= 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            xx += N[ii]*xNode[ii];
            yy += N[ii]*yNode[ii];
          }

          //if(axsy)
            //dvol *= 2.0*PI*yy;

          F[0] = computeValueCur(0, dN_dx) + 1.0;
          F[2] = computeValueCur(0, dN_dy);
          F[1] = computeValueCur(1, dN_dx);
          F[3] = computeValueCur(1, dN_dy) + 1.0;

          detF =  F[0]*F[3] - F[1]*F[2];

          if(detF < 0.0)   return 1;

          veloCur[0] = computeValueDotCur(0, N);
          veloCur[1] = computeValueDotCur(1, N);

          acceCur[0] = computeValueDotDotCur(0, N);
          acceCur[1] = computeValueDotDotCur(1, N);

          if(finite)
          {
            GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

            dvol = gaussweights[gp]*(Jac*thick);
          }

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
        {
          F33 = 1.0;
        }

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        if(finite)
          dvol *= F33;

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        //force[0] = analy.forceX(xx, yy);
        //force[1] = analy.forceY(xx, yy);

        //force[0] -= rho0*acceCur[0];
        //force[1] -= rho0*acceCur[1];

        // Calculate Stiffness and Residual
        //==============================================

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          bc[0][0] = bb1 * cc[0][0] + bb2 * cc[3][0];
          bc[0][1] = bb1 * cc[0][1] + bb2 * cc[3][1];
          bc[0][2] = bb1 * cc[0][3] + bb2 * cc[3][3];

          bc[1][0] = bb2 * cc[1][0] + bb1 * cc[3][0];
          bc[1][1] = bb2 * cc[1][1] + bb1 * cc[3][1];
          bc[1][2] = bb2 * cc[1][3] + bb1 * cc[3][3];

          TI   = 2*ii;
          TIp1 = TI+1;

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[3] + bb2*stre[1];

          Flocal(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal(TIp1) += (bb5*force[1] - sig[1]) ;

          for(jj=0; jj<nlbfU; jj++)
          {
              cc1 = dN_dx[jj];
              cc2 = dN_dy[jj];
              cc3 = N[jj];

              TJ   = 2*jj;
              TJp1 = TJ+1;

              // acceleration term
              //acceFact2 = acceFact1*cc3*rho0;
              acceFact2 = cc3*rho0;

              fact  = bb5*acceFact2;

              //consistent mass matrix
              //Klocal(TI,   TJ)    += fact ;
              //Klocal(TIp1, TJp1)  += fact ;

              Mlocal(TI,   TJ)    += fact ;
              Mlocal(TIp1, TJp1)  += fact ;

              // material Stiffness contribution
              fact  = af*(sig[0]*cc1+sig[1]*cc2)*FiniteFact;

              Klocal(TI,   TJ)    += fact ;
              Klocal(TIp1, TJp1)  += fact ;

              Klocal(TI,   TJ)    +=  af*(bc[0][0] * cc1 + bc[0][2] * cc2) ;
              Klocal(TI,   TJp1)  +=  af*(bc[0][1] * cc2 + bc[0][2] * cc1) ;
              Klocal(TIp1, TJ)    +=  af*(bc[1][0] * cc1 + bc[1][2] * cc2) ;
              Klocal(TIp1, TJp1)  +=  af*(bc[1][1] * cc2 + bc[1][2] * cc1) ;
          }
        }
    }//gp
    //cout << " elemVol = " << elemVol << endl;
    //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);
    //printVector(Flocal);

    //bool MassLumping=false;
    bool MassLumping=true;
    // mass lumping - row-wise sum
    if(MassLumping)
    {
      for(ii=0; ii<nsize; ii++)
      {
        fact=0.0;
        for(jj=0; jj<nsize; jj++)
        {
          fact += Mlocal(ii,jj);
          Mlocal(ii,jj) = 0.0;
        }
        Mlocal(ii,ii) = fact;
      }
    }

    Flocal -= Mlocal*accC;
    Klocal += acceFact1*Mlocal;

    return 0;
}



int BernsteinElem2DEdge3Node::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
    int  ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  fact, dvol, Jac, bb1, cc1, param[2];
    double rho0 = elmDat[5] ;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    if(Mlocal.rows() != nsize)
      Mlocal.resize(nsize, nsize);

    Mlocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    elemVol=0.0;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol = gaussweights[gp]*(thick*Jac);

        xx = yy= 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
        }

        if(axsy)
          dvol *= 2.0*PI*yy;

        elemVol += dvol;

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = (dvol*rho0)*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          for(jj=0; jj<nlbfU; jj++)
          {
              TJ   = 2*jj;
              TJp1 = TJ+1;

              fact  = bb1*N[jj];

              Mlocal(TI,   TJ)    += fact ;
              Mlocal(TIp1, TJp1)  += fact ;
          }
        }
    } //gp
    //cout << " elemVol = " << elemVol << endl;
    //printMatrix(Klocal);  printf("\n\n\n");

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      for(ii=0; ii<nsize; ii++)
      {
        fact=0.0;
        for(jj=0; jj<nsize; jj++)
        {
          fact += Mlocal(ii,jj);
          Mlocal(ii,jj) = 0.0;
        }
        Mlocal(ii,ii) = fact;
      }
    }

    return 0;
}



int BernsteinElem2DEdge3Node::calcResidual(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1, TJ, TJp1;
    int  kk, ind1, ind2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, F33, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  stre[4], cc[4][4], bc[2][4], param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];

    double rho  = elmDat[5] ;
    double rho0 = rho ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
      //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }


    if(Flocal.rows() != nsize)
      Flocal.resize(nsize);

    Flocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    elemVol=0.0;
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol0 = gaussweights[gp]*(Jac*thick);
          dvol  = dvol0;
          elemVol += dvol;

          xx = yy= 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            xx += N[ii]*xNode[ii];
            yy += N[ii]*yNode[ii];
          }

          //if(axsy)
            //dvol *= 2.0*PI*yy;

          F[0] = computeValueCur(0, dN_dx) + 1.0;
          F[2] = computeValueCur(0, dN_dy);
          F[1] = computeValueCur(1, dN_dx);
          F[3] = computeValueCur(1, dN_dy) + 1.0;

          detF =  F[0]*F[3] - F[1]*F[2];

          if(detF < 0.0)   return 1;

          veloCur[0] = computeValueDotCur(0, N);
          veloCur[1] = computeValueDotCur(1, N);

          acceCur[0] = computeValueDotDotCur(0, N);
          acceCur[1] = computeValueDotDotCur(1, N);

          if(finite)
          {
            GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

            dvol = gaussweights[gp]*(Jac*thick);
          }

          // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
        {
          F33 = 1.0;
        }

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        if(finite)
          dvol *= F33;

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        // Calculate Residual
        //==============================================
        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[3] + bb2*stre[1];

          Flocal(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal(TIp1) += (bb5*force[1] - sig[1]) ;
        }
    }//gp


    return 0;
}






int BernsteinElem2DEdge3Node::calcLoadVector(VectorXd& Flocal)
{
    int  ii, jj, kk, gp, TI, TIp1;
    double  param[2], trac[2], normal[2], tang[2], N[3], dN_dx[3], dN_dy[3];
    double  fact, dvol, dvol0, Jac, curvature, nx, ny, bb1, bb2, bb3;

    int  config = (int) elmDat[0];
    double tx = elmDat[1];
    double ty = elmDat[2];
    double gamma = elmDat[4];


    double xNode[6], yNode[6], xx, yy;
    if(config)
    {
      for(ii=0;ii<npElem;ii++)
      {
        yNode[ii] = GeomData->NodePosNew[nodeNums[ii]][1];
        xNode[ii] = GeomData->NodePosNew[nodeNums[ii]][0];
      }
    }
    else
    {
      for(ii=0;ii<npElem;ii++)
      {
        xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
        yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      }
    }

    if(Flocal.rows() != nsize)
      Flocal.resize(nsize);
    Flocal.setZero();


    vector<double>  gausspoints1, gaussweights;
    nGP = 3;
    getGaussPoints1D(nGP, gausspoints1, gaussweights);


    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];

        BernsteinBasisFunsEdge2D(degree, param, xNode, yNode, N, dN_dx, dN_dy, normal, curvature, Jac);
        nx = normal[0];
        ny = normal[1];

        dN_dy[0] = 0.0;  dN_dy[1] = 0.0;  dN_dy[2] = 0.0;

        dvol = gaussweights[gp]*Jac;
        elemVol += dvol;

        //cout << " normal = " << normal[0] << '\t' << normal[1] << '\t' << dvol << '\t' << curvature << endl;
        //cout << " Funcs = " << dN_dx[0] << '\t' << dN_dx[1] << '\t' << dN_dx[2] << endl;
        //cout << " Funcs = " << dN_dy[0] << '\t' << dN_dy[1] << '\t' << dN_dy[2] << endl;

        xx = yy= 0.0;
        for(ii=0;ii<3;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
        }

        if(axsy)
          dvol *= 2.0*PI*xx;

        tang[0] = -normal[1];
        tang[1] =  normal[0];

        //  specified traction values
        trac[0] = tx*normal[0]+ty*tang[0];
        trac[1] = tx*normal[1]+ty*tang[1];

        // specified surface tension
        trac[0] += (2.0*curvature*gamma*normal[0]);
        trac[1] += (2.0*curvature*gamma*normal[1]);

        for(ii=0;ii<3;ii++)
        {
          bb3 = N[ii]*dvol;

          TI   = ii*ndof;
          TIp1 = TI+1;

          Flocal(TI)   += bb3*trac[0];
          Flocal(TIp1) += bb3*trac[1];

          //Flocal(TI)   -= ( (1.0-nx*nx)*dN_dx[ii] -      nx*ny *dN_dy[ii] )*gamma;
          //Flocal(TIp1) -= (     -nx*ny *dN_dx[ii] + (1.0-ny*ny)*dN_dy[ii] )*gamma;
        }
    } // for(gp=0; gp<nGP; gp++)
    //cout << " Volume = " << elemVol << endl;

    //printVector(Flocal);

    return 0;
}




int BernsteinElem2DEdge3Node::calcInternalForces()
{
  return 0;
}



void BernsteinElem2DEdge3Node::elementContourplot(int vartype, int varindex, int index)
{
  return;
}



void BernsteinElem2DEdge3Node::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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

                projectStress(extrapolateFlag, vartype, varindex, index, outval);

              break;

       case 5:  // plot element internal variables

                projectInternalVariable(extrapolateFlag, vartype, varindex, index, outval);

              break;

       default:

              cout  << " Invalid Variable Type to project in 'BernsteinElem2DEdge3Node::projectToNodes'" << endl;
              break;
    }


    assert(vals2project.size() == npElem);

    if(extrapolateFlag)
    {
      stressrecovery_extrapolate_Triangle(degree, outval, &vals2project[0]);
    }
    else
    {
       for(int ii=0; ii<npElem; ii++)
         vals2project[ii] = outval[ii];
    }

    //printVector(vals2project);

    return;
}


void BernsteinElem2DEdge3Node::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    double F[4], detF, F33, dvol, Jac, param[2];
    double  stre[4], cc[4][4];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    int   err,  isw,  count,  count1, ll = 0, ii, jj, gp, nn;

    double dt   = myTime.dt;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    if(extrapolateFlag)
    {
      nn = nGP;
      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);
    }
    else
    {
      nn = npElem;
      gausspoints1.resize(npElem);
      gausspoints2.resize(npElem);
      gaussweights.resize(npElem);

      gausspoints1[0] =  0.0;    gausspoints2[0] =  0.0;
      gausspoints1[1] =  1.0;    gausspoints2[1] =  0.0;
      gausspoints1[2] =  0.0;    gausspoints2[2] =  1.0;
      gausspoints1[3] =  0.5;    gausspoints2[3] =  0.0;
      gausspoints1[4] =  0.5;    gausspoints2[4] =  0.5;
      gausspoints1[5] =  0.0;    gausspoints2[5] =  0.5;
    }

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nn; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[2] = computeValueCur(0, dN_dy);
        F[1] = computeValueCur(1, dN_dx);
        F[3] = computeValueCur(1, dN_dy) + 1.0;

        detF = F[0]*F[3] - F[1]*F[2];

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
        {
          F33 = 1.0;
        }

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(varindex < 4)
           outval[gp] = stre[varindex];
        else if(varindex == 6)
           outval[gp] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2);
        else if(varindex == 7)
           outval[gp] = (stre[0]+stre[1]+stre[2])/3.0;
    }//gp

    return;
}



void BernsteinElem2DEdge3Node::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void BernsteinElem2DEdge3Node::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}




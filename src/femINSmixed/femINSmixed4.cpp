
#include "femINSmixed.h"
#include "headersVTK.h"
#include "BasisFunctionsBernstein.h"
#include "BasisFunctionsLagrange.h"
#include "MyTime.h"
#include "TimeFunction.h"

extern   std::vector<unique_ptr<TimeFunction> > timeFunction;
extern MyTime                 myTime;



int  femINSmixed::applyInterfaceTerms2D()
{
    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the immersed boundary points
    // using Lagrange multipliers
    //
    ////////////////////////////////////////////////////////

    myPoint coords_global, coords_local;
    vector<int>  elemforAssyVecVelo, lmenodeNums, lmeforAssyVec;
    vector<double>  Nf(npElemVelo), dNf_du1(npElemVelo), dNf_du2(npElemVelo);
    VectorXd  Nb, dNb, xNodesIIE, yNodesIIE, specValx, specValy;
    MatrixXd  Khorz;
    VectorXd  Flocal, Flocal2;
    myPoint geom, lagmults, velFluid, velSolid, velDiff, res;

    ImmersedIntegrationElement  *lme;
    ElementBase *elemCur;

    int  ii, jj, ee, row, col, ind, ind1, ind2, r, c, c1, c2, nlb, gp, nGP, elnum;
    double  fact, fact1, fact2, dvol, af = 1.0, detJ;

    nlb = ImmersedBodyObjects[0]->ImmIntgElems[0]->nodeNums.size();

    Nb.resize(nlb);
    dNb.resize(nlb);
    xNodesIIE.resize(nlb);
    yNodesIIE.resize(nlb);
    specValx.resize(nlb);
    specValy.resize(nlb);

    for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      for(int iie=0; iie<ImmersedBodyObjects[bb]->getNumberOfElements(); iie++)
      {
        lme = ImmersedBodyObjects[bb]->ImmIntgElems[iie];
        lmenodeNums  = lme->nodeNums;
        lmeforAssyVec = lme->forAssyVec;

        for(ii=0;ii<nlb;ii++)
        {
          ind = lmenodeNums[ii];

          xNodesIIE[ii] = ImmersedBodyObjects[bb]->nodePosCur[ind][0];
          yNodesIIE[ii] = ImmersedBodyObjects[bb]->nodePosCur[ind][1];

          specValx[ii] = ImmersedBodyObjects[bb]->specVelocity[ind][0];
          specValy[ii] = ImmersedBodyObjects[bb]->specVelocity[ind][1];
        }

        
        for(gp=0; gp<lme->gausspoints.size(); gp++)
        {
           //TODO
            //computeLagrangeBFsLine2D(nlb-1, lme->gausspoints[gp], &xNodesIIE(0), &yNodesIIE(0), &Nb(0), &dNb(0), detJ);

           dvol = lme->gaussweights[gp] * detJ;
           
           elnum   = lme->elemNums[gp];
           elemCur = elems[elnum];


              lagmults.setZero();
              geom.setZero();
              velSolid.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0]      +=  Nb[ii] * xNodesIIE[ii];
                geom[1]      +=  Nb[ii] * yNodesIIE[ii];

                velSolid[0]  +=  Nb[ii] * specValx[ii];
                velSolid[1]  +=  Nb[ii] * specValy[ii];

                lagmults[0]  +=  Nb[ii] * lambdasCur(lmenodeNums[ii]*ndim);
                lagmults[1]  +=  Nb[ii] * lambdasCur(lmenodeNums[ii]*ndim+1);
              }

              ImmersedBodyObjects[bb]->computePointAtGP(iie, lme->gausspoints[gp], geom);

              //cout << " Nb            " << Nb[0] << '\t' << Nb[1] << endl;
              //cout << " coords_global " << geom[0] << '\t' << geom[1] << endl;

              //TODO
              //elemCur->findLocalCoordinates(nodeCoords, geom, coords_local);

              //cout << " coords_local " << coords_local[0] << '\t' << coords_local[1] << endl;

              if(npElemVelo == 6)
                BernsteinBasisFunsTria(2,  coords_local(0), coords_local(1), &Nf[0], &dNf_du1[0], &dNf_du2[0]);
              else if(npElemVelo == 7)
                BernsteinBasisFunsTria(21, coords_local(0), coords_local(1), &Nf[0], &dNf_du1[0], &dNf_du2[0]);
              else if(npElemVelo == 9)
                BernsteinBasisFunsQuad(2,  coords_local(0), coords_local(1), &Nf[0], &dNf_du1[0], &dNf_du2[0]);

              velFluid.setZero();
              for(ii=0;ii<npElemVelo;ii++)
              {
                jj = ndim*elemCur->nodeNums[ii];

                velFluid(0) += (Nf[ii]*velo[jj]);
                velFluid(1) += (Nf[ii]*velo[jj+1]);
              }

              //cout << " Lambdas " << ime << '\t' << lagmults[0] << '\t' << lagmults[1] << endl;
              //cout << " Velocity " << ime << '\t' << velF(0) << '\t' << velF(1) << endl;
              //cout << endl;  cout << endl;

              velDiff = velFluid - velSolid;

              ind1 = lme->nodeNums.size();
              ind2 = elemCur->nodeNums.size();

              //cout << " ime " << ime << '\t' << ind1 << '\t' << ind2 << endl;
              //cout << " Lambdas " << lagmults[0] << '\t' << lagmults[1] << endl;
              //cout << " npElem = " << npElem << endl;

              Khorz.resize(ind1*ndim, ind2*ndof);    Khorz.setZero();

              Flocal.resize(ind2*ndim);   Flocal.setZero();
              Flocal2.resize(ind1*ndim);   Flocal2.setZero();

              for(ii=0;ii<ind2;ii++)
              {
                r = ndim*ii;

                fact = Nf[ii]*dvol;

                //cout << ii << '\t' << r << '\t' << fact << endl;

                Flocal(r)   -= fact*lagmults[0];
                Flocal(r+1) -= fact*lagmults[1];
              }

              //cout << " aaaaaaaaaaa " << endl;
              
              for(ii=0;ii<ind1;ii++)
              {
                r = ndim*ii;
                fact1 = Nb[ii]*dvol;

                Flocal2(r)   -= fact1*velDiff(0);
                Flocal2(r+1) -= fact1*velDiff(1);

                fact1 = af*fact1;
                //cout << ii << '\t' << r << '\t' << fact1 << endl;

                for(jj=0;jj<ind2;jj++)
                {
                  c = ndof*jj;

                  fact2 = fact1*Nf[jj];
                  //cout << jj << '\t' << c << '\t' << fact2 << endl;
                  Khorz(r,   c)   += fact2;
                  Khorz(r+1, c+1) += fact2;
                }
              }

              //printMatrix(Khorz);
              //printf("\n\n");
              //printVector(Flocal);
              //printf("\n\n");
              //printVector(Flocal2);
              //printf("\n\n");

              ind1 = lme->forAssyVec.size();

              elemforAssyVecVelo = elemCur->forAssyVecVelo;
              ind2 = elemforAssyVecVelo.size();

              //printVector(lme->forAssyVec);
              //printf("\n\n");

              for(ii=0;ii<ind2;ii++)
              {
                //c1 = ndim*elemCur->nodeNums[ii];
                //c2 = ndim*ii;
                r = elemCur->forAssyVecVelo[ii];

                //cout << ii << '\t' << "c " << '\t' << c << endl;

                if( r != -1)
                  rhsVec[r]   += Flocal(ii) ;
              }

              //printVector(elemCur->forAssyVecVelo);
              //printf("\n\n");

              for(ii=0;ii<ind1;ii++)
              {
                r = totalDOF_Velo + totalDOF_Pres + lme->forAssyVec[ii];

                //cout << ii << '\t' << "r " << '\t' << r << endl;

                rhsVec[r] += Flocal2(ii);

                for(jj=0;jj<ind2;jj++)
                {
                  c = elemforAssyVecVelo[jj];

                  fact = Khorz(ii, jj);

                  if(c != -1)
                  {
                    matK.coeffRef(r, c) += fact;
                    matK.coeffRef(c, r) += fact;
                  }
                }
              }
        }

    } //for(ime=0;)
    } //for(bb=0;)

    return 0;
}



/*
for semi-implicit scheme
int  femINSmixed::applyInterfaceTerms2D()
{
    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the immersed boundary points
    // using Lagrange multipliers
    //
    ////////////////////////////////////////////////////////

    myPoint coords_global, coords_local;
    vector<double>  Nf(9), dNf_du1(9), dNf_du2(9), Nb(9);
    VectorXd  Khorz(18);
    VectorXd  Flocal, Flocal2;
    myPoint geom, lagmults, velF, velS;

    int  ii, jj, ee, row, col, ind1, ind2, r, c, c1, c2, nlb = 1;
    double  fact, fact1, fact2, dvol, af = 1.0;

    for(int ime=0; ime<nImmersedElems; ime++)
    {
        int elnum = ImmersedElements[ime]->elemNums[0];

        elems[elnum]->findLocalCoordinates(nodeCoords, ImmersedElements[ime]->coords, coords_local);

        if(npElemVelo == 6)
          BernsteinBasisFunsTria(2, coords_local(0), coords_local(1), &Nf[0], &dNf_du1[0], &dNf_du2[0]);
        else if(npElemVelo == 7)
          BernsteinBasisFunsTria(21, coords_local(0), coords_local(1), &Nf[0], &dNf_du1[0], &dNf_du2[0]);
        else if(npElemVelo == 9)
          BernsteinBasisFunsQuad(2, coords_local(0), coords_local(1), &Nf[0], &dNf_du1[0], &dNf_du2[0]);

        Nb[0] = 1.0;
        Nb[1] = 1.0;
        dvol  = 1.0;

              lagmults.setZero();
              geom.setZero();

              lagmults[0] = lambdas(ime*ndim);
              lagmults[1] = lambdas(ime*ndim+1);

              ind1 = ImmersedElements[ime]->nodeNums.size();
              ind2 = elems[elnum]->nodeNums.size();

              //cout << " ime " << ime << '\t' << ind1 << '\t' << ind2 << endl;
              //cout << " Lambdas " << lagmults[0] << '\t' << lagmults[1] << endl;
              //cout << " npElem = " << npElem << endl;

              Khorz.resize(ind1*ndim, ind2*ndof);    Khorz.setZero();

              Flocal.resize(ind2*ndim);   Flocal.setZero();
              Flocal2.resize(ind1*ndim);   Flocal2.setZero();

              velF.setZero();
              velS.setZero();

              for(ii=0;ii<npElem;ii++)
              {
                jj = elems[elnum]->nodeNums[ii];

                velF(0) += (Nf[ii]*velo[jj*ndim]);
                velF(1) += (Nf[ii]*velo[jj*ndim+1]);
              }

              //cout << " Lambdas " << ime << '\t' << lagmults[0] << '\t' << lagmults[1] << endl;
              //cout << " Velocity " << ime << '\t' << velF(0) << '\t' << velF(1) << endl;
              //cout << endl;  cout << endl;

              velF -= velS;

              for(ii=0;ii<ind1;ii++)
              {
                r = ndim*ii;
                fact1 = Nb[ii]*dvol;

                Flocal2(r)   -= fact1*velF(0);
                Flocal2(r+1) -= fact1*velF(1);

                //fact1 = af*fact1;
                //cout << ii << '\t' << r << '\t' << fact1 << endl;

                //
                for(jj=0;jj<ind2;jj++)
                {
                  c = ndof*jj;

                  fact2 = fact1*Nf[jj];
                  //cout << jj << '\t' << c << '\t' << fact2 << endl;
                  Khorz(r,   c)   += fact2;
                  Khorz(r+1, c+1) += fact2;
                }
                //
              }

              for(ii=0;ii<ind2;ii++)
              {
                r = ndim*ii;

                fact = Nf[ii]*dvol;

                //cout << ii << '\t' << r << '\t' << fact << endl;

                Flocal(r)   -= fact*lagmults[0];
                Flocal(r+1) -= fact*lagmults[1];
              }

              //printMatrix(Khorz);
              //printf("\n\n");
              //printVector(Flocal);
              //printf("\n\n");
              //printVector(Flocal3);
              //printf("\n\n");

              ind1 = ImmersedElements[ime]->forAssyVec.size();
              ind2 = elems[elnum]->forAssyVecVelo.size();
              ind2 = elems[elnum]->nodeNums.size();

              //printVector(lme->forAssyVec);
              //printf("\n\n");
              //printVector(nd->forAssyVec);
              //printf("\n\n");

              for(ii=0;ii<ind2;ii++)
              {
                c1 = ndim*elems[elnum]->nodeNums[ii];
                c2 = ndim*ii;

                //cout << ii << '\t' << "c " << '\t' << c << endl;

                rhsVecVeloTemp[c1]   += Flocal(c2) ;
                rhsVecVeloTemp[c1+1] += Flocal(c2+1) ;
              }

              for(ii=0;ii<ind1;ii++)
              {
                r = totalDOF_Pres + ImmersedElements[ime]->forAssyVec[ii];

                //cout << ii << '\t' << "r " << '\t' << r << endl;

                rhsVecPresTemp[r] += Flocal2(ii);
              }

              for(ii=0;ii<ind1;ii++)
              {
                r = fluidDOF + lme->forAssyVec[ii];

                solverEigen->rhsVec[r] += Flocal2(ii);

                for(jj=0;jj<ind2;jj++)
                {
                  c = nd->forAssyVec[jj];

                  fact = Khorz(ii, jj);

                  solverEigen->mtx.coeffRef(r, c) += fact;
                  solverEigen->mtx.coeffRef(c, r) += fact;
                }
              }

    } //for(ee=0;)
  return 0;
}
*/


/*
    int ii, jj, aa, bb, c, r, gp, nr;

    MatrixXd  Klocal, Klocal2;
    VectorXd  Flocal;
    node *nd;

    double PENALTY;
    ImmersedIntegrationElement  *lme;

    int nlf=(degree[0]+1)*(degree[0]+1);

    int nlb, ind1, ind2, nGauss;

      VectorXd  NN(nlf), dNN_dx(nlf), dNN_dy(nlf), dN_dx, dN_dy, Nf;
      VectorXd  Flocal2, vel(DIM), vel2(DIM), lagmults(DIM), Nb, dNb, xx, yy,  specValx, specValy, res(DIM);
      MatrixXd  Khorz;
      myPoint  knotIncr, knotBegin, knotEnd;

      double  detJ, af, dvol, fact, fact1, fact2;

      af = SolnData.td(2);
      
      bool axsy = (GeomData.FluidProps[2] == 1);
      double  rho = GeomData.FluidProps[3];
      double  mu = GeomData.FluidProps[4];


      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        nlb = ImmersedBodyObjects[0]->ImmIntgElems[0]->nodeNums.size();

        Nb.resize(nlb);
        dNb.resize(nlb);
        xx.resize(nlb);
        yy.resize(nlb);
        specValx.resize(nlb);
        specValy.resize(nlb);

        if(ImmersedBodyObjects[bb]->isBoundaryConditionTypeLagrange())
        {
          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->nodeNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->nodeNums[ii]][1];

              specValx[ii] = lme->GeomDataLag->specValCur[lme->nodeNums[ii]][0];
              specValy[ii] = lme->GeomDataLag->specValCur[lme->nodeNums[ii]][1];
            }

            //cout << specValx[0] << '\t' << specValy[0] << endl;
            //cout << specValx[1] << '\t' << specValy[1] << endl;

            //cout << " lme->gausspoints.size() = " << lme->gausspoints.size() << endl;

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine2D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << lme->gaussweights[gp] << '\t' << Nb[0] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              lagmults.setZero();
              geom.setZero();
              vel2.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0]     += Nb[ii] * xx[ii];
                geom[1]     += Nb[ii] * yy[ii];

                vel2[0]     += Nb[ii] * specValx[ii];
                vel2[1]     += Nb[ii] * specValy[ii];

                lagmults[0] += Nb[ii] * SolnData.var3Cur(lme->nodeNums[ii]*DIM);
                lagmults[1] += Nb[ii] * SolnData.var3Cur(lme->nodeNums[ii]*DIM+1);
              }

              if(axsy)
              {
                dvol *= 2.0*PI*geom[0];
              }

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel2[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d \n", Nb[0], nlb);
              //printf("vel2[0] = %12.6f, vel2[1] = %12.6f \n", vel2[0], vel2[1]);
              //printf("lagmults[0] = %12.6f, lagmults[1] = %12.6f \n", lagmults[0], lagmults[1]);
              //printf("detJ = %12.6f, gw = %12.6f \n", detJ, lme->gaussweights[gp]);

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              knotBegin = nd->getKnotBegin();
              knotEnd   = nd->getKnotEnd();
              knotIncr  = nd->getKnotIncrement();

              GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

              if(nd->getParent() == NULL )
              {
                Nf = NN;
                dN_dx = dNN_dx;
                dN_dy = dNN_dy;
              }
              else
              {
                Nf    = nd->SubDivMat * NN;
                dN_dx = nd->SubDivMat * dNN_dx;
                dN_dy = nd->SubDivMat * dNN_dy;
              }

              //printVector(Nf);
              //printf("\n\n");

              ind1 = lme->nodeNums.size();
              ind2 = nd->GlobalBasisFuncs.size();

              Khorz.resize(ind1*DIM, ind2*ndof);    Khorz.setZero();

              Flocal.resize(ind2*ndof);   Flocal.setZero();
              Flocal2.resize(ind1*DIM);   Flocal2.setZero();

              vel(0) = vel2(0) - nd->computeValueCur(0, Nf);
              vel(1) = vel2(1) - nd->computeValueCur(1, Nf);

              for(ii=0;ii<ind1;ii++)
              {
                r = DIM*ii;
                fact1 = Nb[ii]*dvol;

                Flocal2(r)   += fact1*vel(0);
                Flocal2(r+1) += fact1*vel(1);

                fact1 = af*fact1;
                //cout << ii << '\t' << r << '\t' << fact1 << endl;

                for(jj=0;jj<ind2;jj++)
                {
                  c = ndof*jj;

                  fact2 = fact1*Nf[jj];
                  //cout << jj << '\t' << c << '\t' << fact2 << endl;
                  Khorz(r,   c)   += fact2;
                  Khorz(r+1, c+1) += fact2;
                }
              }

              for(ii=0;ii<ind2;ii++)
              {
                r = ndof*ii;

                fact = Nf[ii]*dvol;

                Flocal(r)   -= fact*lagmults[0];
                Flocal(r+1) -= fact*lagmults[1];
              }

              //printMatrix(Khorz);
              //printf("\n\n");
              //printVector(Flocal);
              //printf("\n\n");
              //printVector(Flocal3);
              //printf("\n\n");

              ind1 = lme->forAssyVec.size();
              ind2 = nd->forAssyVec.size();

              //printVector(lme->forAssyVec);
              //printf("\n\n");
              //printVector(nd->forAssyVec);
              //printf("\n\n");

              for(ii=0;ii<ind2;ii++)
              {
                c = nd->forAssyVec[ii];

                solverEigen->rhsVec[c] += Flocal(ii) ;
              }

              for(ii=0;ii<ind1;ii++)
              {
                r = fluidDOF + lme->forAssyVec[ii];

                solverEigen->rhsVec[r] += Flocal2(ii);

                for(jj=0;jj<ind2;jj++)
                {
                  c = nd->forAssyVec[jj];

                  fact = Khorz(ii, jj);

                  solverEigen->mtx.coeffRef(r, c) += fact;
                  solverEigen->mtx.coeffRef(c, r) += fact;
                }
              }

            }//for(gp=0...
          }//for(aa=0...
        }//if(
        else
        {
          myDataIntegrateCutFEM  myData;

          PENALTY = ImmersedBodyObjects[bb]->getPenaltyParameter();

          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->nodeNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->nodeNums[ii]][1];
            }

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine2D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << gaussweights[gp] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              vel.setZero();

              geom.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0] += Nb[ii] * xx[ii];
                geom[1] += Nb[ii] * yy[ii];
              }

              if(axsy)
              {
                dvol *= 2.0*PI*geom[0];
              }

              //cout << " uuuuuuuuuuu " << endl;

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              //cout << " uuuuuuuuuuu " << endl;
              //cout << " dvol " << dvol << endl;

              nr = nd->forAssyVec.size();

              myData.K1 = MatrixXd::Zero(nr, nr);
              myData.F1 = VectorXd::Zero(nr);

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //if(geom[1] == 1.0 )
                //vel[0] = 1.0;

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", param[0], param[1], param[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d, ee = %5d \n", Nb[0], nlb, lme->elemNums[gp]);

              //cout << " uuuuuuuuuuu " << endl;

              myData.geom  = geom;
              myData.param = param;
              myData.dvol  = PENALTY * dvol;

              for(jj=0;jj<DIM;jj++)
              {
                myData.dir = jj;
                myData.specVal[0] = vel[jj];

                //nd->applyBoundaryConditionsAtApoint(jj, param, vel[jj], fact, myData.K1, myData.F1);
                nd->applyBoundaryConditionsAtApoint(myData);
              }

              //solverEigen->assembleMatrixAndVector(velDOF, 0, nd->forAssyVec, nd->forAssyVec2, Klocal, Flocal);
              solverEigen->assembleMatrixAndVector(velDOF, 0, nd->forAssyVec, myData.K1, myData.F1);
              //cout << " uuuuuuuuuuu " << endl;
            }//for(gp=0...
          }//for(aa=0...
        }//else
      }//for(bb=0;...
      //
*/






int  femINSmixed::writeReadResult(int index, string& filename, int stride)
{
    if( (fileCount % outputFreq) !=  0)
        return 0;

    cout << "  femINSmixed::writeReadResult  " << endl;

    int  ii, jj, ind, ee, bb, count;
    double  val1, val2;

    if(index == 1)                                          // write result to a file
    {
        char tmp[500], ff[500];
        sprintf(ff, "%s%s%06d%s", filename.c_str(), "-", fileCount, ".rst");

        ofstream  fout;
        fout.open(filename);

        if(!fout)
        {
            runtime_error("femINSmixed::writeReadResult ... could not open file for writing!");
            return 0;
        }

        fout << "Dimension" << endl;
        fout << ndim << endl;
        fout << "NumberOfNodes" << endl;
        fout << nNode_global << endl;
        fout << "NumberOfElements" << endl;
        fout << nElem_global << endl;
        fout << "Time" << endl;
        fout << timeNow << endl;

        /*
        fout << "Node coordinates" << endl;
        for(ii=0; ii<nNode_global; ii++)
        {
            sprintf(tmp,                " %6d", (ii+1) );
            sprintf(&(tmp[strlen(tmp)])," %16.12f", nodeCoords[ii][0]);
            sprintf(&(tmp[strlen(tmp)])," %16.12f", nodeCoords[ii][1]);
            fout << tmp << "\n";
        }
        */

        fout << "NodalSolutions" << endl;
        if(ndim == 2)
        {
          //fout << '\t' << "Number" << '\t' << "X-Velocity" << '\t' << "Y-Velocity" << '\t' << "Pressure" << '\t' << "X-Acceleration" << '\t' <<  "Y-Acceleration" << endl;
          for(ii=0; ii<nNode_global; ii++)
          {
            ind = ii*ndim;

            sprintf(tmp,                " %6d", (ii+1) );
            sprintf(&(tmp[strlen(tmp)])," %16.12f", velo(ind));
            sprintf(&(tmp[strlen(tmp)])," %16.12f", velo(ind+1));
            sprintf(&(tmp[strlen(tmp)])," %16.12f", pres(ii));

            sprintf(&(tmp[strlen(tmp)])," %16.12f", veloDot(ind));
            sprintf(&(tmp[strlen(tmp)])," %16.12f", veloDot(ind+1));
            fout << tmp << "\n";
          }
        }

        /*
        fout << "Element Node Connectivity" << endl;
        fout << '\t' << "Number" << '\t' << "Node1" << '\t' << "Node2" << '\t' << "Node3" << '\t' << "Node4" << '\t' << "Node5" << '\t' <<  "Node6" << endl;
        for(ee=0; ee<nElem; ee++)
        {
            sprintf(tmp,                " %6d", ee+1);
            for(ii=0; ii<npElem; ii++)
              sprintf(&(tmp[strlen(tmp)])," %6d", (elemConn[ee][ii]+1) );

            fout << tmp << "\n";
        }
        */

        fout.close();
    }
    else                                                    // read result from a file
    {
        cout << "Input file: " << filename.c_str() << endl;
        ifstream  infile(filename);

        if(!infile)
        {
            cout << "femINSmixed::writeReadResult ... could not open file for reading! \n" << endl;
            exit(1);
        }

        string  line, stringVal, stringVec[10];
        int  ii, arrayInt[100], valInt;
        double  tempDbl, arrayDbl[10];
        std::stringstream   linestream(filename);
        vector<string>  stringlist;


        // read the dimension
        infile >> stringVal;
        infile >> valInt;
        assert(ndim ==  valInt);

        // read the number of nodes
        infile >> stringVal;
        infile >> valInt;
        cout << "valInt=" << valInt << endl;
        assert(nNode ==  valInt);

        // read the number of nodes
        infile >> stringVal;
        infile >> valInt;
        assert(nElem ==  valInt);

        // read the Time instant
        infile >> stringVal;
        infile >> tempDbl;
        timeNow = tempDbl;

        // read the nodal solutions
        infile >> stringVal;
        //getline(infile, line);
        //getline(infile, line);

        if(ndim == 2)
        {
          for(ii=0; ii<nNode_global; ii++)
          {
            infile >> arrayDbl[0] >> arrayDbl[1] >> arrayDbl[2] >> arrayDbl[3] >> arrayDbl[4] >> arrayDbl[5];

            ind = ii*ndim;

            veloPrev(ind)      = arrayDbl[1];
            veloPrev(ind+1)    = arrayDbl[2];

            presPrev(ii)       = arrayDbl[3];

            veloDotPrev(ind)   = arrayDbl[4];
            veloDotPrev(ind+1) = arrayDbl[5];
          }
        }
        else   //if(ndim == 2)
        {
          for(ii=0; ii<nNode_global; ii++)
          {
            infile >> arrayDbl[0] >> arrayDbl[1] >> arrayDbl[2] >> arrayDbl[3] >> arrayDbl[4] >> arrayDbl[5] >> arrayDbl[6] >> arrayDbl[7];

            ind = ii*ndim;

            veloPrev(ind)      = arrayDbl[1];
            veloPrev(ind+1)    = arrayDbl[2];
            veloPrev(ind+2)    = arrayDbl[3];

            presPrev(ii)       = arrayDbl[4];

            veloDotPrev(ind)   = arrayDbl[5];
            veloDotPrev(ind+1) = arrayDbl[6];
            veloDotPrev(ind+2) = arrayDbl[7];
          }
        }

        infile.close();

        //fileCount--;
        //timeUpdate();
    }

    return 0;
}




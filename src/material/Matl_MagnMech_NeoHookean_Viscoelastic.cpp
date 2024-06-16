#include "Matl_MagnMech_NeoHookean_Viscoelastic.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "util.h"
#include "TimeFunction.h"
#include "MyTime.h"

extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime                myTime;
extern bool debug;

using namespace std;


Matl_MagnMech_NeoHookean_Viscoelastic::Matl_MagnMech_NeoHookean_Viscoelastic()
{

}


Matl_MagnMech_NeoHookean_Viscoelastic::~Matl_MagnMech_NeoHookean_Viscoelastic()
{

}



int Matl_MagnMech_NeoHookean_Viscoelastic::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = mu/2 (J^{-2/3} trb - 3) + U(J)
    // U(J) is selected based on the input value for 'Utype'

    int err =  computeStressAndTangent_NeoHooke(sss, MIXED_ELEMENT, Utype, Kinv, data_deviatoric, Fn, F, pres, stre, Cmat);


    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0, r2d9 = 2.0/9.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, ii, jj;

    double  mu0   = data_magnfield[0];      // permeability
    double  J     = F.determinant();

    /////////////////////////////////
    // Cauchy stress tensor

    MatrixXd  streMagn(3,3);

    //printVector(ApplMagnfield);
    //printVector(ResiMagnfield);

    //printMatrix(streMagn);

    // similar to Zhao's papers in JMPS
    VectorXd  ResiMagnfieldCur = F*(ResiMagnfield/J);
    VectorXd  ApplMagnfieldCur = timeFunctions[0]->getValue()*ApplMagnfield;
    streMagn = (-1.0/mu0) * (ApplMagnfieldCur*ResiMagnfieldCur.transpose() );


    stre(0) += streMagn(0,0);    stre(3) += streMagn(0,1);    stre(6) += streMagn(0,2);
    stre(1) += streMagn(1,0);    stre(4) += streMagn(1,1);    stre(7) += streMagn(1,2);
    stre(2) += streMagn(2,0);    stre(5) += streMagn(2,1);    stre(8) += streMagn(2,2);

    /////////////////////////////////
    // material tangent tensor
    // this part due to the residual-applied field is zero


    /////////////////////////////////
    // Viscoelastic part

    tis = "BDF1";
    SetTimeParametersFluid(tis, spectralRadius, dt, td);

    printMatrix(Cmat);
    err = computeStressAndTangent_Viscoelasticity_Model1(data_viscoelastic, gp, F, td, ivar, dt, stre, Cmat);
    printMatrix(Cmat);

    return err;
}





int Matl_MagnMech_NeoHookean_Viscoelastic::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre, VectorXd&  magnDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = mu/2 (J^{-2/3} trb - 3) + K/2 (J-1)^2

    double  fact, fact1, fact2, fact3, fact31, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r4d3 = 4.0/3.0,  r2d9 = 2.0/9.0,  r5d3 = 5.0/3.0,  r7d3 = 7.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K      = 1.0/Kinv;  // bulk modulus

    double  mu     = data_deviatoric[0];                    // shear modulus
    double  mu0    = data_magnfield[0];                     // permeability

    double  alpha  = mu0*data_magnfield[1];                 // constant alpha
    double  beta   = mu0*data_magnfield[2];                 // constant beta
    double  eta    = mu0*data_magnfield[3];                 // constant eta



    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose();
    MatrixXd  bSq = b*b;
    VectorXd  bh = b*magnField, bSqh = bSq*magnField;


    double  J     = F.determinant();
    double  Jm1d3 = pow(J, -r1d3);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Ib = b(0,0) + b(1,1) + b(2,2);

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);
    double  dij, dik, dil, djk, djl, dkl;

    Cmat.setZero();
    Bmat.setZero();
    Amat.setZero();
    stre.setZero();
    magnDisp.setZero();

    double  hdoth = magnField[0]*magnField[0] + magnField[1]*magnField[1] + magnField[2]*magnField[2];
    double  bhdotbh = bh[0]*bh[0] + bh[1]*bh[1] + bh[2]*bh[2];

    /////////////////////////////////
    // Cauchy stress tensor

    // Mechanical part
    fact1 = mu*Jm5d3;

    fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    double  factI4_1  = alpha*2.0/J;
    double  factI5_1  = beta*2.0*Jm5d3;
    double  factI5t_1 = eta*2.0*Jm1d3;


    MatrixXd  stressMat(3,3);
    stressMat.setZero();

    for(i=0; i<3; i++)
    {
        // Magnetic displacement
        // 
        // free space
        magnDisp[i] += mu0*magnField[i];
        // I4
        magnDisp(i) -= factI4_1*bh(i);
        // I5
        magnDisp(i) -= factI5_1*bSqh(i);
        // I5tilde
        magnDisp(i) -= factI5t_1*magnField(i);


        for(j=0; j<3; j++)
        {
            dij = KronDelta(i,j);

            // Cauchy stress
            // 
            // Neo-Hooke part
            stressMat(i,j) += ( fact1*b(i,j) + fact2*dij);
            // free space
            stressMat(i,j) += mu0*( magnField(i)*magnField(j) - 0.5*hdoth*dij );
            // I5
            stressMat(i,j) += factI5_1*( bh(i)*bh(j) - r1d3*bhdotbh*dij );
            // I5tilde
            stressMat(i,j) += factI5t_1*( -magnField(i)*magnField(j) + r1d3*hdoth*dij );

            // permeability tensor
            // 
            // free space
            Amat(i,j) += mu0*dij;
            // I4
            Amat(i,j) -= factI4_1*b(i,j);
            // I5
            Amat(i,j) -= factI5_1*bSq(i,j);
            // I5tilde
            Amat(i,j) -= factI5t_1*dij;
        }
    }

    matrix2vector(stressMat, stre);


    /////////////////////////////////
    // material tangent tensor

    // with stress terms
    fact1 = mu*Jm5d3;
    fact2 = mu*Jm5d3*Ib;

    fact7 = J*d2UdJ2 + dUdJ;
    fact8 = -dUdJ;
    fact9 = -dUdJ;


    for(i1=0; i1<9; i1++)
    {
        i = map_tensor4[i1][0];
        j = map_tensor4[i1][1];

        dij = KronDelta(i,j);

        for(j1=0; j1<9; j1++)
        {
            k = map_tensor4[j1][0];
            l = map_tensor4[j1][1];

            dik = KronDelta(i,k);
            dil = KronDelta(i,l);
            djk = KronDelta(j,k);
            djl = KronDelta(j,l);
            dkl = KronDelta(k,l);

            // Devioteric term (without stress terms)
            //Cmat(i1,j1) += fact1*(dil*djk+dik*djl);
            //Cmat(i1,j1) += fact2*(b(i,j)*dkl + dij*b(k,l));
            //Cmat(i1,j1) += fact3*(dij*dkl);

            // Devioteric terms, with stress terms
            Cmat(i1,j1) += fact1*( dik*b(j,l) - r2d3*(b(i,j)*dkl + dij*b(k,l)) );
            Cmat(i1,j1) += fact2*( r1d3*dil*djk + r2d9*dij*dkl );


            if(MIXED_ELEMENT)
            {
                // pressure term (without stress terms)
                //Cmat(i1,j1) += pres*(dij*dkl - djk*dil - djl*dik);

                // pressure term (with stress terms)
                Cmat(i1,j1) += pres*(dij*dkl - djk*dil);
            }
            else
            {
                // volumetric term (without stress terms)
                //Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk)+fact9*(djl*dik);

                // volumetric term (with stress terms)
                Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk);
            }

            // free space
            Cmat(i1,j1) += mu0*magnField(i)*( magnField(j)*dkl - magnField(k)*djl - magnField(l)*djk );
            Cmat(i1,j1) += mu0*magnField(k)*( magnField(l)*dij - magnField(j)*dil );
            Cmat(i1,j1) += mu0*0.5*hdoth*( djk*dil - dij*dkl );

            // I5
            Cmat(i1,j1) += factI5_1*( r2d9*dij*dkl*bhdotbh - r2d3*dkl*bh(i)*bh(j));
            Cmat(i1,j1) += factI5_1*( bh(j)*bh(l)*dik + r1d3*djk*dil*bhdotbh - r2d3*dij*bh(k)*bh(l) );

            // I5tilde
            Cmat(i1,j1) += factI5t_1*( (r2d9*dij*dkl-r1d3*dil*djk)*hdoth);
            Cmat(i1,j1) -= factI5t_1*r2d3*(dij*magnField(k)*magnField(l) + dkl*magnField(i)*magnField(j) );
            Cmat(i1,j1) += factI5t_1*(dil*magnField(k)*magnField(j)+djk*magnField(i)*magnField(l)+djl*magnField(i)*magnField(k));
        }

        for(k=0; k<3; k++)
        {
            dik = KronDelta(i,k);
            djk = KronDelta(j,k);

            // free space
            Bmat(i1,k) -= mu0*( magnField(i)*djk + magnField(j)*dik - magnField(k)*dij );
            // I5
            Bmat(i1,k) -= factI5_1*( b(i,k)*bh(j) + b(j,k)*bh(i) - r2d3*bSqh(k)*dij );
            // I5tilde
            Bmat(i1,k) -= factI5t_1*( r2d3*dij*magnField(k) - dik*magnField(j) - djk*magnField(i) );

        }
    }

    cerr << " Matl_MagnMech_NeoHookean_Viscoelastic::computeStressAndTangent COUPLED needs updating with Viscoelastic model " << endl;
    exit(1);

    /*
    /////////////////////////////////
    // Viscoelastic part

    // Transpose of F
    MatrixXd  Ft = F.transpose();

    // Right Cauchy-Green tensor;
    MatrixXd  C = Ft*F, Cinv = C.inverse(), CbarInv = pow(J,r2d3)*Cinv;

    SetTimeParametersFluid(tis, spectralRadius, dt, td);

    double af    =  td(2);
    double am    =  td(3);
    double gamm  =  td(4);

    //cout << af << '\t' << am << '\t' << gamm << endl;

    int  numSeries = (int) data_viscoelastic[0];  // number of terms in the Prony series

    for(int seriescount=0; seriescount<numSeries; seriescount++)
    {
        int jj = 2*seriescount;

        double  mu   = data_viscoelastic[1+jj];    // modulus
        double  tau  = data_viscoelastic[1+jj+1];  // relaxation time


        VectorXd  AnVec = ivar.varPrev.block<9,1>(9*seriescount,gp), AVec(9), AdotVec(9);
        VectorXd  AdotnVec = ivar.varDotPrev.block<9,1>(9*seriescount,gp);

        MatrixXd  AnMat(3,3), AdotnMat(3,3);
        vector2matrix(AnVec, AnMat);
        vector2matrix(AdotnVec, AdotnMat);

        //MatrixXd  AMat = ((dt/tau)*CbarInv+AnMat)/(1.0+dt/tau);

        MatrixXd  AMat = ( (gamm*dt/am/tau)*CbarInv + (1.0 -(1.0-af)*gamm*dt/am/tau)*AnMat - ((gamm-am)*dt/am/gamm)*AdotnMat)/(1.0+af*gamm*dt/am/tau);

        matrix2vector(AMat, AVec);

        MatrixXd  AdotMat = (1.0/gamm/dt)*(AMat-AnMat) + ((gamm-1.0)/gamm)*AdotnMat;

        matrix2vector(AdotMat, AdotVec);

        // update the internal variables vector
        ivar.var.block<9,1>(9*seriescount,gp) = AVec ;
        ivar.varDot.block<9,1>(9*seriescount,gp) = AdotVec ;


        MatrixXd  Atilde = F*AMat*Ft;
        MatrixXd  Ahat = Atilde;
        //MatrixXd  Ahat = 0.5*(Atilde+Atilde.transpose());
        MatrixXd AC = AMat*C.transpose();
        double  IAC = AC.trace();

        /////////////////////////////////
        // Cauchy stress tensor

        fact1 = mu*Jm5d3;
        fact2 = fact1*IAC*r1d3;

        stre[0] += fact1*Ahat(0,0);    stre[3] += fact1*Ahat(0,1);    stre[6] += fact1*Ahat(0,2);
        stre[1] += fact1*Ahat(1,0);    stre[4] += fact1*Ahat(1,1);    stre[7] += fact1*Ahat(1,2);
        stre[2] += fact1*Ahat(2,0);    stre[5] += fact1*Ahat(2,1);    stre[8] += fact1*Ahat(2,2);
        stre[0] -= fact2 ;             stre[4] -= fact2 ;             stre[8] -= fact2 ;

        /////////////////////////////////
        // material tangent tensor

        fact1 = mu*Jm5d3;
        fact2 = fact1*IAC;
        fact3 = (mu/J)*( gamm*dt/(af*gamm*dt+am*tau) );

        for(i1=0; i1<9; i1++)
        {
            i = map_tensor4[i1][0];
            j = map_tensor4[i1][1];

            dij = KronDelta(i,j);

            for(j1=0; j1<9; j1++)
            {
              k = map_tensor4[j1][0];
              l = map_tensor4[j1][1];

              dik = KronDelta(i,k);
              dil = KronDelta(i,l);
              djk = KronDelta(j,k);
              djl = KronDelta(j,l);
              dkl = KronDelta(k,l);

              // Devioteric term (with stress terms)
              Cmat(i1,j1) += (fact1*dik*Ahat(l,j) - fact1*r2d3* (Ahat(i,j)*dkl+dij*Ahat(l,k)) );
              Cmat(i1,j1) += (fact2* ( r2d9*dij*dkl + r1d3*dil*djk) );
              Cmat(i1,j1) += (fact3*(-dik*djl-dil*djk + r2d3*dij*dkl) );
            }
        }
    }
    */


    return err;
}




int Matl_MagnMech_NeoHookean_Viscoelastic::computeMagneticComponents(VectorXd&  magnField, VectorXd&  magnDisp, MatrixXd&  Amat)
{
    double  eps0  = matData[2]; // permitivity

    /////////////////////////////////
    // Electrical displacement vector

    magnDisp[0] = eps0*magnField[0];
    magnDisp[1] = eps0*magnField[1];
    magnDisp[2] = eps0*magnField[2];


    /////////////////////////////////
    // electric permittivity tensor

    Amat.setZero();
    Amat(0,0) = eps0;  Amat(1,1) = eps0;  Amat(2,2) = eps0;

    return 0;
}



int Matl_MagnMech_NeoHookean_Viscoelastic::computeMagneticDisplacement(VectorXd&  magnField, VectorXd&  magnDisp)
{
    double  eps0  = matData[2]; // permitivity

    magnDisp[0] = eps0*magnField[0];
    magnDisp[1] = eps0*magnField[1];
    magnDisp[2] = eps0*magnField[2];

    return 0;
}



int Matl_MagnMech_NeoHookean_Viscoelastic::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;

    double  K      = 1.0/Kinv;                              // bulk modulus
    double  mu     = data_deviatoric[0];                    // shear modulus

    VectorXd  b(9), FF(9);
    FF(0) = F(0,0);    FF(3) = F(0,1);    FF(6) = F(0,2);
    FF(1) = F(1,0);    FF(4) = F(1,1);    FF(7) = F(1,2);
    FF(2) = F(2,0);    FF(5) = F(2,1);    FF(8) = F(2,2);
    computeLeftCGtensor(FF, b);

    double  J     = determinant3b3(FF);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Ib = b(0) + b(4) + b(8);

    double  fact1 = mu*Jm5d3;
    double  fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += K*(J-1.0);

    stre.setZero();
    for(int i=0; i<9; i++)
      stre[i] = fact1*b[i];

    stre[0] += fact2 ;
    stre[4] += fact2 ;
    stre[8] += fact2 ;

    return 0;
}



int Matl_MagnMech_NeoHookean_Viscoelastic::computeMaxwellStress(VectorXd&  magnField, VectorXd&  stre)
{
    double  eps0   = data_magnfield[0]; // permeability

    stre.setZero();

    addElectricPartToCauchyStress(eps0, magnField, stre);

    return 0;
}



int Matl_MagnMech_NeoHookean_Viscoelastic::computeTotalStress(MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre)
{
    double  eps0   = data_magnfield[0]; // permeability

    // Mechanical part
    computeMechanicalStress(F, pres, stre);

    // Electrical part
    addElectricPartToCauchyStress(eps0, magnField, stre);

    return 0;
}




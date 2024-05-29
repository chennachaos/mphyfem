
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "InternalVariables.h"
#include "util.h"


/*
Volumetric Energy functions

Type 1:
------
U(J) = (K/50)*(J^5 + 1/J^5 - 2)

Type 2:
------
U(J) = (K/32)*(J^2-1/J^2)^2

Type 3:
------
U(J) = 0.5*K*(J-1)^2

Type 4:
------
U(J) = 0.25*K*(J*J-1) - 0.5*K*ln(J)

Type 5:
------
U(J) = 0.25*K*(J-1)^2 + 0.25*K*(ln(J))^2

Type 6:
------
U(J) = K*(J*ln(J) - J + 1)

Type 7:
------
U(J) = K*(J - ln(J) - 1 )

Type 8:
------
U(J) = 0.5*K*(ln(J))^2
*/


double  volumetricfunctions_getFirstDerivative(int type, double K, double J)
{
    double  val=0.0;
    switch(type)
    {
        case -1:
            val = 0.0;
          break;

        case 1:
            val = 0.1*K*( J*J*J*J - 1.0/(J*J*J*J*J*J) );
          break;

        case 2:
            val = 0.125*K*(J*J-1.0/J/J)*(J+1.0/J/J/J);
          break;

        case 3:
            val = K*(J-1.0);
          break;

        case 4:
            val = 0.5*K*(J-1.0/J);
          break;

        case 5:
            val = 0.5*K*(J-1.0+log(J)/J);
          break;

        case 6:
            val = K*log(J);
          break;

        case 7:
            val = K*(1.0-1.0/J);
          break;

        case 8:
            val = K*log(J)/J;
          break;

        default:
          std::cerr << " This type is not implemented for Volumetric energy functions " << std::endl;
          break;
    }
    return val;
}


double  volumetricfunctions_getSecondDerivative(int type, double K, double J)
{
    double  val=0.0;
    switch(type)
    {
        case -1:
            val = 0.0;
          break;

        case 1:
            val = 0.2*K*( 2.0*J*J*J + 3.0/(J*J*J*J*J*J*J) );
          break;

        case 2:
            val = 0.125*K*( 3.0*J*J+5.0/(J*J*J*J*J*J) );
          break;

        case 3:
            val = K;
          break;

        case 4:
            val = 0.5*K*(1.0+1.0/J/J);
          break;

        case 5:
            val = 0.5*K*(J*J+1.0-log(J))/J/J;
          break;

        case 6:
            val = K/log(J);
          break;

        case 7:
            val = K/J/J;
          break;

        case 8:
            val = K*(1.0-log(J))/J/J;
          break;

        default:
          std::cerr << " This type is not implemented for Volumetric energy functions " << std::endl;
          break;
    }
    return val;
}



void  compute_constants_volumetric_functions(bool finite, int Utype, double Jn, double BULK, double& Jhat, double& thetahat)
{
    if(finite)
    {
        if(Utype == -1)
        {
            Jhat     = 1.0;
            thetahat = 0.0;
        }
        else
        {
            double dUdJ   = volumetricfunctions_getFirstDerivative(Utype, BULK, Jn);
            double d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, BULK, Jn);

            Jhat = Jn - dUdJ/d2UdJ2;
            thetahat = 1.0/d2UdJ2;
        }
    }
    else
    {
        Jhat  = 1.0;

        if(Utype == -1)
            thetahat = 0.0;
        else
            thetahat = 1.0/BULK;
    }

    return;
}




// Compute stress and tangent for the Saint Venant-Kirchhoff model
int computeStressAndTangent_SaintVenantKirchhoff(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat)
{
    return 0;
}


// Compute stress and tangent for the Neo-Hookean model
int computeStressAndTangent_NeoHooke(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r5d3 = 5.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K      = 1.0/Kinv;                              // bulk modulus
    double  mu     = matData[0];                            // shear modulus

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose();

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Ib = b.trace();

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);
    double  dij, dik, dil, djk, djl, dkl;

    Cmat.setZero();
    stre.setZero();

    /////////////////////////////////
    // Cauchy stress tensor

    fact1 = mu*Jm5d3;
    fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    stre[0] = fact1*b(0,0);    stre[3] = fact1*b(0,1);    stre[6] = fact1*b(0,2);
    stre[1] = fact1*b(1,0);    stre[4] = fact1*b(1,1);    stre[7] = fact1*b(1,2);
    stre[2] = fact1*b(2,0);    stre[5] = fact1*b(2,1);    stre[8] = fact1*b(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

    /////////////////////////////////
    // material tangent tensor

        // without stress terms
        fact1 =  r1d3*mu*Jm5d3*Ib;
        fact2 = -r2d3*mu*Jm5d3;
        fact3 =  r2d3*fact1;

        // with stress terms
        fact1 =  r1d3*mu*Jm5d3*Ib;
        fact2 = -r2d3*mu*Jm5d3;
        fact3 =  r2d3*r1d3*mu*Jm5d3*Ib;
        fact4 =  mu*Jm5d3;

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
              Cmat(i1,j1) += fact1*(dil*djk);
              Cmat(i1,j1) += fact2*(b(i,j)*dkl + dij*b(k,l));
              Cmat(i1,j1) += fact3*(dij*dkl);
              Cmat(i1,j1) += fact4*(b(j,l)*dik);

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
                Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk)+fact9*(djl*dik);

                // volumetric term (with stress terms)
                //Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk);
              }
            }
        }

    return 0;
}



// Compute stress and tangent for the Mooney-Rivlin model
int computeStressAndTangent_MooneyRivlin(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r4d3 = 4.0/3.0,  r2d9 = 2.0/9.0,  r5d3 = 5.0/3.0,  r7d3 = 5.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K      = 1.0/Kinv;                              // bulk modulus
    double  c1     = matData[0];                           // coefficient 1
    double  c2     = matData[1];                           // coefficient 2

    // Left Cauchy-Green tensor;
    MatrixXd  C = F.transpose()*F;
    MatrixXd  b = F*F.transpose();
    MatrixXd  bSq = b*b;

    double  J     = F.determinant();
    double  Jm1d3 = pow(J, -r1d3);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Jm7d3 = pow(J, -r7d3);
    double  I1 = b.trace();
    double  I2 = 0.5*(I1*I1 - bSq.trace());

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);
    double  dij, dik, dil, djk, djl, dkl;

    Cmat.setZero();
    stre.setZero();

    /////////////////////////////////
    // Cauchy stress tensor

    fact1 = 2.0*c1*Jm5d3 + 2.0*c2*Jm7d3*I1;
    fact2 = -c1*r2d3*Jm5d3*I1 - c2*r4d3*Jm7d3*I2;
    fact3 = -2.0*c2*Jm7d3;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    stre[0]  = fact1*b(0,0);      stre[3]  = fact1*b(0,1);      stre[6]  = fact1*b(0,2);
    stre[1]  = fact1*b(1,0);      stre[4]  = fact1*b(1,1);      stre[7]  = fact1*b(1,2);
    stre[2]  = fact1*b(2,0);      stre[5]  = fact1*b(2,1);      stre[8]  = fact1*b(2,2);

    stre[0] += fact3*bSq(0,0);    stre[3] += fact3*bSq(0,1);    stre[6] += fact3*bSq(0,2);
    stre[1] += fact3*bSq(1,0);    stre[4] += fact3*bSq(1,1);    stre[7] += fact3*bSq(1,2);
    stre[2] += fact3*bSq(2,0);    stre[5] += fact3*bSq(2,1);    stre[8] += fact3*bSq(2,2);

    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

    /////////////////////////////////
    // material tangent tensor

    // without stress terms
    //fact1 =  r1d3*mu*Jm5d3*I1;
    //fact2 = -r2d3*mu*Jm5d3;
    //fact3 =  r2d3*fact1;

    // with stress terms
    fact1 =  2.0*c1*Jm5d3;
    fact2 =  2.0*c2*Jm7d3;
    fact3 =  r4d3*c2*Jm7d3;

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
            Cmat(i1,j1) += fact1*(dik*b(j, l) + r1d3*I1*dil*djk + r2d9*I1*dij*dkl - r2d3*b(i,j)*dkl - r2d3*b(k,l)*dij );
            Cmat(i1,j1) += fact2*(2*b(i,j)*b(k,l) + I1*dik*b(j,l) - b(i,k)*b(j,l) - b(i,l)*b(k,j) - dik*bSq(j,l)  );
            Cmat(i1,j1) += fact3*(I2*dil*djk + r4d3*I2*dkl*dij + 2.0*dij*bSq(k,l) + 2.0*dkl*bSq(i,j) - 2.0*I1*dkl*b(i,j) - 2.0*I1*dij*b(k,l) );

            if(MIXED_ELEMENT)
            {
                // pressure term
                //Cmat(i1,j1) += pres*(dij*dkl - djk*dil - djl*dik);

                Cmat(i1,j1) += pres*(dij*dkl - djk*dil);
            }
            else
            {
                // volumetric term (without stress terms)
                //Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk)+fact9*(djl*dik);

                // volumetric term (with stress terms)
                Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk);
            }
        }
    }

    return 0;
}



// Compute stress and tangent for the Gent model
int computeStressAndTangent_Gent(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r5d3 = 5.0/3.0, r7d3 = 7.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K     = 1.0/Kinv;                               // bulk modulus
    double  mu    = matData[0];                             // shear modulus
    double  Im    = matData[1];

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose(),  bbdev;

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Jm7d3 = pow(J, -r7d3);
    double  Ib = b(0,0) + b(1,1) + b(2,2);
    double  Ibbar = Jm2d3*Ib;
    double  gfact = 1.0/(1.0-(Ibbar-3.0)/Im);

    fact = Ib*r1d3;

    bbdev = b;
    bbdev(0,0) -= fact;
    bbdev(1,1) -= fact;
    bbdev(2,2) -= fact;


    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);
    double  dij, dik, dil, djk, djl, dkl;

    Cmat.setZero();
    stre.setZero();

    /////////////////////////////////
    // Cauchy stress tensor

    // Mechanical part
    fact1 = mu*Jm5d3*gfact;
    fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    stre[0] = fact1*b(0,0);    stre[3] = fact1*b(0,1);    stre[6] = fact1*b(0,2);
    stre[1] = fact1*b(1,0);    stre[4] = fact1*b(1,1);    stre[7] = fact1*b(1,2);
    stre[2] = fact1*b(2,0);    stre[5] = fact1*b(2,1);    stre[8] = fact1*b(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

    /////////////////////////////////
    // material tangent tensor

    // without stress terms
    //fact1 = -r2d3*mu*Jm5d3*gfact;
    //fact2 =   2.0*mu*Jm7d3*gfact*gfact/Im;
    //fact3 =  r1d3*mu*Jm5d3*gfact*Ib;
    //fact4 = -r2d3*mu*Jm5d3*gfact;

    // with stress terms
    fact1 = -r2d3*mu*Jm5d3*gfact;
    fact2 =   2.0*mu*Jm7d3*gfact*gfact/Im;
    fact3 =       mu*Jm5d3*gfact;


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

            // Devioteric terms, without stress terms
            //Cmat(i1,j1) += fact1*(bbdev(i,j)*dkl);
            //Cmat(i1,j1) += fact2*(bbdev(i,j)*bbdev(k,l));
            //Cmat(i1,j1) += fact3*(dik*djl+dil*djk);
            //Cmat(i1,j1) += fact4*(b(k,l)*dij);

            // Devioteric terms, with stress terms
            Cmat(i1,j1) += fact1*(bbdev(i,j)*dkl);
            Cmat(i1,j1) += fact2*(bbdev(i,j)*bbdev(k,l));
            Cmat(i1,j1) += fact3*(dik*b(j,l)-r2d3*b(k,l)*dij+r1d3*Ib*dil*djk);

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
        }
    }

    return 0;
}



// Compute stress and tangent for the Arruda-Boyce model
int computeStressAndTangent_ArrudaBoyce(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;

    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;
    double  r7d3 = 7.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K      = 1.0/Kinv;  // bulk modulus

    double  mu     = matData[0];      // shear modulus
    double  N      = matData[1]; // number of links per polymer chain 
    double  c1     = mu/2.0;
    double  c2     = mu/20.0/N;
    double  c3     = mu*11.0/1050.0/N/N;
    double  c4     = mu*19.0/7000.0/N/N/N;
    double  c5     = mu*519.0/673750.0/N/N/N/N;


    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose(),  bbdev;

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = pow(J, -r5d3);
    double  Jm7d3 = pow(J, -r7d3);
    double  Ib    = b.trace();
    double  Ibbar = Jm2d3*Ib;

    fact = Ib*r1d3;

    bbdev = b;
    bbdev(0,0) -= fact;
    bbdev(1,1) -= fact;
    bbdev(2,2) -= fact;


    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);
    double  dij, dik, dil, djk, djl, dkl;

    Cmat.setZero();
    stre.setZero();

    double a1 = c1 + 2.0*c2*Ibbar + 3.0*c3*Ibbar*Ibbar + 4.0*c4*Ibbar*Ibbar*Ibbar + 5.0*c5*Ibbar*Ibbar*Ibbar*Ibbar;
    double a2 = 2.0*c2 + 6.0*c3*Ibbar + 12.0*c4*Ibbar*Ibbar + 20.0*c5*Ibbar*Ibbar*Ibbar;

    /////////////////////////////////
    // Cauchy stress tensor

    // Mechanical part
    fact1 = 2.0*a1*Jm5d3;
    fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    stre[0] = fact1*b(0,0);    stre[3] = fact1*b(0,1);    stre[6] = fact1*b(0,2);
    stre[1] = fact1*b(1,0);    stre[4] = fact1*b(1,1);    stre[7] = fact1*b(1,2);
    stre[2] = fact1*b(2,0);    stre[5] = fact1*b(2,1);    stre[8] = fact1*b(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

    /////////////////////////////////
    // material tangent tensor

    fact1 =  4.0*Jm7d3*a2;
    fact2 = -4.0*r1d3*a1*Jm5d3;
    fact3 =      r2d3*a1*Jm5d3*Ib;
    fact4 =  4.0*r1d3*r1d3*a1*Jm5d3*Ib;

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


            // Devioteric terms, without stress terms
            Cmat(i1,j1) += fact1*(bbdev(i,j)*bbdev(k,l));
            Cmat(i1,j1) += fact2*(b(i,j)*dkl+b(k,l)*dij);
            Cmat(i1,j1) += fact3*(dik*djl+dil*djk);
            Cmat(i1,j1) += fact4*(dij*dkl);

            if(MIXED_ELEMENT)
            {
                // pressure term
                Cmat(i1,j1) += pres*(dij*dkl - djk*dil - djl*dik);

                //Cmat(i1,j1) += pres*(dij*dkl - djk*dil);
            }
            else
            {
                // volumetric term (without stress terms)
                //Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk)+fact9*(djl*dik);

                // volumetric term (with stress terms)
                Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk);
            }
        }
    }

    return err;
}



int computeStressAndTangent_Longevin8chain(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r5d3 = 5.0/3.0, r7d3 = 7.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K      = 1.0/Kinv;                              // bulk modulus
    double  mu     = matData[0];                            // shear modulus
    double  N      = matData[1];                            // factor N

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose(), bdev;

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Jm7d3 = pow(J, -r7d3);
    double  Ib    = b.trace();
    double  Ibbar = Jm2d3*Ib;

    fact = Ib*r1d3;

    bdev = b;
    bdev(0,0) -= fact;
    bdev(1,1) -= fact;
    bdev(2,2) -= fact;

    double  lrbar = sqrt(Ibbar/3.0/N);

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);
    double  dij, dik, dil, djk, djl, dkl;

    Cmat.setZero();
    stre.setZero();

    /////////////////////////////////
    // Cauchy stress tensor

    fact1 = mu*r1d3*Jm5d3*(3.0-lrbar*lrbar)/(1.0-lrbar*lrbar);
    fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    stre[0] = fact1*b(0,0);    stre[3] = fact1*b(0,1);    stre[6] = fact1*b(0,2);
    stre[1] = fact1*b(1,0);    stre[4] = fact1*b(1,1);    stre[7] = fact1*b(1,2);
    stre[2] = fact1*b(2,0);    stre[5] = fact1*b(2,1);    stre[8] = fact1*b(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

    /////////////////////////////////
    // material tangent tensor

        // without stress terms
        //fact1 =  r1d3*mu*Jm5d3*Ib;
        //fact2 = -r2d3*mu*Jm5d3;
        //fact3 =  r2d3*fact1;

        // with stress terms
        fact1 =  r1d3*mu*Jm7d3*2.0/(3*N*(1.0-lrbar*lrbar)*(1.0-lrbar*lrbar));
        fact2 = -r2d3*r1d3*mu*Jm5d3*(3.0-lrbar*lrbar)/(1.0-lrbar*lrbar);
        fact3 =  r1d3*mu*Jm5d3*(3.0-lrbar*lrbar)/(1.0-lrbar*lrbar);

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
              Cmat(i1,j1) += fact1*(bdev(i,j)*bdev(k,l));
              Cmat(i1,j1) += fact2*(bdev(i,j)*dkl);
              Cmat(i1,j1) += fact3*(dik*b(j,l)+r1d3*djk*dil*Ib-r2d3*dij*b(k,l));

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
            }
        }

    return err;
}




int computeStressAndTangent_Viscoelasticity_Model1(vector<double>& data_viscoelastic, int gp, MatrixXd& F, VectorXd& td, InternalVariables& ivar, double dt, VectorXd&  stre, MatrixXd&  Cmat)
{
    double  fact, fact1, fact2, fact3, fact4, fact5;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0, r2d9 = 2.0/9.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  ii, jj, i1, j1, i, j, k, l, err=0;

    double  J = F.determinant();

    // Transpose of F
    MatrixXd  Ft = F.transpose();

    // Left Cauchy-Green tensor;
    //MatrixXd  b = F*Ft;

    // Right Cauchy-Green tensor;
    MatrixXd  C = Ft*F, Cinv = C.inverse(), CbarInv = pow(J,r2d3)*Cinv;

    double  Jm5d3 = pow(J, -5.0/3.0);
    double  dij, dik, dil, djk, djl, dkl;


    double af    =  td(2);
    double am    =  td(3);
    double gamm  =  td(4);

    //cout << af << '\t' << am << '\t' << gamm << endl;

    int  numSeries = (int) data_viscoelastic[0];  // number of terms in the Prony series

    //cout << "numSeries = " << numSeries << endl;

    for(int seriescount=0; seriescount<numSeries; seriescount++)
    {
        jj = 2*seriescount;

        double  mu   = data_viscoelastic[1+jj];    // modulus
        double  tau  = data_viscoelastic[1+jj+1];  // relaxation time


        VectorXd  AnVec = ivar.varPrev.block<9,1>(9*seriescount,gp), AVec(9), AdotVec(9);
        VectorXd  AdotnVec = ivar.varDotPrev.block<9,1>(9*seriescount,gp);

        MatrixXd  AnMat(3,3), AdotnMat(3,3);
        vector2matrix(AnVec, AnMat);
        vector2matrix(AdotnVec, AdotnMat);

        //MatrixXd  AMat = ((dt/tau)*CbarInv+AnMat)/(1.0+dt/tau);

        MatrixXd  AMat = ( (gamm*dt/am/tau)*CbarInv + (1.0 -(1.0-af)*gamm*dt/am/tau)*AnMat - ((gamm-am)*dt/am/gamm)*AdotnMat)/(1.0+af*gamm*dt/am/tau);

        //printMatrix(AMat);

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

    return err;
}










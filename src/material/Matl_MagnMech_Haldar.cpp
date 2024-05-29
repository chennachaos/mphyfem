#include "Matl_MagnMech_Haldar.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"


using namespace std;


Matl_MagnMech_Haldar::Matl_MagnMech_Haldar()
{

}


Matl_MagnMech_Haldar::~Matl_MagnMech_Haldar()
{

}



int Matl_MagnMech_Haldar::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    //Base material model is Mooney-Rivlin
    // U(J) is selected based on the input value for 'Utype'

    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r4d3 = 4.0/3.0,  r2d9 = 2.0/9.0,  r5d3 = 5.0/3.0,  r7d3 = 7.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    int     Utype  = int(matData[0]); // volumetric energy function
    double  K      = 1.0/matData[1];  // matData[1] is the inverse of bulk modulus
    //double  c1     = matData[2];                           // coefficient 1
    //double  c2     = matData[3];                           // coefficient 2

    double    c1  = 500.0, c2 = 11.9e3;                     // g/mm/s^2
    double    zeta2 = 1.0e6;                                // g/mm/s^2
    double    eta2 = 4.91e-5;                               // mm/A
    double    gamma = -0.0737;                              // -
    double    beta = 2.45e-6;                               // mm^2/A^2
    double    params_nu[] = {1.1, -0.6, 1.2};               // milliTesla
    double     params_k[] = {0.017, 8.3e-3, 0.017};         // mm/A



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

    return err;
}


int Matl_MagnMech_Haldar::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre, VectorXd&  magnDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, InternalVariables& ivar, int gp, double dt)
{
    //Base material model is Mooney-Rivlin
    // U(J) is selected based on the input value for 'Utype'

    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r4d3 = 4.0/3.0,  r2d9 = 2.0/9.0,  r5d3 = 5.0/3.0,  r7d3 = 5.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    int     Utype  = int(matData[0]); // volumetric energy function
    double  K      = 1.0/matData[1];  // matData[1] is the inverse of bulk modulus
    //double  c1     = matData[2];                           // coefficient 1
    //double  c2     = matData[3];                           // coefficient 2

    double    c1  = 500.0, c2 = 11.9e3;                     // g/mm/s^2
    double    zeta2 = 1.0e6;                                // g/mm/s^2
    double    eta2 = 4.91e-5;                               // mm/A
    double    gamma = -0.0737;                              // -
    double    beta = 2.45e-6;                               // mm^2/A^2
    double    params_nu[] = {1.1, -0.6, 1.2};               // milliTesla
    double     params_k[] = {0.017, 8.3e-3, 0.017};         // mm/A
    double    mu0 = 1.25663706143592;                       // g.mm/s^2/A^2



    // Left Cauchy-Green tensor;
    MatrixXd  C = F.transpose()*F;
    MatrixXd  b = F*F.transpose();
    MatrixXd  bSq = b*b;

    double  J     = F.determinant();
    double  Jm1d3 = pow(J, -r1d3);
    double  Jp1d3 = pow(J,  r1d3);
    double  Jp2d3 = pow(J,  r2d3);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Jm7d3 = pow(J, -r7d3);
    double  I1 = b.trace();
    double  I2 = 0.5*(I1*I1 - bSq.trace());

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);
    double  dij, dik, dil, djk, djl, dkl;


    VectorXd  magnFieldOC = F.transpose()*magnField;

    double  I4 = magnFieldOC[0]*magnFieldOC[0] + magnFieldOC[1]*magnFieldOC[1] + magnFieldOC[2]*magnFieldOC[2];

    MatrixXd  CbarInv = Jp2d3*C.inverse();
    MatrixXd  HcrossH = magnField*magnField.transpose();

    double  magnFieldMagnSq = magnField[0]*magnField[0] + magnField[1]*magnField[1] + magnField[2]*magnField[2];
    double  I5tilde = Jp2d3*magnFieldMagnSq;


    double  dPsiMinf_dI4   = zeta2*eta2/(2.0*sqrt(I4)*(1.0+eta2*eta2*I4));
    double  d2PsiMinf_dI42 = -zeta2*eta2*(1.0+3.0*eta2*eta2*I4)/(4.0*pow(I4,1.5)*(1.0+eta2*eta2*I4)*(1.0+eta2*eta2*I4));

    double  dPsiMinf_dI2bar = zeta2*atan(eta2*I4);

    double  dPsiMinf_dI5tilde   = -0.5*mu0*(gamma/(1.0+beta*I5tilde) - 1.0);
    double  d2PsiMinf_dI5tilde2 =  0.5*mu0*gamma*beta/(1.0+beta*I5tilde)/(1.0+beta*I5tilde);


    Cmat.setZero();
    Bmat.setZero();
    Amat.setZero();
    stre.setZero();
    magnDisp.setZero();

    /////////////////////////////////
    // Cauchy stress tensor

    fact1 = 2.0*c1*Jm5d3 + 2.0*c2*Jm7d3*I1 + dPsiMinf_dI2bar*2.0*Jm7d3*I1;
    fact2 = -c1*r2d3*Jm5d3*I1 - c2*r4d3*Jm7d3*I2 - dPsiMinf_dI2bar*2.0*Jm7d3*r2d3*I2;
    fact3 = -2.0*c2*Jm7d3 - dPsiMinf_dI2bar*2.0*Jm7d3;

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


    MatrixXd    hTPh = magnField*magnField.transpose();
    MatrixXd    bh = b*magnField;

    // Magnetic part - Minf
    fact1 = -dPsiMinf_dI5tilde*2.0*Jm1d3;
    fact2 =  dPsiMinf_dI5tilde*2.0*Jm1d3*r1d3;

    stre[0]  += fact1*hTPh(0,0);      stre[3]  += fact1*hTPh(0,1);      stre[6]  += fact1*hTPh(0,2);
    stre[1]  += fact1*hTPh(1,0);      stre[4]  += fact1*hTPh(1,1);      stre[7]  += fact1*hTPh(1,2);
    stre[2]  += fact1*hTPh(2,0);      stre[5]  += fact1*hTPh(2,1);      stre[8]  += fact1*hTPh(2,2);

    stre[0]  += fact2 ;            stre[4]  += fact2 ;            stre[8]  += fact2 ;


    // Magnetic part - free-space
    addElectricPartToCauchyStress(mu0, magnField, stre);




    /////////////////////////////////
    // Magnetic displacement vector

    // free-space contribution
    magnDisp[0] = mu0*magnField[0];
    magnDisp[1] = mu0*magnField[1];
    magnDisp[2] = mu0*magnField[2];


    /////////////////////////////////
    // material tangent tensor

    c2 += dPsiMinf_dI2bar;

    // with stress terms
    fact1 =  2.0*c1*Jm5d3;
    fact2 =  2.0*c2*Jm7d3;
    fact3 =  r4d3*c2*Jm7d3;

    fact4 = d2PsiMinf_dI5tilde2*4.0*Jp1d3;
    fact5 = dPsiMinf_dI5tilde*2.0*Jm1d3;


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

            // Devioteric terms, with stress terms
            Cmat(i1,j1) += fact1*(dik*b(j, l) + r1d3*I1*dil*djk + r2d9*I1*dij*dkl - r2d3*b(i,j)*dkl - r2d3*b(k,l)*dij );
            Cmat(i1,j1) += fact2*(2*b(i,j)*b(k,l) + I1*dik*b(j,l) - b(i,k)*b(j,l) - b(i,l)*b(k,j) - dik*bSq(j,l)  );
            Cmat(i1,j1) += fact3*(I2*dil*djk + r4d3*I2*dkl*dij + 2.0*dij*bSq(k,l) + 2.0*dkl*bSq(i,j) - 2.0*I1*dkl*b(i,j) - 2.0*I1*dij*b(k,l) );

            // Minf
            Cmat(i1,j1) += fact4*( (r1d3*dij*magnFieldMagnSq - magnField(i)*magnField(j)) * (r1d3*dkl*magnFieldMagnSq - magnField(k)*magnField(l)));

            Cmat(i1,j1) += fact5*( -r1d3*djk*dil*magnFieldMagnSq + djk*magnField(i)*magnField(l) + r2d9*dkl*dij*magnFieldMagnSq - r2d3*dkl*magnField(i)*magnField(j) );


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

        for(k=0; k<3; k++)
        {
            dik = KronDelta(i,k);
            djk = KronDelta(j,k);


            // Minf
            Bmat(i1,k) -= fact4*( (r1d3*dij*magnFieldMagnSq - magnField(i)*magnField(j)) * magnField(k) );

            Bmat(i1,k) -= fact5*( r2d3*dij*magnField(k) - dik*magnField(j) - djk*magnField(i) );

        }

    }

    // Magnetic part - free-space
    addMagneticPartToMaterialTensor(mu0, magnField, Cmat);


    /////////////////////////////////
    // coupling tensor - free-space

    addToCouplingTensor(mu0, magnField, Bmat);

    ////////////////////////////////////////////
    // magnetic permeability tensor - free-space

    Amat(0,0) += mu0;  Amat(1,1) += mu0;  Amat(2,2) += mu0;

    fact1 = -d2PsiMinf_dI42*4.0/J;
    fact2 = -dPsiMinf_dI4*2.0/J;
    fact3 = -d2PsiMinf_dI5tilde2*4.0*Jp1d3;
    fact4 = -dPsiMinf_dI5tilde*2.0*Jm1d3;


    for(i=0; i<3; i++)
    {
        magnDisp(i) += fact2*bh(i);
        magnDisp(i) += fact4*magnField(i);

        for(j=0; j<3; j++)
        {
            Amat(i,j) += fact1*bh(i)*bh(j);
            Amat(i,j) += fact2*b(i,j);
            Amat(i,j) += fact3*magnField(i)*magnField(j);
        }

        Amat(i,i) += fact4;
    }


    return err;
}






int Matl_MagnMech_Haldar::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;

    int     Utype  = int(matData[0]); // volumetric energy function
    double  K      = 1.0/matData[1];  // matData[1] is the inverse of bulk modulus
    double  mu     = matData[2];      // shear modulus

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose();

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  I1 = b(0,0) + b(1,1) + b(2,2);

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);

    double  fact1 = mu*Jm5d3;
    double  fact2 = -fact1*r1d3*I1;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    stre.setZero();

    stre[0] = fact1*b(0,0);    stre[3] = fact1*b(0,1);    stre[6] = fact1*b(0,2);
    stre[1] = fact1*b(1,0);    stre[4] = fact1*b(1,1);    stre[7] = fact1*b(1,2);
    stre[2] = fact1*b(2,0);    stre[5] = fact1*b(2,1);    stre[8] = fact1*b(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

    return 0;
}






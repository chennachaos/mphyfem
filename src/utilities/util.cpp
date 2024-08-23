
#include "headersBasic.h"
#include "util.h"
#include "myConstants.h"

using namespace std;
using namespace Eigen;




bool doubleGreater(double left, double right, bool orequal)
{
  if (fabs(left - right) < DBL_EPSILON) 
    return (orequal);

  return (left > right);
}


bool doubleLess(double left, double right, bool orequal)
{
  if (fabs(left - right) < DBL_EPSILON)
    return (orequal);

  return (left < right);
}



void TensorProduct(MatrixXd& A, MatrixXd& B, MatrixXd& C)
{
    int ii, jj, r1, r2, c1, c2;

    r1 = A.rows();
    c1 = A.cols();

    r2 = B.rows();
    c2 = B.cols();

    C.resize(r1*r2, c1*c2);
    for(ii=0;ii<r1;ii++)
    {
      for(jj=0;jj<c1;jj++)
      {
        C.block(r2*ii,c2*jj,r2,c2) = A(ii,jj)*B;
      }
    }

    return;
}



void printMatrix(MatrixXd& AA)
{
    int ii, jj;
    printf("\n\n");
    for(ii=0;ii<AA.rows();ii++)
    {
        for(jj=0;jj<AA.cols();jj++)
           printf("\t%14.8f", AA(ii,jj));
        printf("\n");
    }
    printf("\n\n");

    return;
}


void printSparseMatrixToFile(SparseMatrixXd& mat, char* filename)
{
    FILE * pFile;
    pFile = fopen(filename,"w");

    for(int k=0; k<mat.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(mat,k); it; ++it)
      {
        fprintf(pFile, "%9d \t %9d \t %20.16f \n", it.row(), it.col(), it.value());
      }
    }
    fclose(pFile);

    return;
}


void printVector(VectorXd& AA)
{
    printf("\n\n");
    for(int ii=0;ii<AA.rows();ii++)
        printf("\t%6d\t%12.8f\n", ii, AA(ii));
    printf("\n\n");

   return;
}



void printVector(vector<string>&  vec)
{
    printf("\n\n");
    for(int ii=0;ii<vec.size();ii++)
        printf("\t%s ", vec[ii].c_str());
    printf("\n\n");

   return;
}


void printVector(vector<int>&  vec)
{
    printf("\n\n");
    for(int ii=0;ii<vec.size();ii++)
        printf("\t%6d ", vec[ii]);
    printf("\n\n");

   return;
}


void printVector(vector<double>&  vec)
{
    printf("\n\n");
    for(int ii=0;ii<vec.size();ii++)
        printf("\t%12.6f ", vec[ii]);
    printf("\n\n");

   return;
}



void printVector(double* data, int nn)
{
    printf("\n\n");
    for(int ii=0;ii<nn;ii++)
      printf("\t%6d\t%12.8f\n", ii, data[ii]);
    printf("\n\n");

   return;
}



double  HeavisideFunction(double uu, double ee)
{
   double  val=0.0;
   if(uu < -ee)
     return 0.0;
   else if(abs(uu) <= ee)
     return  (0.5 + 0.5*tanh(uu/ee/ee));
     //return  (0.5 + 0.5*uu/ee - sin(PI*uu/ee)/PI);
   else
     return  1.0;
}




TISFLUID convert2TISFLUID(string& str)
{
    if(str == "STEADY")        return  TISFLUID::STEADY;
    else if(str == "BDF1")
    {
      //cout << " BDF1 " << endl;
      return  TISFLUID::BDF1;
    }
    else if(str == "BDF2")     return  TISFLUID::BDF2;
    else if(str == "BDF3")     return  TISFLUID::BDF3;
    else if(str == "BDF4")     return  TISFLUID::BDF4;
    else if(str == "Midpoint") return  TISFLUID::Midpoint;
    else if(str == "Galpha")   return  TISFLUID::Galpha;
    else
    {
        throw runtime_error("Option not available in convert2TISFLUID");
    }
}



void SetTimeParametersFluid(string& tis, double rho, double dt, VectorXd& td)
{
  td.setZero();

  double  alpf, alpm, beta, gamm;

  td[0] = dt;

  //cout << " convert2TISFLUID(tis) = " << convert2TISFLUID(tis) << endl;

  switch(convert2TISFLUID(tis))
  {
      case  TISFLUID::STEADY: // quasi static

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

       break;

      case  TISFLUID::Midpoint: // generalised midpoint rule

            alpf = 1.0/(1.0 + rho);
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpf*dt; //af*dt
            td[6]  = alpm/gamm;
            td[7]  = 1.0-alpm/gamm;
            td[8]  = 1.0/dt;

            td[9]  = 1.0/dt;  // v_{n+1}
            td[10] = -td[9];  // v_n
            td[11] = 0.0;     // v_{n-1}
            td[12] = 0.0;     // v_{n-2}
            td[15] = 0.0;     // a_n

        break;

      case  TISFLUID::Galpha: // generalised alpha-method

            alpf = 1.0/(1.0 + rho);
            alpm = 0.5*(3.0 - rho)/(1.0 + rho);
            gamm = 0.5 + alpm - alpf;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpf*dt;
            td[6]  = alpm/gamm;
            td[7]  = 1.0-alpm/gamm;
            td[8]  = alpm/gamm/dt;

            td[9]  = 1.0/gamm/dt;     // v_{n+1}
            td[10] = -td[9];          // v_n
            td[11] = 0.0;             // v_{n-1}
            td[12] = 0.0;             // v_{n-2}
            td[15] = 1.0 - 1.0/gamm;  // a_n

         break;

      case  TISFLUID::BDF1: // Backward Euler or BDF1

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  1.0/dt;    // v_{n+1}
            td[10] = -1.0/dt;    // v_n
            td[11] =  0.0;       // v_{n-1}
            td[12] =  0.0;       // v_{n-2}
            td[13] =  0.0;       // v_{n-3}
            td[15] =  0.0;       // a_n

         break;

      case  TISFLUID::BDF2: // BDF2

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  1.5/dt;    // v_{n+1}
            td[10] = -2.0/dt;    // v_n
            td[11] =  0.5/dt;    // v_{n-1}
            td[12] =  0.0;       // v_{n-2}
            td[13] =  0.0;       // v_{n-3}
            td[15] =  0.0;       // a_n

         break;

      case  TISFLUID::BDF3: // BDF3

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  (11.0/6.0)/dt;    // v_{n+1}
            td[10] = -3.0/dt;           // v_n
            td[11] =  1.5/dt;           // v_{n-1}
            td[12] = -(1.0/3.0)/dt;     // v_{n-2}
            td[13] =  0.0;              // v_{n-3}
            td[15] =  0.0;              // a_n

         break;

      case  TISFLUID::BDF4: // BDF4

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  (25.0/12.0)/dt;    // v_{n+1}
            td[10] = -4.0/dt;           // v_n
            td[11] =  3.0/dt;           // v_{n-1}
            td[12] = -(4.0/3.0)/dt;     // v_{n-2}
            td[13] =  (1.0/4.0)/dt;     // v_{n-3}
            td[15] =  0.0;              // a_n

         break;

      default:
            cerr << " SetTimeParametersFluid ... invalid value of tis!" << endl;

         break;
  }

  return;
}




TISSOLID convert2TISSOLID(string& str)
{
    if( (str == "static") || (str == "STATIC") || (str == "STEADY") )        return  TISSOLID::STATIC;
    else if(str == "BDF1")     return  TISSOLID::BDF1;
    else if(str == "Newmark")  return  TISSOLID::Newmark;
    else if(str == "HHTalpha") return  TISSOLID::HHTalpha;
    else if(str == "CHalpha")  return  TISSOLID::CHalpha;
    else if(str == "KDPalpha") return  TISSOLID::KDPalpha;
    else
    {
        throw runtime_error("Option not available in convert2TISSOLID");
    }
}




void SetTimeParametersSolid(string& tis, double rho, double dt, VectorXd& td)
{
  td.setZero();

  double  alpf, alpm, beta, gamm;

  td[0] = dt;

  switch( convert2TISSOLID(tis) )
  {
      case  TISSOLID::STATIC: // quasi static

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;
            td[7]  = alpf;

            // velocity is used as the primary variable for
            // solid dynamics problem
            // td[10] is the multiplication when converting from
            // displacement based formulation to velocity based formulation
            // It is set to ONE for static problem
            td[10] = 1.0;

            // for displacement-based formulation
            td[10] = 0.0;

      break;

      case  TISSOLID::BDF1: // Backward Euler

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = 1.0/dt/dt;
            td[6]  = 1.0/dt;
            td[7]  = alpf;

            //displacement as the primary variable
            //v_{n+1}    = td[10]*d_{n+1} + td[11]*d_n + td[12]*v_n + td[13]*a_n + td[14]*ddot_n;
            //a_{n+1}    = td[15]*d_{n+1} + td[16]*d_n + td[17]*v_n + td[18]*a_n + td[19]*ddot_n;
            //ddot_{n+1} = td[20]*d_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n;

            td[10] = 1.0/dt;    // d_{n+1}
            td[11] = -td[10];   // d_n
            td[12] = 0.0;       // v_n
            td[13] = 0.0;       // a_n
            td[14] = 0.0;       // ddot_n

            td[15] = 1.0/dt/dt; // d_{n+1}
            td[16] = -td[15];   // d_n
            td[17] = -1.0/dt;   // vn
            td[18] = 0.0;       // an
            td[19] = 0.0;       // ddot_n

            td[20] = 1.0/dt;    // d_{n+1}
            td[21] = -td[20];   // d_n
            td[22] = 0.0;       // v_n
            td[23] = 0.0;       // a_n
            td[24] = 0.0;       // ddot_n

            //velocity as the primary variable
            //d_{n+1}    = td[20]*v_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n ;
            //a_{n+1}    = td[25]*v_{n+1} + td[26]*d_n + td[27]*v_n + td[28]*a_n + td[29]*ddot_n ;
            //ddot_{n+1} = td[30]*v_{n+1} + td[31]*d_n + td[32]*v_n + td[33]*a_n + td[34]*ddot_n;

            td[40] = dt;        // v_{n+1}
            td[41] = 1.0;       // d_n
            td[42] = 0.0;       // v_n
            td[43] = 0.0;       // a_n
            td[44] = 0.0;       // ddot_n

            td[45] = 1.0/dt;    // v_{n+1}
            td[46] = 0.0;       // d_n
            td[47] = -td[45];   // v_n
            td[48] = 0.0;       // a_n
            td[49] = 0.0;       // ddot_n

            td[50] = 1.0;       // v_{n+1}
            td[51] = 0.0;       // d_n
            td[52] = 0.0;       // v_n
            td[53] = 0.0;       // a_n
            td[54] = 0.0;       // ddot_n

      break;

      case  TISSOLID::CHalpha: // generalised alpha-method --- CH-alpha one

            alpm = (2.0-rho)/(rho+1.0);
            alpf = 1.0/(rho+1.0);

            gamm = 0.5 + alpm - alpf;
            beta = 0.25*(1.0+alpm-alpf)*(1.0+alpm-alpf);

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpm/beta/dt/dt;
            td[6]  = alpf*gamm/beta/dt;
            td[7]  = alpf;

            //displacement as the primary variable
            //v_{n+1}  = td[10]*d_{n+1} + td[11]*d_n + td[12]*v_n + td[13]*a_n + td[14]*ddot_n;
            //a_{n+1}  = td[15]*d_{n+1} + td[16]*d_n + td[17]*v_n + td[18]*a_n + td[19]*ddot_n;

            td[10] = gamm/beta/dt;           // d_{n+1}
            td[11] = -td[10];                // d_n
            td[12] = 1.0-gamm/beta;          // v_n
            td[13] = dt*(1.0-gamm/2.0/beta); // a_n
            td[14] = 0.0;                    // ddot_n

            td[15] = 1.0/beta/dt/dt;         // d_{n+1}
            td[16] = -td[15];                // d_n
            td[17] = -1.0/beta/dt;           // v_n
            td[18] = 1.0-1.0/2.0/beta;       // a_n
            td[19] = 0.0;                    // ddot_n

            //velocity as the primary variable
            //d_{n+1}  = td[20]*v_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n ;
            //a_{n+1}  = td[25]*v_{n+1} + td[26]*d_n + td[27]*v_n + td[28]*a_n + td[29]*ddot_n ;

            td[40] = dt*beta/gamm;                     // v_{n+1}
            td[41] = 1.0;                              // d_n
            td[42] = dt*(gamm-beta)/gamm;              // v_n
            td[43] = dt*dt*(gamm-2.0*beta)/(2.0*gamm); // a_n
            td[44] = 0.0;                              // ddot_n

            td[45] = 1.0/(gamm*dt);                    // v_{n+1}
            td[46] = 0.0;                              // d_n
            td[47] = -td[45];                          // v_n
            td[48] = (gamm-1.0)/gamm;                  // a_n
            td[49] = 0.0;                              // ddot_n

      break;

      case  TISSOLID::KDPalpha: // generalised alpha-method --- state-space formulation

            alpf = 1.0/(1.0 + rho);
            alpm = 0.5*(3.0 - rho)/(1.0 + rho);
            gamm = 0.5 + alpm - alpf;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = (alpm*alpm)/(alpf*gamm*gamm*dt*dt);
            td[6]  = alpm/gamm/dt;
            td[7]  = alpf;

            //displacement as the primary variable
            //v_{n+1}    = td[10]*d_{n+1} + td[11]*d_n + td[12]*v_n + td[13]*a_n + td[14]*ddot_n;
            //a_{n+1}    = td[15]*d_{n+1} + td[16]*d_n + td[17]*v_n + td[18]*a_n + td[19]*ddot_n;
            //ddot_{n+1} = td[20]*d_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n;

            td[10] = alpm/(alpf*gamm*dt);             // d_{n+1}
            td[11] = -td[10];                         // d_n
            td[12] = (alpf-1.0)/alpf;                 // v_n
            td[13] = 0.0;                             // a_n
            td[14] = (gamm-alpm)/alpf/gamm;           // ddot_n

            td[15] = alpm/(alpf*gamm*gamm*dt*dt);     // d_{n+1}
            td[16] = -td[15];                         // d_n
            td[17] = -1.0/(alpf*gamm*dt);             // v_n
            td[18] = (gamm-1.0)/gamm;                 // a_n
            td[19] = (gamm-alpm)/(alpf*gamm*gamm*dt); // ddot_n

            td[20] = 1.0/(gamm*dt);                   // d_{n+1}
            td[21] = -td[20];                         // d_n
            td[22] = 0.0;                             // v_n
            td[23] = 0.0;                             // a_n
            td[24] = (gamm-1.0)/gamm;                 // ddot_n

            //velocity as the primary variable
            //d_{n+1}    = td[20]*v_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n ;
            //a_{n+1}    = td[25]*v_{n+1} + td[26]*d_n + td[27]*v_n + td[28]*a_n + td[29]*ddot_n ;
            //ddot_{n+1} = td[30]*v_{n+1} + td[31]*d_n + td[32]*v_n + td[33]*a_n + td[34]*ddot_n;

            td[40] = alpf*gamm*dt/alpm;           // v_{n+1}
            td[41] = 1.0;                         // d_n
            td[42] = (1.0-alpf)*gamm*dt/alpm;     // v_n
            td[43] = 0.0;                         // a_n
            td[44] = (alpm-gamm)*dt/alpm;        // ddot_n

            td[45] = 1.0/gamm/dt;                 // v_{n+1}
            td[46] = 0.0;                         // d_n
            td[47] = -td[45];                     // v_n
            td[48] = 1.0-1.0/gamm;                // v_n
            td[49] = 0.0;                         // ddot_n

            td[50] = alpf/alpm;                   // v_{n+1}
            td[51] = 0.0;                         // d_n
            td[52] = (1.0-alpf)/alpm;             // v_n
            td[53] = 0.0;                         // a_n
            td[54] = (alpm-1.0)/alpm;             // ddot_n

      break;

      case TISSOLID::Newmark: // Newmark-beta

            alpm = 1.0;
            alpf = 1.0;

            gamm = 0.5;
            beta = 0.5*gamm;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpm/beta/dt/dt;
            td[6]  = alpf*gamm/beta/dt;
            td[7]  = alpf;

            //displacement as the primary variable
            //v_{n+1}  = td(10)*d_{n+1} + td(11)*d_n + td(12)*v_n + td(13)*a_n + td(14)*ddot_n;
            //a_{n+1}  = td(15)*d_{n+1} + td(16)*d_n + td(17)*v_n + td(18)*a_n + td(19)*ddot_n;

            td[10] = gamm/beta/dt;           // d_{n+1}
            td[11] = -td(10);                // d_n
            td[12] = 1.0-gamm/beta;          // v_n
            td[13] = dt*(1.0-gamm/2.0/beta); // a_n
            td[14] = 0.0;                    // ddot_n

            td[15] = 1.0/beta/dt/dt;         // d_{n+1}
            td[16] = -td[15];                // d_n
            td[17] = -1.0/beta/dt;           // v_n
            td[18] = 1.0-1.0/2.0/beta;       // a_n
            td[19] = 0.0;                    // ddot_n

            //velocity as the primary variable
            //d_{n+1}  = td(20)*v_{n+1} + td(21)*d_n + td(22)*v_n + td(23)*a_n + td(24)*ddot_n ;
            //a_{n+1}  = td(25)*v_{n+1} + td(26)*d_n + td(27)*v_n + td(28)*a_n + td(29)*ddot_n ;

            td[40] = dt*beta/gamm;                     // v_{n+1}
            td[41] = 1.0;                              // d_n
            td[42] = dt*(gamm-beta)/gamm;              // v_n
            td[43] = dt*dt*(gamm-2.0*beta)/(2.0*gamm); // a_n
            td[44] = 0.0;                              // ddot_n

            td[45] = 1.0/(gamm*dt);                    // v_{n+1}
            td[46] = 0.0;                              // d_n
            td[47] = -td[45];                          // v_n
            td[48] = (gamm-1.0)/gamm;                  // a_n
            td[49] = 0.0;                              // ddot_n

      break;

      case TISSOLID::CHalphaExplicit:

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

      break;

      default:

            cerr << "SetTimeParametersSolid ... invalid value of tis!"  << endl;

      break;
  }

  return;
}




void create_vector(double start, double end, double incr, vector<double>& uuu)
{
  int steps = int (((end-start)/incr) + 1);
  uuu.resize(steps);

  uuu[0] = start;
  for(int ii=1;ii<steps;ii++)
    uuu[ii] = uuu[ii-1] + incr;

  return;
}



void map2DPointTo3DPoint(int side, myPoint& ptTemp, double val3)
{
    switch(side)
    {
        case 0:
        case 1:
                ptTemp[2] = ptTemp[1] ;
                ptTemp[1] = ptTemp[0] ;
                ptTemp[0] = val3 ;

        break;

        case 2:
        case 3:

                ptTemp[0] = ptTemp[0] ;
                ptTemp[2] = ptTemp[1] ;
                ptTemp[1] = val3 ;

        break;

        case 4:
        case 5:

                ptTemp[0] = ptTemp[0] ;
                ptTemp[1] = ptTemp[1] ;
                ptTemp[2] = val3 ;

        break;

        default :

            cout << " Invalid 'side' in map2DPointTo3DPoint in myMappings.h " << endl;
        break;
    } //switch(side)
  return;
}





double factorial(unsigned int nn)
{
  if(nn == 0 || nn == 1)
    return 1.0;
  else
  {
    double result=1.0;
    for(unsigned int ii=1;ii<=nn;ii++)
      result *= ii;

    return result;
  }
}


double Bin(unsigned int m, unsigned int n)
{
  if((m == n) || (n == 0))
    return 1;
  else if(n > m)
    return 0;
  else
  {
    double num=1.0;
    int decr = m;
    do
    {
      num = num * decr;
      decr--;
    }while(decr > (m-n));

    return num/factorial(n);
  }
}



//RED, BLUE, GREEN, YELLOW, CYAN, MAGENTA, LIGHTBLUE, WHITE, BLACK

void  getColorValue(int col, double* color)
{
    switch(col)
    {
       case 0 : //RED
           color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.0;
         break;

       case 1 : //BLUE
           color[0] = 0.0;   color[1] = 0.0;   color[2] = 1.0;
         break;

       case 2 : //GREEN
           color[0] = 0.0;   color[1] = 1.0;   color[2] = 0.0;
         break;

       case 3 : //YELLOW
           color[0] = 1.0;   color[1] = 1.0;   color[2] = 0.0;
         break;

       case 4 : //CYAN
           color[0] = 0.0;   color[1] = 1.0;   color[2] = 1.0;
         break;

       case 5 : //MAGENTA
           color[0] = 1.0;   color[1] = 0.0;   color[2] = 1.0;
         break;

       case 6 : //LIGHT BLUE
           color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.0;
         break;

       case 7 : //WHILE
           color[0] = 1.0;   color[1] = 1.0;   color[2] = 1.0;
         break;

       case 8 : //BLACK
           color[0] = 0.0;   color[1] = 0.0;   color[2] = 0.0;
         break;

    }

  return;
}




double dotProductVecs(double* vec1, double* vec2, int N)
{
   double  val=0.0;
   
   for(int ii=0;ii<N;ii++)
     val += (vec1[ii] * vec2[ii]);

   return val;
}



void  split_by_whitespace(const string& input_string, vector<string>& output_string)
{
  string result;

  size_t  len_s = input_string.length();
  //cout << input_string << endl;
  //cout << "length = " << input_string.length() << endl;
  output_string.clear();
  for(size_t count=0; count<=len_s; ++count)
  {
    //cout << count << '\t' << input_string[count] << endl;
    if( (input_string[count] == ' ') || (count == len_s) )
    {
      output_string.push_back(result);
      result.clear();
    }
    else
    {
      result += input_string[count];
    }
  }
  //cout << " done stripping " << endl;

  return;
}






int  compute_transformations_flat_shell(vector<myPoint>& coordsGlobal, MatrixXd& RotMat)
{
  //
  // ylocal
  // ^
  // |
  // 4------3
  // |      |
  // |      |
  // 1------2--->xlocal
  //
  //

  // base vectors in global coordinate system
  myPoint EX, EY, EZ;
  EX << 1.0, 0.0, 0.0;
  EY << 0.0, 1.0, 0.0;
  EZ << 0.0, 0.0, 1.0;

  // vectors in local coordinate system
  // line 1-2
  myPoint  v21 = coordsGlobal[1] - coordsGlobal[0];
  // line 1-4
  myPoint  v41 = coordsGlobal[3] - coordsGlobal[0];

  // normalise
  v21 = v21/v21.norm();
  v41 = v41/v41.norm();

  myPoint exl = v21;
  myPoint ezl = exl.cross(v41);
  myPoint eyl = ezl.cross(exl);

  // transformation matrix from global to local
  RotMat.resize(3,3);
  RotMat.setZero();
  RotMat(0,0) = exl.dot(EX);    RotMat(0,1) = exl.dot(EY);    RotMat(0,2) = exl.dot(EZ);
  RotMat(1,0) = eyl.dot(EX);    RotMat(1,1) = eyl.dot(EY);    RotMat(1,2) = eyl.dot(EZ);
  RotMat(2,0) = ezl.dot(EX);    RotMat(2,1) = ezl.dot(EY);    RotMat(2,2) = ezl.dot(EZ);

  return 0;
}



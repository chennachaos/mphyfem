#ifndef ViscoElasticModels_h
#define ViscoElasticModels_h



inline  double  compute_viscosity_CarreauYasuda(double grad[][2])
{
    //double  n = 0.22, Cu = 1.0, a = 1.25, lambda = 1.902;
    //double  mu0 = 0.056;
    //double  muInf = 0.00345;

    double  mu0 = 0.01;
    double  muInf = 0.01*mu0;
    double  n = 0.5, Cu = 1.0, a = 2.0, lambda = Cu;

    double  gammaDot = grad[0][0]*grad[0][0]+grad[1][1]*grad[1][1]+0.5*(grad[0][1]+grad[1][0])*(grad[0][1]+grad[1][0]);

    gammaDot = sqrt(2.0*gammaDot);

    //double  mu = muInf + (mu0-muInf)*pow(1.0+pow(lambda*gammaDot, a), (n-1)/a);
    double  mu = 0.01;

//    cout << I2 << '\t' << mu << endl;

    return mu;
}


#endif

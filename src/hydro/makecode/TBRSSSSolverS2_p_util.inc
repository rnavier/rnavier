{
    //Inputs: model, and aux unevaluated
    //Outputs: aux evaluated
    using namespace TBRSSSConsts;
    double u[4] ;

    // Evaluate e and n 
    double g33 = TBRSSS::g33(tau) ;
    u[1] = gsl_vector_get(x,0);
    u[2] = gsl_vector_get(x,1);
    u[3] = gsl_vector_get(x,2);
    u[0] = sqrt(1.+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
    double e = (qs[0]*u[0] - qs[1]*u[1] - qs[2]*u[2] - qs[3]*u[3])/u[0]/g33;
    double n =  qs[kNb]/u[0]/g33;

    // Determine the the transport coefficients and entropy
    double vc[20] ; //20 max size
    model->setauxIdeal(e, n, &u[1], aux) ;
    model->getIMVSC(aux, vc) ; 
    double s = model->getS(aux) ;

    // Solve for piIJ and fill up the fill auxiliary
    double C1=1./(1.+aimdt*vc[0]/u[0]);
    double C2=1./(s*g33*u[0]);
    double piIJ[NPi33] ;
    piIJ[0] = qs[kPiIJ + 0]*C1*C2;
    piIJ[1] = qs[kPiIJ + 1]*C1*C2;
    piIJ[2] = qs[kPiIJ + 2]*C1*C2;
    piIJ[3] = qs[kPiIJ + 3]*C1*C2;
    piIJ[4] = qs[kPiIJ + 4]*C1*C2;

    model->setaux(e,n, &u[1], piIJ, aux);
}


{
using namespace TBRSSSConsts;
double u[4] ;
model->getU(aux,u) ;
double e = model->getE(aux) ;
double p = model->getP(aux); 
double cs = model->getCs(aux) ;
double cs2=cs*cs;
double pi[NPi44];
model->getPi(aux, pi);

double deUI[3];
double g33 = TBRSSS::g33(tau) ;
deUI[0] = u[1]*qs[0]/u[0]/u[0]/g33 - u[1]*e/u[0]/u[0] - qs[1]/u[0]/g33;
deUI[1] = u[2]*qs[0]/u[0]/u[0]/g33 - u[2]*e/u[0]/u[0] - qs[2]/u[0]/g33;
deUI[2] = u[3]*qs[0]/u[0]/u[0]/g33 - u[3]*e/u[0]/u[0] - qs[3]/u[0]/g33;

double dpitIuJ[3][3];
dpitIuJ[0][0] = (u[0]*pi[4]-2*u[1]*pi[1])/u[0]/u[0] - pi[1]/(e+p)*deUI[0];
dpitIuJ[0][1] = (u[0]*pi[5]-2*u[2]*pi[1])/u[0]/u[0] - pi[1]/(e+p)*deUI[1];
dpitIuJ[0][2] = (u[0]*pi[6]-2*u[3]*pi[1])/u[0]/u[0] - pi[1]/(e+p)*deUI[2];
dpitIuJ[1][0] = (u[0]*pi[5]-2*u[1]*pi[2])/u[0]/u[0] - pi[2]/(e+p)*deUI[0];
dpitIuJ[1][1] = (u[0]*pi[7]-2*u[2]*pi[2])/u[0]/u[0] - pi[2]/(e+p)*deUI[1];
dpitIuJ[1][2] = (u[0]*pi[8]-2*u[3]*pi[2])/u[0]/u[0] - pi[2]/(e+p)*deUI[2];
dpitIuJ[2][0] = (u[0]*pi[6]-2*u[1]*pi[3])/u[0]/u[0] - pi[3]/(e+p)*deUI[0] + 2*u[3]/u[0]*(u[0]*pi[1]-u[1]*pi[0])/(u[0]*u[0]-u[3]*u[3]);
dpitIuJ[2][1] = (u[0]*pi[8]-2*u[2]*pi[3])/u[0]/u[0] - pi[3]/(e+p)*deUI[1] + 2*u[3]/u[0]*(u[0]*pi[2]-u[2]*pi[0])/(u[0]*u[0]-u[3]*u[3]);
dpitIuJ[2][2] = (u[0]*pi[9]-2*u[3]*pi[3])/u[0]/u[0] - pi[3]/(e+p)*deUI[2] + 2*u[3]/u[0]*(u[0]*pi[3]-u[3]*pi[0])/(u[0]*u[0]-u[3]*u[3]);

JJ[0][0] = ( (1+cs2)*deUI[0]*u[0]*u[1] + (e+p)*(u[1]*u[1]/u[0] + u[0]) + dpitIuJ[0][0] )*g33; 
JJ[0][1] = ( (1+cs2)*deUI[1]*u[0]*u[1] + (e+p)*(u[2]*u[1]/u[0]       ) + dpitIuJ[0][1] )*g33; 
JJ[0][2] = ( (1+cs2)*deUI[2]*u[0]*u[1] + (e+p)*(u[3]*u[1]/u[0]       ) + dpitIuJ[0][2] )*g33; 

JJ[1][0] = ( (1+cs2)*deUI[0]*u[0]*u[2] + (e+p)*(u[1]*u[2]/u[0]       ) + dpitIuJ[1][0] )*g33; 
JJ[1][1] = ( (1+cs2)*deUI[1]*u[0]*u[2] + (e+p)*(u[2]*u[2]/u[0] + u[0]) + dpitIuJ[1][1] )*g33; 
JJ[1][2] = ( (1+cs2)*deUI[2]*u[0]*u[2] + (e+p)*(u[3]*u[2]/u[0]       ) + dpitIuJ[1][2] )*g33; 

JJ[2][0] = ( (1+cs2)*deUI[0]*u[0]*u[3] + (e+p)*(u[1]*u[3]/u[0]       ) + dpitIuJ[2][0] )*g33; 
JJ[2][1] = ( (1+cs2)*deUI[1]*u[0]*u[3] + (e+p)*(u[2]*u[3]/u[0]       ) + dpitIuJ[2][1] )*g33; 
JJ[2][2] = ( (1+cs2)*deUI[2]*u[0]*u[3] + (e+p)*(u[3]*u[3]/u[0] + u[0]) + dpitIuJ[2][2] )*g33; 
}

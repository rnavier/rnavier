/*
 * Copyright (c) 2015, Aleksas Mazeliauskas and Derek Teaney
 * All rights reserved.
 *
 * rnavier is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/rnavier/rnavier/
 */
{
complex<double>  Compile_$1, Compile_$9, Compile_$2, Compile_$3, Compile_$4,
Compile_$5, Compile_$8, Compile_$11, Compile_$29,
Compile_$14, Compile_$17, Compile_$25, Compile_$27,
Compile_$32, Compile_$42, Compile_$44, Compile_$56 ;
Compile_$1=M(0,0);
Compile_$9=sigmar*sigmar;
Compile_$2=M(0,1);
Compile_$3=M(0,2);
Compile_$4=M(0,3);
Compile_$5=M(0,4);
Compile_$8=M(1,0);
Compile_$11=M(1,1);
Compile_$29=pow(sigmar,4.);
Compile_$14=M(1,2);
Compile_$17=M(1,3);
Compile_$25=M(2,0);
Compile_$27=M(2,1);
Compile_$32=M(2,2);
Compile_$42=M(3,0);
Compile_$44=M(3,1);
Compile_$56=M(4,0);
Msm[0][0]=Compile_$1;
Msm[0][1]=Compile_$2;
Msm[0][2]=Compile_$3;
Msm[0][3]=Compile_$4;
Msm[0][4]=Compile_$5;
Msm[0][5]=M(0,5);
Msm[0][6]=M(0,6);
Msm[1][0]=Compile_$8;
Msm[1][1]=Compile_$11+Compile_$1*Compile_$9;
Msm[1][2]=Compile_$14+2.*Compile_$2*Compile_$9;
Msm[1][3]=Compile_$17+3.*Compile_$3*Compile_$9;
Msm[1][4]=4.*Compile_$4*Compile_$9+M(1,4);
Msm[1][5]=5.*Compile_$5*Compile_$9+M(1,5);
Msm[1][6]=0.;
Msm[2][0]=Compile_$25;
Msm[2][1]=Compile_$27+2.*Compile_$8*Compile_$9;
Msm[2][2]=2.*Compile_$1*Compile_$29+Compile_$32+4.*Compile_$11*Compile_$\
9;
Msm[2][3]=6.*Compile_$2*Compile_$29+6.*Compile_$14*Compile_$9+M(2,3);
Msm[2][4]=12.*Compile_$29*Compile_$3+8.*Compile_$17*Compile_$9+M(2,4);
Msm[2][5]=0.;
Msm[2][6]=0.;
Msm[3][0]=Compile_$42;
Msm[3][1]=Compile_$44+3.*Compile_$25*Compile_$9;
Msm[3][2]=6.*Compile_$29*Compile_$8+6.*Compile_$27*Compile_$9+M(3,2);
Msm[3][3]=18.*Compile_$11*Compile_$29+9.*Compile_$32*Compile_$9+M(3,3)\
+6.*Compile_$1*pow(sigmar,6.);
Msm[3][4]=0.;
Msm[3][5]=0.;
Msm[3][6]=0.;
Msm[4][0]=Compile_$56;
Msm[4][1]=4.*Compile_$42*Compile_$9+M(4,1);
Msm[4][2]=12.*Compile_$25*Compile_$29+8.*Compile_$44*Compile_$9+M(4,2)\
;
Msm[4][3]=0.;
Msm[4][4]=0.;
Msm[4][5]=0.;
Msm[4][6]=0.;
Msm[5][0]=M(5,0);
Msm[5][1]=5.*Compile_$56*Compile_$9+M(5,1);
Msm[5][2]=0.;
Msm[5][3]=0.;
Msm[5][4]=0.;
Msm[5][5]=0.;
Msm[5][6]=0.;
Msm[6][0]=M(6,0);
Msm[6][1]=0.;
Msm[6][2]=0.;
Msm[6][3]=0.;
Msm[6][4]=0.;
Msm[6][5]=0.;
Msm[6][6]=0.;
}

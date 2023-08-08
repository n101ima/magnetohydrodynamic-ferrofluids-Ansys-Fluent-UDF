#include "udf.h"


#define PI 3.14159265358979323846264338327950288419716939937510
#define e 2.71828182845904523536028747135266249775724709369995
real Mn=5.92e5;

real csch1(real);
real coth1(real);

DEFINE_SOURCE(sourcey, c, t , dS, eqn)
{

real SourceTermy, a,b, Tem , I, h, dp,alfap, Rof, Rop, Muf,xx, Mu0, alfaff, Mub, Kb, mp;
real a1,a2,a3,fyz,Hyz,Hr, Roff, Muff, M, dHy, dHz, d2Hy, d2Hz, dMy, dMz, y,z;
real L=1.25, ri=0.04957443, ro=0.04957443, tt=0.001;
real x[ND_ND]=0;

a=ro-ri-tt;
b=ro;
h=2*ri;

//Mn=5.33e6;    //it is megnetic number
dp=10e-9;
alfap=0.04;
Rof=1109;
Rop=5200;
Muf=0.0163583;



xx=0.348586;
alfaff=10.5994815e-7;
Mu0=PI*4e-7;
Mub=9.27e-24; 
Kb=1.3806503e-23;
mp=(4*Mub*PI*pow(dp,3)/(6*91.25e-30));
Roff=alfap*Rop+(1+alfap)*Rof;
I=2*PI*b*sqrt( (Mn*Roff*pow(alfaff,2))/(Mu0*xx*pow(h,2)));
//Mn=(Mu0*x*pow(Hr,2)*pow(h,2))/(Roff*pow(alfaff,2));

C_CENTROID(x,c,t);
y=x[1];
z=x[2];
Tem=C_T(c,t);
a1=(6*mp)/(PI*pow(dp,3));
a2=(Mu0*mp)/(Kb*Tem);
a3=1/a2;
fyz=pow((y-a),2)+pow((z-b),2);
Hyz=I/(2*PI*sqrt(fyz));
Hr=I/(2*PI*b);

Muff=(1+2.5*alfap)*Muf;
M=a1*(coth1(a2*Hyz)-(a3/Hyz));

dHy=(I*(a-y))/(2*PI*pow(fyz,1.5));
//d2Hy=(I*(3*pow((y-a),2)*sqrt(fyz)-pow(fyz,1.5)))/(2*PI*pow(fyz,3));
//dMy=a1*((a3*dHy/pow(Hyz,2))-(a2*pow(csch1(a2*Hyz),2)*(a-y)/(pow(fyz,1.5))));
d2Hy=I/(2*PI)*((3*pow((y-a),2)*pow(fyz,-2.5))-pow(fyz,-1.5));
dMy=a1*a2*(a-y)*pow(csch1(a2*Hyz),2)/(2*PI*pow(fyz,1.5))-(2*PI*a1*a3*(y-a))/pow(fyz,0.5);

SourceTermy=Mu0*M*dHy;

//SourceTermy=0.001;
dS[eqn]=(SourceTermy*dMy)+(Mu0*M*d2Hy);
C_UDMI(c,t,0)=Mu0;
C_UDMI(c,t,1)=M;
C_UDMI(c,t,2)=I;
C_UDMI(c,t,3)=dHy;
C_UDMI(c,t,4)=Mu0*M*dHy;
return SourceTermy;

}


DEFINE_SOURCE(sourcez, c, t , dS, eqn)
{
real SourceTermz, a,b, Tem , I, h, dp,alfap, Rof, Rop, Muf,xx, Mu0, alfaff, Mub, Kb, mp;
real a1,a2,a3,fyz,Hyz,Hr, Roff, Muff, M, dHy, dHz, d2Hy, d2Hz, dMy, dMz, y,z;
real L=1, ri=0.025, ro=0.05, tt=0.005;
real x[ND_ND]=0;



a=ro-ri-tt;
b=ro;
h=2*ri;

//Mn=5.33e6;    //it is megnetic number
dp=10e-9;
alfap=0.04;
Rof=1024;
Rop=5200;
Muf=0.00108;



xx=0.348586;
alfaff=10.5994815e-7;
Mu0=PI*4e-7;
Mub=9.27e-24; 
Kb=1.3806503e-23;
mp=(4*Mub*PI*pow(dp,3)/(6*91.25e-30));
Roff=alfap*Rop+(1+alfap)*Rof;
I=2*PI*b*sqrt( (Mn*Roff*pow(alfaff,2))/(Mu0*xx*pow(h,2)));
//Mn=(Mu0*x*pow(Hr,2)*pow(h,2))/(Roff*pow(alfaff,2));

C_CENTROID(x,c,t);
y=x[1];
z=x[2];
Tem=C_T(c,t);
a1=(6*mp)/(PI*pow(dp,3));
a2=(Mu0*mp)/(Kb*Tem);
a3=1/a2;
fyz=pow((y-a),2)+pow((z-b),2);
Hyz=I/(2*PI*sqrt(fyz));
Hr=I/(2*PI*b);

Muff=(1+2.5*alfap)*Muf;

M=a1*(coth1(a2*Hyz)-(a3/Hyz));

dHz=(I*(b-z))/(2*PI*pow(fyz,1.5));
//d2Hz=(I*(3*pow(z-b,2)*sqrt(fyz)-pow(fyz,1.5)))/(2*PI*pow(fyz,3));
//dMz=a1*((a3*dHz/pow(Hyz,2))-(a2*pow(csch1(a2*Hyz),2)*(b-z)/(pow(fyz,1.5))));

d2Hz=I/(2*PI)*((3*pow((z-b),2)*pow(fyz,-2.5))-pow(fyz,-1.5));

dMz=a1*a2*(b-z)*pow(csch1(a2*Hyz),2)/(2*PI*pow(fyz,1.5))-(2*PI*a1*a3*(z-b))/pow(fyz,0.5);

SourceTermz=Mu0*M*dHz;
//SourceTermz=0.001;
dS[eqn]=(SourceTermz*dMz)+(Mu0*M*d2Hz);

return SourceTermz;

}

real csch1(real u)
{
	real a;
	
	a=2/(exp(u)-exp(-u));
	
	return a;
}

real coth1(real u)
{
	real a;
	
	a=(exp(u)+exp(-u))/(exp(u)-exp(-u));
	
	return a;
}

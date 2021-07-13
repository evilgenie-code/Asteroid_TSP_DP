
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>

#include "astro.h" 

const double pai =3.1415926535897932/180.0;
const double pi2 =3.1415926535897932*2.0;
const double HToR =3.1415926535897932 / 12;
const double RToH =12 / 3.1415926535897932;

/* common.c */

/*******************************************************************************/
/*                      ÉÜÉäÉEÉXì˙ÇÃåvéZ                                       */
/*                                                                             */
/*               égópéûä‘Ç™ínï˚ïWèÄéû(LST)Ç»ÇÁÅF LST=1                         */
/*               égópéûä‘Ç™ínï˚ïΩãœéû(LMT)Ç»ÇÁÅF LST=0                         */
/*******************************************************************************/

double JD(double year,double month,double day,double hour,int LST,double longitude)
{
	double s,juday,a,flag;

	flag =1;
	if (longitude <0){ 
		longitude = -longitude;
		flag=-1;
		}

	if ( month < 3 ){
		year = year -1;
		month = month +12;
		};
	s = floor( 365.25 * year) + floor( 30.59 * (month-2.0));
	if (LST==1)
	{juday = s + day + hour/24.0 -  flag*floor(longitude /15.0+0.5)/24.0	+ 1721086.5;}
		else
		{juday = s + day + hour/24.0 -  flag*longitude /15.0/24.0	+ 1721086.5;}


	if ( juday > 2299150 ) {
		a = floor ( year / 400 ) - floor( year / 100)+2;
		juday = juday + a;
		};

	return(juday);

}
/*******************************************************************************/
/*                      óÔì˙ÇÃåvéZ                                             */
/*                                                                             */
/*               égópéûä‘Ç™ínï˚ïWèÄéû(LST)Ç»ÇÁÅF LST=1                         */
/*               égópéûä‘Ç™ínï˚ïΩãœéû(LMT)Ç»ÇÁÅF LST=0                         */
/*******************************************************************************/

void Calendar( double juday, double *year,double *month,double *day,double *hour,int LST,double longitude)
{
	double flag,Z,Alfa,A,B,C,D,E,F,Day0;
	
	flag = 1;

	if (longitude<0 ){
		longitude = - longitude;
		flag=-1.0;
	}
	
	if (LST==1) juday =  juday +  flag*(double)(int)(longitude /15.0+0.5)/24.0;
		else    juday =  juday +  flag*longitude /15.0/24.0	;
		
	juday = juday + 0.5;
	Z = (double)(int)(juday);
	F = juday - Z;

	if (Z < 2299161.0) A = Z;
		else { 
			Alfa = (double)(int) ((Z - 1867216.25)/36524.25);
			A = Z + 1 + Alfa - (double)(int) (Alfa/4.0);}

	B = A + 1524;
	C = (double)(int) ((B-122.1)/365.25);
	D = (double)(int) (365.25 * C);
	E = (double)(int) ((B-D)/30.6001);

	Day0 = B - D - (double)(int)(30.60001*E) + F;
	*day = (double)(int)(Day0);
	*hour= F*24.0;
	if(E<14.0) *month = E-1.0;
	if(E==14.0 || E==15.0)	 *month = E- 13.0;
	if(*month>2.0) *year=C-4716;
	if(*month==1.0 || *month==2.0) *year = C - 4715;


}

/*******************************************************************************/
/*                      ëÂêîÇÃÇ‹ÇÈÇﬂåvéZ                                       */
/*******************************************************************************/

double Compact (double y, double x)
{
	y = x * ( y / x - floor (y / x));
	if((int)x==360 && y <0) y=y+360.0;
	return(y);
}


/*******************************************************************************/
/*                      ï\ÇÃílÇÃãﬂéóÅiíºê¸Åj                                   */
/*******************************************************************************/

double Table_Value(double a, double b, double c, double d, double Diff_ba, double Diff_ca)
{

	double ac, bd, acbd;
	
/*	
               Diff_ca
           a      ac       c
  Diff_ba        abcd 
           b      bd       d
*/

	ac = a + (c-a)*Diff_ca/10.0;
	bd = b + (d-b)*Diff_ca/10.0;
	acbd = ac +(bd-ac)*Diff_ba/10.0;

 return acbd;
}


/*******************************************************************************/
/*                      ÇkÇrÇsÇÃåvéZ                                           */
/*******************************************************************************/

double LST( double jd, double longitude)
{
	double lst,GST0 = 18.6461, JD0=2415020.0, R,b ;

	R = 366.2422/365.2422;
	b = jd - JD0;
	lst = GST0 + 24 * b * R + 3.24E-14 * b * b +longitude/15.0;
	lst = Compact(lst,24);
	if (lst<0.0){
		lst = lst +24;
		};
	return lst;
}

/*******************************************************************************/
/*                      ÇQÇOÇOÇOîNäÓèÄÇÃê¢ãIêîÇÃåvéZ                           */
/*******************************************************************************/

double T_2000 (double jd)
{
	double JD0=2451545.0 , t ;
	t=(jd-JD0)/(100.0*365.25);
	return t ;
}

/*******************************************************************************/
/*                      ÇPÇXÇOÇOîNäÓèÄÇÃê¢ãIêîÇÃåvéZ                           */
/*******************************************************************************/

double C_1900 (double jd)
{
	double JD0=2415021.0 , c ;
	c=(jd-JD0)/(100.0*365.25);
	return c ;
}


/*******************************************************************************/
/*                      ÇPÇWÇOÇOîNäÓèÄÇÃê¢ãIêîÇÃåvéZ                           */
/*******************************************************************************/

double J_1800 (double jd)
{
	double j;
	j = ( jd - 2378496.0 ) / 36525.0;
	return (j);
}


/*******************************************************************************/
/*                      â©ìπåXéŒäpÇÃåvéZ                                        */
/*******************************************************************************/
/*double Obl(double jd)
{
	double c,obl;
	
	c=C_1900(jd);
	obl = 23.4523 - 0.013*c -0.00000016388*c*c ; 

	return obl;
}
*/
/*******************************************************************************/
/*                      â©ìπåXéŒäpÇÃåvéZ                                        */
/*******************************************************************************/
double Obl0_2000(double jd)
{
	double t,ob,/*dF,dE,*/obl;

/*	Nutation_data( jd,&dF, &dE);*/
	
	t=T_2000(jd);
	obl = 23.439291 - 0.0130042*t - 0.00000016*t*t + 0.000000504*t*t*t; 

/*	obl  = obl + dE/pai;*/
	return obl;
}
double Obl_2000(double jd)
{
	double t,ob,dF,dE,obl;

	Nutation_data( jd,&dF, &dE);
	
	t=T_2000(jd);
//	obl = 23.439291 - 0.0130042*t - 0.00000016*t*t + 0.000000504*t*t*t; 

    t   = t/100;
    obl = 23.44484666666667 + (-4680.93*t - 1.55*t*t + 1999.25*t*t*t - 51.38*t*t*t*t 
           - 249.67*t*t*t*t*t - 39.05*t*t*t*t*t*t + 7.12*t*t*t*t*t*t*t 
            + 27.87*t*t*t*t*t*t*t*t + 5.79*t*t*t*t*t*t*t*t*t + 2.45*t*t*t*t*t*t*t*t*t*t)/3600;

	obl  = obl + dE/pai;
   
	return obl;
}

/*******************************************************************************/
/*                      â©åoâ©à‹Ç©ÇÁê‘åoê‘à‹Ç÷ÇÃïœä∑                           */
/*******************************************************************************/

void LamToRa(double jd,double Lam,double Bet,double *RA,double *DEC)
{

	double obl;
	
	obl = Obl_2000(jd);
	
	*RA  = atan2((cos(obl*pai)*cos(Bet*pai)*sin(Lam*pai)-sin(obl*pai)*sin(Bet*pai)),cos(Bet*pai)*cos(Lam*pai))/pai;
	if (*RA <0) *RA = *RA +360.0;
	*DEC = asin (cos(obl*pai)*sin(Bet*pai)+sin(obl*pai)*cos(Bet*pai)*sin(Lam*pai))/pai;

}

/*******************************************************************************/
/*                      ê‘åoê‘à‹Ç©ÇÁâ©åoâ©à‹Ç÷ÇÃïœä∑                           */
/*******************************************************************************/

void RaToLam(double jd,double RA,double DEC,double *Lam,double *Bet)
{

	double obl;
	
	obl = Obl_2000(jd);
	
	*Lam  = atan2((sin(obl*pai)*sin(DEC*pai)+cos(obl*pai)*cos(DEC*pai)*sin(RA*pai)),cos(DEC*pai)*cos(RA*pai))/pai;
	if (*Lam <0) *Lam = *Lam +360.0;
	*Bet = asin (cos(obl*pai)*sin(DEC*pai)-sin(obl*pai)*cos(DEC*pai)*sin(RA*pai))/pai;

}

/*******************************************************************************/
/*                      ê‘åoê‘à‹Ç©ÇÁï˚à çÇìxÇÃåvéZ                             */
/*          íAÇµÅAÇ±ÇÃåvéZÇÕÅuå√ìVï∂äwÅvãLç⁄ÇÃåvéZéÆÇ≈ÇÕÇ†ÇËÇ‹ÇπÇÒÅB           */
/*******************************************************************************/
 
void Az_El(double jd,double Ra, double Dec, double latitude,double longitude,  double *Az, double *El )
{
	double lst,lat,Ha,a,b;

    lst = LST(jd,longitude)*15.0;
	lat = latitude;
	Ha  = lst - Ra;

	a   = sin (Ha * pai);
	b   = cos (Ha * pai) * sin (lat * pai) - tan(Dec*pai)*cos(lat*pai);
	*Az = atan2 (a,b)/pai;
	*El =asin(sin (lat *pai) * sin ( Dec * pai) + cos (lat *pai)* cos ( Dec *pai)*cos(Ha*pai))/pai;

}


/*******************************************************************************/
/*                      ì˙ïtïœêîÇ©ÇÁîNåéì˙ÇÃï™ó£                               */
/*******************************************************************************/

void DateToYMD(double Date,double *year,double *month,double *day)

{
	double flag,Year,Month;
	
	flag=0;
	
	if (Date<0) {
		Date=-Date;
	    flag = 1;
	    };
	Year  = floor(Date);
	Month = floor((Date-Year)*100);
	*day  = (Date*10000-Year*10000-Month*100);
	*year = Year;
	*month= Month;
	if(flag==1) *year = -*year;

}

/*******************************************************************************/
/*                      éûä‘ï™ïbï\é¶Ç©ÇÁéûä‘Ç÷ÇÃïœä∑                           */
/*******************************************************************************/
double HMSToHour(double HMS)

{
	int flag=0;
	double H,Hour,Min,Sec;
	
	if (HMS<0){
		flag=1;
		HMS=-HMS;
		};
	Hour  = floor(HMS);
	Min   = floor((HMS - Hour)*100);
	Sec   = (HMS*10000-Hour*10000-Min*100);
	H     = Hour + Min/60.0 +  Sec/3600.0;
    if(flag==1)H=-H;
	return H;
}

/*******************************************************************************/
/*                      éûä‘ï™ïbï\é¶Ç©ÇÁäpìxÇ÷ÇÃïœä∑  ( 1999/07/01 í«â¡ )      */
/*******************************************************************************/
double HMSToDeg(double HMS)

{
	double H,Deg;
	
	H=HMSToHour( HMS);
	
	Deg = H*15.0;
	return Deg;
}

/*******************************************************************************/
/*                      ìxï™ïbï\é¶Ç©ÇÁìxÇ÷ÇÃïœä∑                               */
/*******************************************************************************/

double DMSToDeg(double DMS)

{
	int flag=0;
	double D,Deg,Min,Sec;
	
	
	if (DMS<0){
		flag=1;
		DMS=-DMS;
		};
	Deg   = floor(DMS);
	Min   = floor((DMS - Deg)*100);
	if(Min>99.0){Min=Min-100;Deg=Deg+1;}
	Sec   = DMS*10000-Deg*10000-Min*100;
	if(Sec>99.0){Sec=Sec-100;Min=Min+1;}
	D     = Deg +  Min/60.0 +  Sec/3600.0;
	if(flag==1)D=-D;
    return D;

}

/*******************************************************************************/
/*                      ìxï\é¶Ç©ÇÁéûï™ïbï\é¶Ç÷ÇÃïœä∑                           */
/*******************************************************************************/

double DegToHMS(double D)

{
	int    flag=0;
	double Hour,Min, Sec,H;
	
	if(D<0) {
		D=-D;
		flag=1;
		};

	D     = D/15.0;
	Hour  = floor(D);
	Min   = floor((D - Hour)*60);
	Sec   = D*3600-Hour*3600-Min*60;
	H     = Hour +  Min/100.0 +  Sec/10000.0;
	if(flag==1)H=-H;

    return H;
}

/*******************************************************************************/
/*                      ìxï\é¶Ç©ÇÁìxï™ïbï\é¶Ç÷ÇÃïœä∑                           */
/*******************************************************************************/

double DegToDMS(double D)

{

	int    flag=0;
	double Deg,Min,Sec,H;
	
	if(D<0) {
		D=-D;
		flag=1;
		};
	Deg   = floor(D);
	Min   = floor((D - Deg)*60.0);
	Sec   = (D*3600-Deg*3600-Min*60);
	H     = Deg +  Min/100.0 + Sec/10000.0;
	if(flag==1)H=-H;
    return H;
}

/*******************************************************************************/
/*                      è¨êîì_ëÊìÒà Ç≈ÇÃéléÃå‹ì¸                               */
/*******************************************************************************/
double Digit2( double x)
{
	int flag=0;

	if(x<0){
		x=-x;
		flag=1;
	}
	x=(double)((int)(x*100.0+0.5));
	x=x/100.0;
	if(flag==1)x=-x;
	return x;

}

/*******************************************************************************/
/*                      è¨êîì_ëÊéOà Ç≈ÇÃéléÃå‹ì¸                               */
/*******************************************************************************/

double Digit3( double x)
{
	int flag=0;

	if(x<0){
		x=-x;
		flag=1;
	}
	x=(double)((int)(x*1000.0+0.5));
	x=x/1000.0;
	if(flag==1)x=-x;
	return x;

}

/*******************************************************************************/
/*                     çŒç∑ÇÃåvéZ                                              */
/*              (Astronomical algorithmsÇ…ÇÊÇÈ)                                */
/*******************************************************************************/

void Precession (double jd, double Year0, double RA0,double DEC0,double V1,double V2,double *RA,double *DEC)
{
	double jd0,t,T, Zeta,Zeta1,Zeta2,Zeta3,Theta,Theta1,Theta2,Theta3,Zed,Zed1,Zed2,Zed3;
	double CosTheta,SinTheta ,CosDec0,SinDec0,CosRaZeta ,SinRaZeta;

	jd0=JD( Year0,1,1,12,0,0);/*1999/07/02*/

	T=(jd0-2451545.0)/36525.0;/*1999/07/02*/
	t=(jd-jd0)/36525.0;

	Zeta1  = 2306.2181 + 1.39656*T - 0.000139*T*T;
	Zeta2  = 0.30188 - 0.000344*T;
	Zeta3  = 0.017998;
	Theta1 = 2004.3109 - 0.85330*T - 0.000217*T*T;
	Theta2 = - ( 0.42665 + 0.000217*T);
	Theta3 = - 0.041833;
	Zed1   = 2306.2181 + 1.39656*T - 0.000139*T*T;
	Zed2   = 1.09468 + 0.000066*T;
	Zed3   = 0.018203;

/* for B1900 system Åuå√ìVï∂äwÇÃíËêîÅvÅ@*/
/*
	jd0=JD( Year0,1,0,19.31264,0,0);
	t=(jd-jd0)/36525.0;

	Zeta1  = 2304.250;
	Zeta2  =    0.302;
	Zeta3  =    0.018;
	Theta1 = 2004.68 ;
	Theta2 =   -0.426;
	Theta3 =   -0.042;
	Zed1   = 2304.250;
	Zed2   =    1.093;
	Zed3   =    0.018;
*/

	Zeta  = (  Zeta1*t +  Zeta2*t*t + Zeta3*t*t*t )/3600.0;
	Zed   = (   Zed1*t +   Zed2*t*t +  Zed3*t*t*t )/3600.0;
	Theta = ( Theta1*t + Theta2*t*t +Theta3*t*t*t )/3600.0;

/*printf("Zeta=%10.3f Zed =%10.3f Theta=%10.3f \n",Zeta,Zed,Theta);*/

	Zeta = Zeta*pai;
	Zed  = Zed*pai;
	Theta= Theta*pai;
	RA0  = (RA0 + V1*t*100/3600.0/cos(DEC0*pai))*pai;
	DEC0 = (DEC0+ V2*t*100/3600.0)*pai;

	CosTheta  = cos (Theta);
	SinTheta  = sin (Theta);
	CosDec0   = cos (DEC0);
	SinDec0   = sin (DEC0);
	CosRaZeta = cos(RA0+Zeta);
	SinRaZeta = sin(RA0+Zeta);

	*DEC   = asin (CosTheta*SinDec0 + SinTheta*CosDec0*CosRaZeta)/pai;
	*RA =atan2(CosDec0*SinRaZeta ,CosTheta*CosDec0*CosRaZeta-SinTheta*SinDec0)/pai;
	*RA =*RA + Zed/pai;
	if (*RA<0) *RA += 360.0;

}






/* Nutation.c */

static double Nut[106][10]={

{ 0, 0, 0, 0, 1, 6798.4,-171996, -174.2, 92025, 8.9},
{ 0, 0, 0, 0, 2, 3399.2,   2062,    0.2,  -895, 0.5},
{-2, 0, 2, 0, 1, 1305.5,     46,      0,   -24, 0  },
{ 2, 0,-2, 0, 0, 1095.2,     11,      0,     0, 0  },
{-2, 0, 2, 0, 2, 1615.7,     -3,      0,     1, 0  },
{ 1,-1, 0,-1, 0, 3232.9,     -3,      0,     0, 0  },
{ 0,-2, 2,-2, 1, 6786.3,     -2,      0,     1, 0  },
{ 2, 0,-2, 0, 1,  943.2,      1,      0,     0, 0  },
{ 0, 0, 2,-2, 2,  182.6, -13187,   -1.6,  5736,-3.1},
{ 0, 1, 0, 0, 0,  365.3,   1426,   -3.4,    54,-0.1},
{ 0, 1, 2,-2, 2,  121.7,   -517,    1.2,   224,-0.6},
{ 0,-1, 2,-2, 2,  365.2,    217,   -0.5,   -95, 0.3},
{ 0, 0, 2,-2, 1,  177.8,    129,    0.1,   -70, 0  },
{ 2, 0, 0,-2, 0,  205.9,     48,      0,     1, 0  },
{ 0, 0, 2,-2, 0,  173.3,    -22,      0,     0, 0  },
{ 0, 2, 0, 0, 0,  182.6,     17,   -0.1,     0, 0  },
{ 0, 1, 0, 0, 1,  386.0,    -15,      0,     9, 0  },
{ 0, 2, 2,-2, 2,   91.3,    -16,    0.1,     7, 0  },
{ 0,-1, 0, 0, 1,  346.6,    -12,      0,     6, 0  },
{-2, 0, 0, 2, 1,  199.8,     -6,      0,     3, 0  },
{ 0,-1, 2,-2, 1,  346.6,     -5,      0,     3, 0  },
{ 2, 0, 0,-2, 1,  212.3,      4,      0,    -2, 0  },
{ 0, 1, 2,-2, 1,  119.6,      4,      0,    -2, 0  },
{ 1, 0, 0,-1, 0,  411.8,     -4,      0,     0, 0  },
{ 2, 1, 0,-2, 0,  131.7,      1,      0,     0, 0  },
{ 0, 0,-2, 2, 1,  169.0,      1,      0,     0, 0  },
{ 0, 1,-2, 2, 0,  329.8,     -1,      0,     0, 0  },
{ 0, 1, 0, 0, 2,  409.2,      1,      0,     0, 0  },
{-1, 0, 0, 1, 1,  388.3,      1,      0,     0, 0  },
{ 0, 1, 2,-2, 0,  117.5,     -1,      0,     0, 0  },
{ 0, 0, 2, 0, 2,   13.7,  -2274,   -0.2,   977,-0.5},
{ 1, 0, 0, 0, 0,   27.6,    712,    0.1,    -7, 0  },
{ 0, 0, 2, 0, 1,   13.6,   -386,   -0.4,   200, 0  },
{ 1, 0, 2, 0, 2,    9.1,   -301,      0,   129,-0.1},
{ 1, 0, 0,-2, 0,   31.8,   -158,      0,    -1, 0  },
{-1, 0, 2, 0, 2,   27.1,    123,      0,   -53, 0  },
{ 0, 0, 0, 2, 0,   14.8,     63,      0,    -2, 0  },
{ 1, 0, 0, 0, 1,   27.7,     63,    0.1,   -33, 0  },
{-1, 0, 0, 0, 1,   27.4,    -58,   -0.1,    32, 0  },
{-1, 0, 2, 2, 2,    9.6,    -59,      0,    26, 0  },
{ 1, 0, 2, 0, 1,    9.1,    -51,      0,    27, 0  },
{ 0, 0, 2, 2, 2,    7.1,    -38,      0,    16, 0  },
{ 2, 0, 0, 0, 0,   13.8,     29,      0,    -1, 0  },
{ 1, 0, 2,-2, 2,   23.9,     29,      0,   -12, 0  },
{ 2, 0, 2, 0, 2,    6.9,    -31,      0,    13, 0  },
{ 0, 0, 2, 0, 0,   13.6,     26,      0,    -1, 0  },
{-1, 0, 2, 0, 1,   27.0,     21,      0,   -10, 0  },
{-1, 0, 0, 2, 1,   32.0,     16,      0,    -8, 0  },
{ 1, 0, 0,-2, 1,   31.7,    -13,      0,     7, 0  },
{-1, 0, 2, 2, 1,    9.5,    -10,      0,     5, 0  },
{ 1, 1, 0,-2, 0,   34.8,     -7,      0,     0, 0  },
{ 0, 1, 2, 0, 2,   13.2,      7,      0,    -3, 0  },
{ 0,-1, 2, 0, 2,   14.2,     -7,      0,     3, 0  },
{ 1, 0, 2, 2, 2,    5.6,     -8,      0,     3, 0  },
{ 1, 0, 0, 2, 0,    9.6,      6,      0,     0, 0  },
{ 2, 0, 2,-2, 2,   12.8,      6,      0,    -3, 0  },
{ 0, 0, 0, 2, 1,   14.8,     -6,      0,     3, 0  },
{ 0, 0, 2, 2, 1,    7.1,     -7,      0,     3, 0  },
{ 1, 0, 2,-2, 1,   23.9,      6,      0,    -3, 0  },
{ 0, 0, 0,-2, 1,   14.7,     -5,      0,     3, 0  },
{ 1,-1, 0, 0, 0,   29.8,      5,      0,     0, 0  },
{ 2, 0, 2, 0, 1,    6.9,     -5,      0,     3, 0  },
{ 0, 1, 0,-2, 0,   15.4,     -4,      0,     0, 0  },
{ 1, 0,-2, 0, 0,   26.9,      4,      0,     0, 0  },
{ 0, 0, 0, 1, 0,   29.5,     -4,      0,     0, 0  },
{ 1, 1, 0, 0, 0,   25.6,     -3,      0,     0, 0  },
{ 1, 0, 2, 0, 0,    9.1,      3,      0,     0, 0  },
{ 1,-1, 2, 0, 2,    9.4,     -3,      0,     1, 0  },
{-1,-1, 2, 2, 2,    9.8,     -3,      0,     1, 0  },
{-2, 0, 0, 0, 1,   13.7,     -2,      0,     1, 0  },
{ 3, 0, 2, 0, 2,    5.5,     -3,      0,     1, 0  },
{ 0,-1, 2, 2, 2,    7.2,     -3,      0,     1, 0  },
{ 1, 1, 2, 0, 2,    8.9,      2,      0,    -1, 0  },
{-1, 0, 2,-2, 1,   32.6,     -2,      0,     1, 0  },
{ 2, 0, 0, 0, 1,   13.8,      2,      0,    -1, 0  },
{ 1, 0, 0, 0, 2,   27.8,     -2,      0,     1, 0  },
{ 3, 0, 0, 0, 0,    9.2,      2,      0,     0, 0  },
{ 0, 0, 2, 1, 2,    9.3,      2,      0,    -1, 0  },
{-1, 0, 0, 0, 2,   27.3,      1,      0,    -1, 0  },
{ 1, 0, 0,-4, 0,   10.1,     -1,      0,     0, 0  },
{-2, 0, 2, 2, 2,   14.6,      1,      0,    -1, 0  },
{-1, 0, 2, 4, 2,    5.8,     -2,      0,     1, 0  },
{ 2, 0, 0,-4, 0,   15.9,     -1,      0,     0, 0  },
{ 1, 1, 2,-2, 2,   22.5,      1,      0,    -1, 0  },
{ 1, 0, 2, 2, 1,    5.6,     -1,      0,     1, 0  },
{-2, 0, 2, 4, 2,    7.3,     -1,      0,     1, 0  },
{-1, 0, 4, 0, 2,    9.1,      1,      0,     0, 0  },
{ 1,-1, 0,-2, 0,   29.3,      1,      0,     0, 0  },
{ 2, 0, 2,-2, 1,   12.8,      1,      0,    -1, 0  },
{ 2, 0, 2, 2, 2,    4.7,     -1,      0,     0, 0  },
{ 1, 0, 0, 2, 1,    9.6,     -1,      0,     0, 0  },
{ 0, 0, 4,-2, 2,   12.7,      1,      0,     0, 0  },
{ 3, 0, 2,-2, 2,    8.7,      1,      0,     0, 0  },
{ 1, 0, 2,-2, 0,   23.8,     -1,      0,     0, 0  },
{ 0, 1, 2, 0, 1,   13.1,      1,      0,     0, 0  },
{-1,-1, 0, 2, 1,   35.0,      1,      0,     0, 0  },
{ 0, 0,-2, 0, 1,   13.6,     -1,      0,     0, 0  },
{ 0, 0, 2,-1, 2,   25.4,     -1,      0,     0, 0  },
{ 0, 1, 0, 2, 0,   14.2,     -1,      0,     0, 0  },
{ 1, 0,-2,-2, 0,    9.5,     -1,      0,     0, 0  },
{ 0,-1, 2, 0, 1,   14.2,     -1,      0,     0, 0  },
{ 1, 1, 0,-2, 1,   34.7,     -1,      0,     0, 0  },
{ 1, 0,-2, 2, 0,   32.8,     -1,      0,     0, 0  },
{ 2, 0, 0, 2, 0,    7.2,      1,      0,     0, 0  },
{ 0, 0, 2, 4, 2,    4.8,     -1,      0,     0, 0  },
{ 0, 1, 0, 1, 0,   27.3,      1,      0,     0, 0  }};

/*
static double dNut[4][10]={

{ 0, 0, 0, 0, 1, 6798.4,   -725,    417,   213, 224},
{ 0, 1, 0, 0, 0,  365.3,    523,     61,   208, -24},
{ 0, 0, 2,-2, 2,  182.6,    102,   -118,   -41, -47},
{ 0, 0, 2, 0, 2,   13.7,    -81,      0,    32, 0  }};
*/

void Nutation_data(double jd,double *dF, double *dE)
{
	double L,Ld,F,D,Omg,Arg,t;
	int i;

	t=T_2000 (jd);
	L=    134.962981389 + 477198.867398056*t + 0.008697222*t*t + 1.0/ 56250.0*t*t*t;
	Ld=   357.527723333 +  35999.050340000*t - 0.000160278*t*t - 1.0/300000.0*t*t*t;
	F=     93.271910278 + 483202.017538056*t - 0.003682500*t*t + 1.0/ 56250.0*t*t*t;
	D=    297.850363056 + 445267.111480000*t - 0.001914167*t*t + 1.0/189474.0*t*t*t;
	Omg=  125.044522222 -   1934.136260833*t + 0.002070833*t*t + 1.0/450000.0*t*t*t;
	*dF=0.0;
	*dE=0.0;
	for(i=0;i<106;i++){
		Arg=(Nut[i][0]*L+Nut[i][1]*Ld+Nut[i][2]*F+Nut[i][3]*D+Nut[i][4]*Omg)*pai;
		*dF=*dF+(Nut[i][6]+Nut[i][7]*t)/3600.0*1e-4*sin(Arg);
		*dE=*dE+(Nut[i][8]+Nut[i][9]*t)/3600.0*1e-4*cos(Arg);
		}
/*
	for(i=0;i<4;i++){
		Arg=(dNut[i][0]*L+dNut[i][1]*Ld+dNut[i][2]*F+dNut[i][3]*D+dNut[i][4]*Omg)*pai;
		*dF=*dF+(dNut[i][6]*sin(Arg)+dNut[i][7]*cos(Arg))/3600.0*1e-5;
		*dE=*dE+(dNut[i][8]*cos(Arg)+dNut[i][9]*sin(Arg))/3600.0*1e-5;
		}
*/
	*dF=*dF*pai;
	*dE=*dE*pai;

}
void Nutation (FILE *fp_Bin ,int DE, double jd, Matrix_3 *Nu)
{

    double E,E0,dE,dE0,dF,dF0,r[6];
	int nctr,ntarg,inside;

	nctr  = 12;
	/* Where does the ephemeris think things should be? */
	stcomm.km=FALSE_;/* AU & AU/Day*/
	stcomm.bary=TRUE_;/*Barycentre (Solar-sys. Barycentre)*/

	if(DE=406)	Nutation_data( jd,&dF, &dE);
	else{
	ntarg =  14;
	PLEPH(fp_Bin,jd,ntarg,nctr,r,&inside);
	dF = r[0];
	dE = r[1];
	}	
	dF0=dF;
	dE0=dE;

	E0 = Obl0_2000(jd)*pai;
	E  = E0 + dE;

	Matrix_R1(Nu1, - E, &Nu1);
	Matrix_R3(Nu2, -dF, &Nu2);
	Matrix_R1(Nu3, +E0, &Nu3);
	Multi_Matrix_3(Nu1,Nu2,&Nu4);
	Multi_Matrix_3(Nu4,Nu3,Nu);

}
void Precession_Data (double jd,Matrix_3 *Pr)
{
	double jd0,t,T, Zeta1,Zeta2,Zeta3,Theta1,Theta2,Theta3,Zed1,Zed2,Zed3;
	double Zeta,Zed,Theta;

	jd0=JD( 2000,1,1,12,0,0);

	T=(jd0-2451545.0)/36525.0;
	t=(jd-jd0)/36525.0;
	Zeta1  = 2306.2181 + 1.39656*T - 0.000139*T*T;
	Zeta2  = 0.30188 - 0.000344*T;
	Zeta3  = 0.017998;
	Theta1 = 2004.3109 - 0.85330*T - 0.000217*T*T;
	Theta2 = - ( 0.42665 + 0.000217*T);
	Theta3 = - 0.041833;
	Zed1   = 2306.2181 + 1.39656*T - 0.000139*T*T;
	Zed2   = 1.09468 + 0.000066*T;
	Zed3   = 0.018203;

	Zeta  = (  Zeta1*t +  Zeta2*t*t + Zeta3*t*t*t )/3600.0;
	Zed   = (   Zed1*t +   Zed2*t*t +  Zed3*t*t*t )/3600.0;
	Theta = ( Theta1*t + Theta2*t*t +Theta3*t*t*t )/3600.0;

	Matrix_R3(Pr1, - Zed*pai, &Pr1);
	Matrix_R2(Pr2, Theta*pai, &Pr2);
	Matrix_R3(Pr3, -Zeta*pai, &Pr3);
	Multi_Matrix_3(Pr1,Pr2,&Pr4);
	Multi_Matrix_3(Pr4,Pr3,Pr);

}






/* Matrix.c */

void Sum_Matrix_3(Matrix_3 abc,Matrix_3 def,Matrix_3 *ghi)
{
	int i,j;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++)
		ghi->max[j][i]=abc.max[j][i]+def.max[j][i];
		}
}
void Multi_Matrix_3(Matrix_3 abc,Matrix_3 def,Matrix_3 *ghi)
{
	ghi->max[0][0]=abc.max[0][0]*def.max[0][0]+abc.max[0][1]*def.max[1][0]+abc.max[0][2]*def.max[2][0];
	ghi->max[0][1]=abc.max[0][0]*def.max[0][1]+abc.max[0][1]*def.max[1][1]+abc.max[0][2]*def.max[2][1];
	ghi->max[0][2]=abc.max[0][0]*def.max[0][2]+abc.max[0][1]*def.max[1][2]+abc.max[0][2]*def.max[2][2];
	ghi->max[1][0]=abc.max[1][0]*def.max[0][0]+abc.max[1][1]*def.max[1][0]+abc.max[1][2]*def.max[2][0];
	ghi->max[1][1]=abc.max[1][0]*def.max[0][1]+abc.max[1][1]*def.max[1][1]+abc.max[1][2]*def.max[2][1];
	ghi->max[1][2]=abc.max[1][0]*def.max[0][2]+abc.max[1][1]*def.max[1][2]+abc.max[1][2]*def.max[2][2];
	ghi->max[2][0]=abc.max[2][0]*def.max[0][0]+abc.max[2][1]*def.max[1][0]+abc.max[2][2]*def.max[2][0];
	ghi->max[2][1]=abc.max[2][0]*def.max[0][1]+abc.max[2][1]*def.max[1][1]+abc.max[2][2]*def.max[2][1];
	ghi->max[2][2]=abc.max[2][0]*def.max[0][2]+abc.max[2][1]*def.max[1][2]+abc.max[2][2]*def.max[2][2];
}
void Sub_Matrix_3(Matrix_3 abc,Matrix_3 def,Matrix_3 *ghi)
{
	int i,j;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++)
		ghi->max[j][i]=abc.max[j][i]-def.max[j][i];
		}
}
void Multi_Const_Matrix_3(Matrix_3 abc,double a,Matrix_3 *ghi)
{
	int i,j;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++)
		ghi->max[j][i]=a*abc.max[j][i];
		}
}

void Sum_Matrix_1(Matrix_1 ab,Matrix_1 de,Matrix_1 *gh)
{
	int i;

	for(i=0;i<3;i++){
		gh->max[i]=ab.max[i]+de.max[i];
		}
}
void Multi_Matrix_1(Matrix_3 abc,Matrix_1 de,Matrix_1 *gh)
{
	gh->max[0]=abc.max[0][0]*de.max[0]+abc.max[0][1]*de.max[1]+abc.max[0][2]*de.max[2];
	gh->max[1]=abc.max[1][0]*de.max[0]+abc.max[1][1]*de.max[1]+abc.max[1][2]*de.max[2];
	gh->max[2]=abc.max[2][0]*de.max[0]+abc.max[2][1]*de.max[1]+abc.max[2][2]*de.max[2];
}
void Sub_Matrix_1(Matrix_1 ab,Matrix_1 de,Matrix_1 *gh)
{
	int i;

	for(i=0;i<3;i++){
		gh->max[i]=ab.max[i]-de.max[i];
		}
}
void Multi_Const_Matrix_1(Matrix_1 ab,double a,Matrix_1 *gh)
{
	int i;

	for(i=0;i<3;i++){
		gh->max[i]=a*ab.max[i];
		}
}
void Eql_Matrix_1(Matrix_1 ab, Matrix_1 *gh)
{
	int i;

	for(i=0;i<3;i++){
		gh->max[i]=ab.max[i];
		}
}
void Scalar (Matrix_1 ab, double *a)
{
	double sum=0.0;
	int i;

	for(i=0;i<3;i++){
		sum=sum+ab.max[i]*ab.max[i];
		}
	*a= sqrt(sum);
}
void Scalar_Dot (Matrix_1 ab,Matrix_1 de, double *a)
{
	double sum=0.0;
	int i;

	for(i=0;i<3;i++){
		sum=sum+ab.max[i]*de.max[i];
		}
	*a= sum;
}

void Matrix_R1(Matrix_3 abc,double a, Matrix_3 *ghi)
{

			ghi->max[0][0]=  1.0;
			ghi->max[1][0]=  0.0;
			ghi->max[2][0]=  0.0;
			ghi->max[0][1]=  0.0;
			ghi->max[1][1]=  cos(a);
			ghi->max[2][1]= -sin(a);
			ghi->max[0][2]=  0.0;
			ghi->max[1][2]=  sin(a);
			ghi->max[2][2]=  cos(a);
}
void Matrix_R2(Matrix_3 abc,double a, Matrix_3 *ghi)
{

			ghi->max[0][0]=  cos(a);
			ghi->max[1][0]=  0.0;
			ghi->max[2][0]=  sin(a);
			ghi->max[0][1]=  0.0;
			ghi->max[1][1]=  1.0;
			ghi->max[2][1]=  0.0;
			ghi->max[0][2]=  -sin(a);
			ghi->max[1][2]=  0.0;
			ghi->max[2][2]=  cos(a);
}
void Matrix_R3(Matrix_3 abc,double a, Matrix_3 *ghi)
{

			ghi->max[0][0]=  cos(a);
			ghi->max[1][0]=  -sin(a);
			ghi->max[2][0]=  0.0;
			ghi->max[0][1]=  sin(a);
			ghi->max[1][1]=  cos(a);
			ghi->max[2][1]=  0.0;
			ghi->max[0][2]=  0.0;
			ghi->max[1][2]=  0.0;
			ghi->max[2][2]=  1.0;
}




/*jpl_eph.c*/

/* The following portion is moved to astro.h by S.Takesako 1999/11/10 */
/**************************************************************/
/*
#define WAIT getch()
#define TRUE_ (1)
#define FALSE_ (0)
*/
/* Common Block Declarations */
/*struct { char **cnam;
									double *cval,ss[3],au,emrat;
									int denum, ncon, ipt[36],lpt[3],ksize;
							} HEADER;

struct { double *buf;}epib;
*/
/*     COMMON AREA STCOMM:                                               */
/*          KM   LOGICAL FLAG DEFINING PHYSICAL UNITS OF THE OUTPUT      */
/*               STATES. KM = .TRUE., KM AND KM/SEC                      */
/*                          = .FALSE., AU AND AU/DAY                     */
/*               DEFAULT VALUE = .FALSE.  (KM DETERMINES TIME UNIT       */
/*              FOR NUTATIONS AND LIBRATIONS. ANGLE UNIT IS ALWAYS RAD.) */
/*        BARY   LOGICAL FLAG DEFINING OUTPUT CENTER.                    */
/*               ONLY THE 9 PLANETS ARE AFFECTED.                        */
/*                       BARY = .TRUE. =\ CENTER IS SOLAR-SYS. BARYCEN.  */
/*                             = .FALSE. =\ CENTER IS SUN                */
/*               DEFAULT VALUE = .FALSE.                                 */
/*       PVSUN   DP 6-WORD ARRAY CONTAINING THE BARYCENTRIC POSITION AND */
/*               VELOCITY OF THE SUN.                                    */
/*struct { int km, bary;
									double pvsun[6];
							} stcomm;

*/
/**************************************************************/


int const2 = 2;
int const3 = 3;

/*========================================================================+
|                                 VOID ERROR                              |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| Print out the error string.                                             |
+-------------------------------------------------------------------------+
| INPUTS : error_text = String for the error text                         |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
+-------------------------------------------------------------------------+
| REFERENCE: %                                                            |
+========================================================================*/
void ERROR(int i,char *msg)
{
		fprintf(stderr,"\n\n...Run-time error...");
		fprintf(stderr,"\n%d %s",i,msg);
		fprintf(stderr,"\n...now exiting to system...");
		exit(1);
}

/*========================================================================+
|                                 DVECTOR                                 |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| Allocate a double vector with subscript range v[nl..nh].                |
+-------------------------------------------------------------------------+
| INPUTS : nl = Start dimension (subscript range) of the vector           |
|          nh = End dimension (subscript range) of the vector             |
+-------------------------------------------------------------------------+
| OUTPUTS: dvector = Allocated vector                                     |
+-------------------------------------------------------------------------+
| REFERENCE: NUMERICAL RECEPIES IN C                                      |
+========================================================================*/
double *dvector(int nl,int nh)
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
	if (!v)
				ERROR(0,"allocation failure in dvector()");
	return v-nl+1;
}

/*========================================================================+
|                          FREE_DVECTOR(DOUBLE)                           |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| Free a double vector allocated with dvector().                          |
+-------------------------------------------------------------------------+
| INPUTS : v  = the vector, which should be free.                         |
|          nl = The starting number of row.                               |
+-------------------------------------------------------------------------+
| OUTPUTS: %                                                              |
+-------------------------------------------------------------------------+
| REFERENCE: NUMERICAL RECEPIES IN C                                      |
+========================================================================*/
void free_dvector(double *v,int nl)
{
		free((char *) (v+nl-1));
}

/*========================================================================+
|                              CHAR CMATRIX                               |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| Allocate a double matrix with subscript range m[nrl..nrh][ncl..nch].    |
+-------------------------------------------------------------------------+
| INPUTS : nrl = Starting row number (subscript range) of the matrix      |
|          nrh = Ending row number (subscript range) of the matrix        |
|          ncl = Starting column number (subscript range) of the matrix   |
|          nch = Ending column number (subscript range) of the matrix     |
+-------------------------------------------------------------------------+
| OUTPUTS: dmatrix = Allocated matrix                                     |
+-------------------------------------------------------------------------+
| REFERENCE: %                                                            |
+========================================================================*/
char **Cmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	char **m;

	/* Allocate pointers to rows: */
	m=(char **) malloc((size_t)((nrow+1)*sizeof(char*)));
	if (!m)
				ERROR(0,"allocation failure 1 in matrix()");
	m += 1;
	m -= nrl;

	/* Allocate rows and set pointers to them: */
	m[nrl]=(char *) malloc((size_t)((nrow*ncol+1)*sizeof(char)));
	if (!m[nrl])
				ERROR(0,"allocation failure 2 in matrix()");
	m[nrl] += 1;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++)
					m[i]=m[i-1]+ncol;

	/* Return pointer to array of pointers to rows : */
	return m;
}

/*========================================================================+
|                             FREE_MATRIX(CHAR)                           |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| Free a double matrix allocated with dmatrix().                          |
+-------------------------------------------------------------------------+
| INPUTS : m  = the matrix, which should be free.                         |
|          nrl= The starting number of row.                               |
|          ncl= The starting number of column.                            |
+-------------------------------------------------------------------------+
| OUTPUTS: %                                                              |
+-------------------------------------------------------------------------+
| REFERENCE: %                                                            |
+========================================================================*/
void free_Cmatrix(char **m, int nrl, int ncl)
/* Free a double matrix allocated by dmatrix() : */
{
		free((char*) (m[nrl]+ncl-1));
		free((char*) (m+nrl-1));
}

/*========================================================================+
|                                GO2POS                                   |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| Position the cursor about x characters from the actual position of the  |
| cursor in a specified file to the right side.                           |
+-------------------------------------------------------------------------+
| INPUTS : fp = File identification                                       |
|          nr = Number of characters to be read.                          |
+-------------------------------------------------------------------------+
| OUTPUTS: %                                                              |
+-------------------------------------------------------------------------+
| REFERENCE: %                                                              |
+========================================================================*/
void GO2POS(FILE *fp,int nr)
{
	int i;
	char ch;
	for(i=1;i<=nr;++i)fscanf_s(fp,"%c",&ch);

}

/*========================================================================+
|                               INITIAL                                   |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| INITIALIZE THE CONSTANTS FROM THE JPL BINARY FILE.                      |
+-------------------------------------------------------------------------+
| INPUTS : %                                                              |
+-------------------------------------------------------------------------+
| OUTPUTS: FP_BIN = POINTER TO THE JPL BINARY FILE.                       |
+-------------------------------------------------------------------------+
| REFERENCE: %                                                            |
+========================================================================*/
FILE *INITIAL(int DE)
//void INITIAL(int DE,FILE *fp_Bin)
{
/* The following portion is modified by S.Takesako 1999/11/10 */
/**************************************************************/
		int  I,J,K;
		char *Name="        ";
		long n;
		FILE *fp_Bin;

		if (DE==200){
		Name="dos.200";
		HEADER.ksize=1652;
		}
		else
			if (DE==405){
				Name="dos.405";
				HEADER.ksize=2036;
				}
		else
				{
				Name="dos.406";
				HEADER.ksize=1456;
				}
/**************************************************************/


		/* READ THE SIZE AND NUMBER OF MAIN EPHEMERIS: */
		if((fopen_s(&fp_Bin,Name,"rb")) != NULL)
				{ fprintf(stderr, "Cannot open %s file.\n",Name);
						exit(1);
				};
/* The following portion is commented out by S.Takesako 1999/11/10 */
/*
		fseek(fp_Bin,-1L*sizeof(int),SEEK_END);
		fread(&HEADER.ncon,sizeof(int),1, fp_Bin);
		SIZE=-42L*sizeof(int)-2L*sizeof(double)-3L*sizeof(double)-
							(long)HEADER.ncon*(sizeof(double)+9L*sizeof(char));
		fseek(fp_Bin,SIZE,SEEK_END);
*/
		/* READ THE SIZE AND NUMBER OF MAIN EPHEMERIS: */

/* The following portion is commented out by S.Takesako 1999/11/10 */
/*
		fread(&HEADER.ksize,sizeof(int),1, fp_Bin);
*/

		/* READ THE NUMBER AND NAMES OF CONSTANTS: (GROUP 1040/4):*/
/* The following portion is commented out by S.Takesako 1999/11/10 */
/*
		HEADER.cnam = Cmatrix(1,HEADER.ncon,0,8);
		for(I=1;I<=HEADER.ncon;++I)
					fread(HEADER.cnam[I],9,1,fp_Bin);
*/
/* The following portion is added by S.Takesako 1999/11/10 */
/***************************************************/   
		for (I=0; I<3 ; I++)
		fread  (TTL[I] ,sizeof(char),84,fp_Bin);

		for (I=1; I<=400 ; I++){
		fread  (CNAM[I] ,sizeof(char),6,fp_Bin);
		HEADER.cnam[I]=*CNAM[I];
		}
/***************************************************/   

		/* READ NUMBER OF VALUES AND VALUES (GROUP 1041/4): */
/* The following portion is commented out by S.Takesako 1999/11/10 */
/*
		HEADER.cval = dvector(1,HEADER.ncon);
		for(I = 1;I <= HEADER.ncon; ++I)
					fread(&HEADER.cval[I],sizeof(double),1, fp_Bin);
*/
		for(I=0;I<=2;++I)
					fread(&HEADER.ss[I],sizeof(double),1, fp_Bin);

    	fread(&HEADER.ncon,sizeof(long),1, fp_Bin);
		fread(&HEADER.au,sizeof(double),1, fp_Bin);
		fread(&HEADER.emrat,sizeof(double),1, fp_Bin);

		/* READ POINTERS NEEDED BY INTERP (GROUP 1050): */

/* The following portion is modified by S.Takesako 1999/11/10 */
/***************************************************/   
		ipt=&IPT[0][0];
		for(I = 0;I <= 35; ++I)	fread(ipt+I,sizeof(long),1, fp_Bin);

		I=0;
		for(J=0;J<3;++J){
			for(K=0;K<12;++K){
			HEADER.ipt[I]=IPT[K][J];
			I=I+1;
			}
		}

		fread(&HEADER.denum,sizeof(long),1, fp_Bin);

		for(I=0;I<=2;++I)
					fread(&HEADER.lpt[I],sizeof(long),1, fp_Bin);

		n=1;
		fseek(fp_Bin,HEADER.ksize*n*4,0);

		HEADER.cval = dvector(1,400/*HEADER.ncon*/);
		for(I = 1;I <= 400/*HEADER.ncon*/; ++I)
					fread(&HEADER.cval[I],sizeof(double),1, fp_Bin);

		epib.buf=dvector(1,HEADER.ksize/2);
/***************************************************/   

		return(fp_Bin);
}

/*========================================================================+
|                               PRINTOUT                                  |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+========================================================================*/
void PRINTOUT(void)
{
	int I;
	printf("\n ksize = %d",HEADER.ksize);
	printf("\n ncon  = %d\n\n Please enter a key...\n\n",HEADER.ncon);
//	WAIT;
	for(I=1;I<=HEADER.ncon;++I)
				{ printf("\n %s = %24.15g",HEADER.cnam[I],HEADER.cval[I]);
						if((I%20)==0)
								{ printf("\n\n Please enter a key...\n\n"); // WAIT;
						}
				}

	for(I=0;I<=2;++I)
				printf("\n ss[%d] = %24.15g",I,HEADER.ss[I]);

	printf("\n AU    = %24.15g",HEADER.au);
	printf("\n EMRAT = %24.15g",HEADER.emrat);
	printf("\n DENUM = %d",HEADER.denum);
	printf("\n\n Please enter a key...\n\n"); // WAIT;

	for(I = 0;I <= 35; ++I)
				printf("\n ipt[%d]=%d",I,HEADER.ipt[I]);
	printf("\n\n Please enter a key...\n\n"); // WAIT;
	for(I=0;I<=2;++I)
				printf("\nlpt[%d]=%d",I,HEADER.lpt[I]);
	printf("\n\n Please wait the program is searching for");
	printf(" the correct time span...\n\n");

}

/*========================================================================+
|                           SPLIT & DINT                                  |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| This subroutine breaks a D.P. number into a D.P. integer and a D.P.     |
| fractional part.                                                        |
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|         TT = D.P. input number                                          |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|          FR = D.P. 2-word output array.                                 |
|               FR(1) contains integer part                               |
|               FR(2) contains fractional part                            |
|               For negative input numbers, FR(1) contains the next       |
|               more negative integer; FR(2) contains a positive fraction.|
+-------------------------------------------------------------------------+
| REFERENCE: JPL (FORTRAN CODE)                                           |
+========================================================================*/
double Dint(double x)
{ return( (x>0) ? floor(x) : -floor(-x) );}

/*=======================================================================*/

int SPLIT(double tt,double *fr)
{
		/* Function Body */
		fr[0] = Dint(tt);
		fr[1] = tt - fr[0];
		if(tt < 0. && fr[1] != 0.)
				{ /* Make adjustments for negative input number. */
						fr[0] += -1.;
						fr[1] += 1.;
				}
		return 0;
}


/*========================================================================+
|                               DMOD                                      |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| DOUBLE PRECISION MODULA                                                 |
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|          x = NUMERATOR                                                  |
|          y = DENOMINATOR                                                |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|         DMOD = D.P. MODULA                                              |
+-------------------------------------------------------------------------+
| REFERENCE: %                                                            |
+========================================================================*/
double Dmod(double x,double y)
{
		double quotient;

		quotient = x / y;
		if(quotient >= 0)
					quotient = floor(quotient);
		else
					quotient = -floor(-quotient);
		return(x - y * quotient );
}



/*========================================================================+
|                               INTERP                                    |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| This subroutine differentiates and interpolates a set of Chebyshev      |
| coefficients to give position and velocity.                             |
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|  BUF   1st location of array of D.P. Chebyshev coefficients of position |      Ó‰ÌÓÏÂÌ˚È Ï‡ÒÒË‚ ÍÓ˝ÙËˆËÂÌÚÓ‚ ÔÓÎËÌÓÏ‡ ◊Â·˚¯Ó‚‡ 
|    T   T(1) IS DP FRACTIONAL TIME IN INTERVAL COVERED BY                |      ‚ÂÏˇ
|             COEFFICIENTS AT WHICH INTERPOLATION IS WANTED               |      T(1) - ‚ÂÏˇ ‚ ËÌÚÂ‚‡ÎÂ (0 .LE. T(1) .LE. 1)
|             (0 .LE. T(1) .LE. 1).  T(2) IS DP LENGTH OF WHOLE           |      T(2) - ‰ÎËÌÌ‡ ËÚÂ‚‡Î‡
|             INTERVAL IN INPUT TIME UNITS.                               |
|  NCF   # OF COEFFICIENTS PER COMPONENT                                  |      ˜ËÒÎÓ ÍÓÏÔÓÌÂÌÚ ‚ÂÍÚÓ‡ ‚ÒÂ„‰‡ 3
|  NCM   # OF COMPONENTS PER SET OF COEFFICIENTS                          |      ˜ËÒÎÓ ˜ÎÂÌÓ‚ ˇ‰‡
|   NA   # OF SETS OF COEFFICIENTS IN FULL ARRAY                          |      ÍÓÎË˜ÂÒÚ‚Ó ÒÛ·ËÌÚÂ‚‡ÎÓ‚
|          (I.E., # OF SUB-INTERVALS IN FULL INTERVAL)                    |
|   FL   INTEGER FLAG: =1 FOR POSITIONS ONLY                              |      ÙÎ‡„: =1 - ÚÓÎ¸ÍÓ ÔÓÎÓÊÂÌËÂ
|                      =2 FOR POS AND VEL                                 |            =2 - ÔÓÎÓÊÂÌËÂ Ë ÒÍÓÓÒÚ¸
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|  PV   INTERPOLATED QUANTITIES REQUESTED.  DIMENSION  EXPECTED IS        |      Ï‡ÒÒË‚ ˜ÎÂÌÓ‚ ˇ‰‡ ‡ÁÏÂ‡ PV(NCM,FL)
|       PV(NCM,FL), DP.                                                   |
+-------------------------------------------------------------------------+
| REFERENCE: JPL (FORTRAN CODE)                                           |
+========================================================================*/
int INTERP(double *buf,double *t,int ncf,int ncm,int na,int fl,double *pv)
{
	/* Initialized data */
	static int np = 2;
	static int nv = 3;
	static double twot = 0.;
	static double pc[18] = { 1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };
	static double vc[18] = { 0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };

	/* Local variables */
	double vfac, temp;
	int    i, j, l;
	double tc, dt1, dna;

	/* ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET */
	/* OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME    */
	/* WITHIN THAT SUBINTERVAL.                                  */
	dna = (double)na;
	dt1 = Dint(t[0]);
	temp = dna * t[0];
	l = (Dint) (temp - dt1)+1;

	/* TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1) */
	tc = (Dmod(temp, 1.) + dt1) * 2. - 1.;

	/* CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,    */
	/* AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.        */
	/* (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE */
	/* CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)     */
	if(tc != pc[1])
			{       np = 2;
					nv = 3;
					pc[1] = tc;
					twot = tc + tc;
			}

	/* BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED */
	/* AND ARE STORED IN THE ARRAY 'PC'.                           */
	if(np < ncf)
			{ for(i = np + 1; i <= ncf; ++i)
								pc[i - 1] = twot * pc[i - 2] - pc[i - 3];                 //  œÓÎÂÌÓÏ ◊Â·˚¯Ó‚‡ 1 Ó‰‡
					np = ncf;
			}

	/* INTERPOLATE TO GET POSITION FOR EACH COMPONENT */
	for(i = 1; i <= ncm; ++i)
				{ pv[i-1] = 0.;
						for(j = ncf; j >= 1; --j)
									pv[i-1]+=pc[j-1]*buf[j-1+(i-1+(l-1)*ncm)*ncf];         //  ÒÛÏÏËÓ‚‡ÌËÂ ˇ‰‡
				}
	if (fl <= 1)
				return 0;

	/* IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH    */
	/* DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED. */
	vfac = (dna + dna) / t[1];
	vc[2] = twot + twot;
	if(nv < ncf)
			{ for (i = nv + 1; i <= ncf; ++i)
								vc[i-1]=twot*vc[i-2]+pc[i-2]+pc[i-2]-vc[i-3];
					nv = ncf;
			}

	/* INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT */
	for(i = 1; i <= ncm; ++i)
				{ pv[i-1+ncm] = 0.;
						for (j =ncf; j >= 2; --j)
										pv[i-1+ncm]+=vc[j-1]*buf[j-1+(i-1+(l-1)*ncm)*ncf];
						pv[i-1+ncm] *= vfac;
				}
	return 0;
}

/*========================================================================+
|                               STATE                                     |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| THIS SUBROUTINE READS AND INTERPOLATES THE JPL PLANETARY EPHEMERIS FILE.|
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|         JED   DP 2-WORD JULIAN EPHEMERIS EPOCH AT WHICH INTERPOLATION   |
|               IS WANTED.  ANY COMBINATION OF JED(1)+JED(2) WHICH FALLS  |
|               WITHIN THE TIME SPAN ON THE FILE IS A PERMISSIBLE EPOCH.  |
|                A. FOR EASE IN PROGRAMMING, THE USER MAY PUT THE         |
|                   ENTIRE EPOCH IN JED(1) AND SET JED(2)=0.              |
|                B. FOR MAXIMUM INTERPOLATION ACCURACY, SET JED(1) =      |
|                   THE MOST RECENT MIDNIGHT AT OR BEFORE INTERPOLATION   |
|                   EPOCH AND SET JED(2) = FRACTIONAL PART OF A DAY       |
|                   ELAPSED BETWEEN JED(1) AND EPOCH.                     |
|                C. AS AN ALTERNATIVE, IT MAY PROVE CONVENIENT TO SET     |
|                   JED(1) = SOME FIXED EPOCH,SUCH AS START OF INTEGRATION|
|                   AND JED(2) = ELAPSED INTERVAL BETWEEN THEN AND EPOCH. |
|        LIST   12-WORD INTEGER ARRAY SPECIFYING WHAT INTERPOLATION       |
|               IS WANTED FOR EACH OF THE BODIES ON THE FILE.             |
|                         LIST(I)=0, NO INTERPOLATION FOR BODY I          |
|                                =1, POSITION ONLY                        |
|                                =2, POSITION AND VELOCITY                |
|               THE DESIGNATION OF THE ASTRONOMICAL BODIES BY I IS:       |
|                         I = 1: MERCURY                                  |
|                           = 2: VENUS                                    |
|                           = 3: EARTH-MOON BARYCENTER                    |
|                           = 4: MARS                                     |
|                           = 5: JUPITER                                  |
|                           = 6: SATURN                                   |
|                           = 7: URANUS                                   |
|                           = 8: NEPTUNE                                  |
|                           = 9: PLUTO                                    |
|                           =10: GEOCENTRIC MOON                          |
|                           =11: NUTATIONS IN LONGITUDE AND OBLIQUITY     |
|                           =12: LUNAR LIBRATIONS (IF ON FILE)            |
|                                                                         |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|          PV   DP 6 X 11 ARRAY THAT WILL CONTAIN REQUESTED INTERPOLATED  |
|               QUANTITIES.  THE BODY SPECIFIED BY LIST(I) WILL HAVE ITS  |
|               STATE IN THE ARRAY STARTING AT PV(1,I).  (ON ANY GIVEN    |
|               CALL, ONLY THOSE WORDS IN 'PV' WHICH ARE AFFECTED BY THE  |
|              FIRST 10 'LIST' ENTRIES (AND BY LIST(12) IF LIBRATIONS ARE |
|               ON THE FILE) ARE SET.  THE REST OF THE 'PV' ARRAY         |
|               IS UNTOUCHED.)  THE ORDER OF COMPONENTS STARTING IN       |
|               PV(1,I) IS: X,Y,Z,DX,DY,DZ.                               |
|               ALL OUTPUT VECTORS ARE REFERENCED TO THE EARTH MEAN       |
|               EQUATOR AND EQUINOX OF EPOCH. THE MOON STATE IS ALWAYS    |
|               GEOCENTRIC; THE OTHER NINE STATES ARE EITHER HELIOCENTRIC |
|               OR SOLAR-SYSTEM BARYCENTRIC, DEPENDING ON THE SETTING OF  |
|               COMMON FLAGS (SEE BELOW).                                 |
|               LUNAR LIBRATIONS, IF ON 12, ARE PUT INTO PV(K,11) IF      |
|               LIST(12) IS 1 OR 2.                                       |
|         NUT   DP 4-WORD ARRAY THAT WILL CONTAIN NUTATIONS AND RATES,    |
|               DEPENDING ON THE SETTING OF LIST(11).  THE ORDER OF       |
|               QUANTITIES IN NUT IS:                                     |
|                        D PSI  (NUTATION IN LONGITUDE)                   |
|                        D EPSILON (NUTATION IN OBLIQUITY)                |
|                        D PSI DOT                                        |
|                        D EPSILON DOT                                    |
+-------------------------------------------------------------------------+
| REFERENCE: JPL (FORTRAN CODE)                                           |
+========================================================================*/
int STATE(FILE *fp_Bin,double *jed,int *list,double *pv,double *nut)
{
		/* Initialized data */
		static double aufac  = 1.;
		static int FIRST = TRUE_;
		static long nrl  = 0L;
		static double t[2];

		/* Local variables */
		int i, j, k;
		double s,jd[4];
		long nr;
		/* 1ST TIME IN, GET POINTER DATA, ETC., FROM EPH FILE */
		if(FIRST)
			{ FIRST = FALSE_;
				if(stcomm.km)
					{ t[1] = HEADER.ss[2] * 86400.;}
				else
					{ t[1] = HEADER.ss[2];
						aufac = 1. / HEADER.au;
					}
			}

		/* MAIN ENTRY POINT -- CHECK EPOCH AND READ RIGHT RECORD */
		s = jed[0] - .5;
		SPLIT(s, &jd[0]);
		SPLIT(jed[1], &jd[2]);
		jd[0] = jd[0] + jd[2] + .5;
		jd[1] += jd[3];
		SPLIT(jd[1], &jd[2]);
		jd[0] += jd[2];

		/* ERROR RETURN OF EPOCH OUT OF RANGE */
		if((jd[0] < HEADER.ss[0]) || ((jd[0] + jd[3]) > HEADER.ss[1]))
				{ printf("\n STATE: Epoch out of range.");
						exit(1);
				}

		/* CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL */
		nr  = (long) ((jd[0] - HEADER.ss[0]) / HEADER.ss[2])+1;
		if(jd[0]==HEADER.ss[1])
			--nr;

		t[0]= ((jd[0]-((double)(nr-1)*HEADER.ss[2]+HEADER.ss[0]))+jd[3])/
					HEADER.ss[2];

		/* READ CORRECT RECORD IF NOT IN CORE */
		if(nr != nrl)
				{ fseek(fp_Bin,(long)(nr+1)*HEADER.ksize*4/*sizeof(double)*/,SEEK_SET);
						nrl = nr;
						for(k=1;k<=HEADER.ksize/2;++k)
									fread(&epib.buf[k],sizeof(double),1, fp_Bin);
				}

		/* INTERPOLATE SSBARY SUN */
		INTERP(&epib.buf[HEADER.ipt[10]],t,HEADER.ipt[22],const3,
										HEADER.ipt[34],const2,stcomm.pvsun);

		for(i = 1; i <= 6; ++i)
						stcomm.pvsun[i - 1] *= aufac;

		/* CHECK AND INTERPOLATE WHICHEVER BODIES ARE REQUESTED */
		for(i = 1; i <= 10; ++i)
				{ if(list[i-1] > 0)
									{ if(HEADER.ipt[11+i] <= 0)
														ERROR(i, "th body requested - not on file");
											INTERP(&epib.buf[HEADER.ipt[i-1]],t,HEADER.ipt[11+i],
																			const3,HEADER.ipt[23+i],list[i-1],&pv[(i-1)*6]);

											for(j = 1; j <= list[i-1] * 3; ++j)
														{ if(i <= 9 && !stcomm.bary)
																			pv[6*i-7+j]=pv[6*i-7+j]*aufac-stcomm.pvsun[j-1];
																else
																			pv[j +6*i -7] *= aufac;

														}
									}
					}

		/* DO NUTATIONS IF REQUESTED (AND IF ON FILE) */
		if(list[10] > 0 && HEADER.ipt[23] > 0)
					INTERP(&epib.buf[HEADER.ipt[11]],t,HEADER.ipt[23],const2,
													HEADER.ipt[35],list[10],&nut[0]);

		/* GET LIBRATIONS IF REQUESTED (AND IF ON FILE) */
		if(HEADER.lpt[1] > 0 && list[11] > 0)
					INTERP(&epib.buf[HEADER.lpt[0]],t,HEADER.lpt[1],const3,
													HEADER.lpt[2],list[11],&pv[60]);

		/* THAT'S ALL */
		return 0;

} /* state_ */

/*========================================================================+
|                                PLEPH0                                   |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| THIS SUBROUTINE READS THE JPL PLANETARY EPHEMERIS & GIVES THE POSITION  |
| AND VELOCITY OF THE POINT 'TARG'WITH RESPECT TO 'CENT'.                 |
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|       JD = D.P. JULIAN EPHEMERIS DATE AT WHICH INTERPOLATION IS WANTED. |
|       ** NOTE THE ENTRY DPLEPH FOR A DOUBLY-DIMENSIONED TIME **         |
|          THE REASON FOR THIS OPTION IS DISCUSSED IN THE SUB. STATE      |
|      TARG = int NUMBER OF 'TARGET' POINT.                               |
|      CENT = int NUMBER OF CENTER POINT.                                 |
|             THE NUMBERING CONVENTION FOR 'TARG' AND 'CENT' IS:          |
|                 1 = MERCURY         8 = NEPTUNE                         |
|                 2 = VENUS           9 = PLUTO                           |
|                 3 = EARTH          10 = MOON                            |
|                 4 = MARS           11 = SUN                             |
|                 5 = JUPITER        12 = SOLAR-SYSTEM BARYCENTER         |
|                 6 = SATURN         13 = EARTH-MOON BARYCENTER           |
|                 7 = URANUS         14 = NUTATIONS (intITUDE AND OBLIQ)  |
|                15 = LIBRATIONS, IF ON EPH FILE                          |
|               (IF NUTATIONS ARE WANTED, SET TARG = 14. FOR LIBRATIONS,  |
|               SET TARG = 15. 'CENT' WILL BE IGNORED ON EITHER CALL.)    |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|       RRD = OUTPUT 6-WORD D.P. ARRAY CONTAINING POSITION AND VELOCITY   |
|             OF POINT 'TARG' RELATIVE TO 'CENT'. THE UNITS ARE AU AND    |
|             AU/DAY. FOR LIBRATIONS THE UNITS ARE RADIANS AND RADIANS    |
|             PER DAY. IN THE CASE OF NUTATIONS THE FIRST FOUR WORDS OF   |
|             RRD WILL BE SET TO NUTATIONS AND RATES, HAVING UNITS OF     |
|             RADIANS AND RADIANS/DAY.                                    |
|             NOTE: IN MANY CASES THE USER WILL NEED ONLY POSITION        |
|                   VALUES FOR EPHEMERIDES OR NUTATIONS. FOR              |
|                   POSITION-ONLY OUTPUT, THE int VARIABLE 'IPV'          |
|                   IN THE COMMON AREA /PLECOM/ SHOULD BE SET = 1         |
|                   BEFORE THE NEXT CALL TO PLEPH. (ITS DEFAULT           |
|                   VALUE IS 2, WHICH RETURNS BOTH POSITIONS AND          |
|                   RATES.)                                               |
|                                                                         |
|      INSIDE is .TRUE. if the input Julian Ephemeris Date (JD) is within |
|             the ephemeris time span.  If not, INSIDE is set to .FALSE.  |
+-------------------------------------------------------------------------+
| REFERENCE: JPL (Fortran code)                                           |
+========================================================================*/
int PLEPH0(FILE *fp_Bin,int n,double jd,int targ,int cent,double *rrd,
											int *inside,double *jd2)
{
		/* INITIALIZED DATA: */
		double     fac =   0.;
		double embf[2] = { -1.,1. };
		double  pv[78] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
		static int  nemb = 1;
		static int FIRST = TRUE_;
		static int bsave;
		static double ve[2];

		/* DECLARATION: */
		int list[12] = { 0,0,0,0,0,0,0,0,0,0,0,0 };
		int llst[13] = { 1,2,10,4,5,6,7,8,9,10,11,11,3 };

		int l[2] = {0}, tc[2], i, lme;
		double jed[2],jdtot;


/*		if(targ==14) cent= 0; *//* 1999/07/25 */
		
		/* 1ST TIME IN, BE SURE EPHEMERIS IS INITIALIZED */
		if(FIRST)
				{ FIRST=FALSE_;
						ve[0] = 1. / (HEADER.emrat + 1.);
						ve[1] = HEADER.emrat * ve[0];
				}

		if(n==1)
				{ /* ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT */
						/* (SEE THE DISCUSSION IN THE SUBROUTINE STATE)              */
						jed[0] = jd2[0];
						jed[1] = jd2[1];
				}
		else
				{ /* INITIALIZE JED FOR 'STATE' AND SET UP COMPONENT COUNT */
						jed[0] = jd;
						jed[1] = 0.;
				}

		jdtot = jed[0] + jed[1];

		if(jdtot < HEADER.ss[0] || jdtot > HEADER.ss[1])
				{ *inside = FALSE_;
						return(0);
				}

		*inside = TRUE_;

		/* CHECK FOR NUTATION CALL */
		if(targ == 14)
				{ if(HEADER.ipt[34] > 0)
								{ list[10] = 2;
										STATE(fp_Bin,jed, list, pv, rrd);
										list[10] = 0;
										return(0);
								}
						else
								{ printf("\n *** NO NUTATION ON THE EPHEMERIS FILE ***");
										exit(1);
								}
				}

		/* CHECK FOR LIBRATIONS */
		if(targ == 15)
				{ if(HEADER.lpt[1] > 0)
								{ list[11] = 2;
										STATE(fp_Bin,jed, list, pv,rrd);
										list[11] = 0;
										for(i = 1; i <= 6; ++i)
														rrd[i-1] = pv[i + 59];
										return(0);
								}
						else
								{ printf("\n *** NO LIBRATIONS ON THE EPHEMERIS FILE ***");
										exit(1);
								}
				}

		/* CHECK FOR TARGET POINT = CENTER POINT */
		if(targ == cent)
				{ for(i = 1; i <= 6; ++i)
									rrd[i-1] = 0.;
						return(0);
				}

		/* FORCE BARYCENTRIC OUTPUT BY 'STATE' */
		bsave = stcomm.bary;
		stcomm.bary = TRUE_;

		/* SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL */
		tc[0] = targ;
		tc[1] = cent;
		lme = 0;

		for(i = 1; i <= 2; ++i)
					{ 
			l[i - 1] = llst[tc[i - 1] - 1];
							if( l[i - 1] < 11)
								list[l[i - 1] - 1] = 2;
							if(tc[i - 1] == 3)
									{ lme = 3;
											fac = -ve[0];
									}
							else if( tc[i - 1] == 10)
														{ lme = 10; fac = ve[1];}
							else if( tc[i - 1] == 13)
														{ nemb = i; }
					}

		if((list[9] == 2) && (l[0] != l[1]))
		list[2] = 2 - list[2];

		/* MAKE CALL TO STATE */
		STATE(fp_Bin,jed, list, pv, rrd);



		/* CASE: EARTH-TO-MOON */
		if((targ == 10) && (cent == 3))
				{ for(i = 1; i <= 6; ++i)
									rrd[i-1] = pv[i + 53];
				}

		/* CASE: MOON-TO-EARTH */
		else if((targ == 3) && (cent == 10))
									{ for(i = 1; i <= 6; ++i)
															rrd[i-1] = -pv[i + 53];
									}
		/* CASE: EMBARY TO MOON OR EARTH */
		else if((targ == 13 || cent == 13) && list[9] == 2)
									{ for(i = 1; i <= 6; ++i)
															rrd[i-1] = pv[i + 53] * fac * embf[nemb - 1];
									}
		/* OTHERWISE, GET EARTH OR MOON VECTOR AND THEN GET OUTPUT VECTOR */
		else
				{ for(i = 1; i <= 6; ++i)
									{ pv[i + 59] = stcomm.pvsun[i - 1];
											pv[i + 71] = pv[i + 11];
											if(lme > 0)
													   pv[i + lme * 6 - 7] = pv[i + 11] + fac * pv[i + 53];
											rrd[i-1] = pv[i + targ*6 - 7] - pv[i + cent * 6 - 7];
									}
				}

		/* CLEAR 'STATE' BODY ARRAY AND RESTORE BARYCENTER FLAG */
		list[2] = 0;
		list[l[0] - 1] = 0;
		list[l[1] - 1] = 0;
		stcomm.bary = bsave;

		/*  THAT'S ALL */
		return(0);
} /* pleph0 */

/*========================================================================+
|                                PLEPH                                    |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+========================================================================*/
int PLEPH(FILE *fp_Bin,double jd,int targ,int cent,double *rrd,int *inside)
{
		return PLEPH0(fp_Bin,0,jd,targ,cent,rrd,inside,(double *)0);
}

/*========================================================================+
|                               DPLEPH                                    |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+========================================================================*/
int DPLEPH(FILE *fp_Bin,double *jd2,int targ,int cent,double *rrd,int *inside)
{
//		return PLEPH0(fp_Bin,1,(double)0,targ,cent,rrd,(int *)0,jd2);
		return PLEPH0(fp_Bin,1,(double)0,targ,cent,rrd,inside,jd2);
}




/*========================================================================+
|                               INTERP                                    |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| This subroutine differentiates and interpolates a set of Chebyshev      |
| coefficients to give position and velocity.                             |
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|  BUF   1st location of array of D.P. Chebyshev coefficients of position |
|    T   T(1) IS DP FRACTIONAL TIME IN INTERVAL COVERED BY                |
|             COEFFICIENTS AT WHICH INTERPOLATION IS WANTED               |
|             (0 .LE. T(1) .LE. 1).  T(2) IS DP LENGTH OF WHOLE           |
|             INTERVAL IN INPUT TIME UNITS.                               |
|  NCF   # OF COEFFICIENTS PER COMPONENT                                  |
|  NCM   # OF COMPONENTS PER SET OF COEFFICIENTS                          |
|   NA   # OF SETS OF COEFFICIENTS IN FULL ARRAY                          |
|          (I.E., # OF SUB-INTERVALS IN FULL INTERVAL)                    |
|   FL   INTEGER FLAG: =1 FOR POSITIONS ONLY                              |
|                      =2 FOR POS AND VEL                                 |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|  PV   INTERPOLATED QUANTITIES REQUESTED.  DIMENSION  EXPECTED IS        |
|       PV(NCM,FL), DP.                                                   |
+-------------------------------------------------------------------------+
| REFERENCE: JPL (FORTRAN CODE)                                           |
+========================================================================*/
int INTERPx9(double *buf,double *t,int ncf,int ncm,int na,int fl,double *pv)
{
	/* Initialized data */
	static int np = 2;
	static int nv = 3;
	static int nA = 4;
	static double twot = 0.;
	static double pc[18] = { 1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };
	static double vc[18] = { 0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };
	static double ac[18] = { 0.,0.,4.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };

	/* Local variables */
	double vfac, afac, temp;
	int    i, j, l;
	double tc, dt1, dna;

	/* ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET */
	/* OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME    */
	/* WITHIN THAT SUBINTERVAL.                                  */
	dna = (double)na;
	dt1 = Dint(t[0]);
	temp = dna * t[0];
	l = (Dint) (temp - dt1)+1;

	/* TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1) */
	tc = (Dmod(temp, 1.) + dt1) * 2. - 1.;

	/* CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,    */
	/* AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.        */
	/* (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE */
	/* CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)     */
	if(tc != pc[1])
			{ np = 2;
					nv = 3;
					pc[1] = tc;
					twot = tc + tc;
			}

	/* BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED */
	/* AND ARE STORED IN THE ARRAY 'PC'.                           */
	if(np < ncf)
			{ for(i = np + 1; i <= ncf; ++i)
								pc[i - 1] = twot * pc[i - 2] - pc[i - 3];
					np = ncf;
			}

	/* INTERPOLATE TO GET POSITION FOR EACH COMPONENT */
	for(i = 1; i <= ncm; ++i)
				{ pv[i-1] = 0.;
						for(j = ncf; j >= 1; --j)
									pv[i-1]+=pc[j-1]*buf[j-1+(i-1+(l-1)*ncm)*ncf];
				}
	if (fl <= 1)
				return 0;

	/* IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH    */
	/* DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED. */
	vfac = (dna + dna) / t[1];
	vc[2] = twot + twot;
	if(nv < ncf)
			{ for (i = nv + 1; i <= ncf; ++i)
								vc[i-1]=twot*vc[i-2]+pc[i-2]+pc[i-2]-vc[i-3];
//                              vc[i-1]=(twot*vc[i-2]-vc[i-3])*(i-1);
					nv = ncf;
			}

	/* INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT */
	for(i = 1; i <= ncm; ++i)
				{ pv[i-1+ncm] = 0.;
						for (j =ncf; j >= 2; --j)
										pv[i-1+ncm]+=vc[j-1]*buf[j-1+(i-1+(l-1)*ncm)*ncf];
						pv[i-1+ncm] *= vfac;
				}

	/* IF ACCELERATION INTERPOLATION IS WANTED, BE SURE ENOUGH  */
	/* DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.   */
	afac = vfac*vfac;
//	afac = vfac / t[1];
	ac[3] = 24*tc;
	if(nA < ncf)
			{ for (i = nA + 1; i <= ncf; ++i)
//								ac[i-1]=twot*ac[i-2]+pc[i-2]+pc[i-2]-ac[i-3];
//                              ac[i-1]=(i*pc[i-1]-vc[i-1])*(i-1)/(tc*tc-1);
//                              ac[i-1]=(twot*ac[i-2]-ac[i-3])*(i-1)*(i-2);
//                              ac[i-1]=(tc*vc[i-1]-(i-1)*(i-1)*pc[i-1])/(1-tc*tc);
                                ac[i-1]=twot*ac[i-2]+4*vc[i-2]-ac[i-3];
					nA = ncf;
			}

	/* INTERPOLATE TO GET ACCELERATION FOR EACH COMPONENT */
	for(i = 1; i <= ncm; ++i)
				{ pv[i-1+2*ncm] = 0.;
						for (j =ncf; j >= 3; --j)
										pv[i-1+2*ncm]+=ac[j-1]*buf[j-1+(i-1+(l-1)*ncm)*ncf];
						pv[i-1+2*ncm] = pv[i-1+2*ncm]*afac;
				}

	return 0;
}

/*========================================================================+
|                               STATE                                     |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| THIS SUBROUTINE READS AND INTERPOLATES THE JPL PLANETARY EPHEMERIS FILE.|
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|         JED   DP 2-WORD JULIAN EPHEMERIS EPOCH AT WHICH INTERPOLATION   |
|               IS WANTED.  ANY COMBINATION OF JED(1)+JED(2) WHICH FALLS  |
|               WITHIN THE TIME SPAN ON THE FILE IS A PERMISSIBLE EPOCH.  |
|                A. FOR EASE IN PROGRAMMING, THE USER MAY PUT THE         |
|                   ENTIRE EPOCH IN JED(1) AND SET JED(2)=0.              |
|                B. FOR MAXIMUM INTERPOLATION ACCURACY, SET JED(1) =      |
|                   THE MOST RECENT MIDNIGHT AT OR BEFORE INTERPOLATION   |
|                   EPOCH AND SET JED(2) = FRACTIONAL PART OF A DAY       |
|                   ELAPSED BETWEEN JED(1) AND EPOCH.                     |
|                C. AS AN ALTERNATIVE, IT MAY PROVE CONVENIENT TO SET     |
|                   JED(1) = SOME FIXED EPOCH,SUCH AS START OF INTEGRATION|
|                   AND JED(2) = ELAPSED INTERVAL BETWEEN THEN AND EPOCH. |
|        LIST   12-WORD INTEGER ARRAY SPECIFYING WHAT INTERPOLATION       |
|               IS WANTED FOR EACH OF THE BODIES ON THE FILE.             |
|                         LIST(I)=0, NO INTERPOLATION FOR BODY I          |
|                                =1, POSITION ONLY                        |
|                                =2, POSITION AND VELOCITY                |
|               THE DESIGNATION OF THE ASTRONOMICAL BODIES BY I IS:       |
|                         I = 1: MERCURY                                  |
|                           = 2: VENUS                                    |
|                           = 3: EARTH-MOON BARYCENTER                    |
|                           = 4: MARS                                     |
|                           = 5: JUPITER                                  |
|                           = 6: SATURN                                   |
|                           = 7: URANUS                                   |
|                           = 8: NEPTUNE                                  |
|                           = 9: PLUTO                                    |
|                           =10: GEOCENTRIC MOON                          |
|                           =11: NUTATIONS IN LONGITUDE AND OBLIQUITY     |
|                           =12: LUNAR LIBRATIONS (IF ON FILE)            |
|                                                                         |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|          PV   DP 6 X 11 ARRAY THAT WILL CONTAIN REQUESTED INTERPOLATED  |
|               QUANTITIES.  THE BODY SPECIFIED BY LIST(I) WILL HAVE ITS  |
|               STATE IN THE ARRAY STARTING AT PV(1,I).  (ON ANY GIVEN    |
|               CALL, ONLY THOSE WORDS IN 'PV' WHICH ARE AFFECTED BY THE  |
|              FIRST 10 'LIST' ENTRIES (AND BY LIST(12) IF LIBRATIONS ARE |
|               ON THE FILE) ARE SET.  THE REST OF THE 'PV' ARRAY         |
|               IS UNTOUCHED.)  THE ORDER OF COMPONENTS STARTING IN       |
|               PV(1,I) IS: X,Y,Z,DX,DY,DZ.                               |
|               ALL OUTPUT VECTORS ARE REFERENCED TO THE EARTH MEAN       |
|               EQUATOR AND EQUINOX OF EPOCH. THE MOON STATE IS ALWAYS    |
|               GEOCENTRIC; THE OTHER NINE STATES ARE EITHER HELIOCENTRIC |
|               OR SOLAR-SYSTEM BARYCENTRIC, DEPENDING ON THE SETTING OF  |
|               COMMON FLAGS (SEE BELOW).                                 |
|               LUNAR LIBRATIONS, IF ON 12, ARE PUT INTO PV(K,11) IF      |
|               LIST(12) IS 1 OR 2.                                       |
|         NUT   DP 4-WORD ARRAY THAT WILL CONTAIN NUTATIONS AND RATES,    |
|               DEPENDING ON THE SETTING OF LIST(11).  THE ORDER OF       |
|               QUANTITIES IN NUT IS:                                     |
|                        D PSI  (NUTATION IN LONGITUDE)                   |
|                        D EPSILON (NUTATION IN OBLIQUITY)                |
|                        D PSI DOT                                        |
|                        D EPSILON DOT                                    |
+-------------------------------------------------------------------------+
| REFERENCE: JPL (FORTRAN CODE)                                           |
+========================================================================*/
int STATEx9(FILE *fp_Bin,double *jed,int *list,double *pv,double *nut)
{
		/* Initialized data */
		static double aufac  = 1.;
		static int FIRST = TRUE_;
		static long nrl  = 0L;
		static double t[2]                     , PV[9], hh = 1.e-13;

		/* Local variables */
		int i, j, k;
		double s,jd[4];
		long nr;
		/* 1ST TIME IN, GET POINTER DATA, ETC., FROM EPH FILE */
		if(FIRST)
			{ FIRST = FALSE_;
				if(stcomm.km)
					{ t[1] = HEADER.ss[2] * 86400.;}
				else
					{ t[1] = HEADER.ss[2];
						aufac = 1. / HEADER.au;
					}
			}

		/* MAIN ENTRY POINT -- CHECK EPOCH AND READ RIGHT RECORD */
		s = jed[0] - .5;
		SPLIT(s, &jd[0]);
		SPLIT(jed[1], &jd[2]);
		jd[0] = jd[0] + jd[2] + .5;
		jd[1] += jd[3];
		SPLIT(jd[1], &jd[2]);
		jd[0] += jd[2];

		/* ERROR RETURN OF EPOCH OUT OF RANGE */
		if((jd[0] < HEADER.ss[0]) || ((jd[0] + jd[3]) > HEADER.ss[1]))
				{ printf("\n STATE: Epoch out of range.");
						exit(1);
				}

		/* CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL */
		nr  = (long) ((jd[0] - HEADER.ss[0]) / HEADER.ss[2])+1;
		if(jd[0]==HEADER.ss[1])
			--nr;

		t[0]= ((jd[0]-((double)(nr-1)*HEADER.ss[2]+HEADER.ss[0]))+jd[3])/
					HEADER.ss[2];

		/* READ CORRECT RECORD IF NOT IN CORE */
		if(nr != nrl)
				{ fseek(fp_Bin,(long)(nr+1)*HEADER.ksize*4/*sizeof(double)*/,SEEK_SET);
						nrl = nr;
						for(k=1;k<=HEADER.ksize/2;++k)
									fread(&epib.buf[k],sizeof(double),1, fp_Bin);
				}

		/* INTERPOLATE SSBARY SUN */
		INTERPx9(&epib.buf[HEADER.ipt[10]],t,HEADER.ipt[22],const3,HEADER.ipt[34],const2,stcomm.pvsun);

//		for(i = 1; i <= 6; ++i)
		for(i = 1; i <= 9; ++i)
						stcomm.pvsun[i - 1] *= aufac;

		/* CHECK AND INTERPOLATE WHICHEVER BODIES ARE REQUESTED */
		for(i = 1; i <= 10; ++i)
				{ if(list[i-1] > 0)
									{ if(HEADER.ipt[11+i] <= 0)
														ERROR(i, "th body requested - not on file");
//											INTERP(&epib.buf[HEADER.ipt[i-1]],t,HEADER.ipt[11+i],const3,HEADER.ipt[23+i],list[i-1],&pv[(i-1)*6]);
											INTERPx9(&epib.buf[HEADER.ipt[i-1]],t,HEADER.ipt[11+i],const3,HEADER.ipt[23+i],list[i-1],&pv[(i-1)*9]);

											t[0] += hh;
											INTERPx9(&epib.buf[HEADER.ipt[i-1]],t,HEADER.ipt[11+i],const3,HEADER.ipt[23+i],list[i-1],PV);

											for(j = 1; j <= list[i-1] * 3; ++j)
														{ 
															    if(i <= 9 && !stcomm.bary)
//																			pv[6*i-7 +j]=pv[6*i-7 +j]*aufac-stcomm.pvsun[j-1];
																			pv[9*i-10+j]=pv[9*i-10+j]*aufac-stcomm.pvsun[j-1];
																else
//																			pv[j +6*i -7 ] *= aufac;
																			pv[j +9*i -10] *= aufac;
														}
									}
					}

		/* DO NUTATIONS IF REQUESTED (AND IF ON FILE) */
		if(list[10] > 0 && HEADER.ipt[23] > 0)
					INTERPx9(&epib.buf[HEADER.ipt[11]],t,HEADER.ipt[23],const2,HEADER.ipt[35],list[10],&nut[0]);

		/* GET LIBRATIONS IF REQUESTED (AND IF ON FILE) */
		if(HEADER.lpt[1] > 0 && list[11] > 0)
//					INTERP(&epib.buf[HEADER.lpt[0]],t,HEADER.lpt[1],const3,HEADER.lpt[2],list[11],&pv[60]);
					INTERPx9(&epib.buf[HEADER.lpt[0]],t,HEADER.lpt[1],const3,HEADER.lpt[2],list[11],&pv[90]);

		/* THAT'S ALL */
		return 0;

} /* state_ */

/*========================================================================+
|                                PLEPH0                                   |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| THIS SUBROUTINE READS THE JPL PLANETARY EPHEMERIS & GIVES THE POSITION  |
| AND VELOCITY OF THE POINT 'TARG'WITH RESPECT TO 'CENT'.                 |
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|       JD = D.P. JULIAN EPHEMERIS DATE AT WHICH INTERPOLATION IS WANTED. |
|       ** NOTE THE ENTRY DPLEPH FOR A DOUBLY-DIMENSIONED TIME **         |
|          THE REASON FOR THIS OPTION IS DISCUSSED IN THE SUB. STATE      |
|      TARG = int NUMBER OF 'TARGET' POINT.                               |
|      CENT = int NUMBER OF CENTER POINT.                                 |
|             THE NUMBERING CONVENTION FOR 'TARG' AND 'CENT' IS:          |
|                 1 = MERCURY         8 = NEPTUNE                         |
|                 2 = VENUS           9 = PLUTO                           |
|                 3 = EARTH          10 = MOON                            |
|                 4 = MARS           11 = SUN                             |
|                 5 = JUPITER        12 = SOLAR-SYSTEM BARYCENTER         |
|                 6 = SATURN         13 = EARTH-MOON BARYCENTER           |
|                 7 = URANUS         14 = NUTATIONS (intITUDE AND OBLIQ)  |
|                15 = LIBRATIONS, IF ON EPH FILE                          |
|               (IF NUTATIONS ARE WANTED, SET TARG = 14. FOR LIBRATIONS,  |
|               SET TARG = 15. 'CENT' WILL BE IGNORED ON EITHER CALL.)    |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|       RRD = OUTPUT 6-WORD D.P. ARRAY CONTAINING POSITION AND VELOCITY   |
|             OF POINT 'TARG' RELATIVE TO 'CENT'. THE UNITS ARE AU AND    |
|             AU/DAY. FOR LIBRATIONS THE UNITS ARE RADIANS AND RADIANS    |
|             PER DAY. IN THE CASE OF NUTATIONS THE FIRST FOUR WORDS OF   |
|             RRD WILL BE SET TO NUTATIONS AND RATES, HAVING UNITS OF     |
|             RADIANS AND RADIANS/DAY.                                    |
|             NOTE: IN MANY CASES THE USER WILL NEED ONLY POSITION        |
|                   VALUES FOR EPHEMERIDES OR NUTATIONS. FOR              |
|                   POSITION-ONLY OUTPUT, THE int VARIABLE 'IPV'          |
|                   IN THE COMMON AREA /PLECOM/ SHOULD BE SET = 1         |
|                   BEFORE THE NEXT CALL TO PLEPH. (ITS DEFAULT           |
|                   VALUE IS 2, WHICH RETURNS BOTH POSITIONS AND          |
|                   RATES.)                                               |
|                                                                         |
|      INSIDE is .TRUE. if the input Julian Ephemeris Date (JD) is within |
|             the ephemeris time span.  If not, INSIDE is set to .FALSE.  |
+-------------------------------------------------------------------------+
| REFERENCE: JPL (Fortran code)                                           |
+========================================================================*/
int PLEPH0x9(FILE *fp_Bin,int n,double jd,int targ,int cent,double *rrd,int *inside,double *jd2)
{
		/* INITIALIZED DATA: */
		double     fac =   0.;
		double embf[2] = { -1.,1. };
//		double  pv[78] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
		double  pv[13*9] = { 0. };
		static int  nemb = 1;
		static int FIRST = TRUE_;
		static int bsave;
		static double ve[2];

		/* DECLARATION: */
		int list[12] = { 0,0,0,0,0,0,0,0,0,0,0,0 };
		int llst[13] = { 1,2,10,4,5,6,7,8,9,10,11,11,3 };

		int l[2],tc[2],i,lme;
		double jed[2],jdtot;

/*		if(targ==14) cent= 0; *//* 1999/07/25 */
		
		/* 1ST TIME IN, BE SURE EPHEMERIS IS INITIALIZED */
		if(FIRST)
				{ FIRST=FALSE_;
						ve[0] = 1. / (HEADER.emrat + 1.);
						ve[1] = HEADER.emrat * ve[0];
				}

		if(n==1)
				{ /* ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT */
						/* (SEE THE DISCUSSION IN THE SUBROUTINE STATE)              */
						jed[0] = jd2[0];
						jed[1] = jd2[1];
				}
		else
				{ /* INITIALIZE JED FOR 'STATE' AND SET UP COMPONENT COUNT */
						jed[0] = jd;
						jed[1] = 0.;
				}

		jdtot = jed[0] + jed[1];

		if(jdtot < HEADER.ss[0] || jdtot > HEADER.ss[1])
				{ *inside = FALSE_;
						return(0);
				}

		*inside = TRUE_;

		/* CHECK FOR NUTATION CALL */
		if(targ == 14)
				{ if(HEADER.ipt[34] > 0)
								{ list[10] = 3;
										STATEx9(fp_Bin,jed, list, pv, rrd);
										list[10] = 0;
										return(0);
								}
						else
								{ printf("\n *** NO NUTATION ON THE EPHEMERIS FILE ***");
										exit(1);
								}
				}

		/* CHECK FOR LIBRATIONS */
		if(targ == 15)
				{ if(HEADER.lpt[1] > 0)
								{ list[11] = 3;
										STATEx9(fp_Bin,jed, list, pv,rrd);
										list[11] = 0;
//										for(i = 1; i <= 6; ++i)
//														rrd[i-1] = pv[i + 59];
										for(i = 1; i <= 9; ++i)
										                rrd[i-1] = pv[i + 10*9-1];
										return(0);
								}
						else
								{ printf("\n *** NO LIBRATIONS ON THE EPHEMERIS FILE ***");
										exit(1);
								}
				}

		/* CHECK FOR TARGET POINT = CENTER POINT */
		if(targ == cent)
				{ 
//					for(i = 1; i <= 6; ++i)
					for(i = 1; i <= 9; ++i)
									rrd[i-1] = 0.;
					return(0);
				}

		/* FORCE BARYCENTRIC OUTPUT BY 'STATE' */
		bsave = stcomm.bary;
		stcomm.bary = TRUE_;

		/* SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL */
		tc[0] = targ;
		tc[1] = cent;
		lme = 0;

		for(i = 1; i <= 2; ++i)
					{ l[i - 1] = llst[tc[i - 1] - 1];
							if( l[i - 1] < 11)list[l[i - 1] - 1] = 3;
							if(tc[i - 1] == 3)
									{ lme = 3;
											fac = -ve[0];
									}
							else if( tc[i - 1] == 10)
														{ lme = 10; fac = ve[1];}
							else if( tc[i - 1] == 13)
														{ nemb = i; }
					}

		if((list[9] == 3) && (l[0] != l[1]))
		list[2] = 3 - list[2];

		/* MAKE CALL TO STATE */
		STATEx9(fp_Bin,jed, list, pv, rrd);



		/* CASE: EARTH-TO-MOON */
		if((targ == 10) && (cent == 3))
				{ 
//					for(i = 1; i <= 6; ++i) rrd[i-1] = pv[i + 53];
					for(i = 1; i <= 9; ++i)	rrd[i-1] = pv[i + 9*9-1];
				}

		/* CASE: MOON-TO-EARTH */
		else if((targ == 3) && (cent == 10))
									{ 
//										for(i = 1; i <= 6; ++i) rrd[i-1] = -pv[i + 53];
										for(i = 1; i <= 9; ++i) rrd[i-1] = -pv[i + 9*9-1];
									}
		/* CASE: EMBARY TO MOON OR EARTH */
		else if((targ == 13 || cent == 13) && list[9] == 3)
									{ 
//										for(i = 1; i <= 6; ++i) rrd[i-1] = pv[i + 53   ] * fac * embf[nemb - 1];
										for(i = 1; i <= 9; ++i)	rrd[i-1] = pv[i + 9*9-1] * fac * embf[nemb - 1];
									}
		/* OTHERWISE, GET EARTH OR MOON VECTOR AND THEN GET OUTPUT VECTOR */
		else
				{ 
//					for(i = 1; i <= 6; ++i)
					for(i = 1; i <= 9; ++i)
									{ 
//										    pv[i + 59] = stcomm.pvsun[i - 1];
//											pv[i + 71] = pv[i + 11];
										    pv[i + 9*10-1] = stcomm.pvsun[i - 1];
											pv[i + 9*12-1] = pv[i + 9*2-1];
											if(lme > 0)
//													   pv[i + lme *6 - 7] = pv[i + 11] + fac * pv[i + 53];
//											rrd[i-1] = pv[i + targ*6 - 7] - pv[i + cent * 6 - 7];
													   pv[i + lme *9 -10] = pv[i + 17] + fac * pv[i + 9*9-1];
											rrd[i-1] = pv[i + targ*9 -10] - pv[i + cent * 9 - 10];
									}
				}

		/* CLEAR 'STATE' BODY ARRAY AND RESTORE BARYCENTER FLAG */
		list[2] = 0;
		list[l[0] - 1] = 0;
		list[l[1] - 1] = 0;
		stcomm.bary = bsave;

		/*  THAT'S ALL */
		return(0);
} /* pleph0 */

/*========================================================================+
|                                PLEPH                                    |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+========================================================================*/
int PLEPHx9(FILE *fp_Bin,double jd,int targ,int cent,double *rrd,int *inside)
{
		return PLEPH0x9(fp_Bin,0,jd,targ,cent,rrd,inside,(double *)0);
}

/*========================================================================+
|                               DPLEPH                                    |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+========================================================================*/
int DPLEPHx9(FILE *fp_Bin,double *jd2,int targ,int cent,double *rrd,int *inside)
{
//		return PLEPH0(fp_Bin,1,(double)0,targ,cent,rrd,(int *)0,jd2);
		return PLEPH0x9(fp_Bin,1,(double)0,targ,cent,rrd,inside,jd2);
}



/*========================================================================+
|                               INTERP                                    |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| This subroutine differentiates and interpolates a set of Chebyshev      |
| coefficients to give position and velocity.                             |
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|  BUF   1st location of array of D.P. Chebyshev coefficients of position |
|    T   T(1) IS DP FRACTIONAL TIME IN INTERVAL COVERED BY                |
|             COEFFICIENTS AT WHICH INTERPOLATION IS WANTED               |
|             (0 .LE. T(1) .LE. 1).  T(2) IS DP LENGTH OF WHOLE           |
|             INTERVAL IN INPUT TIME UNITS.                               |
|  NCF   # OF COEFFICIENTS PER COMPONENT                                  |
|  NCM   # OF COMPONENTS PER SET OF COEFFICIENTS                          |
|   NA   # OF SETS OF COEFFICIENTS IN FULL ARRAY                          |
|          (I.E., # OF SUB-INTERVALS IN FULL INTERVAL)                    |
|   FL   INTEGER FLAG: =1 FOR POSITIONS ONLY                              |
|                      =2 FOR POS AND VEL                                 |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|  PV   INTERPOLATED QUANTITIES REQUESTED.  DIMENSION  EXPECTED IS        |
|       PV(NCM,FL), DP.                                                   |
+-------------------------------------------------------------------------+
| REFERENCE: JPL (FORTRAN CODE)                                           |
+========================================================================*/
int INTERPx12(double *buf,double *t,int ncf,int ncm,int na,int fl,double *pv)
{
	/* Initialized data */
	static int np = 2;
	static int nv = 3;
	static int nA = 4;
	static int nJ = 5;
	static double twot = 0.;
	static double pc[18] = { 1.,0.,0., 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };
	static double vc[18] = { 0.,1.,0., 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };
	static double ac[18] = { 0.,0.,4., 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };
	static double jc[18] = { 0.,0.,0.,24.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. };

	/* Local variables */
	double vfac, afac, jfac, temp;
	int    i, j, l;
	double tc, dt1, dna;

	/* ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET */
	/* OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME    */
	/* WITHIN THAT SUBINTERVAL.                                  */
	dna = (double)na;
	dt1 = Dint(t[0]);
	temp = dna * t[0];
	l = (Dint) (temp - dt1)+1;

	/* TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1) */
	tc = (Dmod(temp, 1.) + dt1) * 2. - 1.;

	/* CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,    */
	/* AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.        */
	/* (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE */
	/* CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)     */
	if(tc != pc[1])
			{		np = 2;
					nv = 3;
					pc[1] = tc;
					twot = tc + tc;
			}

	/* BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED */
	/* AND ARE STORED IN THE ARRAY 'PC'.                           */
	if(np < ncf)
			{ for(i = np + 1; i <= ncf; ++i)
								pc[i - 1] = twot * pc[i - 2] - pc[i - 3];
					np = ncf;
			}

	/* INTERPOLATE TO GET POSITION FOR EACH COMPONENT */
	for(i = 1; i <= ncm; ++i)
				{ pv[i-1] = 0.;
						for(j = ncf; j >= 1; --j)
									pv[i-1]+=pc[j-1]*buf[j-1+(i-1+(l-1)*ncm)*ncf];
				}
	if (fl <= 1)
				return 0;

	/* IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH    */
	/* DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED. */
	vfac = (dna + dna) / t[1];
	vc[2] = twot + twot;
	if(nv < ncf)
			{ for (i = nv + 1; i <= ncf; ++i)
								vc[i-1]=twot*vc[i-2]+pc[i-2]+pc[i-2]-vc[i-3];
					nv = ncf;
			}

	/* INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT */
	for(i = 1; i <= ncm; ++i)
				{ pv[i-1+ncm] = 0.;
						for (j =ncf; j >= 2; --j)
										pv[i-1+ncm]+=vc[j-1]*buf[j-1+(i-1+(l-1)*ncm)*ncf];
						pv[i-1+ncm] *= vfac;
				}

	/* IF ACCELERATION INTERPOLATION IS WANTED, BE SURE ENOUGH  */
	/* DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.   */
	afac = vfac*vfac;
	ac[3] = 24*tc;
	if(nA < ncf)
			{ for (i = nA + 1; i <= ncf; ++i)
                                ac[i-1]=twot*ac[i-2]+4*vc[i-2]-ac[i-3];
					nA = ncf;
			}

	/* INTERPOLATE TO GET ACCELERATION FOR EACH COMPONENT */
	for(i = 1; i <= ncm; ++i)
				{ pv[i-1+2*ncm] = 0.;
						for (j =ncf; j >= 3; --j)
										pv[i-1+2*ncm]+=ac[j-1]*buf[j-1+(i-1+(l-1)*ncm)*ncf];
						pv[i-1+2*ncm] = pv[i-1+2*ncm]*afac;
				}

	/* IF JERK INTERPOLATION IS WANTED, BE SURE ENOUGH			*/
	/* DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.   */
	jfac = vfac*vfac*vfac;
	jc[4] = 184*tc;
	if(nJ < ncf)
			{ for (i = nJ + 1; i <= ncf; ++i)
                                jc[i-1]=twot*jc[i-2]+6*ac[i-2]-jc[i-3];
					nJ = ncf;
			}

	/* INTERPOLATE TO GET JERK FOR EACH COMPONENT */
	for(i = 1; i <= ncm; ++i)
				{ pv[i-1+3*ncm] = 0.;
						for (j =ncf; j >= 4; --j)
										pv[i-1+3*ncm]+=jc[j-1]*buf[j-1+(i-1+(l-1)*ncm)*ncf];
						pv[i-1+3*ncm] *= jfac;
				}

	return 0;
}

/*========================================================================+
|                               STATE                                     |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| THIS SUBROUTINE READS AND INTERPOLATES THE JPL PLANETARY EPHEMERIS FILE.|
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|         JED   DP 2-WORD JULIAN EPHEMERIS EPOCH AT WHICH INTERPOLATION   |
|               IS WANTED.  ANY COMBINATION OF JED(1)+JED(2) WHICH FALLS  |
|               WITHIN THE TIME SPAN ON THE FILE IS A PERMISSIBLE EPOCH.  |
|                A. FOR EASE IN PROGRAMMING, THE USER MAY PUT THE         |
|                   ENTIRE EPOCH IN JED(1) AND SET JED(2)=0.              |
|                B. FOR MAXIMUM INTERPOLATION ACCURACY, SET JED(1) =      |
|                   THE MOST RECENT MIDNIGHT AT OR BEFORE INTERPOLATION   |
|                   EPOCH AND SET JED(2) = FRACTIONAL PART OF A DAY       |
|                   ELAPSED BETWEEN JED(1) AND EPOCH.                     |
|                C. AS AN ALTERNATIVE, IT MAY PROVE CONVENIENT TO SET     |
|                   JED(1) = SOME FIXED EPOCH,SUCH AS START OF INTEGRATION|
|                   AND JED(2) = ELAPSED INTERVAL BETWEEN THEN AND EPOCH. |
|        LIST   12-WORD INTEGER ARRAY SPECIFYING WHAT INTERPOLATION       |
|               IS WANTED FOR EACH OF THE BODIES ON THE FILE.             |
|                         LIST(I)=0, NO INTERPOLATION FOR BODY I          |
|                                =1, POSITION ONLY                        |
|                                =2, POSITION AND VELOCITY                |
|               THE DESIGNATION OF THE ASTRONOMICAL BODIES BY I IS:       |
|                         I = 1: MERCURY                                  |
|                           = 2: VENUS                                    |
|                           = 3: EARTH-MOON BARYCENTER                    |
|                           = 4: MARS                                     |
|                           = 5: JUPITER                                  |
|                           = 6: SATURN                                   |
|                           = 7: URANUS                                   |
|                           = 8: NEPTUNE                                  |
|                           = 9: PLUTO                                    |
|                           =10: GEOCENTRIC MOON                          |
|                           =11: NUTATIONS IN LONGITUDE AND OBLIQUITY     |
|                           =12: LUNAR LIBRATIONS (IF ON FILE)            |
|                                                                         |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|          PV   DP 6 X 11 ARRAY THAT WILL CONTAIN REQUESTED INTERPOLATED  |
|               QUANTITIES.  THE BODY SPECIFIED BY LIST(I) WILL HAVE ITS  |
|               STATE IN THE ARRAY STARTING AT PV(1,I).  (ON ANY GIVEN    |
|               CALL, ONLY THOSE WORDS IN 'PV' WHICH ARE AFFECTED BY THE  |
|              FIRST 10 'LIST' ENTRIES (AND BY LIST(12) IF LIBRATIONS ARE |
|               ON THE FILE) ARE SET.  THE REST OF THE 'PV' ARRAY         |
|               IS UNTOUCHED.)  THE ORDER OF COMPONENTS STARTING IN       |
|               PV(1,I) IS: X,Y,Z,DX,DY,DZ.                               |
|               ALL OUTPUT VECTORS ARE REFERENCED TO THE EARTH MEAN       |
|               EQUATOR AND EQUINOX OF EPOCH. THE MOON STATE IS ALWAYS    |
|               GEOCENTRIC; THE OTHER NINE STATES ARE EITHER HELIOCENTRIC |
|               OR SOLAR-SYSTEM BARYCENTRIC, DEPENDING ON THE SETTING OF  |
|               COMMON FLAGS (SEE BELOW).                                 |
|               LUNAR LIBRATIONS, IF ON 12, ARE PUT INTO PV(K,11) IF      |
|               LIST(12) IS 1 OR 2.                                       |
|         NUT   DP 4-WORD ARRAY THAT WILL CONTAIN NUTATIONS AND RATES,    |
|               DEPENDING ON THE SETTING OF LIST(11).  THE ORDER OF       |
|               QUANTITIES IN NUT IS:                                     |
|                        D PSI  (NUTATION IN LONGITUDE)                   |
|                        D EPSILON (NUTATION IN OBLIQUITY)                |
|                        D PSI DOT                                        |
|                        D EPSILON DOT                                    |
+-------------------------------------------------------------------------+
| REFERENCE: JPL (FORTRAN CODE)                                           |
+========================================================================*/
int STATEx12(FILE *fp_Bin,double *jed,int *list,double *pv,double *nut)
{
		/* Initialized data */
		static double aufac  = 1.;
		static int FIRST = TRUE_;
		static long nrl  = 0L;
		static double t[2]; //                     , PV[12], hh = 1.e-13;

		/* Local variables */
		int i, j, k;
		double s,jd[4];
		long nr;
		/* 1ST TIME IN, GET POINTER DATA, ETC., FROM EPH FILE */
		if(FIRST)
			{ FIRST = FALSE_;
				if(stcomm.km)
					{ t[1] = HEADER.ss[2] * 86400.;}
				else
					{ t[1] = HEADER.ss[2];
						aufac = 1. / HEADER.au;
					}
			}

		/* MAIN ENTRY POINT -- CHECK EPOCH AND READ RIGHT RECORD */
		s = jed[0] - .5;
		SPLIT(s, &jd[0]);
		SPLIT(jed[1], &jd[2]);
		jd[0] = jd[0] + jd[2] + .5;
		jd[1] += jd[3];
		SPLIT(jd[1], &jd[2]);
		jd[0] += jd[2];

		/* ERROR RETURN OF EPOCH OUT OF RANGE */
		if((jd[0] < HEADER.ss[0]) || ((jd[0] + jd[3]) > HEADER.ss[1]))
				{ printf("\n STATE: Epoch out of range.");
						exit(1);
				}

		/* CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL */
		nr  = (long) ((jd[0] - HEADER.ss[0]) / HEADER.ss[2])+1;
		if(jd[0]==HEADER.ss[1])
			--nr;

		t[0]= ((jd[0]-((double)(nr-1)*HEADER.ss[2]+HEADER.ss[0]))+jd[3])/
					HEADER.ss[2];

		/* READ CORRECT RECORD IF NOT IN CORE */
		if(nr != nrl)
				{ fseek(fp_Bin,(long)(nr+1)*HEADER.ksize*4/*sizeof(double)*/,SEEK_SET);
						nrl = nr;
						for(k=1;k<=HEADER.ksize/2;++k)
									fread(&epib.buf[k],sizeof(double),1, fp_Bin);
				}

		/* INTERPOLATE SSBARY SUN */
		INTERPx12(&epib.buf[HEADER.ipt[10]],t,HEADER.ipt[22],const3,HEADER.ipt[34],const2,stcomm.pvsun);

//		for(i = 1; i <=  6; ++i)
		for(i = 1; i <= 12; ++i)
						stcomm.pvsun[i - 1] *= aufac;

		/* CHECK AND INTERPOLATE WHICHEVER BODIES ARE REQUESTED */
		for(i = 1; i <= 10; ++i)
				{ if(list[i-1] > 0)
									{ if(HEADER.ipt[11+i] <= 0)
														ERROR(i, "th body requested - not on file");
//											INTERP   (&epib.buf[HEADER.ipt[i-1]],t,HEADER.ipt[11+i],const3,HEADER.ipt[23+i],list[i-1],&pv[(i-1)*6 ]);
											INTERPx12(&epib.buf[HEADER.ipt[i-1]],t,HEADER.ipt[11+i],const3,HEADER.ipt[23+i],list[i-1],&pv[(i-1)*12]);

//											t[0] += hh;
//											INTERPx12(&epib.buf[HEADER.ipt[i-1]],t,HEADER.ipt[11+i],const3,HEADER.ipt[23+i],list[i-1],PV);

											for(j = 1; j <= list[i-1] * 3; ++j)
														{ 
															    if(i <= 9 && !stcomm.bary)
//																			pv[ 6*i-7 +j]=pv[ 6*i-7 +j]*aufac-stcomm.pvsun[j-1];
																			pv[12*i-13+j]=pv[12*i-13+j]*aufac-stcomm.pvsun[j-1];
																else
//																			pv[j + 6*i -7 ] *= aufac;
																			pv[j +12*i -13] *= aufac;
														}
									}
					}

		/* DO NUTATIONS IF REQUESTED (AND IF ON FILE) */
		if(list[10] > 0 && HEADER.ipt[23] > 0)
					INTERPx12(&epib.buf[HEADER.ipt[11]],t,HEADER.ipt[23],const2,HEADER.ipt[35],list[10],&nut[0]);

		/* GET LIBRATIONS IF REQUESTED (AND IF ON FILE) */
		if(HEADER.lpt[1] > 0 && list[11] > 0)
//					INTERP   (&epib.buf[HEADER.lpt[0]],t,HEADER.lpt[1],const3,HEADER.lpt[2],list[11],&pv[ 60]);
					INTERPx12(&epib.buf[HEADER.lpt[0]],t,HEADER.lpt[1],const3,HEADER.lpt[2],list[11],&pv[120]);

		/* THAT'S ALL */
		return 0;

} /* state_ */

/*========================================================================+
|                                PLEPH0                                   |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+-------------------------------------------------------------------------+
| DESCRIPTION:                                                            |
| THIS SUBROUTINE READS THE JPL PLANETARY EPHEMERIS & GIVES THE POSITION  |
| AND VELOCITY OF THE POINT 'TARG'WITH RESPECT TO 'CENT'.                 |
+-------------------------------------------------------------------------+
| INPUTS :                                                                |
|       JD = D.P. JULIAN EPHEMERIS DATE AT WHICH INTERPOLATION IS WANTED. |
|       ** NOTE THE ENTRY DPLEPH FOR A DOUBLY-DIMENSIONED TIME **         |
|          THE REASON FOR THIS OPTION IS DISCUSSED IN THE SUB. STATE      |
|      TARG = int NUMBER OF 'TARGET' POINT.                               |
|      CENT = int NUMBER OF CENTER POINT.                                 |
|             THE NUMBERING CONVENTION FOR 'TARG' AND 'CENT' IS:          |
|                 1 = MERCURY         8 = NEPTUNE                         |
|                 2 = VENUS           9 = PLUTO                           |
|                 3 = EARTH          10 = MOON                            |
|                 4 = MARS           11 = SUN                             |
|                 5 = JUPITER        12 = SOLAR-SYSTEM BARYCENTER         |
|                 6 = SATURN         13 = EARTH-MOON BARYCENTER           |
|                 7 = URANUS         14 = NUTATIONS (intITUDE AND OBLIQ)  |
|                15 = LIBRATIONS, IF ON EPH FILE                          |
|               (IF NUTATIONS ARE WANTED, SET TARG = 14. FOR LIBRATIONS,  |
|               SET TARG = 15. 'CENT' WILL BE IGNORED ON EITHER CALL.)    |
+-------------------------------------------------------------------------+
| OUTPUTS:                                                                |
|       RRD = OUTPUT 6-WORD D.P. ARRAY CONTAINING POSITION AND VELOCITY   |
|             OF POINT 'TARG' RELATIVE TO 'CENT'. THE UNITS ARE AU AND    |
|             AU/DAY. FOR LIBRATIONS THE UNITS ARE RADIANS AND RADIANS    |
|             PER DAY. IN THE CASE OF NUTATIONS THE FIRST FOUR WORDS OF   |
|             RRD WILL BE SET TO NUTATIONS AND RATES, HAVING UNITS OF     |
|             RADIANS AND RADIANS/DAY.                                    |
|             NOTE: IN MANY CASES THE USER WILL NEED ONLY POSITION        |
|                   VALUES FOR EPHEMERIDES OR NUTATIONS. FOR              |
|                   POSITION-ONLY OUTPUT, THE int VARIABLE 'IPV'          |
|                   IN THE COMMON AREA /PLECOM/ SHOULD BE SET = 1         |
|                   BEFORE THE NEXT CALL TO PLEPH. (ITS DEFAULT           |
|                   VALUE IS 2, WHICH RETURNS BOTH POSITIONS AND          |
|                   RATES.)                                               |
|                                                                         |
|      INSIDE is .TRUE. if the input Julian Ephemeris Date (JD) is within |
|             the ephemeris time span.  If not, INSIDE is set to .FALSE.  |
+-------------------------------------------------------------------------+
| REFERENCE: JPL (Fortran code)                                           |
+========================================================================*/
int PLEPH0x12(FILE *fp_Bin,int n,double jd,int targ,int cent,double *rrd,int *inside,double *jd2)
{
		/* INITIALIZED DATA: */
		double     fac =   0.;
		double embf[2] = { -1.,1. };
//		double  pv[78] = { 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
//							0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
		double  pv[13*12] = { 0. };
		static int  nemb = 1;
		static int FIRST = TRUE_;
		static int bsave;
		static double ve[2];

		/* DECLARATION: */
		int list[12] = { 0,0,0,0,0,0,0,0,0,0,0,0 };
		int llst[13] = { 1,2,10,4,5,6,7,8,9,10,11,11,3 };

		int l[2],tc[2],i,lme;
		double jed[2],jdtot;

/*		if(targ==14) cent= 0; *//* 1999/07/25 */
		
		/* 1ST TIME IN, BE SURE EPHEMERIS IS INITIALIZED */
		if(FIRST)
				{ FIRST=FALSE_;
						ve[0] = 1. / (HEADER.emrat + 1.);
						ve[1] = HEADER.emrat * ve[0];
				}

		if(n==1)
				{ /* ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT */
				  /* (SEE THE DISCUSSION IN THE SUBROUTINE STATE)              */
						jed[0] = jd2[0];
						jed[1] = jd2[1];
				}
		else
				{ /* INITIALIZE JED FOR 'STATE' AND SET UP COMPONENT COUNT */
						jed[0] = jd;
						jed[1] = 0.;
				}

		jdtot = jed[0] + jed[1];

		if(jdtot < HEADER.ss[0] || jdtot > HEADER.ss[1])
				{ *inside = FALSE_;
						return(0);
				}

		*inside = TRUE_;

		/* CHECK FOR NUTATION CALL */
		if(targ == 14)
				{ if(HEADER.ipt[34] > 0)
								{ list[10] = 4;
										STATEx12(fp_Bin,jed, list, pv, rrd);
										list[10] = 0;
										return(0);
								}
						else
								{ printf("\n *** NO NUTATION ON THE EPHEMERIS FILE ***");
										exit(1);
								}
				}

		/* CHECK FOR LIBRATIONS */
		if(targ == 15)
				{ if(HEADER.lpt[1] > 0)
								{ list[11] = 4;
										STATEx12(fp_Bin,jed, list, pv,rrd);
										list[11] = 0;
//										for(i = 1; i <= 6; ++i)
//														rrd[i-1] = pv[i + 59];
										for(i = 1; i <= 12; ++i)
										                rrd[i-1] = pv[i + 10*12-1];
										return(0);
								}
						else
								{ printf("\n *** NO LIBRATIONS ON THE EPHEMERIS FILE ***");
										exit(1);
								}
				}

		/* CHECK FOR TARGET POINT = CENTER POINT */
		if(targ == cent)
				{ 
//					for(i = 1; i <=  6; ++i)
					for(i = 1; i <= 12; ++i)
									rrd[i-1] = 0.;
					return(0);
				}

		/* FORCE BARYCENTRIC OUTPUT BY 'STATE' */
		bsave = stcomm.bary;
		stcomm.bary = TRUE_;

		/* SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL */
		tc[0] = targ;
		tc[1] = cent;
		lme = 0;

		for(i = 1; i <= 2; ++i)
					{ l[i - 1] = llst[tc[i - 1] - 1];
							if( l[i - 1] < 11)list[l[i - 1] - 1] = 4;
							if(tc[i - 1] == 3)
									{ lme = 3;
											fac = -ve[0];
									}
							else if( tc[i - 1] == 10)
														{ lme = 10; fac = ve[1];}
							else if( tc[i - 1] == 13)
														{ nemb = i; }
					}

		if((list[9] == 4) && (l[0] != l[1]))
		list[2] = 4 - list[2];

		/* MAKE CALL TO STATE */
		STATEx12(fp_Bin,jed, list, pv, rrd);



		/* CASE: EARTH-TO-MOON */
		if((targ == 10) && (cent == 3))
				{ 
//					for(i = 1; i <=  6; ++i) rrd[i-1] = pv[i + 53];
					for(i = 1; i <= 12; ++i) rrd[i-1] = pv[i + 12*12-1];
				}

		/* CASE: MOON-TO-EARTH */
		else if((targ == 3) && (cent == 10))
									{ 
//										for(i = 1; i <=  6; ++i) rrd[i-1] = -pv[i + 53];
										for(i = 1; i <= 12; ++i) rrd[i-1] = -pv[i + 12*12-1];
									}
		/* CASE: EMBARY TO MOON OR EARTH */
		else if((targ == 13 || cent == 13) && list[9] == 4)
									{ 
//										for(i = 1; i <=  6; ++i) rrd[i-1] = pv[i + 53     ] * fac * embf[nemb - 1];
										for(i = 1; i <= 12; ++i) rrd[i-1] = pv[i + 12*12-1] * fac * embf[nemb - 1];
									}
		/* OTHERWISE, GET EARTH OR MOON VECTOR AND THEN GET OUTPUT VECTOR */
		else
				{ 
//					for(i = 1; i <=  6; ++i)
					for(i = 1; i <= 12; ++i)
									{ 
//										    pv[i + 59] = stcomm.pvsun[i - 1];
//											pv[i + 71] = pv[i + 11];
										    pv[i + 12*10-1] = stcomm.pvsun[i - 1];
											pv[i + 12*12-1] = pv[i + 12*2-1];
											if(lme > 0)
//													   pv[i + lme *6 -  7] = pv[i + 11] + fac * pv[i + 53];
//											rrd[i-1] = pv[i + targ*6 -  7] - pv[i + cent * 6 - 7];
													   pv[i + lme *12 -13] = pv[i + 23] + fac * pv[i + 9*12-1];
											rrd[i-1] = pv[i + targ*12 -13] - pv[i + cent * 12 - 13];
									}
				}

		/* CLEAR 'STATE' BODY ARRAY AND RESTORE BARYCENTER FLAG */
		list[2] = 0;
		list[l[0] - 1] = 0;
		list[l[1] - 1] = 0;
		stcomm.bary = bsave;

		/*  THAT'S ALL */
		return(0);
} /* pleph0 */

/*========================================================================+
|                                PLEPH                                    |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+========================================================================*/
int PLEPHx12(FILE *fp_Bin,double jd,int targ,int cent,double *rrd,int *inside)
{
		return PLEPH0x12(fp_Bin,0,jd,targ,cent,rrd,inside,(double *)0);
}

/*========================================================================+
|                               DPLEPH                                    |
+=========================================================================+
| AUTOR: K. ARFA-KABOODVAND (Aero-Space Engineer)                         |
+========================================================================*/
int DPLEPHx12(FILE *fp_Bin,double *jd2,int targ,int cent,double *rrd,int *inside)
{
//		return PLEPH0(fp_Bin,1,(double)0,targ,cent,rrd,(int *)0,jd2);
		return PLEPH0x12(fp_Bin,1,(double)0,targ,cent,rrd,inside,jd2);
}



void Ephemerides(FILE *fp_Bin,int DE, double jd,int Planet)
{

	
	/* The detailed calculation method is described in the
	   "The Astromical Almanac"
	   Page B36 to B39 in case of 1999 version */
	
	
		/* VARIABLES FOR JPL COMPUTATION:*/

		/* VARIABLE DECLARATION: */
		int nctr,ntarg,inside;
		double r[6],et,et0;
		double Tau;
		double sd;

/**********************
			step. 1
***********************/
						et0=jd;
						nctr  = 12;
						/* Where does the ephemeris think things should be? */
						stcomm.km=FALSE_;/* AU & AU/Day*/
						stcomm.bary=TRUE_;/*Barycentre (Solar-sys. Barycentre)*/

/**********************
			step. 2
***********************/
						ntarg =  3;
						PLEPH(fp_Bin,et0,ntarg,nctr,r,&inside);
						Eb.max[0]=r[0],Eb.max[1]=r[1],Eb.max[2]=r[2];
						Ed.max[0]=r[3],Ed.max[1]=r[4],Ed.max[2]=r[5];
						Scalar ( Eb, &S_Eb);
						ntarg =  11;
						PLEPH(fp_Bin,et0,ntarg,nctr,r,&inside);
						Sb.max[0]=r[0],Sb.max[1]=r[1],Sb.max[2]=r[2];
						Scalar ( Sb, &S_Sb);
						if(Planet==11){
						Sub_Matrix_1(Sb,Eb,&P);
						Scalar ( P, &S_P);
						Eph.d=S_P;
						S_P0=S_P;
						Tau = S_P/173.1446;
						}
						else{
						ntarg =  Planet;
						PLEPH(fp_Bin,et0,ntarg,nctr,r,&inside);
						Qb.max[0]=r[0],Qb.max[1]=r[1],Qb.max[2]=r[2];
						Scalar ( Qb, &S_Qb);
						Sub_Matrix_1(Qb,Eb,&P);
						Scalar ( P, &S_P);
						Eph.d=S_P;
						S_P0=S_P;
						Sub_Matrix_1(Eb,Sb,&E);
						Scalar ( E, &S_E);
						Sub_Matrix_1(Qb,Sb,&Q);
						Scalar ( Q, &S_Q);
						Tau = (S_P+2.0*9.87e-9*log((S_E+S_P+S_Q)/(S_E-S_P+S_Q)))/173.1446;

						}

/******************/


						/* Where does the ephemeris think things should be? */
						stcomm.km=FALSE_;/* AU & AU/Day*/
						stcomm.bary=TRUE_;/*Barycentre (Solar-sys. Barycentre)*/

						if(Planet==11){

							while(1){	
							et    = et0-Tau;
							nctr  = 12;
						
							ntarg =  11;
							PLEPH(fp_Bin,et,ntarg,nctr,r,&inside);
							Sb.max[0]=r[0],Sb.max[1]=r[1],Sb.max[2]=r[2];
							Scalar ( Sb, &S_Sb);

							Sub_Matrix_1(Sb,Eb,&P);
							Scalar ( P, &S_P);
							Tau = S_P/173.1446;
							if(fabs(S_P0-S_P)<1e-9) break;
							S_P0 = S_P;
								}
							}

						else{ 
							while(1){	
							et    = et0-Tau;
							nctr  = 12;
							ntarg =  Planet;
							PLEPH(fp_Bin,et,ntarg,nctr,r,&inside);
							Qb.max[0]=r[0],Qb.max[1]=r[1],Qb.max[2]=r[2];
							Scalar ( Qb, &S_Qb);
						
							ntarg =  11;
							PLEPH(fp_Bin,et,ntarg,nctr,r,&inside);
							Sb.max[0]=r[0],Sb.max[1]=r[1],Sb.max[2]=r[2];
							Scalar ( Sb, &S_Sb);

							Sub_Matrix_1(Qb,Eb,&P);
							Scalar ( P, &S_P);
							Sub_Matrix_1(Qb,Sb,&Q);
							Scalar ( Q, &S_Q);
							Tau = (S_P+2.0*9.87e-9*log((S_E+S_P+S_Q)/(S_E-S_P+S_Q)))/173.1446;
							if(fabs(S_P0-S_P)<1e-9) break;
							S_P0 = S_P;
							}
						}

/********************/

			if(Planet==11){
			Multi_Const_Matrix_1(P,1/S_P,&p);
			Scalar ( p, &S_p);}
			else{

			Multi_Const_Matrix_1(P,1/S_P,&p);
			Scalar ( p, &S_p);
			Multi_Const_Matrix_1(Q,1/S_Q,&q);
			Scalar ( q, &S_q);
			Multi_Const_Matrix_1(E,1/S_E,&e);
			Scalar ( e, &S_e);
			}
/**********************
			step. 3
***********************/
			if(Planet==11){
				Eql_Matrix_1(p, &p1);
			}
			else{
			Scalar_Dot ( p,q, &S_pq);
			Scalar_Dot ( e,p, &S_ep);
			Scalar_Dot ( q,e, &S_qe);
			Multi_Const_Matrix_1(e,S_pq,&p01);
			Multi_Const_Matrix_1(q,S_ep,&p02);
			Sub_Matrix_1(p01,p02,&p00);
			Multi_Const_Matrix_1(p00,2.0*9.87e-9/S_E/(1+S_qe),&p00);
			Sum_Matrix_1(p,p00,&p1);
			}

/**********************
			step. 4
***********************/

			Multi_Const_Matrix_1(Ed,0.0057755,&V);
			Scalar ( V, &S_V);
			Beta=1.0/sqrt(1-S_V*S_V);
			Scalar_Dot ( p1,V, &S_p1V);
			Multi_Const_Matrix_1(p1,1.0/(Beta*(1.0+S_p1V)),&p20);
			Multi_Const_Matrix_1(V,(1.0+S_p1V/(1.0+1.0/Beta))/(1.0+S_p1V),&p21);
			Sum_Matrix_1(p20,p21,&p2);

/**********************
			step. 5
***********************/

			Precession_Data ( jd,&Pr);
			Nutation ( fp_Bin ,DE, jd, &Nu);
			Multi_Matrix_3(Pr,Nu,&R);

			Multi_Matrix_1(R,p2,&p3);
			Scalar ( p3, &S_p3);

/**********************
			step. 6
***********************/
			Eph.ra=atan2(p3.max[1],p3.max[0])/pai;
			if(Eph.ra<0)Eph.ra=Eph.ra+360.0;
			Eph.dec=asin(p3.max[2])/pai;

/**********************
			step. 7
***********************/

			Matrix_R1(Nu3, Obl_2000(jd)*pai, &Nu3);
			Multi_Matrix_1(Nu3,p3,&p3);
 			Eph.lam=atan2(p3.max[1],p3.max[0])/pai;
			if(Eph.lam<0)Eph.lam=Eph.lam+360.0;
			Eph.bet=asin(p3.max[2])/pai;

			switch(Planet){
			case 1	: sd=3.36;break;
			case 2	: sd=8.34;break;
			case 3	: sd=0.0;break;
			case 4	: sd=4.64;break;
			case 5	: sd=98.44;break;
			case 6	: sd=82.73;break;
			case 7	: sd=35.02;break;
			case 8	: sd=33.05;break;
			case 9	: sd=2.07;break;
			case 10	: sd=6378.14;break;
			case 11	: sd=959.63;break;
			default	: sd=0.0;break;
				}

			if(Planet==10)
				{
				Eph.ss=asin(sd*0.272481/(Eph.d*1.4959787066e8))/pai;
				}
				
				
				else
				{Eph.ss=sd/Eph.d;}

/*			printf  ("\nSemi_Diameter=%15.9f\n",Eph.ss);*/

}

void Conversion       ( const double *state1,  double *state3, const double MUSUN, double jd)
{

   double x, y, z, Vx, Vy, Vz, ratio;

   x  = state1[0];
   y  = state1[1];
   z  = state1[2];
   Vx = state1[3];
   Vy = state1[4];
   Vz = state1[5];

  
    ratio  =     cos(pai*Obl_2000(jd)) * y  +  sin(pai*Obl_2000(jd)) * z;
    state3[2] = -sin(pai*Obl_2000(jd)) * y  +  cos(pai*Obl_2000(jd)) * z;
	state3[1] = ratio;
	state3[0] = x;	

// radius distance 
	ratio = Vx*Vx + Vy*Vy + Vz*Vz;
// Convert from equatorial to ecliptic coordinates 
	ratio  =      cos(pai*Obl_2000(jd)) * Vy  +  sin(pai*Obl_2000(jd)) * Vz;
	state3[5]  = -sin(pai*Obl_2000(jd)) * Vy  +  cos(pai*Obl_2000(jd)) * Vz;
	state3[4] = ratio;
	state3[3] = Vx;

	return;
}
void Conversion_J2000 ( const double *state1,  double *state3, int n)
{

   double Obl_J2000, cosO, sinO;

//	Obl_J2000=pai*23.4392911;
	Obl_J2000=pai*23.439279444444444;
	cosO = cos(Obl_J2000);
	sinO = sin(Obl_J2000);

//  Convert from equatorial to ecliptic coordinates   

	state3[0] =  state1[0];	
	state3[1] =  cosO * state1[1]  +  sinO * state1[2];
    state3[2] = -sinO * state1[1]  +  cosO * state1[2];

	if (n==1) return;

 	state3[3] =  state1[3]; 
	state3[4] =  cosO * state1[4]  +  sinO * state1[5];
	state3[5] = -sinO * state1[4]  +  cosO * state1[5];

	if (n==2) return;

	state3[6] =  state1[6];	
	state3[7] =  cosO * state1[7]  +  sinO * state1[8];
    state3[8] = -sinO * state1[7]  +  cosO * state1[8];

	if (n==3) return;

	state3[ 9] =  state1[9];	
	state3[10] =  cosO * state1[10]  +  sinO * state1[11];
    state3[11] = -sinO * state1[10]  +  cosO * state1[11];

	return;
}





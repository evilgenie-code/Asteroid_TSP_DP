//-----------------------------------------------------------------------------------------------------------------------------------------------------
//  —истемы уравнений  и  Ёфемериды
//-----------------------------------------------------------------------------------------------------------------------------------------------------
#include <cmath>
#include <complex>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h> 
#include "Constants.h"
#include "Ephemeris.h"

extern "C" {
#include "astro.h" 
} 

extern "C" {
	__declspec(dllimport) void INITDE405(void);
	__declspec(dllimport) void FREEDE405(void);
	__declspec(dllimport) void PLNPOSVEL(double *, int *, double *);
} 


double Kep_elements[100][6], jd0_elements[100];

extern FILE *fp_Bin;
extern double UnitM, UnitP, Unitb;

void  Unit_MPb(double UnitM, double &UnitP, double &Unitb){
 
	    UnitP   = UnitAcc*UnitM;               // newton
	    Unitb   = UnitP/(1000*UnitV);          // kg/sec	  
  
  }

namespace Ephemeris{

	void Calendar(double juday, double* year, double* month, double* day, double* hour, double* min, double* sec, int LST, double longitude){

		double flag, Z, Alfa, A, B, C, D, E, F, Day0, hour0, min0;

		flag = 1;

		if (longitude < 0) {
			longitude = -longitude;
			flag = -1.0;
		}

		if (LST == 1) juday = juday + flag * (double)(int)(longitude / 15.0 + 0.5) / 24.0;
		else    juday = juday + flag * longitude / 15.0 / 24.0;

		juday = juday + 0.5;
		Z = (double)(int)(juday);
		F = juday - Z;

		if (Z < 2299161.0) A = Z;
		else {
			Alfa = (double)(int)((Z - 1867216.25) / 36524.25);
			A = Z + 1 + Alfa - (double)(int)(Alfa / 4.0);
		}

		B = A + 1524;
		C = (double)(int)((B - 122.1) / 365.25);
		D = (double)(int)(365.25 * C);
		E = (double)(int)((B - D) / 30.6001);

		Day0 = B - D - (double)(int)(30.60001 * E) + F;
		*day = (double)(int)(Day0);
		hour0 = F * 24.0;
		*hour = floor(hour0);
		min0 = (hour0 - *hour) * 60;
		*min = floor(min0);
		*sec = (min0 - *min) * 60;
		if (E < 14.0) *month = E - 1.0;
		if (E == 14.0 || E == 15.0)	 *month = E - 13.0;
		if (*month > 2.0) *year = C - 4716;
		if (*month == 1.0 || *month == 2.0) *year = C - 4715;

	}

	double JD				(double year,double month,double day,double hour,int LST,double longitude)
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
	double JD_epf           (double year,double month,double day,double hour){

		int LST = 0;
		double longitude = 0;
		double jd = JD(year,month,day,hour,LST,longitude);
	
    	return jd;
	}
	FILE* INITIAL_DE        (int DE){
		FILE* fp_Bin;
		fp_Bin=INITIAL(DE);
		return fp_Bin;
	}
	void Ephemeris          (FILE *fp_Bin,double  jd,int targ,int cent,double *rrd,int *inside){

		double rrd0[6];
		PLEPH(fp_Bin,jd,targ,cent,rrd0,inside);	

	    rrd0[0]=rrd0[0]*HEADER.au;          //  km
	    rrd0[1]=rrd0[1]*HEADER.au;          //  km
	    rrd0[2]=rrd0[2]*HEADER.au;          //  km
	    rrd0[3]=rrd0[3]*(HEADER.au/86400);  //  km/sec
	    rrd0[4]=rrd0[4]*(HEADER.au/86400);  //  km/sec
	    rrd0[5]=rrd0[5]*(HEADER.au/86400);  //  km/sec

		Conversion_J2000 ( rrd0,  rrd, 2);
	}

	void Cartesian_coordinate     (double* Kep, double JD0, double JD, double* Dec){

		// Kep[0] -  a, больша€ полуось            AU
		// Kep[1] -  e, эксцентриситет             -
		// Kep[2] -  i, наклонение                rad
		// Kep[3] -  w, аргумент перигели€        rad
		// Kep[4] - OM, долгота восход€щего узла  rad
		// Kep[5] - M0, средн€€ аномали€          rad
		// JD0    -     эпоха                     day
		// JD     -     момент времени            day
		// Dec[]  - x,y,x,Vx,Vy,Vz

		double dM, dE, E, M, tet, gam, R, V, tetw, tetwgam;

		M = Kep[5] + sqrt(1/(Kep[0]*Kep[0]*Kep[0]))*(JD-JD0)/UnitT;
		M = fmod(M,PI2);
		E = M;
		dE = 1;
		while (abs(dE) > 1.e-14){
			dE = M + Kep[1]*sin(E) - E;
			E = E + dE;
		}

		tet = 2*atan2(sqrt(1+Kep[1])*tan(E*0.5),sqrt(1-Kep[1]));
		gam = atan2(Kep[1]*sin(tet),(1+Kep[1]*cos(tet)));

		R = Kep[0]*(1-Kep[1]*cos(E));
		V = sqrt(2/R-1/Kep[0]);

		tetw = tet+Kep[3];
		tetwgam = tet+Kep[3] - gam;

		Dec[0] = R*(cos(tetw)*cos(Kep[4]) - sin(tetw)*cos(Kep[2])*sin(Kep[4]));
		Dec[1] = R*(cos(tetw)*sin(Kep[4]) + sin(tetw)*cos(Kep[2])*cos(Kep[4]));
		Dec[2] = R*(sin(tetw)*sin(Kep[2]));

		Dec[3] = V*(-sin(tetwgam)*cos(Kep[4]) - cos(tetwgam)*cos(Kep[2])*sin(Kep[4]));
		Dec[4] = V*(-sin(tetwgam)*sin(Kep[4]) + cos(tetwgam)*cos(Kep[2])*cos(Kep[4]));
		Dec[5] = V*(cos(tetwgam)*sin(Kep[2]));
	}
    void Cartesian_coordinate     (double* Kep, double JD0, double JD, double* Dec, double *workd){

	    // workd[9]

		// Kep[0] -  a, больша€ полуось            AU
		// Kep[1] -  e, эксцентриситет             -
		// Kep[2] -  i, наклонение                rad
		// Kep[3] -  w, аргумент перигели€        rad
		// Kep[4] - om, долгота восход€щего узла  rad
		// Kep[5] - M0, средн€€ аномали€          rad
		// JD0    -     эпоха                     day
		// JD     -     момент времени            day
		// Dec[]  - x,y,x,Vx,Vy,Vz

		// workd[0] =      dE, - 
	    // workd[1] =       E, - эксцентрическа€ анамали€
	    // workd[2] =       M, - средн€€ аномали€
	    // workd[3] =     tet, - истинна€ аномали€
	    // workd[4] =     gam, 
	    // workd[5] =       R, 
	    // workd[6] =       V, 
	    // workd[7] =    tetw, 
	    // workd[8] = tetwgam

		workd[2] = Kep[5] + sqrt(1/(Kep[0]*Kep[0]*Kep[0]))*(JD-JD0)/UnitT;
		workd[2] = fmod(workd[2],PI2);
		workd[1] = workd[2];
		workd[0] = 1;
		while (abs(workd[0]) > 1.e-14){
			workd[0] = workd[2] + Kep[1]*sin(workd[1]) - workd[1];
			workd[1] = workd[1] + workd[0];
		}

		workd[3] = 2*atan2(sqrt(1+Kep[1])*tan(workd[1]*0.5),sqrt(1-Kep[1]));
		workd[4] = atan2(Kep[1]*sin(workd[3]),(1+Kep[1]*cos(workd[3])));

		workd[5] = Kep[0]*(1-Kep[1]*cos(workd[1]));
		workd[6] = sqrt(2/workd[5]-1/Kep[0]);

		workd[7] = workd[3]+Kep[3];
		workd[8] = workd[7] - workd[4];

		Dec[0] = workd[5]*(cos(workd[7])*cos(Kep[4]) - sin(workd[7])*cos(Kep[2])*sin(Kep[4]));
		Dec[1] = workd[5]*(cos(workd[7])*sin(Kep[4]) + sin(workd[7])*cos(Kep[2])*cos(Kep[4]));
		Dec[2] = workd[5]*(sin(workd[7])*sin(Kep[2]));

		Dec[3] = workd[6]*(-sin(workd[8])*cos(Kep[4]) - cos(workd[8])*cos(Kep[2])*sin(Kep[4]));
		Dec[4] = workd[6]*(-sin(workd[8])*sin(Kep[4]) + cos(workd[8])*cos(Kep[2])*cos(Kep[4]));
		Dec[5] = workd[6]*(cos(workd[8])*sin(Kep[2]));
	}
    void Cartesian_coordinate     (double* Kep, double   R,            double* Dec, double *workd){

	    // workd[7]

		// Kep[0] -  hp, высота перицентра          -
		// Kep[1] -  ha, высота апоцентра           -
		// Kep[2] -   i, наклонение                rad
		// Kep[3] -   w, аргумент перигели€        rad
		// Kep[4] -  om, долгота восход€щего узла  rad
		// Kep[5] - tet, истинна€ аномали€         rad
		// Dec[]  - x,y,x,Vx,Vy,Vz

		// workd[0] =       a, больша€ полуось      -
	    // workd[1] =       e, эксцентриситет       -
	    // workd[2] =       R,                      -
	    // workd[3] =       V, 
	    // workd[4] =     gam, 
	    // workd[5] =    tetw, 
	    // workd[6] = tetwgam

		workd[0] = (Kep[0] + Kep[1] + 2.*R);
		workd[1] = (Kep[1] - Kep[0])/workd[0];
		workd[0] = 0.5*workd[0];

		workd[2] = workd[0]*(1.-workd[1]*workd[1])/(1.+workd[1]*cos(Kep[5]));
		workd[3] = sqrt(2./workd[2]-1./workd[0]);

		workd[4] = atan2(workd[1]*sin(Kep[5]),(1.+workd[1]*cos(Kep[5])));

		workd[5] = Kep[5]+Kep[3];
		workd[6] = workd[5] - workd[4];

		Dec[0] = workd[2]*(cos(workd[5])*cos(Kep[4]) - sin(workd[5])*cos(Kep[2])*sin(Kep[4]));
		Dec[1] = workd[2]*(cos(workd[5])*sin(Kep[4]) + sin(workd[5])*cos(Kep[2])*cos(Kep[4]));
		Dec[2] = workd[2]*(sin(workd[5])*sin(Kep[2]));

		Dec[3] = workd[3]*(-sin(workd[6])*cos(Kep[4]) - cos(workd[6])*cos(Kep[2])*sin(Kep[4]));
		Dec[4] = workd[3]*(-sin(workd[6])*sin(Kep[4]) + cos(workd[6])*cos(Kep[2])*cos(Kep[4]));
		Dec[5] = workd[3]*( cos(workd[6])*sin(Kep[2]));
	}

	void Ephemerisx6        (FILE *fp_Bin,double  jd,int targ,int cent,double *rrd,int *inside, double* Kep, double jd0, int dimension){
/*		+-------------------------------------------------------------------------+
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
		|                100 < asteroid                                           |
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
		+-------------------------------------------------------------------------+*/

		double rrd0[9], JD[2]; JD[0] = jd; JD[1] = 0;

		if (targ < 100) {

			PLEPH(fp_Bin,jd,targ,cent,rrd0,inside);	

			rrd0[0] = rrd0[0] * HEADER.au;                                   //  km
			rrd0[1] = rrd0[1] * HEADER.au;                                   //  km
			rrd0[2] = rrd0[2] * HEADER.au;                                   //  km
			rrd0[3] = rrd0[3] * (HEADER.au / 86400);                         //  km/sec
			rrd0[4] = rrd0[4] * (HEADER.au / 86400);                         //  km/sec
			rrd0[5] = rrd0[5] * (HEADER.au / 86400);                         //  km/sec

			Conversion_J2000 ( rrd0,  rrd, 2);

			if (!dimension){
				rrd[0] = rrd[0] / UnitR;
				rrd[1] = rrd[1] / UnitR;
				rrd[2] = rrd[2] / UnitR;

				rrd[3] = rrd[3] / UnitV;
				rrd[4] = rrd[4] / UnitV;
				rrd[5] = rrd[5] / UnitV;

			}
		}

		else if (targ >= 100) {

			Cartesian_coordinate (Kep, jd0, jd, rrd);

			if (dimension){
				rrd[0] = rrd[0] * UnitR;                                       //  km
				rrd[1] = rrd[1] * UnitR;                                       //  km
				rrd[2] = rrd[2] * UnitR;                                       //  km

				rrd[3] = rrd[3] * UnitV;                                       //  km/sec
				rrd[4] = rrd[4] * UnitV;                                       //  km/sec
				rrd[5] = rrd[5] * UnitV;                                       //  km/sec

			}
		}

		else { rrd[0] = 0; rrd[1] = 0; rrd[2] = 0; rrd[3] = 0; rrd[4] = 0; rrd[5] = 0; *inside=-1; }
       
	}
	
	void NAME                (char *Name, int &n, int i){

		int k=0;

		while (_strnicmp(&Name[k], " ",1)==0 & k<15) k++;
	
		int size = strlen(&Name[k]);

		if (_strnicmp(&Name[k], "1",size)==0 || _strnicmp(&Name[k],"MERCURY"   ,size)==0)   { n =  1; return; }
		if (_strnicmp(&Name[k], "2",size)==0 || _strnicmp(&Name[k],"VENUS"     ,size)==0)   { n =  2; return; }
		if (_strnicmp(&Name[k], "3",size)==0 || _strnicmp(&Name[k],"EARTH"     ,size)==0)   { n =  3; return; }
		if (_strnicmp(&Name[k], "4",size)==0 || _strnicmp(&Name[k],"MARS"      ,size)==0)   { n =  4; return; }
		if (_strnicmp(&Name[k], "5",size)==0 || _strnicmp(&Name[k],"JUPITER"   ,size)==0)   { n =  5; return; }
		if (_strnicmp(&Name[k], "6",size)==0 || _strnicmp(&Name[k],"SATURN"    ,size)==0)   { n =  6; return; }
		if (_strnicmp(&Name[k], "7",size)==0 || _strnicmp(&Name[k],"URANUS"    ,size)==0)   { n =  7; return; }
		if (_strnicmp(&Name[k], "8",size)==0 || _strnicmp(&Name[k],"NEPTUNE"   ,size)==0)   { n =  8; return; }
		if (_strnicmp(&Name[k], "9",size)==0 || _strnicmp(&Name[k],"PLUTO"     ,size)==0)   { n =  9; return; }
		if (_strnicmp(&Name[k],"13",size)==0 || _strnicmp(&Name[k],"EARTH-MOON",size)==0)   { n = 13; return; }

		n = 100+i; 	
		Ephemeris::Kep_elements_MPCORB (&Name[k], &Kep_elements[i][0], jd0_elements[i]);
		Kep_elements[i][2] = Kep_elements[i][2] / raddeg;
		Kep_elements[i][3] = Kep_elements[i][3] / raddeg;
		Kep_elements[i][4] = Kep_elements[i][4] / raddeg;
		Kep_elements[i][5] = Kep_elements[i][5] / raddeg;

	}
	void Epoch_MPCORB_JD     (char * Epoch, double &jd0){

		double year, month, day, hour = 0;
		char Year[2];
		double year34;
		char ALF[] = {'1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V'};
		int i;

		Year[0] = Epoch[1];
		Year[1] = Epoch[2];
		sscanf_s(Year, "%lf", &year34);

		switch (Epoch[0]) {
			case 'I':  year = 1800+year34; break; 
			case 'J':  year = 1900+year34; break; 
			case 'K':  year = 2000+year34; break; 
		}

		for (i = 0; i < 12; i++) if (Epoch[3] == ALF[i]) { month = i + 1; break; }
		for (i = 0; i < 31; i++) if (Epoch[4] == ALF[i]) {   day = i + 1; break; }

		jd0 = JD_epf (year,month,day,hour);
	
	}
	void Kep_elements_MPCORB (char * Name, double* Kep, double &jd0){
	
		// Kep[0] -  a, больша€ полуось            AU
		// Kep[1] -  e, эксцентриситет             -
		// Kep[2] -  i, наклонение                deg
		// Kep[3] -  w, аргумент перигели€        deg
		// Kep[4] - OM, долгота восход€щего узла  deg
		// Kep[5] - M0, средн€€ аномали€          deg
		// jd0    -     эпоха                     day

	    ifstream fin("mpcorb.dat");	
		char line[1000], Epoch[5], word[50], * tok;
		int size = 0, j; 
		char *next_token = NULL;

		size = strlen(Name);

		while( fin.getline(line, 1000)) {

			int da = _strnicmp(&line[175],Name,size);

			if (da==0) { 

				j=0;
				for(char * tok = strtok_s( line," ", &next_token); tok; tok = strtok_s( NULL," ", &next_token))	{

					switch (j) {
						case  3: Epoch_MPCORB_JD (tok, jd0);	break; 					
						case  4: Kep[5] = atof( tok ) / raddeg; break;
						case  5: Kep[3] = atof( tok ) / raddeg; break;
						case  6: Kep[4] = atof( tok ) / raddeg; break;
						case  7: Kep[2] = atof( tok ) / raddeg; break;
						case  8: Kep[1] = atof( tok );			break; 
						case 10: Kep[0] = atof( tok );			break; 
					}
						j++;
						if (j>10) break; 
				}
				break; 
			}
		}
		fin.close();
	}

	void Kep_elements_MPCORB_data (char Name[][19], double Kep[][6], double* jd0, double a[2], double e[2], double i[2], int nn[2], int n) {

		// Kep[0] -  a, больша€ полуось            AU
		// Kep[1] -  e, эксцентриситет             -
		// Kep[2] -  i, наклонение                deg
		// Kep[3] -  w, аргумент перигели€        deg
		// Kep[4] - OM, долгота восход€щего узла  deg
		// Kep[5] - M0, средн€€ аномали€          deg
		// jd0    -     эпоха                     day

		ifstream fin("mpcorb.dat");
		char line[1000], Epoch[5], word[50], *tok;
		int size = 0, j, k;
		char *next_token = NULL;
		bool aa, ee, ii, in;

		size = 18;

		for (k = 0; k < 41 + nn[0]; k++) 
			fin.getline(line, 1000);

		for (in = nn[0], k = 0; in < nn[1] && fin.getline(line, 1000) && k < n; in++) {

			strncpy(Name[k], &line[175], size); Name[k][size] = '\0';

			j = 0;
			for (char * tok = strtok_s(line, " ", &next_token); tok; tok = strtok_s(NULL, " ", &next_token)) {

				switch (j) {
					case  3: Epoch_MPCORB_JD(tok, jd0[k]);		break;
					case  4: Kep[k][5] = atof(tok) / raddeg;	break;
					case  5: Kep[k][3] = atof(tok) / raddeg;	break;
					case  6: Kep[k][4] = atof(tok) / raddeg;	break;
					case  7: Kep[k][2] = atof(tok) / raddeg;	break;
					case  8: Kep[k][1] = atof(tok);				break;
					case 10: Kep[k][0] = atof(tok);				break;
					}
					j++;
					if (j>10) break;
				}

			aa = a[0] < Kep[k][0] && Kep[k][0] < a[1];
			ee = e[0] < Kep[k][1] && Kep[k][1] < e[1];
			ii = i[0] < Kep[k][2] && Kep[k][2] < i[1];

			if ( aa && ee && ii ) k++;

		}
		fin.close();
	}
	void Kep_elements_MPCORB_dataName(char Name[][19], double Kep[][6], double* jd0, int n) {

		// Kep[0] -  a, больша€ полуось            AU
		// Kep[1] -  e, эксцентриситет             -
		// Kep[2] -  i, наклонение                deg
		// Kep[3] -  w, аргумент перигели€        deg
		// Kep[4] - OM, долгота восход€щего узла  deg
		// Kep[5] - M0, средн€€ аномали€          deg
		// jd0    -     эпоха                     day

		ifstream fin("mpcorb.dat");
		char line[1000], Epoch[5], word[50], * tok;
		int size = 0, i, j, k;
		char* next_token = NULL;

		for (i = 0; i < n; i++) {

			size = strlen(Name[i]);

			fin.seekg(0, ios::beg);
			for (k = 0; k < 41; k++)
				fin.getline(line, 1000);

			while (fin.getline(line, 1000)) {

				int da = _strnicmp(&line[175], Name[i], size);
				
				j = 0;
				if (da == 0) {

					for (char* tok = strtok_s(line, " ", &next_token); tok; tok = strtok_s(NULL, " ", &next_token)) {

						switch (j) {
						case  3: Epoch_MPCORB_JD(tok, jd0[i]);		break;
						case  4: Kep[i][5] = atof(tok) / raddeg;	break;
						case  5: Kep[i][3] = atof(tok) / raddeg;	break;
						case  6: Kep[i][4] = atof(tok) / raddeg;	break;
						case  7: Kep[i][2] = atof(tok) / raddeg;	break;
						case  8: Kep[i][1] = atof(tok);				break;
						case 10: Kep[i][0] = atof(tok);				break;
						}
						j++;
						if (j > 10)
							break;
					}
				}
				if (j > 10)
					break;
			}

		}
		fin.close();
	}

	void Keplerian_elements(double* Dec, double  t, double mu, double* Kep, int par) {

		// Dec[] - x,y,x,Vx,Vy,Vz UnitR, UnitV
		// t     - момент времени UnitT

		// Kep[ 0] -  a, больша€ полуось				UnitR
		// Kep[ 1] -  p, фокальный параметр				UnitR
		// Kep[ 2] -  e, эксцентриситет					-
		// Kep[ 3] - rp, радиус перицентра				UnitR
		// Kep[ 4] - ra, радиус апоцентра				UnitR

		// Kep[ 5] - Incl, наклонение					rad
		// Kep[ 6] - Node, долгота восход€щего узла		rad
		// Kep[ 7] - Peri, аргумент перигели€			rad

		// Kep[ 8] - Nu, истинна€ аномали€				rad		
		// Kep[ 9] -  E, эксцентрическа€ аномали€		rad		
		// Kep[10] -  M, средн€€ аномали€				rad		

		// Kep[11] - tp, врем€ прохождени€ перицентра	UnitT		

		static double r, V, h, c[3], Cxy, C, p, f[3], F, e, a, ra, rp, Incl, Node, Peri, l[3], cosPeri, cosNu, frc, Nu, E, n, M, tp, tanNu2;

		r = sqrt(Dec[0] * Dec[0] + Dec[1] * Dec[1] + Dec[2] * Dec[2]);
		V = sqrt(Dec[3] * Dec[3] + Dec[4] * Dec[4] + Dec[5] * Dec[5]);

		h = V * V - 2. * mu / r;
		a = -mu / h;

		c[0] = Dec[1] * Dec[5] - Dec[2] * Dec[4];
		c[1] = Dec[2] * Dec[3] - Dec[0] * Dec[5];
		c[2] = Dec[0] * Dec[4] - Dec[1] * Dec[3];

		Cxy = sqrt(c[0] * c[0] + c[1] * c[1]);
		C = sqrt(Cxy * Cxy + c[2] * c[2]);

		p = C * C / mu;

		f[0] = Dec[4] * c[2] - Dec[5] * c[1] - mu / r * Dec[0];
		f[1] = Dec[5] * c[0] - Dec[3] * c[2] - mu / r * Dec[1];
		f[2] = Dec[3] * c[1] - Dec[4] * c[0] - mu / r * Dec[2];

		F = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);

		e = F / mu;
		double e2 = sqrt(1. + C * C / (mu * mu) * h);

		rp = p / (1. + e);
		ra = p / (1. - e);

		Incl = acos(c[2] / C);

		if (abs(c[0]) < 1.e-10 && abs(c[1]) < 1.e-10)	Node = 0;
		else if (c[0] > 0)								Node = acos(-c[1] / Cxy);
		else if (c[0] < 0)								Node = PI2 - acos(-c[1] / Cxy);
		else if (c[0] == 0 && c[1] < 0)						Node = 0;
		else if (c[0] == 0 && c[1] > 0)						Node = PI;

		l[0] = cos(Node);
		l[1] = sin(Node);
		l[2] = 0.;

		cosPeri = (l[0] * f[0] + l[1] * f[1] + l[2] * f[2]) / F;

		if (f[2] / F > 1.e-8) Peri = acos(cosPeri);
		else if (f[2] / F < 1.e-8) Peri = PI2 - acos(cosPeri);
		else if (-1.e-8 < f[2] / F < 1.e-8 && l[0] * f[0]>0) Peri = 0.;
		else if (-1.e-8 < f[2] / F < 1.e-8 && l[0] * f[0] < 0) Peri = PI;

		cosNu = (Dec[0] * f[0] + Dec[1] * f[1] + Dec[2] * f[2]) / (r * F);
		frc = (f[1] * Dec[0] - f[2] * Dec[1]) * c[0] + (f[2] * Dec[0] - f[0] * Dec[2]) * c[1] + (f[0] * Dec[1] - f[1] * Dec[0]) * c[2];

		if (frc > 0) Nu = acos(cosNu);
		else       Nu = PI2 - acos(cosNu);

		tanNu2 = tan(Nu / 2.);

		if (h < 0) {

			E = 2. * atan(sqrt((1. - e) / (1. + e)) * tanNu2);
			n = sqrt(mu / (a * a * a));
			M = E - e * sin(E);

		}
		else if (h == 0) {

			E = tanNu2;
			n = 2. * sqrt(mu / p * p * p);
			M = E + 1. / 3. * E * E * E;

		}
		else if (h > 0) {

			E = 2. * atanh(complex <double>(sqrt((e - 1.) / (e + 1.)) * tanNu2, 0.)).real();
			n = sqrt(mu / (-a * a * a));
			M = e * sinh(E) - E;

		}
		tp = t - M / n;

		if (par) {

			Kep[0] = a; Kep[1] = p; Kep[2] = e;
			Kep[3] = rp; Kep[4] = ra;

			Kep[5] = Incl; Kep[6] = Node; Kep[7] = Peri;

			Kep[8] = Nu;   Kep[9] = E;    Kep[10] = M;

			Kep[11] = tp;

		}
		else {

			Kep[0] = a;
			Kep[1] = e;
			Kep[2] = Incl;
			Kep[3] = Peri;
			Kep[4] = Node;
			Kep[5] = M;

		}
	}

}



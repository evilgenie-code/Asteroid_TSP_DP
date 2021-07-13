/***************************************************************************/
/*                                                                         */
/*   astro.h      :  「古天文学」記載の方法による天文計算ヘッダファイル    */
/*      Ver. 1.2    1999/07/02          By   takesako@mrj.biglobe.ne.jp    */
/*                                           (c)1999 Shinobu Takesako      */
/***************************************************************************/
/*   1999/06/29 彗星位置計算追加に伴う追加                                 */
/*   1999/07/01 恒星位置計算追加に伴う追加                                 */
/***************************************************************************/

#ifndef ASTRO_H
#define ASTRO_H

//#define pai (3.1415926535897932/180.0)
//#define pi2 (3.1415926535897932*2.0)
//#define HToR (3.1415926535897932 / 12)
//#define RToH (12 / 3.1415926535897932)

/*#define DToR (pai / 180)
#define HToR (pai / 12)
#define RToD (180 / pai)
#define RToH (12 / pai)
#define SToR (DToR / 3600)
#define EarthRadius 6378.14
#define SolarParallax (8.794 * SToR)
*/
#define WAIT getch()
#define TRUE_ (1)
#define FALSE_ (0)

typedef struct tagPosition {
     double lam,  /* 黄経   */
            bet,  /* 黄緯   */
            ss,   /* 視半径 */
            x,    /* 地心直交座標成分 */
            y,    /* 地心直交座標成分 */
            z,    /* 地心直交座標成分 */
            d,    /* 地球との距離 */
            r,    /* 地球との距離 */
            br,   /* 光度  */
            ra,   /* 赤経  */
            dec,  /* 赤緯  */
			elong,/* 彗星の太陽からの離角 */
			phase;/* 太陽－彗星－地球の位相角 */

} Position;

/*
Position Sun;   
Position Planet;
Position Moon;
Position Comet;
Position Star; */
Position Eph;

typedef struct tagParameter {
  char *CometName;
  double	T  , /* 近日点通過日時*/
			Epo, /* 軌道要素epoch */
			Mu,  /* 日日平均運動  */
			Qq,  /* 近日点距離    */
			Omg, /* 昇交点黄経    */
            Somg,/* 近日点引数    */
            Inc, /* 軌道傾斜角    */
			Ax,  /* 軌道半長径    */
			Ec,  /* 軌道離心率    */
			Equ, /* 軌道要素基準年*/
			G0,  /* 彗星絶対光度  */
			Ne;  /* 彗星光度係数  */
			
} Parameter;

Parameter Comet_Para;


typedef struct tagTDATE {
  int yy, mm;
  double dd;
  int h, m;
  double s;
  double H;
  } TDATE;


/* Common Block Declarations */
struct { char/* **cnam*/ *cnam[400];
									double *cval,ss[3],au,emrat;
									int denum, ncon, ipt[36],lpt[3],ksize;
							} HEADER;

struct { double *buf;}epib;
char  *CNAM[401][6];
char  *TTL[3][84];
long IPT[12][3],*ipt;
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
struct { int km, bary;
									double pvsun[12];
							} stcomm;

/******************************************************************************/
/*　　　行列式の計算　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　 */
/******************************************************************************/

typedef struct tagMatrix_3{
			double max[3][3];
				} Matrix_3;

Matrix_3 abc,def,ghi,R,Pr,Pr1,Pr2,Pr3,Pr4,Nu,Nu1,Nu2,Nu3,Nu4;

typedef struct tagMatrix_1{
			double max[3];
				} Matrix_1;

Matrix_1 Eb,Ed,E,Sb,Qb,P,Q,p,q,e,pq,ep,qe,p00,p01,p02,p1,p2,p3,V,p20,p21;
double   S_Eb,S_Ed,S_E,S_Sb,S_Qb,S_Q,S_P,S_P0,S_p,S_q,S_e,S_pq,S_ep,S_qe,S_V,S_p1V,S_p3,Beta,Ra,Dec;

/* common.c */
double JD(double year,double month,double day,double hour,int LST,double longitude);
void   Calendar( double juday, double *year,double *month,double *day,double *hour,int LST,double longitude);
double Compact (double y, double x);
double Table_Value(double a, double b, double c, double d, double Diff_ba, double Diff_ca);
double LST( double jd, double longitude);
double T_2000 (double jd);
double C_1900 (double jd);
double J_1800 (double jd);
/*double Obl(double jd);*/
double Obl0_2000(double jd);
double Obl_2000(double jd);
void   LamToRa(double jd,double Lam,double Bet,double *Ra,double *Dec);
void   RaToLam(double jd,double RA,double DEC,double *Lam,double *Bet);
void   Az_El(double jd,double Ra, double Dec, double latitude,double longitude,  double *Az, double *El );
void   DateToYMD(double Date,double *Year,double *Month,double *Day);
double HMSToHour(double HMS);
double HMSToDeg(double HMS);/* 1999/07/01 added*/
double DMSToDeg(double DMS);
double DegToHMS(double D);
double DegToDMS(double D);
double Digit3( double x);
double Digit2( double x);
void   Precession (double jd, double Year0, double RA0,double DEC0,double V1,double V2,double *RA,double *DEC);

/* Nutation.c */
void Nutation_data(double jd,double *dF, double *dE);
void Precession_Data (double jd,Matrix_3 *Pr);
void Nutation (FILE *fp_Bin ,int DE,double jd,Matrix_3 *Nu);

/* Matrix.c */
void Sum_Matrix_3(Matrix_3 abc,Matrix_3 def,Matrix_3 *ghi);
void Multi_Matrix_3(Matrix_3 abc,Matrix_3 def,Matrix_3 *ghi);
void Sub_Matrix_3(Matrix_3 abc,Matrix_3 def,Matrix_3 *ghi);
void Multi_Const_Matrix_3(Matrix_3 abc,double a,Matrix_3 *ghi);
void Sum_Matrix_1(Matrix_1 ab,Matrix_1 de,Matrix_1 *gh);
void Multi_Matrix_1(Matrix_3 abc,Matrix_1 de,Matrix_1 *gh);
void Sub_Matrix_1(Matrix_1 ab,Matrix_1 de,Matrix_1 *gh);
void Multi_Const_Matrix_1(Matrix_1 ab,double a,Matrix_1 *gh);
void Eql_Matrix_1(Matrix_1 ab, Matrix_1 *gh);
void Scalar (Matrix_1 ab, double *a);
void Scalar_Dot (Matrix_1 ab,Matrix_1 de, double *a);
void Matrix_R1(Matrix_3 abc,double a, Matrix_3 *ghi);
void Matrix_R2(Matrix_3 abc,double a, Matrix_3 *ghi);
void Matrix_R3(Matrix_3 abc,double a, Matrix_3 *ghi);


/*jpl_eph.c*/
void ERROR(int i,char *msg);
double *dvector(int nl,int nh);
void free_dvector(double *v,int nl);
char **Cmatrix(int nrl, int nrh, int ncl, int nch);
void free_Cmatrix(char **m, int nrl, int ncl);
void GO2POS(FILE *fp,int nr);
FILE *INITIAL(int DE);
//void INITIAL(int DE,FILE *fp_Bin);
void PRINTOUT(void);
double Dint(double x);
int SPLIT(double tt,double *fr);
double Dmod(double x,double y);

int INTERP(double *buf,double *t,int ncf,int ncm,int na,int fl,double *pv);
int STATE (FILE *fp_Bin,double *jed,int *list,double *pv,double *nut);
int PLEPH0(FILE *fp_Bin,int n,double jd,int targ,int cent,double *rrd,int *inside,double *jd2);
int PLEPH (FILE *fp_Bin,double   jd,int targ,int cent,double *rrd,int *inside);
int DPLEPH(FILE *fp_Bin,double *jd2,int targ,int cent,double *rrd,int *inside);

int INTERPx9(double *buf,double *t,int ncf,int ncm,int na,int fl,double *pv);
int STATEx9 (FILE *fp_Bin,double *jed,int *list,double *pv,double *nut);
int PLEPH0x9(FILE *fp_Bin,int n,double jd,int targ,int cent,double *rrd,int *inside,double *jd2);
int PLEPHx9 (FILE *fp_Bin,double   jd,int targ,int cent,double *rrd,int *inside);
int DPLEPHx9(FILE *fp_Bin,double *jd2,int targ,int cent,double *rrd,int *inside);

int INTERPx12(double *buf,double *t,int ncf,int ncm,int na,int fl,double *pv);
int STATEx12 (FILE *fp_Bin,double *jed,int *list,double *pv,double *nut);
int PLEPH0x12(FILE *fp_Bin,int n,double jd,int targ,int cent,double *rrd,int *inside,double *jd2);
int PLEPHx12 (FILE *fp_Bin,double   jd,int targ,int cent,double *rrd,int *inside);
int DPLEPHx12(FILE *fp_Bin,double *jd2,int targ,int cent,double *rrd,int *inside);

void Ephemerides(FILE *fp_Bin,int DE, double jd,int Planet);

void Conversion       ( const double *state1,  double *state3, const double MUSUN, double jd);
void Conversion_J2000 ( const double *state1,  double *state3, int n);


#endif	//ASTRO_H
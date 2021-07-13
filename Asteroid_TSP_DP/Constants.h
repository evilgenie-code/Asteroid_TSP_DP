#ifndef Constants_H
#define Constants_H

#include <cmath>

extern double UnitM, UnitP, Unitb;

const double M_E = 2.71828182845904523536;
const double PI  = 3.14159265358979323846264338327950288419716939937510;  
const double PI2 = 2.*PI;
const double raddeg = 180/PI;

const double AU = 1.49597870691e8;                  // km             астрономическа€ еденица
const double fM_Sun = 1.32712440018e11;             // km^3/sec^2     гравитационный параметр —олнца
const double fM_E = 398600.4415;                    // km^3/sec^2     гравитационный параметр «емли
const double R_E = 6371;                            // km             радиус «емли
const double Gacc_E = 9.80665*1.e-3;                // km/sec^2       ускорение свободного падени€ у «емли

const double UnitR   = AU;					          
const double UnitV   = sqrt(fM_Sun/UnitR);          // km/sec
const double UnitT   = (UnitR/UnitV)/86400;         // day  
const double UnitA   = fM_Sun/(UnitR*UnitR);        // km/sec^2
const double UnitJ   = UnitA/(UnitT*86400);         // km/sec^3
const double UnitAcc = 1000*fM_Sun/(UnitR*UnitR);   // m/sec^2 
//    double UnitM  ;// = 1000;                        // kg
//    double UnitP  ;// = UnitAcc*UnitM;               // newton
//    double Unitb  ;// = UnitP/(1000*UnitV);          // kg/sec
//    double UnitN      = UnitP*UnitV;			       // kW (к¬т) (newton*km/s)

const double UnitR_E   = R_E;			
const double UnitV_E   = sqrt(fM_E/UnitR_E);          // km/sec
const double UnitT_E   = (UnitR_E/UnitV_E)/86400;     // day  
const double UnitA_E   = fM_E/(UnitR_E*UnitR_E);      // km/sec^2 
const double UnitAcc_E = 1000*fM_E/(UnitR_E*UnitR_E); // m/sec^2 
const double UnitP_E   = UnitAcc_E*UnitM;             // newton

#endif	// Constants_H
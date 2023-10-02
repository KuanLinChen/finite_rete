#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <vector>
#include <map>
#include "nasa_poly.hpp"
#include "read_reaction.hpp"

using namespace std;
class jacobian
{
  public:
  // CThermo *thermo;
  // reaction *react;
  double **A, **A2, *B, *BB,**XLAM, *GMOLCC, *F, *FM3, *XKN, *XKF, *XN1, *XN2, *XN3, *WDT;
  double TMOLN;
  void init( reaction *react, CThermo *thermo)
  { 
    // thermo = tt;
    // react  = rr;

    A = new double *[react->number_of_species];
    for(int i=0; i < react->number_of_species ; i++ ){
      A[i] = new double[ react->number_of_species ];
    }
    A2 = new double *[react->number_of_species];
    for(int i=0; i < react->number_of_species ; i++ ){
      A2[i] = new double[ react->number_of_species ];
    }

    B = new double [react->number_of_reactions];
    BB = new double [react->number_of_reactions];
    FM3 = new double [react->number_of_reactions];
    XKN = new double [react->number_of_reactions];
    XKF = new double [react->number_of_reactions];

    GMOLCC = new double [react->number_of_species];
    F      = new double [react->number_of_species];//gibb's free energy
    XN1      = new double [react->number_of_species];
    XN2      = new double [react->number_of_species];
    XN3      = new double [react->number_of_species];
    WDT      = new double [react->number_of_species];
    XLAM = new double *[react->number_of_reactions];
    for(int i=0; i < react->number_of_reactions ; i++ ){
      XLAM[i] = new double[ react->number_of_species ];
    }

  } ;
  void compute_jac_matrix(double density, double pressure, double temperature, double *alpha, double DT, reaction *react, CThermo *thermo);

};

void jacobian::compute_jac_matrix(double density, double pressure, double temperature, double *alpha, double DT,  reaction *react, CThermo *thermo)
{

  //double ESP   = 0.20   ;  // time step increment control 
  //double DTMIN = 1.0E-10;  // minimum time step size
  //double ATOL  = 1.0E-8 ;  // species sensitivity
  //double NSTEP = 500    ;  // maximum integration steps
  //double XDT = DTMIN;

  //-----START CALCULATIONS
  for (int i=0;i<react->number_of_species;i++){//100
    if(alpha[i] <= 1.0E-29) alpha[i] = 0.0;
    WDT[i] = 0.0;
    XN2[i] = 0.0;
    XN1[i] = alpha[i];
  }


  int NORDER = 1;
  //------ 1st order implicit
  double C1JAC=-1.0;
  double C2JAC=0.0;
  //------ 2nd order parasol 
  if(NORDER == 2){
    C1JAC=-0.5;
  } else if(NORDER==3){
  //------ 4th order parasol 
    C1JAC=-0.5;
    C2JAC=1./12;
  }

  TMOLN=0.0;

  for (int i=0;i<react->number_of_species;i++){//4010
    XN2[i]=alpha[i];
    XN3[i]=XN1[i]/thermo->thermo_data[i].get_molecular_wt();
    XN3[i]=alpha[i]/thermo->thermo_data[i].get_molecular_wt();
    TMOLN=TMOLN+XN3[i];
  }//4010


  for (int i=0;i<react->number_of_species;i++){//4020
    XN3[i]=XN3[i]/TMOLN;
  }

  //cout<<endl;
  double FACT  = 1.0E-03;
  for (int i=0;i<react->number_of_reactions;i++){
    for(int j=0;j<react->number_of_species;j++){
      XLAM[i][j]=0.0;
    }
  }//4100
  double SUM3=0.0;

  for(int j=0;j<react->number_of_species;j++){

    GMOLCC[j]=density*alpha[j]/thermo->thermo_data[j].get_molecular_wt()*FACT;
    SUM3 = SUM3 + GMOLCC[j];
    F[j]=thermo->thermo_data[j].get_g(temperature, pressure)*thermo->thermo_data[j].get_molecular_wt()/thermo->thermo_data[j].get_R_univ()/temperature;//convert to non-dimensional
    //cout<<"GMOLCC: "<<GMOLCC[j]<<endl;
    //cout<<"F: "<<F[j]<<endl;
  }
  //exit(1);
  //cout<<"SUM3: "<<SUM3<<endl<<endl;

  double RGGCGS= 82.047163;
  double TQZ   = temperature;
  double TQX   = 1.0/(RGGCGS*TQZ);
  double TQY   = 1.0/TQZ;

  double RNT=0.0, RRF=0.0;
  for (int i=0;i<react->number_of_reactions;i++){//4200

    RNT = 0.0;
    RRF = 0.0;

    FM3[i]=1.0;
    if( react->channels[i].Third_body == "yes" ){
      FM3[i] = SUM3;
    }


    for(int j=0;j<react->number_of_species;j++){//4250
      RNT = RNT+react->STCOEF[i][j];
      RRF = RRF+react->STCOEF[i][j]*F[j];
    }//4250

    //XKN(I) = TQX**RNT*EXP(AMAX1(-65.0,AMIN1(65.0,-RRF)))
    XKN[i]= pow( TQX, RNT )*exp(max(-65.0,min(65.0,-RRF)));
    //cout<<"i: "<<i<<endl;
    //cout<<"type: "<<react->channels[i].type <<endl;
    //cout<<"XKN: "<<XKN[i]<<endl;

    //XKF(I) = ARRHA(I)*TQZ**(-ARRHN(I))*EXP(-ARRHB(I)*TQY)
    XKF[i]=  react->channels[i].A*pow(TQZ, (react->channels[i].B) )*exp(-react->channels[i].E_R*TQY) ;
    //cout<<"XKF: "<<XKF[i]<<endl<<endl;;

  }//4200


//
//     CALCULATE RATES FOR EACH REACTION
//
  double RATEF=0.0, RATEB=0.0;

  for (int i=0;i<react->number_of_reactions;i++){//4300

    RATEF = FM3[i];
    RATEB = FM3[i];
  //--REGULAR
    //if( react->channels[i].type == "Elementary" or react->channels[i].type == "3rd-body" ){
    if( react->channels[i].type == "Elementary" ){

      for(int j=0;j<react->number_of_species;j++){//4400

        if(react->STCOEF[i][j] < 0.0){
          if(GMOLCC[j] > 1.0E-30){
            RATEF = RATEF*pow( GMOLCC[j], (-react->STCOEF[i][j]) );
          }else{
            RATEF = 0.0;
          }
        }else if(react->STCOEF[i][j] > 0.0){
          if( GMOLCC[j]>1.0E-30){
            RATEB = RATEB*pow( GMOLCC[j], ( react->STCOEF[i][j]) );
          }else{
            RATEB = 0.0;
          }
        }
      }//4400

      RATEF=XKF[i]*RATEF;
      RATEB=XKF[i]*RATEB/XKN[i];
      // cout<<"Elementary"<<endl;
      // cout<<"RATEF: "<<RATEF<<endl;
      // cout<<"RATEB: "<<RATEB<<endl;


      for(int j=0;j<react->number_of_species;j++){//4450
        if( react->STCOEF[i][j] < 0.0){
          if( GMOLCC[j] > 1.0E-30){
            XLAM[i][j]=-react->STCOEF[i][j]*RATEF/GMOLCC[j];
          }
        }else if( react->STCOEF[i][j] > 0.0){
          if(GMOLCC[j] > 1.0E-30){
            XLAM[i][j]=-react->STCOEF[i][j]*RATEB/GMOLCC[j];
          } 
        } 
      }//4450

      BB[i]=RATEF-RATEB;
      //cout<<"B: "<<BB[i]<<endl;
//--GLOBAL-1      
    // }else if( react->channels[i].type == "tmp" ) {

    //   RATEB = 0.0;
    //   for(int j=0;j<react->number_of_species;j++){
    //     if(react->STCOEF[i][j] < 0.0){
    //       if( GMOLCC[j] > 1.0E-30){
    //         RATEF = RATEF*pow(GMOLCC[j],(-react->STCOEF[i][j]));
    //       }else{
    //         RATEF = 0.0;
    //       }
    //     } 
    //   }//4500

    //   RATEF=XKF[i]*RATEF;

    //   for(int j=0;j<react->number_of_species;j++){
    //     if(react->STCOEF[i][j] < 0.0){
    //       if( GMOLCC[j] > 1.0E-30 ) {
    //         XLAM[i][j]=-react->STCOEF[i][j]*RATEF/GMOLCC[j];
    //       }
    //     }
    //   }//4550
    //   BB[i]=RATEF;
  //--GLOBAL-2
    }else if( react->channels[i].type == "Global" ){

      RATEB = 0.0;
      for(int j=0;j<react->number_of_species;j++){//4600
        if(react->STCOEG[i][j] < 0.0){
          if (GMOLCC[j] > 1.0E-30){
            RATEF = RATEF*pow(GMOLCC[j],(-react->STCOEG[i][j]));
          }else{
            RATEF = 0.0;
          }
        } 
      }//4600

      RATEF=XKF[i]*RATEF;
      for(int j=0;j<react->number_of_species;j++){
        if(react->STCOEG[i][j] < 0.0){
          if(GMOLCC[j] > 1.0E-30){
            XLAM[i][j]=-react->STCOEG[i][j]*RATEF/GMOLCC[j];
          } 
        } 
      }//4650
      BB[i]=RATEF;
    }//endif
  }//4300


//
//     Fill arrays for Gauss elimination
//
  for(int k=0;k<react->number_of_species;k++){//4700
    B[k]=0.0;

    for (int i=0;i<react->number_of_reactions;i++){//4705
      B[k]=B[k]+react->STCOEF[i][k]*BB[i];
    }//4705

    for(int j=0;j<react->number_of_species;j++){//4710
       A[k][j]=0.0;
      A2[k][j]=0.0;
      for (int i=0;i<react->number_of_reactions;i++){//4720
       A[k][j]=A[k][j]+react->STCOEF[i][k]*XLAM[i][j];
      }//4720
    }//4710
  }//4700


  double SQJAC=0.0; 
  if( NORDER == 3){
    for(int i=0;i<react->number_of_species;i++){
      for(int j=0;j<react->number_of_species;j++){//4740
        SQJAC=0.0;
        for(int k=0;k<react->number_of_species;k++){//4750
          SQJAC=SQJAC+A[i][k]*A[k][j];
        }
        A2[i][j]=SQJAC;
      }
    }
  }  


  for(int i=0;i<react->number_of_species;i++){
    for(int j=0;j<react->number_of_species;j++){
      A[i][j]=C1JAC*A[i][j]*DT+C2JAC*A2[i][j]*DT*DT;
    }
      A[i][i]=A[i][i]+1.0;
      B[i]=B[i]*DT;
  }
    // for(int i=0;i<react->number_of_reactions;i++){
    //   for(int j=0;j<react->number_of_species;j++){
    //     cout<<react->STCOEF[i][j]<<setw(5);
    //   }
    //   cout<<endl;
    // }

  // for(int i=0;i<react->number_of_species;i++){
  //   for(int j=0;j<react->number_of_species;j++){
  //       //cout<<std::fixed <<scientific<< std::setprecision(3)<<A[i][j]<<setw(10);
  //       cout<<scientific<< std::setprecision(2)<<A[i][j]<<setw(10);
  //     }
  //     cout<<endl;
  //   }
  //   cout<<endl;
  // for(int i=0;i<react->number_of_species;i++){
  //   cout<<B[i]<<endl;
  // }


};
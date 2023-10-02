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

using namespace std;
class CFormula
{
  public:
    string name;
    CFormula( string ){}
    double get_cp_0( double, double ) { return 0.0; }
    double get_h_0 ( double, double ) { return 0.0; }
    double get_s_0 ( double, double ) { return 0.0; }
    double get_g_0 ( double, double ) { return 0.0; }
    double h_f_0;
    const double R_univ = 8.314510;   // J/mol-K 
    //const double R_univ = 8.314510E3;   
    double molecular_wt;
    double get_R(){
      return R_univ/molecular_wt;
    }
    double get_molecular_wt(){
      return molecular_wt;
    }
};
class CRealFluid
{
  public:
    string name;
    CRealFluid( string ){}
    double get_cp_0( double, double ) { return 0.0; }
    double get_h_0 ( double, double ) { return 0.0; }
    double get_s_0 ( double, double ) { return 0.0; }
    double get_g_0 ( double, double ) { return 0.0; }
    double h_f_0;
    const double R_univ = 8.314510;   // J/mol-K 
    //const double R_univ = 8.314510E3;   
    double molecular_wt; 
    double get_R(){
      return R_univ/molecular_wt;
    }
    double get_molecular_wt(){
      return molecular_wt;
    }

};

class CCEA
{
  public:
    CCEA( string name ){
      get_species(name) ;
    }
    void get_species( string species_name ) {
      //To uppercase
      for (long unsigned int x=0; x<std::strlen(species_name.c_str()); x++){
        species_name[x]=toupper(species_name[x]);
      }

      ifstream infile("nasa-9.dat");
      string line, tmp_line ;
      while( infile.peek() != EOF ){
        getline( infile, line);
        tmp_line = line.substr( 0,16).c_str();
        tmp_line.erase( remove( tmp_line.begin(), tmp_line.end(), ' ' ), tmp_line.end() );
        //To uppercase
        for (long unsigned int x=0; x<std::strlen(tmp_line.c_str()); x++){
          tmp_line[x]=toupper(tmp_line[x]);
        }
        if ( strcmp(species_name.c_str(), tmp_line.c_str()) == 0 ){
          break;
        }else if(infile.peek() == EOF){
          cerr<<endl;
          cerr<<"species not find!! pls check"<<endl;
          cerr<<"species name: "<<species_name<<endl;
          exit(1);
        }
      }//Find out the species location


      //Line 1:
        name   = line.substr( 0, 15);
        for (long unsigned int x=0; x<std::strlen(name.c_str()); x++)
          putchar(toupper(name[x]));

        source = line.substr(18, 79);

      //Line 2:
        getline( infile, line);
        num_T_inteval = atoi(line.substr( 0,  2).c_str());
        ref_data_code =      line.substr( 3,  9         );
        chem_formula  =      line.substr(10, 50         );
        phase         = atoi(line.substr(50, 52).c_str());
        molecular_wt  = atof(line.substr(52, 65).c_str());
        h_f_0         = atof(line.substr(65, 80).c_str());

        for( int i=0 ; i < num_T_inteval ; i++ ){
        //Line 3:
          getline( infile, line);
          //temperature range
          tmp_line = line.substr(0, 22) ;
          T_range[i][0] = atof( tmp_line.substr( 1, 11).c_str() ) ;
          T_range[i][1] = atof( tmp_line.substr(11, 22).c_str() ) ;
          delta_h_0 = atof(line.substr(65, 80).c_str()) ;
        //Line 4:
          getline( infile, line);
          tmp_line = line.substr(0, 80) ;
          replace(tmp_line.begin(), tmp_line.end(), 'D', 'e'); 
          A[i][0] = atof( tmp_line.substr( 0, 16).c_str() ) ;
          A[i][1] = atof( tmp_line.substr(16, 32).c_str() ) ;          
          A[i][2] = atof( tmp_line.substr(32, 48).c_str() ) ;          
          A[i][3] = atof( tmp_line.substr(48, 64).c_str() ) ;          
          A[i][4] = atof( tmp_line.substr(64, 80).c_str() ) ;          
        //Line 5:
          getline( infile, line);
          tmp_line = line.substr(0, 80) ;
          replace(tmp_line.begin(), tmp_line.end(), 'D', 'e'); 
          A[i][5] = atof( tmp_line.substr( 0, 16).c_str() ) ;  
          A[i][6] = atof( tmp_line.substr(16, 32).c_str() ) ;  
          B[i][0] = atof( tmp_line.substr(48, 64).c_str() ) ;  
          B[i][1] = atof( tmp_line.substr(64, 80).c_str() ) ;  
        }
        infile.clear();
        infile.close();
      //print_info();
    }
    //exit(1);
    void get_coefficient( double T ){
      int inteval ;
      if(num_T_inteval>0){    
        for( int i=0 ; i < num_T_inteval ; i++ ){
          if( T >= T_range[i][0] and T <=  T_range[i][1]){
            inteval = i;
          } else if( T < T_range[0][0] ){
            inteval=0;
            cerr<<name<<" reach lower-bound"<<", T: "<<T<<endl;

          } else if(T>T_range[num_T_inteval-1][1] ){
            inteval=num_T_inteval-1;
            cerr<<name<<" reach upper-bound"<<", T: "<<T<<endl;
          }//end-if
        }//end-for
        a1 = A[inteval][0] ;
        a2 = A[inteval][1] ;
        a3 = A[inteval][2] ;
        a4 = A[inteval][3] ;
        a5 = A[inteval][4] ;
        a6 = A[inteval][5] ;
        a7 = A[inteval][6] ;
        b1 = B[inteval][0] ;
        b2 = B[inteval][1] ;
      }else{
        a1 = A[ 0][0] ;
        a2 = A[0][1] ;
        a3 = A[0][2] ;
        a4 = A[0][3] ;
        a5 = A[0][4] ;
        a6 = A[0][5] ;
        a7 = A[0][6] ;
        b1 = B[0][0] ;
        b2 = B[0][1] ;        
      }
    }
    //heat capacity
    double get_cp_0( double T ){ // J/mol-K
      get_coefficient(T);
      return R_univ*(a1*pow(T,-2)+a2*pow(T,-1)+a3+a4*pow(T, 1)+a5*pow(T, 2)+a6*pow(T, 3)+a7*pow(T, 4))/molecular_wt;//[kJ/kg/k]
    }
    //enthalpy
    double get_h_0( double T ){ // J/mol
      get_coefficient(T);
      return R_univ*T*( -a1*pow(T,-2)+a2*log(T)/T+a3+a4*T/2.0+a5*(T*T)/3.0+a6*(T*T*T)/4.0+a7*(T*T*T*T)/5.0+b1/T )/molecular_wt ; //[kJ/kg]
    }
    //entropy
    double get_s_0( double T ){ // J/mol-K
      get_coefficient(T);
      return R_univ*( -a1*pow(T,-2)/2.0-
                  a2*pow(T,-1)+a3*log(T)+a4*T+a5*(T*T)/2.0+a6*(T*T*T)/3.0+a7*(T*T*T*T)/4.0+b2 )/molecular_wt ; //[kJ/kg/k]
    }
    //Gibb's free energy
    double get_g_0( double T ){ // J/mol-K
      return ( get_h_0(T)-T*get_s_0(T) );
    }
    double get_R(){
      return R_univ/molecular_wt;
    }
    // double get_R_univ(){
    //   return R_univ;
    // }
    string name, source, ref_data_code, chem_formula, T_range_str;
    double a1,a2,a3,a4,a5,a6,a7,b1,b2;
    const double R_univ = 8.314510;   // J/mol-K 
    //const double R_univ = 8.314510E3;   //????
    double molecular_wt, /// g/mol
           T_range[3][2],// three inteval, lower&higer bound
                 A[3][7],// three inteval, 7 coefficient
                 B[3][2],// three inteval, 7 coefficient
           h_f_0,        // Heat of formation at 298.15 K, J/mol                     
           delta_h_0;    // h_0(298.15)–h_0(0) J/mol, if available  

    int num_T_inteval, 
        //coefficients,  // always 7
        phase;           //0: 'gas', 1: 'solid'
        
    void print_info(){
      cout<<"name         : "<<name <<endl;
      cout<<"source       : "<<source <<endl;
      cout<<"num_T_inteval: "<<num_T_inteval<<endl;
      cout<<"ref_data_code: "<<ref_data_code<<endl;
      cout<<"chem_formula : "<<chem_formula<<endl;
      cout<<"phase        : "<<phase<<endl;
      cout<<"molecular_wt : "<<molecular_wt<<endl;
      cout<<"h_f_0        : "<<h_f_0<<endl;
      for( int i=0 ; i < num_T_inteval ; i++ ){
        cout<<"T_range    : "<<T_range[i][0]<<"\t"<<T_range[i][1]<<"\t"<<endl;
        cout<<"A          : "<<A[i][0]<<"\t"<<A[i][1]<<"\t"<<A[i][2]<<"\t"<<A[i][3]<<"\t"<<A[i][4]<<"\t"<<A[i][5]<<"\t"<<A[i][6]<<"\t"<<endl;
        cout<<"B          : "<<B[i][0]<<"\t"<<B[i][1]<<"\t"<<endl;
        //cout<<endl;
      }
      cout<<endl;
    }
};
class CThermoData
{
  public:
  int model;
  string name;
  CCEA *CEA=NULL;
  CFormula *Formula=NULL;
  CRealFluid *RealFluid=NULL;

  CThermoData(){}

  void allocate(int m, string n){
    //cout<<"CThermoData->allocate, m: "<<m<<", n: "<<n<<endl;
    model = m;
    name  = n;
    //To uppercase
    for (long unsigned int x=0; x<std::strlen(name.c_str()); x++){
      name[x]=toupper(name[x]);
    }
    //cout<<"CThermoData->model = "<<model<<endl;
    switch(model){
      case 0:       CEA = new CCEA(name);       break;
      case 1: RealFluid = new CRealFluid(name); break;
      case 2:   Formula = new CFormula(name);   break;
      default: cout<<"Create model error"<<endl;exit(1);break;
    }
  }


  double get_cp(double T, double P){
    double value=0.0;
    switch(model) { 
      case 0: value =       CEA->get_cp_0(T  ); break;
      case 1: value = RealFluid->get_cp_0(T,P); break;
      case 2: value =   Formula->get_cp_0(T,P); break;
      default: cout<<"get_cp_0 error"<<endl;exit(1);break;
    }
    return value;
  }
  double get_h0(double T, double P){
    double value=0.0;
    switch(model) { 
      case 0: value =       CEA->get_h_0(T  ); break;
      case 1: value = RealFluid->get_h_0(T,P); break;
      case 2: value =   Formula->get_h_0(T,P); break;
      default: cout<<"get_h_0 error"<<endl;exit(1);break;
    }
    return value;
  }
  double get_s(double T, double P){
    double value=0.0;
    switch(model) { 
      case 0: value =       CEA->get_s_0(T  ); break;
      case 1: value = RealFluid->get_s_0(T,P); break;
      case 2: value =   Formula->get_s_0(T,P); break;
      default: cout<<"get_s_0 error"<<endl;exit(1);break;
    }
    return value;
  }
  double get_g(double T, double P){
    double value=0.0;
    switch(model) { 
      case 0: value =       CEA->get_g_0(T  ); break;
      case 1: value = RealFluid->get_g_0(T,P); break;
      case 2: value =   Formula->get_g_0(T,P); break;
      default: cout<<"get_g_0 error"<<endl;exit(1);break;
    }
    return value;
  }
  double get_h_f_0(double T, double P){
    double value=0.0;
    switch(model) { 
      case 0: value =       CEA->h_f_0; break;
      case 1: value = RealFluid->h_f_0; break;
      case 2: value =   Formula->h_f_0; break;
      default: cout<<"get_h_f_0 error"<<endl;exit(1);break;
    }
    return value;
  }
  double get_R(){
    double value=0.0;
    switch(model) { 
      case 0: value =       CEA->get_R(); break;
      case 1: value = RealFluid->get_R(); break;
      case 2: value =   Formula->get_R(); break;
      default: cout<<"get_R error"<<endl;exit(1);break;
    }
    return value;
  }
  double get_R_univ(){
    double value=0.0;
    switch(model) { 
      case 0: value =       CEA->R_univ; break;
      case 1: value = RealFluid->R_univ; break;
      case 2: value =   Formula->R_univ; break;
      default: cout<<"get_R error"<<endl;exit(1);break;
    }
    return value;
  }
  double get_molecular_wt(){
    double value=0.0;
    //cout<<"model: "<<model<<endl;
    switch(model) { 
      case 0: value =       CEA->molecular_wt; break;
      case 1: value = RealFluid->molecular_wt; break;
      case 2: value =   Formula->molecular_wt; break;
      default: cout<<"molecular_wt error"<<endl;exit(1);break;
    }
    return value;
  }
  string get_name(){
    string value;
    //cout<<"model: "<<model<<endl;
    switch(model) { 
      case 0: value =       CEA->name; break;
      case 1: value = RealFluid->name; break;
      case 2: value =   Formula->name; break;
      default: cout<<"get name error"<<endl;exit(1);break;
    }
    return value;
  }
};
class CThermo
{
  public:
  int nSpecies;
  vector<CThermoData> thermo_data ;

  CThermo( int n, string *name, int *m ){
    nSpecies = n;
    CThermoData init_vector;
    for(int i=0 ; i<nSpecies ; i++ ){
      thermo_data.push_back(init_vector);
        //to upper
        for (long unsigned int x=0; x<std::strlen(name[i].c_str()); x++){
          name[i][x]=toupper(name[i][x]);
        }
      thermo_data[i].allocate(m[i], name[i]);
      //cout<<", model: "<<m[i]<<endl;
    }
  }
  double get_cp(double T, double P, double *Y){
    double mixtue_cp=0.0;
    for (int i=0 ; i < nSpecies ; i++ )
      mixtue_cp += Y[i]*thermo_data[i].get_cp(T,P);
    return mixtue_cp;
  }
  double get_h0(double T, double P, double *Y){
    double mixtue_h0=0.0;
    for (int i=0 ; i < nSpecies ; i++ ){
      mixtue_h0 += Y[i]*thermo_data[i].get_h0(T,P);
    }
    return mixtue_h0;
  }
  double get_s(double T, double P, double *Y){
    double mixtue_s=0.0;
    for (int i=0 ; i < nSpecies ; i++ )
      mixtue_s += Y[i]*thermo_data[i].get_s(T,P);
    return mixtue_s;
  } 
  double get_g(double T, double P, double *Y){
    double mixtue_g=0.0;
    for (int i=0 ; i < nSpecies ; i++ )
      mixtue_g += Y[i]*thermo_data[i].get_g(T,P);
    return mixtue_g;
  }   
  double get_h_f_0(double T, double P, double *Y){
    double mixtue_h_f_0=0.0;
    for (int i=0 ; i < nSpecies ; i++ )
      mixtue_h_f_0 += Y[i]*thermo_data[i].get_h_f_0(T,P);
    return mixtue_h_f_0;
  } 
  double get_R(double T, double P, double *Y){
    double mixtue_R=0.0;
    for (int i=0 ; i < nSpecies ; i++ )
      mixtue_R += Y[i]*thermo_data[i].get_R();
    return mixtue_R;
  } 

  double get_molecular_wt(double T, double P, double *Y){
    double mixtue_R=0.0;
    for (int i=0 ; i < nSpecies ; i++ )
      mixtue_R += Y[i]*thermo_data[i].get_molecular_wt();
    return mixtue_R;
  } 

  void list_name(){
    //cout<<"AA"<<endl;
    for (int i=0 ; i < nSpecies ; i++ ) cout<<"species["<<i<<"]: "<< thermo_data[i].get_name() <<endl;
  }
  double get_Tm(double mixtue_h0_in, double T, double P, double *Y ){
    double T_old, T_new, P_old, f, df, mixtue_h0_test, mixtue_cp_test;
    int kk=1;
    T_new=T;
    P_old=P;
    /* newton method */
    for( int k=1 ; k <= kk ; k++ ){
      if( kk >= 50 ){
        cout<<"Conversion from enthalpy to temperature does not converge within the set limit(>50)"<<endl;
        exit(1);
      }
      T_old=T_new;
      mixtue_h0_test=get_h0(T_old, P_old, Y);
      mixtue_cp_test=get_cp(T_old, P_old, Y);
      f = mixtue_h0_test - mixtue_h0_in;
      df= mixtue_cp_test;
      T_new = T_old - f/df;
      if ( fabs(f/df) > 1E-6 ) kk++ ;
    }
    return T_new;
  }

  double compute_total_energy( double T, double P, double u, double v, double w, double *Y){
    double Et=0.0; 
    double mixture_Enth = get_h0(T, P, Y);
    double mixture_h_f_0 = get_h_f_0(T, P, Y);
    double mixture_R =  get_R(T, P, Y);
    double rho = P/mixture_R/T;
    double KE=0.5*(u*u+v*v+w*w);
    return Et = mixture_Enth-mixture_h_f_0-(P/rho)+KE; 
  }

  double compute_enthalpy( double Et, double T, double P, double u, double v, double w, double *Y){
    double enthalpy=0.0;
    double mixture_h_f_0 = get_h_f_0(T, P, Y);
    double mixture_R=  get_R(T, P, Y);
    double rho = P/mixture_R/T;
    double KE=0.5*(u*u+v*v+w*w);
    return enthalpy = Et-KE+(P/rho)+mixture_h_f_0; 
  }

};

// int main()
// {
//   /* initialize some variables */

//     int opt;
//     //cout << "To compute mixture enthalpy from temperature enter 1:\n";
//     cout << "To compute mixture enthalpy from temperature enter 1, and ";
//     cout << "to compute mixture temperature from enthalpy enter 2: ";
//     cin >> opt;
//     cout << endl;
//     if(opt < 1 || opt > 2) {
//       cout << "@@Error: Wrong option!!";
//       return 0;
//     }

//     //int nSpecies = 3;
//     int nSpecies;
//     cout << "Enter number of chemical species = ";
//     cin >> nSpecies;
//     cout << endl;
//     double Y[nSpecies];
//     double Y_sum = 0.0;
//     int  Model[nSpecies];
//     string species_name[nSpecies];

//     for(int i = 0; i < nSpecies; ++i) {
//       int k = i+1;
//       cout << "Enter name of species " << k << ": ";
//       cin >> species_name[i];
//       cout << endl;
//       cout << "Enter mass fraction of species (<= 1) " << k << " = ";
//       cin >> Y[i];
//       cout << endl;
//       Y_sum = Y_sum + Y[i];
//       Model[i] = 0;
//     }
//     if(fabs(Y_sum - 1.0) > 1.0e-10) {
//       cout << "@@Error: Sum of species mass fractions = " << Y_sum << " does not equal to 1\n";
//       return 0;
//     }

//     //Y[0] = 0.6; species_name[0]="N2"; Model[0]=0;//0 for CEA, 1 for RealFluid and 2 for Formula
//     //Y[1] = 0.2; species_name[1]="O2"; Model[1]=0;
//     //Y[2] = 0.2; species_name[2]="NO"; Model[2]=0;

//     double P_test=101325.0;
//     if(opt == 1) {
//       //double T_test=365.0;
//       double T_test;
//       cout << "Enter temperature of species mixture(deg K) = ";
//       cin >> T_test;
//       cout << endl;

//       CThermo thermo(nSpecies, species_name, Model);
//       //thermo.list_name();
//       double mixture_Cp;
//       double mixture_Enth;
//       mixture_Cp = thermo.get_cp(T_test, P_test, Y);
//       mixture_Enth = thermo.get_h0(T_test, P_test, Y);
//       cout << "Mixture Cp(kJ/kg-K) = " << mixture_Cp << ", ";
//       cout << "enthalpy(kJ/kg) = " << mixture_Enth;
//       cout << endl;
//     }
//     else {
//       double Enth_test;
//       cout << "Enter enthalpy of species mixture(kJ/kg) = ";
//       cin >> Enth_test;
//       cout << endl;

//       CThermo thermo(nSpecies, species_name, Model);
//       //thermo.list_name();
//       double mixture_Tm;
//       double mixture_Cp;
//       double Tm_guess = 300.0;
//       mixture_Tm = thermo.get_Tm(Enth_test, Tm_guess, P_test, Y);
//       mixture_Cp = thermo.get_cp(mixture_Tm, P_test, Y);
//       cout << "Mixture Cp(kJ/kg-K) = " << mixture_Cp << ", ";
//       cout << "temperature(K) = " << mixture_Tm;
//       cout << endl;
//     }

//     //cout<<"Test temperature: "<<T_test<<endl;

//   /*---------------------------*/
//     //CThermo thermo( nSpecies, species_name, Model);
//     //thermo.list_name();

//   /* compute mixture h_0 for test */
//     //double mixitue_h_0=0.0;
//     //mixitue_h_0 = thermo.get_h0(T_test, P_test, Y);       //dummy
//     //cout<<"mixture temperature: "<<thermo.get_Tm( mixitue_h_0, 298.15, P_test, Y )<<endl;

//   return 0;
// }

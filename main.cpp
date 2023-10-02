#include <iostream>
#include "nasa_poly.hpp"
#include "read_reaction.hpp"
#include "jacobian.hpp"
using namespace std;
// Function to perform Gaussian elimination
// Perform Gaussian elimination
void gaussianElimination(double** A, double* B, int n) {
    for (int i = 0; i < n; i++) {
        // Find pivot row
        int pivotRow = i;
        for (int j = i + 1; j < n; j++) {
            if (abs(A[j][i]) > abs(A[pivotRow][i]))
                pivotRow = j;
        }

        // Swap rows
        swap(A[i], A[pivotRow]);
        swap(B[i], B[pivotRow]);

        // Elimination
        for (int j = i + 1; j < n; j++) {
            double ratio = A[j][i] / A[i][i];
            for (int k = i; k < n; k++)
                A[j][k] -= ratio * A[i][k];
            B[j] -= ratio * B[i];
        }
    }
}

// Perform back substitution and solve for x
void backSubstitution(double** A, double* B, double* x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = B[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
}
int main() {

    reaction react;
    react.read_file( "./react.in" );
    react.bulid_STCOEFG();
    int species_size = react.species_list.size();

    string species_name[species_size];
    int    Model[species_size];
    for (int i=0; i < species_size ; i++){
      species_name[i] = react.species_list[i];
      Model[i]=0;//CEA
    }

    CThermo thermo(species_size, species_name, Model);
    
    int number_of_cell=1;
    jacobian *jac;
    jac = new jacobian[number_of_cell];

    for(int i=0;i< number_of_cell;i++ )
      jac[i].init( &react, &thermo);

    double alpha[species_size];


    //variables="Time","Tmperature","H","H2","H2O","O","O2","OH"
    alpha[0]=0.000; //H  
    alpha[1]=0.111; //H2 
    alpha[2]=0.000; //H2O
    alpha[3]=0.000; //O  
    alpha[4]=0.889; //O2 
    alpha[5]=0.000; //OH 

    // alpha[0]=0.0; //C2H4
    // alpha[1]=0.0; //CO  
    // alpha[2]=0.0; //CO2 
    // alpha[3]=0.0; //H   
    // alpha[4]=0.0; //H2  
    // alpha[5]=0.0; //H2O 
    // alpha[6]=0.0; //O   
    // alpha[7]=0.65; //O2  
    // alpha[8]=0.0; //OH  
    // alpha[9]=0.35; //RP-1

    //variables="Time","Tmperature","CH4","CO","CO2","H","H2","H2O","O","O2","OH"
    // alpha[0]=0.2; //CH4
    // alpha[1]=0.0; //CO 
    // alpha[2]=0.0; //CO2
    // alpha[3]=0.0; //H  
    // alpha[4]=0.0; //H2 
    // alpha[5]=0.0; //H2O
    // alpha[6]=0.0; //O  
    // alpha[7]=0.8; //O2 
    // alpha[8]=0.0; //OH 



    double* solution = new double[species_size];
    double *WDT=new double[species_size];
    double pressure=101325;
    double temperature=300.0;
    double ignition_temperature=1200;
    double initial_temperature=300.0;
    double FACTR = 1.0E+03;
    double dt=1.0e-12;
    double density=pressure/(8314/thermo.get_molecular_wt(initial_temperature, pressure, alpha))/initial_temperature;
    double mixtue_h0_in = thermo.get_h0(initial_temperature, pressure, alpha);
    

    cout<<"enthalpy = "<<mixtue_h0_in<<" "<<"density "<<" "<<density<<endl;
    double time=0;
    int Iter_max=500000000;
    //for (int i=0;i<Iter_max;i++){
    int i=0;
    while( time < 1.5 ){
      i++;
      time = time+dt;
      temperature = max(ignition_temperature, initial_temperature);
      jac[0].compute_jac_matrix( density, pressure, temperature, alpha, dt, &react, &thermo);

      gaussianElimination(jac[0].A, jac[0].B, species_size);
      backSubstitution(jac[0].A, jac[0].B, solution, species_size);

  
      for(int j=0;j<species_size;j++){
        alpha[j]=alpha[j]+solution[j]/density*thermo.thermo_data[j].get_molecular_wt()*FACTR;
      }

      for(int i=0;i<species_size;i++){
        if(alpha[i]<0.0)alpha[i]=0.0;
      }
      initial_temperature = thermo.get_Tm( mixtue_h0_in, initial_temperature, pressure, alpha );


#if false
      if(i%100==0){
      cout<<std::fixed <<scientific<< std::setprecision(5)<<" Time "<<setw(15)
          <<"  T   "<<setw(15)
          <<"  H   "<<setw(15)
          <<"  H2  "<<setw(15)
          <<"  H2O "<<setw(15)
          <<"  O   "<<setw(15)
          <<"  O2  "<<setw(15)
          <<"  OH  "<<setw(15)
          //<<"  SUM "<<setw(15)
          //<<" SUM_O"<<setw(15)
          //<<" SUM_H"<<setw(15)
          <<endl;
      }
#endif

      if(i%1000==0){
      cout<<std::fixed <<scientific<< std::setprecision(5)<<time<<setw(15)
          <<initial_temperature<<setw(15)
          <<alpha[0]<<setw(15)
          <<alpha[1]<<setw(15)
          <<alpha[2]<<setw(15)
          <<alpha[3]<<setw(15)
          <<alpha[4]<<setw(15)
          <<alpha[5]<<setw(15)
          // <<alpha[6]<<setw(15)
          // <<alpha[7]<<setw(15)
          // <<alpha[8]<<setw(15)
          //<<alpha[9]<<setw(15)
          //<<sum<<setw(15)
          //<<sum_O<<setw(15)
          //<<sum_H<<setw(15)
          //<<thermo.get_h0(temperature, pressure, alpha)<<setw(15)
          <<endl;
        //cout<<thermo.get_h0(temperature, pressure, alpha)<<endl;
      }

//C2H4
//CO  
//CO2 
//H   
//H2  
//H2O 
//O   
//O2  
//OH  
//RP-1


          
    }

  // 
  //   double *WDT=new double[species_size];




    return 0;
}
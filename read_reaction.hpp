#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cctype>
#include <iomanip>
using namespace std;
void remove_spaces(string& str) {
    str.erase(remove(str.begin(), str.end(), ' '), str.end());
}
struct species{
  string name;
  double coeff;
};
class channel{
public:
  vector<double> reactants_coeff;
  vector<string> reactants_name;

  vector<double> products_coeff;
  vector<string> products_name;
  double A;
  double B;
  double E_R;
  string type;
  string Third_body;
  // vector<string> dependency_name;
  vector<double> dependency_power;
};
class reaction
{
public:
  int number_of_reactions, number_of_species;
  channel *channels;
  vector<string> species_list;
  void read_file(string);
  void check_number_of_species();
  void bulid_STCOEFG();

  void test(){
    cout<<"test"<<endl;
  };

  double ** STCOEF, **STCOEG ;
};
void reaction::read_file( string file_name )
{

  ifstream reaction_file(file_name); 
    if (!reaction_file){
        cerr << "Failed to open the file." << endl;
        exit(1);
    }

  string line, LHS, RHS;

  //first get the number of reactions
    getline(reaction_file, line,'=');
    getline(reaction_file, line);
    remove_spaces(line);
    number_of_reactions=atoi(line.c_str());
    cout<<"number of reactions: "<< number_of_reactions<<endl;
    channels = new channel[number_of_reactions];



  string channel_str, tmp_str, tmp1;
  string products_str, reactants_str;
  int j = 0 ;
  while ( getline(reaction_file, line) ){
    if( line[0]!='/' ){
      channel_str="";
      channel_str = line.substr(0, line.find(';')); 
      reactants_str= channel_str.substr(0, channel_str.find('='));
      products_str= channel_str.substr(channel_str.find('=')+1);


      string token;
      /* reactants */
      stringstream react(reactants_str);
      bool species_coeff = false ;
      while ( getline(react, token, ' ')) {

          if(token=="+" or token.empty() ){
            //do nothing
          }else{

            if( isdigit( token[0] ) ){//isdigit

              channels[j].reactants_coeff.push_back( atof( token.c_str()  ) );
              species_coeff=true;
            } else if( isalpha( token[0]) ){

              channels[j].reactants_name.push_back(token);

              if( species_coeff ){
                species_coeff=false;
              } else {
                channels[j].reactants_coeff.push_back( 1.0 );
                species_coeff=false;
              }
            }
          }
      }

      /* products */
      stringstream prod(products_str);
      species_coeff = false ;
      while ( getline(prod, token, ' ')) {

          if(token=="+" or token.empty() ){
            //do nothing
          }else{

            if( isdigit( token[0] ) ){//isdigit

              channels[j].products_coeff.push_back( atof( token.c_str()  ) );
              species_coeff=true;
            } else if( isalpha( token[0]) ){

              channels[j].products_name.push_back(token);

              if( species_coeff ){
                species_coeff=false;
              } else {
                channels[j].products_coeff.push_back( 1.0 );
                species_coeff=false;
              }
            }


          }
      }


      tmp_str = line.substr(line.find(';')+1);
      line = tmp_str;

      vector<string> vec_token;
      stringstream ssss(line);
      while ( getline(ssss, token, ' ')) 
      {
        if(!token.empty()){
          vec_token.push_back(token);
        }
      }
      // for (const auto& aaa : vec_token) {
      //    //cout <<aaa<< " ";
      // }
       channels[j].A          = atof( vec_token[0].c_str() );
       channels[j].B          = atof( vec_token[1].c_str() );
       channels[j].E_R        = atof( vec_token[2].c_str() );
       channels[j].type       =       vec_token[3] ;
       channels[j].Third_body =       vec_token[4] ;

#if true
      cout<<"j: "<<j<<" "<<channels[j].A<<" "<<channels[j].B<<" "<<channels[j].E_R<<" "<<channels[j].type <<" "<<channels[j].Third_body<<endl;
#endif


      string buffer;
      if( vec_token.size() > 5  ){ 
        for(long unsigned int i=5 ; i < vec_token.size() ; i++ ) {
          buffer = vec_token[i].substr(vec_token[i].find('^')+1) ;
          //cout<<buffer<<endl;
          channels[j].dependency_power.push_back(atof(buffer.c_str()));
        }
      }

      j++;
    }
  }//end while
  reaction_file.close();

};
void reaction::bulid_STCOEFG()
{

  for(int i=0 ; i < number_of_reactions ; i++ )
  {
    for ( long unsigned  j=0; j < channels[i].reactants_name.size() ; j++ ){
      species_list.push_back(channels[i].reactants_name[j]);
    }
    for ( long unsigned  j=0; j < channels[i].products_name.size() ; j++ ){
      species_list.push_back(channels[i].products_name[j]);
    }  
  }
  /* sort and unique */
  sort( species_list.begin(), species_list.end() );
  species_list.erase( unique( species_list.begin(), species_list.end() ), species_list.end() );
  number_of_species = species_list.size();
  cout<<"number of species: "<< number_of_species <<endl;
  // for( auto & name : species_list ){
  //   cout<<name<<setw(5);
  // }
  // cout<<endl;
  


  /* bulid stoichiometric coefficient */
    STCOEF = new double *[number_of_reactions];
    for(int i=0; i < number_of_reactions ; i++ ){
      STCOEF[i] = new double[ number_of_species ];
    }

    STCOEG = new double *[number_of_reactions];
    for(int i=0; i < number_of_reactions ; i++ ){
      STCOEG[i] = new double[ number_of_species ];
    }




    for(int i=0;i<number_of_reactions;i++){
      for(int j=0;j<number_of_species;j++){
        STCOEF[i][j]=0.0;
        STCOEG[i][j]=0.0;
      }
    }

    for(int i=0;i<number_of_reactions;i++){
      for(int j=0;j<number_of_species;j++){
        for(long unsigned int k=0;k<channels[i].reactants_name.size();k++){
          if(channels[i].reactants_name[k]==species_list[j]){

            STCOEF[i][j]-=channels[i].reactants_coeff[k];

            if(channels[i].dependency_power.size() > 0 ){
              STCOEG[i][j]-=channels[i].dependency_power[k];
            }

          }
        }
      }
    }


    for(int i=0;i<number_of_reactions;i++){
      for(int j=0;j<number_of_species;j++){
        for(long unsigned int  k=0;k<channels[i].products_name.size();k++){
          if(channels[i].products_name[k]==species_list[j]){
           STCOEF[i][j]+=channels[i].products_coeff[k];
          }
        }
      }
    }

    // for(int i=0;i<number_of_reactions;i++){
    //   for(int j=0;j<number_of_species;j++){
    //     cout<<STCOEF[i][j]<<setw(5);
    //   }
    //   cout<<endl;
    // }
    // cout<<endl;
    // for(int i=0;i<number_of_reactions;i++){
    //   for(int j=0;j<number_of_species;j++){
    //     cout<<STCOEG[i][j]<<setw(5);
    //   }
    //   cout<<endl;
    // }
    // exit(1);
}

// int main()
// {
//     reaction react;
//     react.read_file( "./react.in" );
//     react.bulid_STCOEFG();
//     cout<<"end"<<endl;


//     return 0;
// }

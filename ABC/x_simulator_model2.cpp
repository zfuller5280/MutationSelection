//Written by Zach Fuller. Simulation code translated from Python script by Ipsita Agarwal
//Demographic model incorporated using code from Yuval Simons for Simons et al. (2018)
#include <chrono>
#include <random>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <map>
#include <time.h>
#include <string>
#include <iostream>
  using std::cout;
  using std::endl;
#include <iomanip>
  using std::setprecision;
#include <cstdlib>
  using std::atoi;
  using std::atof;
#include <ctime>
  using std::time;
#include <vector>
  using std::vector;
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
  using boost::poisson_distribution;
#include <boost/random/variate_generator.hpp>
  using boost::variate_generator;
#include "population_sex_diff_v2.h"
#include "BRand.hpp"
#define taujump 56000
#include <boost/math/distributions/binomial.hpp>


int gens,RUNS,length;
int popsize(int demographic_model, int dem_uncert, int gen, int i);
int popNe;
double sel,dom,mutU,Uvar;
int idx=0;
int site_class;
int exp_lof_num;
double mut_rate,expec_lof,obs_num;
int initgen;
int demographic_model;
int mut_uncert;
int dem_uncert;
int last_pop_size, pop_size;
int pid_seed;
int epsilon;
int obs;
int distance_function;
int pop_sample_size;

//Population size changes from Schiffles and Durbin (2015) inferred by MSMC of Europeans
int Ne[56]= {14448,14068,14068,14464,14464,15208,15208,16256,16256,17618,17618,19347,19347,21534,21534,24236,24236,27367,27367,30416,30416,32060,32060,31284,29404,26686,23261,18990,16490,16490,12958,12958,9827,9827,7477,7477,5791,5791,4670,4670,3841,3841,3372,3372,3287,3359,3570,4095,4713,5661,7540,11375,14310,13292,14522,613285};
int T[56]={55940,51395,47457,43984,40877,38067,35501,33141,30956,28922,27018,25231,23545,21951,20439,19000,17628,16318,15063,13859,12702,11590,10517,9482,8483,7516,6580,5672,5520,5156,4817,4500,4203,3922,3656,3404,3165,2936,2718,2509,2308,2116,1930,1752,1579,1413,1252,1096,945,798,656,517,383,252,124,0};

int randflag=0;
using namespace std;

//Functions for boost random poisson and binomial variables

int poicond1(double l);
int poi(double l);

boost::mt19937 gent;

double binomial_pmf(int N, double p, int x);
double randnum (double a, double b);
int rand_binom(int n, double p);

int main(int argc, char *argv[]){
  // Show options and order of parameters if no arguments are given
  if (argc == 1)
  {
    std::cerr << "Usage: " << argv[0] << " U obs dom runs mut_uncert dem_uncert demographic_model popNe process_seed distance_function sample_size\n"
        "Simulate an allele frequency "
        "on the X chromosome\n"
        "The selection coefficient and epsilon (if needed) "
        "are read from stdin\n";
    return 1;
  }
  if (argc==11)
     {
       Uvar=atof(argv[1]);
       //length=atof(argv[2]);
       obs=atof(argv[2]);
       //dom=atof(argv[3]);
       RUNS=atof(argv[3]);
       mut_uncert=atof(argv[4]);
       dem_uncert=atof(argv[5]);
       //Input 1 if Schiffles-Durbin demographic model should be used. Else, constant size is used
       demographic_model=atof(argv[6]);
       //Only useful if constant size population model is used
       popNe=atof(argv[7]);
       pid_seed=atof(argv[8]);
       distance_function=atof(argv[9]);
       pop_sample_size=atof(argv[10]);
     }
     else{
       cout<<"Error: Not enough parameters enetered" << endl;
     }

     // For use in pakman, epsilon followed by the selection coefficent must be read from stdin
     std::cin >> epsilon;
     std::cin >> sel >> dom;

     if (std::cin.fail())
     {
         std::cerr << "Error: could not read epsilon or selection coefficient from stdin\n";
         return 1;
     }
     // If running multiple threads or processors, pid_seed allows for different seeds to run instead of using the same CPU clock time
     gent.seed(time(NULL)+pid_seed);
     BRand::Controller.seed(time(NULL));
     std::uniform_real_distribution<double> mt_rand{0.0, 1.0};

     int gen;
     int m_gen, f_gen;

     int alpha = 3.5; //From Gao et al. 2019 Fig 2B

     //pops[0] will hold females, pops[1] will hold males
     population::initialize(sel, dom);
     population_male::initialize(sel, dom);
     population* pops[1];
     population_male* male_pops[1];
     pops[0]= new population((popsize(demographic_model,dem_uncert,0,0)/2));
     male_pops[0]= new population_male((popsize(demographic_model,dem_uncert,0,0)/2));

     //double halftheta=2*popsize(demographic_model,dem_uncert,initgen,0)*mutU;

     for(int run=0;run<RUNS;run++){
       if (mut_uncert==1){
         std::normal_distribution<double> norm_distribution(Uvar, Uvar/10);
         mutU = norm_distribution(gent);
       }
       else{
         mutU =  Uvar;
       }
       //double malemutU = ((1.*alpha)/(1+alpha)) * 2 * mutU;
       //double femalemutU = (1./(1+alpha)) * 2 * mutU;
       double malemutU = (3.*mutU*alpha)/(2.+alpha);
       double femalemutU = (3.*mutU)/(2.+alpha);

       randflag = 0;
       if (demographic_model==1){
         //Burn-in period, only used under Schiffels-Durbin model
         initgen=150000;
       }
       else{
         //Burn-in period, only used under constant population size model
         initgen=popNe*10;
       }
       idx=0;
       gen=initgen;
       pops[0]->size=int(popsize(demographic_model,dem_uncert,gen,0))/2;
       pops[0]->clear();

       male_pops[0]->size=int(popsize(demographic_model,dem_uncert,gen,0))/2;
       male_pops[0]->clear();

       while (gen>1){
         //cout<<mutU<<"\t"<<malemutU<<"\n";
          int total_allele_count = (pops[0]->allelenum()) + (male_pops[0]->allelenum());
          //cout<<gen<<"\t"<<total_allele_count<<"\n";
          if ((total_allele_count > 0)||(gen<taujump)){
            pops[0]->mutateup(poi( femalemutU*(2*pops[0]->size-pops[0]->allelenum())));
            male_pops[0]->mutateup(poi( malemutU*(male_pops[0]->size-male_pops[0]->allelenum())));

            pops[0]->mutatedown(poi( femalemutU*pops[0]->allelenum()));
            male_pops[0]->mutatedown(poi( malemutU*male_pops[0]->allelenum()));
          }
          else if ((total_allele_count==0) ){
            int f_gen = int(-log( BRand::Controller.nextOpened())* (1./(femalemutU * 2 * pops[0]->size)));
            int m_gen = int(-log( BRand::Controller.nextOpened())* (1./(malemutU * male_pops[0]->size)));
            //cout<<f_gen<<" "<<m_gen<<"\n";
            gen-=min(f_gen, m_gen);
            if (m_gen<= f_gen){
              male_pops[0]->size=popsize(demographic_model,dem_uncert,gen-1,0);
              male_pops[0]->clear();
              male_pops[0]->mutateup(poicond1(malemutU*male_pops[0]->size));
              //male_pops[0]->mutatedown(poicond1(malemutU*male_pops[0]->size));
            }
            else if (m_gen > f_gen){
              pops[0]->size=popsize(demographic_model,dem_uncert,gen-1,0);
              pops[0]->clear();
              pops[0]->mutateup(poicond1(femalemutU*2*pops[0]->size));
              //pops[0]->mutatedown(poicond1(femalemutU*2*pops[0]->size));
            }


          }
          else if (total_allele_count==(1.5*popsize(demographic_model,dem_uncert,gen-1,0))){
            int f_gen = int(-log( BRand::Controller.nextOpened())* (1./(femalemutU * 2 * pops[0]->size)));
            int m_gen = int(-log( BRand::Controller.nextOpened())* (1./(malemutU * male_pops[0]->size)));
            //cout<<f_gen<<" "<<m_gen<<"\n";
            gen-=min(f_gen, m_gen);
            if (m_gen<= f_gen){
              male_pops[0]->size=popsize(demographic_model,dem_uncert,gen-1,0);
              male_pops[0]->fix();
              male_pops[0]->mutatedown(poicond1(malemutU*male_pops[0]->size));
              //male_pops[0]->mutatedown(poicond1(malemutU*male_pops[0]->size));
            }
            else if (m_gen > f_gen){
              pops[0]->size=popsize(demographic_model,dem_uncert,gen-1,0);
              pops[0]->fix();
              pops[0]->mutatedown(poicond1(femalemutU*2*pops[0]->size));
              //pops[0]->mutatedown(poicond1(femalemutU*2*pops[0]->size));
            }
          }
          int popsize_male = int(popsize(demographic_model,dem_uncert,gen-1,0))/2;
          int popsize_female = max((popsize(demographic_model,dem_uncert,gen-1,0))-popsize_male,1);
          double f_A_female =  pops[0]->prob();
          double f_A_male = male_pops[0]->prob();
          double f_Y_male = 1.;
          //cout<<popsize_male<<"\t"<<popsize_female<<"\n";
          pops[0]->populate_from(f_A_female, f_A_male, popsize_female);
          male_pops[0]->populate_from(f_A_female, f_Y_male, popsize_male);
          //cout<<gen<<"\t"<<idx<<"\t"<<total_allele_count<<"\t"<<(pops[0]->size + male_pops[0]->size)<<"\n";
          gen--;
       }

       int total_allele_count = (pops[0]->allelenum()) + (male_pops[0]->allelenum());
       int total_pop_size = (pops[0]->size*2 + male_pops[0]->size);

       double raw_allele_freq = (total_allele_count*1.)/total_pop_size; // Get the derived allele frequency in the total population
       double sample_allele_count = rand_binom(pop_sample_size, raw_allele_freq); // Get an allele count with binomial sampling equal to the sample size in observed data

       //cout<<gen<<"\t"<<total_allele_count<<"\t"<<total_pop_size<<"\t"<<(total_allele_count*1.)/(1.5*total_pop_size)<<"\t"<<mutU<<"\t"<<femalemutU<<"\t"<<malemutU<<"\n";
       // Absolute distance function
       if (distance_function==1){
         if (abs(sample_allele_count - obs) <= epsilon){
           std::cout << "accept\n";
         }
         else{
           std::cout << "reject\n";
         }
       }

       // Binomial distance kernel
       if (distance_function==2){
         double obs_freq = (1.*obs)/pop_sample_size;
         double binom_prob = binomial_pmf(pop_sample_size, obs_freq, sample_allele_count);
         double rand_num = randnum(0, 1);
         if (rand_num < binom_prob){
           std::cout << "accept\n";
         }
         else{
           std::cout << "reject\n";
         }
       }
     }
     return EXIT_SUCCESS;
}

int poicond1(double l) //Poisson random variate conditional on the result being at least 1
{double expl=exp(-l);
double p=expl+BRand::Controller.nextClosed()*(1-expl);
int k=1;
while (p>=expl)
{k++;
p=p*BRand::Controller.nextClosed();
}
return k-1;
}


int popsize(int demographic_model, int dem_uncert, int gen, int i) //Calculates the size of population at generation gen according to Schiffles & Durbin's model
{
  if (demographic_model==1){
    if (gen<=T[idx]) idx+=1;
      if (dem_uncert==0){
        return Ne[idx];
      }
      else{
        pop_size = Ne[idx];
        if ((idx==55) && (randflag==0))
          {
          unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

          gent.seed(seed);
          std::uniform_real_distribution<double> distribution(log10(600000),log10(6000000));
          last_pop_size = pow(10,distribution(gent));
          pop_size = last_pop_size;

          randflag = 1;
          }
        else if ((idx>=55) && (randflag==1))
          {
            pop_size = last_pop_size;
          }
        //cout << gen << " " << idx << " " << Ne[idx] << " " << pop_size << "\n" ;
        return pop_size;
      }
  }
  else{
    return popNe;
  }
}

int poi(double l) //Regualr Poisson random variate
{double p=1,expl=exp(-l);
int k=0;
while (p>=expl)
{k++;
p=p*BRand::Controller.nextClosed();
}
return k-1;
}

double randnum (double a, double b) // Get a random uniform number between a and b
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  gent.seed(seed);
  std::uniform_real_distribution<double> distribution (a,b);
  return distribution(gent);
}

int rand_binom(int n, double p) // Random binomial variable
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  gent.seed(seed);
  std::binomial_distribution<> dist(n, p);
  return dist(gent);
}

// Calcualte the probability mass function of binomial distribution
double binomial_pmf(int N, double p, int x){
  using boost::math::binomial_distribution;
  binomial_distribution<>my_dist(N, p);
  double prob = pdf(my_dist, x);
  return prob;
}

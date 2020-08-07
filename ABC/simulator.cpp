//Modified by Eduardo Amorim (guerraamorim@gmail.com) from the original of Yuval Simons (Simons et al. 2014).
//Modified by Zach Fuller from Eduardo Amorim (Amorim et al. 2017) and Yuval Simons (Simons et al. 2014)

#include <chrono>
#include <random>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "population.h"
#include "BRand.hpp"
#include <stdint.h>
#include <map>
#include <time.h>
#include <sstream>
#include <fstream>
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
#include <boost/math/distributions/binomial.hpp>

#define SPLIT 1000
#define Nzero 10000
#define tau 2040
#define taujump 56000
//#define SEED 2710
#define TAU 100

#define ratio 1

boost::mt19937 gent;

double lognormal();
double randnum (double a, double b);
int rand_binom(int n, double p);

int nfinal,RUNS,Ntau=14448;
int randflag=0;

using namespace std;

struct freqs
{
  double freq0;
  double freq1;

  freqs(const double a=0,const double b=0) :
    freq0(a),freq1(b) {}
};

class valueComp
{
public:
  bool operator()(const freqs& A,
          const freqs& B)
    const
    { if (A.freq0!=B.freq0)
         return A.freq0<B.freq0;
      else
         return A.freq1<B.freq1; }
};

char* filename(char *);
char* seriesfilename(char * str);
int popsize(int demographic_model, int dem_uncert, int gen, int i);
int poi(double l);
int poicond1(double l);
int boost_poi(double l);
int boost_poicond1(double l);
double binomial_pmf(int N, double p, int x);
long long int NcR(int n, int r);

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

double lognormal(double mutU);

//

void Print(const vector<int>& v);

void Print (const vector<int>& v){
  //vector<int> v;
  for (int i=0; i<v.size();i++){
    cout << v[i] << endl;
  }
}

int Ne[56]= {14448,14068,14068,14464,14464,15208,15208,16256,16256,17618,17618,19347,19347,21534,21534,24236,24236,27367,27367,30416,30416,32060,32060,31284,29404,26686,23261,18990,16490,16490,12958,12958,9827,9827,7477,7477,5791,5791,4670,4670,3841,3841,3372,3372,3287,3359,3570,4095,4713,5661,7540,11375,14310,13292,14522,613285};
int T[56]={55940,51395,47457,43984,40877,38067,35501,33141,30956,28922,27018,25231,23545,21951,20439,19000,17628,16318,15063,13859,12702,11590,10517,9482,8483,7516,6580,5672,5520,5156,4817,4500,4203,3922,3656,3404,3165,2936,2718,2509,2308,2116,1930,1752,1579,1413,1252,1096,945,798,656,517,383,252,124,0};

int main(int argc, char *argv[])
{

    int initgen=150000,stopover;

    // Show options and order of parameters if no arguments are given
    if (argc == 1)
    {
      std::cerr << "Usage: " << argv[0] << " U obs dom runs mut_uncert dem_uncert demographic_model popNe process_seed distance_function sample_size\n"
          "Simulate an allele frequency "
          "on the autosomes\n"
          "The selection coefficient and epsilon (if needed) "
          "are read from stdin\n";
      return 1;
    }
    if (argc==12) // User may run the script with "./a.out RUNS sel DOM", where RUNS = number of runs, sel = selective coefficient (in this work always set to 1), and DOM = heterozygote effect (in this work set to 0 or 1%)
       {
         Uvar=atof(argv[1]);
         obs=atof(argv[2]);
         dom=atof(argv[3]);
         RUNS=atof(argv[4]);
         mut_uncert=atof(argv[5]);
         dem_uncert=atof(argv[6]);
         //Input 1 if Schiffles-Durbin demographic model should be used. Else, constant size is used
         demographic_model=atof(argv[7]);
         //Only useful if constant size population model is used
         popNe=atof(argv[8]);
         pid_seed=atof(argv[9]);
         distance_function=atof(argv[10]);
         pop_sample_size=atof(argv[11]);

       }
    else
        {
            cout<<"Error: Not enough parameters enetered" << endl;
        }

    // For use in pakman, epsilon followed by the selection coefficent must be read from stdin
    std::cin >> epsilon;
    std::cin >> sel;

    if (std::cin.fail())
    {
        std::cerr << "Error: could not read epsilon or selection coefficient from stdin\n";
        return 1;
    }

    // If running multiple threads or processors, pid_seed allows for different seeds to run instead of using the same CPU clock time
    gent.seed(time(NULL)+pid_seed);
    BRand::Controller.seed(time(NULL)+pid_seed);

    stopover=RUNS;

    std::uniform_real_distribution<double> mt_rand{0.0, 1.0};

    population::initialize(sel,dom);
    population* pops[2];
    pops[0]= new population(popsize(demographic_model,dem_uncert,0,0));
    pops[1]= new population(popsize(demographic_model,dem_uncert,0,0));

    int euroflag=0;


    int gen,thispop;

    time_t tt;
    struct tm *tim;

    double totfreq;
    int age=0;
    int tauallele,before, after;
    freqs f;
    ofstream myfile,results,afrDistrib,eurDistrib,mutDistribution;

    for(int run=0;run<RUNS;run++){
            //If uncertainity is included for the mutation rate, draw from a normal distribution with mean=U, and SD=1/10 * U
            if (mut_uncert==1){
              std::normal_distribution<double> norm_distribution(Uvar, Uvar/10);
              mutU = norm_distribution(gent);
            }
            else{
              mutU =  Uvar;
            }
            {
            randflag = 0;
            if (demographic_model==1){
              //Burn-in period, only used under Schiffels-Durbin model
              initgen=150000;
            }
            else{
              //Burn-in period, only used under constant population size model, equal to 10*Ne
              initgen=popNe*10;
            }
            gen=initgen;
            idx=0;
            double halfthetaVAR = 2*popsize(demographic_model,dem_uncert,initgen,0)*mutU;
            double baseupVAR = (1./halfthetaVAR);
            double basedownVAR = (1./halfthetaVAR);

            pops[0]->size=popsize(demographic_model,dem_uncert,gen,0); //Initializing the population size for the beginning of each run with deleterious allele absent from the population

                pops[0]->clear();

            euroflag=0; //A flag indicating if the African-European population split has occured (0 = it hasn't occurred)

            while (gen>1)
                    {
          if (euroflag==1) //Keep the derived allele derived
          {totfreq=0.5*(pops[0]->freq()+pops[1]->freq());
        if (totfreq==1) {pops[0]->clear();pops[1]->clear();totfreq=0.0;}
      }
          else
      {totfreq=pops[0]->freq();
        if (totfreq==1) {pops[0]->clear();totfreq=0.0;}
      }
          int i=0;
		      if (((pops[i]->allelenum()>0)&&(pops[i]->allelenum()<(2*pops[i]->size)))||(gen<=taujump)) //If the allele is segregating or we are in recent history (last 5920 generations), introduce mutations
                        {
                            pops[i]->mutateup(boost_poi( mutU*(2*pops[i]->size-pops[i]->allelenum()) )) ;
                            pops[i]->mutatedown(boost_poi( mutU*pops[i]->allelenum())) ;
                        }
		      else if (pops[i]->allelenum()==0) //If the ancestral allele is fixed calculate when will the next deleterious allele appear and advance population state to that point in time (unless that point in time is in recent history; then just set time to initial history).
                        {
                             gen-=int(-log( mt_rand(gent))*baseupVAR);
                             if (gen<=taujump)
                                gen=taujump+1;
                             else
                                 {pops[i]->size=popsize(demographic_model,dem_uncert,gen,0);
                                 pops[i]->clear();
                                 pops[i]->mutateup(poicond1(2*mutU*pops[i]->size)); }

                        }
                        else if (pops[i]->allelenum()==(2*pops[i]->size))//If the deleterious allele is fixed calculate when will the next deleterious allele appear and advance population state to that point in time (unless that point in time is in recent history and then just set time to initial history). Note that with s = 1 this cannot happen
                        {
                             gen-=int(-log(mt_rand(gent))*basedownVAR);
                             if (gen<=taujump)
                                gen=taujump+1;
                             else
                                 {pops[i]->size=popsize(demographic_model,dem_uncert,gen,0);
                                 pops[i]->fix();
                                 pops[i]->mutatedown(poicond1(2*mutU*pops[i]->size)); }
                        }
                           pops[i]->populate_from(pops[i]->prob(),popsize(demographic_model,dem_uncert,gen-1,0));

		      gen--; //Next generation please
        }
}

// Get the raw allele frequency and then sample allele counts from a binomial distribution with the number of chromosomes equal to the sample size
double raw_allele_freq = pops[0]->freq();
double sample_allele_count = rand_binom(pop_sample_size, raw_allele_freq);

// Distance function for absolute distnance
if (distance_function==1){
  if (abs(sample_allele_count - obs) <= epsilon){
    std::cout << "accept\n";
  }
  else{
    std::cout << "reject\n";
  }
}

// Distance function for binomial distance kernel (epsilon is not used, accepted with binomial probability)
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

        return pop_size;
      }
  }
  else{
    return popNe;
  }

}

double randnum (double a, double b) // Get a random uniform number between a and b
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  gent.seed(seed);
  std::uniform_real_distribution<double> distribution (a,b);
  return distribution(gent);
}

double lognormal(double mutU)
{
  mutU = log10(mutU) - (0.3249/2) * log(10);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  gent.seed(seed);

  std::normal_distribution<double> distribution(mutU,0.57);

  double u = distribution(gent);
  return pow(10,u);
}

int boost_poi(double l) //Regular Poisson random variate

{
    if (l==0.0){
      return 0.0;
    }

    else{
      boost::random::poisson_distribution<> dist(l);
      return dist(gent);
    }

}

int boost_poicond1(double l) //Regular Poisson random variate

{
    if (l==0.0){
      return 1.0;
    }

    else{
      int res=0;
      while (res<1){
          boost::random::poisson_distribution<> dist(l);
          res = dist(gent);
      }
      return(res);
    }

}

int poicond1(double l) //Poisson random variate conditional on the result being at least 1
{

double expl=exp(-l);
double p=expl+BRand::Controller.nextClosed()*(1-expl);
int k=1;
while (p>=expl)
{k++;
p=p*BRand::Controller.nextClosed();
}
return k-1;

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

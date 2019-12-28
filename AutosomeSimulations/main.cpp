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
#define ln10 2.30258509
#define SEL 0.0005
#define SPLIT 1000
#define Nzero 10000
#define tau -1
#define taujump 56000
//#define SEED 2710
#define TAU 100
#define Urate 0.0000000125//0.00000025
#define ratio 1

int nfinal,RUNS,Ntau=14448;
int idx=0;
double sel;

int Ne[56]= {14448,14068,14068,14464,14464,15208,15208,16256,16256,17618,17618,19347,19347,21534,21534,24236,24236,27367,27367,30416,30416,32060,32060,31284,29404,26686,23261,18990,16490,16490,12958,12958,9827,9827,7477,7477,5791,5791,4670,4670,3841,3841,3372,3372,3287,3359,3570,4095,4713,5661,7540,11375,14310,13292,14522,613285};
int T[56]={55940,51395,47457,43984,40877,38067,35501,33141,30956,28922,27018,25231,23545,21951,20439,19000,17628,16318,15063,13859,12702,11590,10517,9482,8483,7516,6580,5672,5520,5156,4817,4500,4203,3922,3656,3404,3165,2936,2718,2509,2308,2116,1930,1752,1579,1413,1252,1096,945,798,656,517,383,252,124,0};
 
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
int popsize(int gen, int i);
int poi(double l);
int poicond1(double l);
int binom(int n, double p);
int binom1(int n, double p);


int main(int argc, char *argv[])
{
    int initgen=150000,stopover;
    double m=0.00015;    
    if (argc==3)
       {RUNS=atoi(argv[1]);
       sel=atof(argv[2]);
       }
    else    
        {
        cout<<"enter RUNS:";
        cin>>RUNS;
        cout<<"enter sel:";
        cin>>sel;
        }
    stopover=RUNS;
    cout<<"RUNS:"<<RUNS<<" tau:"<<tau<<" sel:"<<sel<<" initgen:"
    <<initgen<<" stopover:"<<stopover<<" m:"<<m<<"\n";
 
    BRand::Controller.seed(time(NULL));
    
    population::initialize(sel);
    population* pops[2];
    pops[0]= new population(popsize(0,0)); 
    pops[1]= new population(popsize(0,0));

    int euroflag=0;
    double halftheta=2*popsize(initgen,0)*Urate; 



    //map<int,int> count[2],counts[2][2],countd[2],taucount,taucountd;
    //map<freqs,int,valueComp> joint,joints[2],jointd;

    int gen,thispop;
    double baseup=(1./halftheta);
    double basedown=(1./(ratio*halftheta));
    double oneovertwoU=(1./(2*Urate));
    double oneovertwoUratio=(1./(2*Urate*ratio));
    int skip=0,stretch;
    time_t tt;
    struct tm *tim;
    char library[100],comm[200],totalfile[200],totalres[200];
    char time_stamp [80];
    double Ep,s1,S2N,Ep2,Ex,Exderived,s12,Ex2,Ex4,Dx,Dx2;
    double europ, euros1, eurox, euroS2N;
    double taufreq,totfreq;
    int mutnum;
    int age=0;
    int tauallele,before, after;
    freqs f;
    
    
    ofstream myfile;
    
    time(&tt); //We use a time stamp to differentiate files in parallel computing
    tim=localtime(&tt);
    sprintf(library, "infSchiffles/regular");
    sprintf(comm, "mkdir %s",library);
    system(comm);
    
    sprintf(library, "%s/res%d_%d_%d",library, tim->tm_mday,tim->tm_mon+1,tim->tm_year+1900);
    sprintf(comm, "mkdir %s",library);
    system(comm);

    		//Open results file
    sprintf(totalfile, "%s/Schifflesdata_sel%f_runs%d_%d_%d_%d.csv", library,log10(sel),RUNS, tim->tm_hour,tim->tm_min,tim->tm_sec);
    cout<<"opening "<<totalfile<<"\n";
    myfile.open (totalfile,std::fstream::out | std::fstream::app); //Data
    cout<<"opened "<<totalfile<<"\n";
    
    for(int run=0;run<RUNS;run++)
            {
            gen=initgen;
            pops[0]->size=popsize(gen,0);//Initializing the population size for the beginning of each run
	    pops[0]->clear();

            
            euroflag=0; //A flag indicating if a Af-Eur population split has occured

            while (gen>0)
                    {
		      if (gen==tau) //tau is the generation of the out of Africa exodus.
                      {
                        euroflag=1;
                        *pops[1]=*pops[0]; //Creating the European population from the African population.
                        taufreq=pops[0]->freq(); //Keeping track of deleterious allele frequencies at the split for calculation of the change in frequency since the split
                      }
		      if (gen==920) m=0.000025; //Reduced migration period
		      
		      if (euroflag==1) //Keep the derived allele derived
			{totfreq=0.5*(pops[0]->freq()+pops[1]->freq());
			  if (totfreq==1) {pops[0]->clear();pops[1]->clear();totfreq=0.0;}
			}
		      else
			{totfreq=pops[0]->freq();
			  if (totfreq==1) {pops[0]->clear();totfreq=0.0;}
			}
		     
		      if (totfreq==0.0) age=gen; //keep track of allele age


		      for(int i=0;i<(1+euroflag);i++) //For each population
                      {

		      if (((pops[i]->allelenum()>0)&&(pops[i]->allelenum()<(2*pops[i]->size)))||(gen<=taujump)) //If the allele is segregating or we are in recent history (last 5920 gens)
                        {    //Introduce mutations
                            pops[i]->mutateup(poi( Urate*(2*pops[i]->size-pops[i]->allelenum()) )) ; 
                            pops[i]->mutatedown(poi(Urate*ratio*pops[i]->allelenum()));
                        }
		      else if ((pops[i]->allelenum()==0)||(pops[i]->allelenum()==(2*pops[i]->size))) //If the allele is fixed calculate when will the next  allele appear and advance population state to that point in time (unless that point in time is in recent history and then just set time to initial history).
                        {
                             
			     pops[i]->clear();
                             gen-=int(-log( BRand::Controller.nextOpened())*baseup);
                             if (gen<=taujump) 
                                gen=taujump+1;
                             else
			       { age=gen;
				 pops[i]->size=popsize(gen,i);
                                 pops[i]->clear();
                                 pops[i]->mutateup(poicond1(2*Urate*pops[i]->size)); }

                        }    
                      
		      if (euroflag==0) //If a split between Af and Eur populations has occured apply migration
                           pops[i]->populate_from(pops[i]->prob(),popsize(gen-1,i));
                        else
                           pops[i]->populate_from((1-m)*(pops[i]->prob())+m*(pops[1-i]->prob()),popsize(gen-1,i)); 

                    }
                    
                                    
		      gen--; //Next generation please
                   }

	    if ((pops[1]->freq()+pops[0]->freq())!=0) myfile<<age<<","<<pops[0]->freq()<<","<<0<<","<<0<<"\n"; //write results of run
            
            if (((100*(run+1))%RUNS)==0) cout<<((100.0*(run+1))/RUNS)<<"%\n";              

    }

    delete pops[0];
    delete pops[1];
    cout<<"deleted pops\n";
    
    cout<<"file closing\n";
    myfile.close();
    cout<<"file closed\n";

    return EXIT_SUCCESS;
}


int popsize(int gen, int i) //Calculates the size of population i (0=African, 1=European) at generation gen according to Schiffles & Durbin's model
{

  if (gen<=T[idx]) idx-=1;
  return Ne[idx];



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

int binom(int n, double p) //binomial random variate generator
{ if (p<0.5) return binom1(n,p);
else return  n-binom1(n,1-p);
}



int binom1(int n, double p) //binomial random variate generator
{int x=0;
double y=0,c=log(1-p);
if (c==0)
   {//std::cout<<"c is zero!\n";
   return 0;}
while (1) 
    {
    y+=ceil(log(BRand::Controller.nextOpened())/c);
    if (y>n) return x;    
    x++;    
    } 
}

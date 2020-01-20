
#ifndef POPULATION_H
#define POPULATION_H
#include <vector>

class population
{
	public:
		// class constructor
		population(int N);

		static void initialize(double sel, double dom);//set selection and dominance coefficients
		void populate_from(population &p,double male_prob, int N=0);//create next generation
		void populate_from(double female_prob, double male_prob, int N=0);//create next generation
		int choose_allele();//choose an alleleic state (homo beneficial, hetero or homo deleterious) from the population. I think this function isn't actually used anywhere anymore.
		int binom(int n, double p); //binomial distribution
		int alleleholders[3]; //vector of the three classes of the pop: homo beneficial, hetero and homo deleterious

		~population();  // class destructor

		double freq();//reutrn frequency
		double prob();//Expected frequency of deleterious allele post selection
		void mutateup(int n);// n deleterious mutations occur
		void mutatedown(int n);// n beneficial mutations occur
		void clear(); //Fix for beneficial
		void fix(); //Fix for deleterious
		int size;
		int allelenum(); //number of copies of  deleterious allel
		population& operator= (const population &Source); //copy one pop into the other. used to create Europeans from Africans


	private:

            static double s;
            static double h;
            static double hs;


};

class population_male
{
	public:
		// class constructor
		population_male(int N);

		static void initialize(double sel, double dom);//set selection and dominance coefficients
		void populate_from(population_male &p,double male_prob, int N=0);//create next generation
		void populate_from(double female_prob, double male_prob, int N=0);//create next generation
		int choose_allele();//choose an alleleic state (homo beneficial, hetero or homo deleterious) from the population. I think this function isn't actually used anywhere anymore.
		int binom(int n, double p); //binomial distribution
		int alleleholders[2]; //vector of the three classes of the pop: homo beneficial, hetero and homo deleterious

		~population_male();  // class destructor

		double freq();//reutrn frequency
		double prob();//Expected frequency of deleterious allele post selection
		void mutateup(int n);// n deleterious mutations occur
		void mutatedown(int n);
		void clear(); //Fix for beneficial
		void fix(); //Fix for deleterious
		int size;
		int allelenum(); //number of copies of  deleterious allel



	private:

            static double s;
            static double h;
            static double hs;


};

#endif // POPULATION_H

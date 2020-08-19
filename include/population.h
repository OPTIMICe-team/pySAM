#ifndef POPULATION_H
#define POPULATION_H

#include "aggregate.h"
#include "math_lib.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>

using namespace boost::random;
//using namespace std;

class population
{
	public:
		/** Default constructor */
		population();
		/** Default destructor */
		~population();
		population(int);

		vector<aggregate> aggregates;
		//uniform_real_distribution<double> uniform;
		//mt19937 gen;
		uniform_real_distribution<float> uniform;
		mt11213b gen;

		int N;
		int CAP_idx;
		int NAP_idx;

		void fill_exp(int,int,exponential_distribution<double>&);
		void fill_gamma(int,int,gamma_distribution<double>&);
		void fill_uniform(int,int,uniform_real_distribution<double>&);
		void get_collision_couple();
		void make_collision(double);
		
		vector<double> get_population_parameters();
		vector<vector<double> > get_hystogram(unsigned int);
		vector<vector<double> > get_hystogram(double);
		void print_size_mass(string);

	protected:
	private:
};

#endif // POPULATION_H

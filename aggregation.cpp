/*
This is an evolution of the collection program used for Ori et al (2014) snowflake generation
It implements a collection kernel similar to the works of Westbrook or Leinonen
It still does not have a nice runtime interface, so you will have to recompile after every modification
*/

#include <iostream>
#include <time.h>

#include "aggregate.h"
#include "plate.h"
#include "column.h"
#include "population.h"
#include "math_lib.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>

using namespace boost::random;
using namespace std;

int main()
{
	uniform_real_distribution<double> angle(0.0,2*M_PI);
	uniform_real_distribution<double> cart(-10.0,10.0);
	uniform_real_distribution<double> distro_uni(50,150);
	exponential_distribution<double> distro_exp(1./200.); // la media è 1/k
	gamma_distribution<double> distro_gamma(20,20); // 10 20

	population pop(2000);
	//pop.fill_exp(20000,0,distro_exp);
	pop.fill_gamma(1000,0,distro_gamma);
	pop.fill_gamma(1000,1,distro_gamma);
	//pop.fill_gamma(500,2,distro_gamma);
	//pop.fill_gamma(2000,3,distro_gamma);
	//pop.fill_uniform(20000,0,distro_uni);

    double d=20.; //interdipole spacing
    unsigned int j=500;
    vector<double> par;
    vector<vector<double> > hystogram;
	time_t cycle_start,cycle_end;
	time(&cycle_start);
    for(int i=0;i<2002;i++)
    {
    	j++;
        pop.make_collision(d);
        cout<<j<<" "<<flush;
        if(j>=499)
        {
		time_t start, end;
		time(&start);
        	par=pop.get_population_parameters();
		time(&end);
		double seconds=difftime(end, start);
        	cout<<"Distro gamma k= "<<par[0]<<"\t theta= "<<par[1]<<" ho impiegato "<<seconds<<endl;
        	double bin_size=100.;
        	hystogram=pop.get_hystogram(bin_size);
        	for(int i=0;i<2;i++)
        	{
        		for(int j=0;j<hystogram[i].size();j++)
        		{
        			cout<<hystogram[i][j]<<"\t";
        		}
        		cout<<endl;
        	}
        	//pop.aggregates[0].save_vtk_mesh();
        	j=0;
		time(&cycle_end);
		seconds=difftime(cycle_end,cycle_start);
		cout<<"Per questo ciclo ho impiegato "<<seconds<<" secondi"<<endl;
		time(&cycle_start);
        }
    }
    return 0;
}

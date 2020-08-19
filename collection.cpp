#include <iostream>
#include <algorithm>
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

// Distro parametrizations
double theta = 1.0/0.07;
double kappa(double Dmax) {
	// Assumes Dmax in um
	return Dmax*0.07*0.07/1.0;
}


int main() {
	time_t begin;
	time(&begin);
	mt11213b gen;
	uniform_real_distribution<float> uniform;

	// Create a seed aggregate
	unsigned int type = 0; // 0 for columns
	double Dmin = 100.0;
	double D = 100.0;      // initial column length
	double Dmax = 2000.0;
	unsigned int N = 100;  // number of collisions
	double d = 20.0; // dipole resolution in microns
	double target = 25000;

	// Initialise computations
	gen.seed(time(NULL));
	aggregate CAP(type, D);

	//for(int i=0; i<N; i++) {
	while(CAP.Maximum_dimension<target) {
		gamma_distribution<double> distro_gamma(kappa(CAP.Maximum_dimension), theta);
		D = distro_gamma(gen);
		D = max(D, Dmin);
		D = min(D, Dmax);
		cout<<D<<"\t"<<CAP.Maximum_dimension<<endl;
		aggregate NAP(type, D);

		// Rotate randomly
		CAP.rotate(10*uniform(gen), 10*uniform(gen), 10*uniform(gen), 2*M_PI*uniform(gen));
		NAP.rotate(10*uniform(gen), 10*uniform(gen), 10*uniform(gen), 2*M_PI*uniform(gen));
		CAP.rotate(10*uniform(gen), 10*uniform(gen), 10*uniform(gen), 2*M_PI*uniform(gen));
		NAP.rotate(10*uniform(gen), 10*uniform(gen), 10*uniform(gen), 2*M_PI*uniform(gen));
		CAP.rotate(10*uniform(gen), 10*uniform(gen), 10*uniform(gen), 2*M_PI*uniform(gen));
		NAP.rotate(10*uniform(gen), 10*uniform(gen), 10*uniform(gen), 2*M_PI*uniform(gen));

		// Displace the aggregates
		double Rmax_sum = CAP.Rmax+NAP.Rmax;
		double dx = 2*(uniform(gen)-0.5)*Rmax_sum;
		double dy = 2*(uniform(gen)-0.5)*sqrt(Rmax_sum*Rmax_sum-dx*dx); // dx<Rmax_sum altrimenti è un bordello
		NAP.translate(dx, dy, 2*Rmax_sum);

		time_t timer_start;
		time_t timer_stop;
		double seconds;
		time(&timer_start);

		bool collision=CAP.collide(NAP); // if it does not collide try again
		while(!collision) {
			NAP.translate(-dx,-dy,0.0);
			dx=2*(uniform(gen)-0.5)*Rmax_sum;
			dy=2*(uniform(gen)-0.5)*sqrt(Rmax_sum*Rmax_sum-dx*dx);
			NAP.translate(dx,dy,0.0);
			collision=CAP.collide(NAP);
		};

		CAP.N+=NAP.N; // update number of pristines (always 1 in this scenario)
		CAP.crystal.reserve(CAP.N); // actually adds the shape information of NAP to CAP
		for(int i=0; i<NAP.N; i++) {
			CAP.crystal.push_back(NAP.crystal[i]);
		}

		time(&timer_stop);
		seconds=difftime(timer_stop, timer_start);
		cout<<"collided in "<<seconds<<endl;

		CAP.update_mic();
		CAP.update_maximum_dimension();

		CAP.save_aggregate(d);
		string prefix = "vtk/";
		CAP.save_vtk_mesh(prefix);
		time(&timer_stop);
		seconds=difftime(timer_stop, begin);
		cout<<"After "<<seconds<<" seconds CAP size is "<<CAP.Maximum_dimension/1000.0<<" millimeters"<<endl;
	}

}
#include "population.h"
#include <omp.h>

population::population() {
	//ctor
}

population::~population() {
	//dtor
}

population::population(int NN) {
	N=NN;
	aggregates.reserve(N);
	gen.seed(time(NULL));
}

void population::fill_exp(int n, int type, exponential_distribution<double>& expo) {
	time_t timer_start;
	time_t timer_stop;
	double seconds;
	time(&timer_start);

	double D=0.0;
	#pragma omp parallel for
	for(int i=0;i<n;i++) {
		while((D<100)||(D>1000)) {
			D=expo(gen);
		}
		aggregate temp(type, D);
		#pragma omp critical
		{
			aggregates.push_back(temp);
		}
		D=0;
		cout<<i<<" "<<flush;
	}

	time(&timer_stop);
	seconds=difftime(timer_stop,timer_start);
	cout<<endl<<"To fill the population it took  "<<seconds<<" seconds"<<endl;
}

void population::fill_gamma(int n, int type, gamma_distribution<double>& gamma) {
	time_t timer_start;
	time_t timer_stop;
	double seconds;
	time(&timer_start);

	double D=0.0;
	#pragma omp parallel for
	for(int i=0;i<n;i++) {
		while((D<40)||(D>2000)) {
			D=gamma(gen);
		}
		aggregate temp(type, D);
		#pragma omp critical
		{
			aggregates.push_back(temp);
		}
		D=0;
		cout<<i<<" "<<flush;
	}

	time(&timer_stop);
	seconds=difftime(timer_stop,timer_start);
	cout<<endl<<"To fill the population it took  "<<seconds<<" seconds"<<endl;
}

void population::fill_uniform(int n, int type, uniform_real_distribution<double>& uniform) {
	time_t timer_start;
	time_t timer_stop;
	double seconds;
	time(&timer_start);

	double D=0.0;
	#pragma omp parallel for
	for(int i=0;i<n;i++) {
		D=uniform(gen);
		aggregate temp(type, D);
		#pragma omp critical
		{
			aggregates.push_back(temp);
		}
		D=0;
	}
	time(&timer_stop);
	seconds=difftime(timer_stop,timer_start);
	cout<<endl<<"To fill the population it took  "<<seconds<<" seconds"<<endl;
}

void population::get_collision_couple() {
	CAP_idx=0;
	NAP_idx=0;
	double total_probability=0.;
	#pragma omp parallel for reduction(+:total_probability)
	for(int i=0;i<N-1;i++) {
		for(int j=i+1;j<N;j++) {
			total_probability+=(aggregates[i].Rmax+aggregates[j].Rmax)*(aggregates[i].Rmax+aggregates[j].Rmax)
			*abs((aggregates[i].v-aggregates[j].v));
		}
	}
	
	//double randomness=0.2;
	//total_probability/=(1.-randomness);
	//double number_density=1./(N*(N-1.));
	//double random_probability=randomness*total_probability*number_density;
	double probability=total_probability*uniform(gen);
	total_probability=0;

	for(int i=0;(i<N-1)&&(CAP_idx==NAP_idx);i++) {
		for(int j=i+1;(j<N)&&(CAP_idx==NAP_idx);j++) {
			total_probability+=(aggregates[i].Rmax+aggregates[j].Rmax)*(aggregates[i].Rmax+aggregates[j].Rmax)
			*abs((aggregates[i].v-aggregates[j].v));//+random_probability*uniform(gen);
			if(total_probability>probability) {
				CAP_idx=i;
				NAP_idx=j;
			}
		}
	}
}

void population::make_collision(double d) {
	time_t timer_start;
	time_t timer_stop;
	double seconds;
	time(&timer_start);
	get_collision_couple();
	time(&timer_stop);
	seconds=difftime(timer_stop,timer_start);

	aggregates[CAP_idx].rotate(10*uniform(gen),10*uniform(gen),10*uniform(gen),2*M_PI*uniform(gen));
	aggregates[NAP_idx].rotate(10*uniform(gen),10*uniform(gen),10*uniform(gen),2*M_PI*uniform(gen));

	double Rmax_sum=aggregates[CAP_idx].Rmax+aggregates[NAP_idx].Rmax;
	double dx=2*(uniform(gen)-0.5)*Rmax_sum;
	double dy=2*(uniform(gen)-0.5)*sqrt(Rmax_sum*Rmax_sum-dx*dx); // dx<Rmax_sum altrimenti è un bordello
	aggregates[NAP_idx].translate(dx,dy,2*Rmax_sum);

	time(&timer_start);

	bool collision=aggregates[CAP_idx].collide(aggregates[NAP_idx]);
	while(!collision) {
		aggregates[NAP_idx].translate(-dx,-dy,0.0);
		dx=2*(uniform(gen)-0.5)*Rmax_sum;
		dy=2*(uniform(gen)-0.5)*sqrt(Rmax_sum*Rmax_sum-dx*dx);
		aggregates[NAP_idx].translate(dx,dy,0.0);
		collision=aggregates[CAP_idx].collide(aggregates[NAP_idx]);
	};
	aggregates[CAP_idx].N+=aggregates[NAP_idx].N;
	aggregates[CAP_idx].crystal.reserve(aggregates[CAP_idx].N);
	for(int i=0;i<aggregates[NAP_idx].N;i++) {
		aggregates[CAP_idx].crystal.push_back(aggregates[NAP_idx].crystal[i]);
	}

	time(&timer_stop);
	seconds=difftime(timer_stop,timer_start);

	aggregates[CAP_idx].update_mic();
	aggregates[CAP_idx].update_maximum_dimension();
	//aggregates[CAP_idx].save_aggregate(d);
	aggregates.erase(aggregates.begin()+NAP_idx);
	N--;
}

vector<double> population::get_population_parameters() {
	double mean_x=0;
	double s=0;
	for(int i=0;i<N;i++) {
		mean_x+=aggregates[i].Maximum_dimension;
		s-=log(aggregates[i].Maximum_dimension);
	}
	mean_x/=N;
	s/=N;
	s+=log(mean_x);
	cout<<"La media di max_dim è "<<mean_x<<endl;	
	vector<double> parameters;
	parameters.push_back((3.-s+sqrt((s-3.)*(s-3.)+24.*s))/(12*s));
	parameters.push_back(mean_x/parameters[0]);
	cout<<parameters[0]<<"\t"<<parameters[1]<<endl;
	double k_new=parameters[0]-(log(parameters[0])-digammaR(parameters[0])-s)/((1./parameters[0])-trigammaR(parameters[0]));
	cout<<"k_new "<<k_new<<endl;
	double err=abs((k_new-parameters[0])/parameters[0]);
	double toll=1.e-5;
	parameters[0]=k_new;
	while(err<toll) {
		k_new=parameters[0]-(log(parameters[0])-digammaR(parameters[0])-s)/((1./parameters[0])-trigammaR(parameters[0]));
	}
	parameters[1]=mean_x/parameters[0];
	return parameters;
}

vector<vector<double> > population::get_hystogram(unsigned int nbins) {
	vector<vector<double> > hyst;
	hyst.resize(2);
	hyst[0].resize(nbins);
	hyst[1].resize(nbins);
	double dmax=0.;
	double dmin=99e99;
	for(int i=0;i<N;i++) {
		if(aggregates[i].Maximum_dimension>dmax) {
			dmax=aggregates[i].Maximum_dimension;
		}
		else if(aggregates[i].Maximum_dimension<dmin) {
			dmin=aggregates[i].Maximum_dimension;
		}
	}
	double bin_size=(dmax-dmin)/((double)nbins);
	for(int i=0;i<nbins;i++) {
		hyst[0][i]=dmin+i*bin_size;
		hyst[1][i]=0.;
	}
	unsigned int index=0;
	for(int i=0;i<N;i++) {
		index=floor((aggregates[i].Maximum_dimension-dmin)/bin_size);
		hyst[1][index]++;
	}
	return hyst;
}

vector<vector<double> > population::get_hystogram(double bin_size) {
	vector<vector<double> > hyst;
	hyst.resize(2);
	double dmax=0.;
	double dmin=0.;
	for(int i=0;i<N;i++) {
		if(aggregates[i].Maximum_dimension>dmax) {
			dmax=aggregates[i].Maximum_dimension;
		}
	}
	unsigned int nbins=ceil((dmax-dmin)/((double)bin_size));
	hyst[0].resize(nbins);
	hyst[1].resize(nbins);
	for(int i=0;i<nbins;i++) {
		hyst[0][i]=dmin+i*bin_size;
		hyst[1][i]=0;
	}
	unsigned int index=0;
	ostringstream counter,bin;
	counter<<N+1;
	string pre=counter.str();
	string prefix;
	for(int i=0;i<N;i++) {
		index=floor((aggregates[i].Maximum_dimension-dmin)/bin_size);
		if(hyst[1][index]<0.90) {
			prefix=pre;
			prefix.append("_");
			bin.str("");
			bin.clear();
			bin<<index;
			prefix.append(bin.str());
			///
			cout<<"Salvo bin"<<index<<" con prefisso"<<prefix<<endl;
			aggregates[i].save_vtk_mesh(prefix); //
			///
		}
		hyst[1][index]++;
	}
	///
	prefix=pre;
	prefix.append("_");
	print_size_mass(prefix);
	///
	return hyst;
}

void population::print_size_mass(string prefix) {
	ostringstream N_aggregate;
	N_aggregate<<N;
	prefix.append("_");
	prefix.append(N_aggregate.str());
	prefix.append("Mass_Size.dat");
	ofstream OUT(prefix.c_str());
	for(int i=0; i<N; i++) {
		OUT<<aggregates[i].Maximum_dimension/1000.<<"\t"<<aggregates[i].mass*1000.<<endl;
	}
	OUT.close();
}

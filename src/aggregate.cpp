#include "aggregate.h"
//#include <omp.h>

aggregate::aggregate() {
	//default ctor do nothing
}

aggregate::aggregate(int type, double D) {
	N=1;
	CM=Vector3d::Zero();
	crystal.resize(1);
	if(type==0) {
		crystal[0]=new column(D);
	}
	else if(type==1) {
		crystal[0]=new plate(D);
	}
	else if(type==2) {
		crystal[0]=new rosette6b(D);
	}
	else if(type==3) {
		crystal[0]=new dendrite(D);
	}
	else {
		cout<<endl<<"aggregate:I don't know the pristine type "<<type<<" please check your code"<<endl;
		exit(0);
	}
	Rgyr=crystal[0]->Rgyr;
	set_Rmax();
	Maximum_dimension=D;//2*Rmax;
	set_mass();
	set_v();
}

aggregate::~aggregate() {
	//dtor
}

void aggregate::set_mass() {
	double aggr_vol=0;
	for(int i=0;i<N;i++)
	{
		aggr_vol+=crystal[i]->vol;   // micrometers
	}
	mass=aggr_vol*rho_ice;          // kg
}

void aggregate::set_v() {
	v=(ni/(Rgyr/1000000))*pow((mass*g)/(rho_air*ni*ni),0.7);    /// I decided to set alpha=0.7
}

void aggregate::set_Rmax() { // we can speed this up with convex hull if total number of vertexes goes up
	double squared_distance=0;
	double provv;
	for(int i=0; i<N; i++) {
		for(int j=0; j<crystal[i]->N_vertexes; j++) {
			provv=crystal[i]->vertex.col(j).dot(crystal[i]->vertex.col(j));
			if(provv>squared_distance) {
				squared_distance=provv;
			}
		}
	}
	Rmax=sqrt(squared_distance);
}

void aggregate::set_Rgyr() {
	/// Suppose that the CM of the aggregate is in the origin and the mass is updated
	double rad=0;
	for(int i=0;i<N;i++) {
		rad+=crystal[i]->vol*(crystal[i]->Rgyr*crystal[i]->Rgyr+crystal[i]->CM.dot(crystal[i]->CM));
	}
	Rgyr=sqrt(rho_ice*rad/mass);
}

void aggregate::set_CM() {
	double x=0.,y=0.,z=0.;
	double vol=0.;
	for(int i=0;i<N;i++) {
		x+=crystal[i]->CM(0)*crystal[i]->vol;
		y+=crystal[i]->CM(1)*crystal[i]->vol;
		z+=crystal[i]->CM(2)*crystal[i]->vol;
		vol+=crystal[i]->vol;
	}
	CM(0)=x/vol;
	CM(1)=y/vol;
	CM(2)=z/vol;
}

void aggregate::rotate(double x, double y, double z, double alpha) {
	for(int i=0;i<N;i++) {
		crystal[i]->rotate(x,y,z,alpha);
	}
}

void aggregate::rotate(Matrix3d& R) {
	for(int i=0;i<N;i++) {
		crystal[i]->rotate(R);
	}
}

void aggregate::translate(double dx, double dy, double dz) {
	CM(0)+=dx;
	CM(1)+=dy;
	CM(2)+=dz;
	for(int i=0;i<N;i++) {
		crystal[i]->translate(dx, dy, dz);
	}
}

void aggregate::translate(Vector3d& T) {
	CM+=T;
	for(int i=0;i<N;i++) {
		crystal[i]->translate(T);
	}
}

bool aggregate::collide(aggregate& NAP) {
//    vector<double> distance(64,9e100);
	double distance=9e100;
	double provvisional;

//#pragma omp parallel for private(provvisional)
	for(int i=0;i<N;i++) {
		for(int j=0;j<NAP.N;j++) {
			provvisional=crystal[i]->distance(NAP.crystal[j]);
//            if(provvisional<distance[omp_get_thread_num()])
//            {
//                distance[omp_get_thread_num()]=provvisional;
//            }
			if(provvisional<distance) {
				distance=provvisional;
			}
		}
	}
//	provvisional=9e100;
//	for(int i=0;i<distance.size();i++)
//	{
//		if(distance[i]<provvisional)
//		{
//			provvisional=distance[i];
//		}
//	}
//    if (provvisional<9e100-1)
	if(distance<9e99) {
//    	NAP.translate(0.,0.,-provvisional);
		NAP.translate(0.0,0.0,-distance);
		return true;
	}
	else {
		return false;
	}
}

void aggregate::update_mic() {   /// CONTROLLARE /// magari non serve
	set_CM();
	reset_CM();
	set_Rmax();
	set_mass();
	set_Rgyr(); // suppose reset_CM and updated mass
	//cout<<"Sto per orientare orizzontalmente"<<endl;
	//orient_horizontally();
	//cout<<"Orizzontaleggiato"<<endl;
	set_v();
	//set_CM();
	//reset_CM();
}

void aggregate::reset_CM() {
	translate(-CM(0),-CM(1),-CM(2));
}

void aggregate::save_aggregate(double d) {
	/////cout<<"Salvo l'aggregato"<<endl;
	time_t timer_start;
	time_t timer_stop;
	double seconds;
	time(&timer_start);

	ostringstream D_max,N_pristine,Mass;
	D_max<<(Maximum_dimension/1000.);   //millimeters
	N_pristine<<N;
	Mass<<(mass*1000.);  // grams
	string Max_Dim=D_max.str();
	string Prist_Num=N_pristine.str();
	string MASS=Mass.str();
	string extension(".dat");
	string underscore("_");
	Prist_Num.append(underscore);
	Prist_Num.append(Max_Dim);
	Prist_Num.append(underscore);
	Prist_Num.append(MASS);
	Prist_Num.append(extension);
	ofstream OUT(Prist_Num.c_str());
	MatrixX3i point;
	//ofstream OUT2("vertexes.dat");
	for(int i=0; i<N; i++) {
		//cout<<"riempo "<<i<<endl;
		crystal[i]->fill(point, d);
		//cout<<" riempito"<<endl;
		OUT<<point<<endl;
		point.resize(0, 3);
		//OUT2<<crystal[i]->vertex;
	}
	OUT.close();
	time(&timer_stop);
	seconds=difftime(timer_stop,timer_start);
	//////cout<<"Per salvare l'aggregato ho impiegato  "<<seconds<<" secondi"<<endl;
	//OUT2.close();
}

void aggregate::update_maximum_dimension() {
	double squared_distance=0.;
	double provv_distance=0.;
	Vector3d diff;
	for(int i=0; i<(N-1); i++) {
		for(int j=(i+1); j<N; j++) {
			for(int ii=0; ii<crystal[i]->N_vertexes; ii++) {
				for(int jj=0; jj<crystal[j]->N_vertexes; jj++) {
					diff=crystal[i]->vertex.col(ii)-crystal[j]->vertex.col(jj);
					provv_distance=diff.dot(diff);
					if(provv_distance>squared_distance) {
						squared_distance=provv_distance;
					}
				}
			}
		}
	}
	double distance=sqrt(squared_distance);
	if(distance>Maximum_dimension) {
		Maximum_dimension=distance;
	}
}

void aggregate::orient_horizontally() {
	unsigned int N_vert=0;
	for(int i=0;i<N;i++) {
		N_vert+=crystal[i]->N_vertexes;
	}
	Matrix3Xd points(3,N_vert);
	N_vert=crystal[0]->N_vertexes;
	//cout<<"vado con la costruzione dei punti"<<endl;
	points.block(0,0,3,N_vert)=crystal[0]->vertex;
	//cout<<"vado iterativamente con la costruzione dei punti"<<endl;
	for(int i=1;i<N;i++) {
		//cout<<"Sto per fare la "<<i<<" iterazione"<<endl;
		points.block(0,N_vert,3,crystal[i]->N_vertexes)=crystal[i]->vertex;
		//cout<<"Aggiorno N_vert "<<N_vert<<endl;
		N_vert+=crystal[i]->N_vertexes;
		//cout<<"Aggiornato N_vert "<<N_vert<<" Ho fatto la "<<i<<" iterazione"<<endl;
	}
	/// estimate tollerance
	
	double tolerance=0.001;
	Vector3d center;
	Vector3d axes;
	Matrix3d rotation;
//	cout<<"invoco min_encl_ellipsoid"<<endl;
	min_encl_ellipsoid(points, tolerance, center, axes, rotation);
//	cout<<"I tre assi sono "<<axes<<endl;
	center=-center;
	translate(center);
	rotation.transposeInPlace();
	rotate(rotation);
	//Maximum_dimension=2*axes(2);
}

void aggregate::save_vtk_mesh(string prefix) {
	ostringstream D_max,N_pristine,Mass;
	D_max<<(Maximum_dimension/1000.);   //millimeters
	N_pristine<<N;
	Mass<<(mass*1000.);  // grams
	string Max_Dim=D_max.str();
	string Prist_Num=N_pristine.str();
	string MASS=Mass.str();
	string extension(".vtk");
	string underscore("_");
	prefix.append(underscore);
	prefix.append(Prist_Num);
	prefix.append(underscore);
	prefix.append(Max_Dim);
	prefix.append(underscore);
	prefix.append(MASS);
	prefix.append(extension);
	ofstream OUT(prefix.c_str());
	OUT<<"# vtk DataFile Version 1.0"<<endl;
	OUT<<Prist_Num.c_str()<<endl;
	unsigned int NV=0;
	unsigned int NF=0;
	unsigned int TOT=0;
	for(unsigned int i=0;i<N;i++) {
		NV+=crystal[i]->N_vertexes;
		NF+=crystal[i]->N_faces;
		for(unsigned int j=0;j<crystal[i]->N_faces;j++) {
			TOT+=(1+crystal[i]->face[j].size());
		}
	}
	OUT<<"ASCII"<<endl<<endl<<"DATASET POLYDATA"<<endl<<"POINTS "<<NV<<" float"<<endl;
	for(unsigned int i=0;i<N;i++) {
		for(unsigned int j=0;j<crystal[i]->N_vertexes;j++) {
			OUT<<crystal[i]->vertex.col(j).transpose()<<endl;
		}
	}	
	OUT<<"POLYGONS "<<NF<<" "<<TOT<<endl;
	NV=0;
	for(unsigned int i=0;i<N;i++) {
		for(unsigned int j=0;j<crystal[i]->N_faces;j++) {
			OUT<<crystal[i]->face[j].size()<<" ";
			for(int k=0;k<crystal[i]->face[j].size();k++) {
				OUT<<NV+crystal[i]->face[j][k]<<" ";
			}
			OUT<<endl;
		}
		NV+=crystal[i]->N_vertexes;
	}
	OUT.close();
}

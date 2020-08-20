#include "pristine.h"
#include "geom_lib.h"

pristine::pristine() {
	//ctor
}

pristine::~pristine() {
	//dtor
}

double pristine::distance(const pristine* NAP) {
	double dist=9.e100;
	double provv;
	double nx,ny,nz;            // coordinates of the normal vector of each face
	double* polyX;              // vectors of coordinates of vertexes for pip algortihm
	double* polyY;
	int N_vert_face;            // number of vertexes for each face
	double d;             // parameters of the equation // per ora mi serve solo d

	bool pip;           //point in polygon (convex)

	for(int i=0;i<NAP->N_faces;i++) {
		// verify that the NAP's face has a NEGATIVE orientation on the z-axis
		plane_normal(NAP->vertex.col(NAP->face[i][0]),NAP->vertex.col(NAP->face[i][1]),NAP->vertex.col(NAP->face[i][2]),nx,ny,nz);
		if(nz<0) {
			N_vert_face=NAP->face[i].size();
			polyX=new double[N_vert_face];
			polyY=new double[N_vert_face];
			for(int k=0;k<N_vert_face;k++) {
				polyX[k]=NAP->vertex(0,NAP->face[i][k]);
				polyY[k]=NAP->vertex(1,NAP->face[i][k]);
			}
			for(int j=0;j<N_vertexes;j++) {
				// verify with a point in polygon algorithm that there will be a collision
				pip=point_in_polygon_convex(N_vert_face,-1,polyX,polyY,vertex(0,j),vertex(1,j));
				if(pip) {
					// find the vertical distance from the CAP_pristine vertex's and the NAP_pristine face's
					d=-(nx*NAP->vertex(0,NAP->face[i][0])+ny*NAP->vertex(1,NAP->face[i][0])+nz*NAP->vertex(2,NAP->face[i][0]));
					provv=-(nx*vertex(0,j)+ny*vertex(1,j)+nz*vertex(2,j)+d)/nz;
					if(provv<dist)dist=provv;
				}
			}
			delete polyX;
			delete polyY;
		}
	}

	for(int i=0;i<N_faces;i++)
	{
		// verify that the CAP's face has a POSITIVE orientation on the z-axis
		plane_normal(vertex.col(face[i][0]),vertex.col(face[i][1]),vertex.col(face[i][2]),nx,ny,nz);
		if(nz>0) {
			N_vert_face=face[i].size();
			polyX=new double[N_vert_face];
			polyY=new double[N_vert_face];
			for(int k=0;k<N_vert_face;k++) {
				polyX[k]=vertex(0,face[i][k]);
				polyY[k]=vertex(1,face[i][k]);
			}
			for(int j=0;j<NAP->N_vertexes;j++) {
				// verify with a point in polygon algorithm that will be a collision
				pip=point_in_polygon_convex(N_vert_face,1,polyX,polyY,NAP->vertex(0,j),NAP->vertex(1,j));
				if(pip) {
					// find the vertical distance from the CAP_pristine vertex's and the NAP_pristine face's
					d=-(nx*vertex(0,face[i][0])+ny*vertex(1,face[i][0])+nz*vertex(2,face[i][0]));
					provv=(nx*NAP->vertex(0,j)+ny*NAP->vertex(1,j)+nz*NAP->vertex(2,j)+d)/nz;
					if(provv<dist)dist=provv;
				}
			}
			delete polyX;
			delete polyY;
		}
	}

	double diff_sidesEF,diff_sidesPQ;
	double px,py,qx,qy,ex,ey,fx,fy;

	double t,r;

	for(int i=0;i<N_edges;i++) {
		for(int j=0;j<NAP->N_edges;j++) {
			// verify that the two segments cross on the horizontal plane
			// it is sufficient to use cross_prod_2D
			ex=vertex(0,edge[i][0]);
			ey=vertex(1,edge[i][0]);
			fx=vertex(0,edge[i][1]);
			fy=vertex(1,edge[i][1]);

			px=NAP->vertex(0,NAP->edge[j][0]);
			py=NAP->vertex(1,NAP->edge[j][0]);
			qx=NAP->vertex(0,NAP->edge[j][1]);
			qy=NAP->vertex(1,NAP->edge[j][1]);

			/// vedere se il calcolo diretto dei parametri r t non sia più veloce del check ///
			diff_sidesEF=(cross_prod_2D(px-ex,py-ey,fx-ex,fy-ey))*(cross_prod_2D(qx-ex,qy-ey,fx-ex,fy-ey));
			if(diff_sidesEF<0) {
				diff_sidesPQ=(cross_prod_2D(ex-qx,ey-qy,px-qx,py-qy))*(cross_prod_2D(fx-qx,fy-qy,px-qx,py-qy));
				if(diff_sidesPQ<0) {
					// find the vertical distance of the two segments
					if((fx-ex)==0) {
						t=(py+(ex-px)*(qy-py)/(qx-px)-ey)/(fy-ey-(fx-ex)*(qy-py)/(qx-px));
						r=(ex+t*(fx-ex)-px)/(qx-px);
					}
					else {
						r=(ey+(px-ex)*(fy-ey)/(fx-ex)-py)/(qy-py-(qx-px)*(fy-ey)/(fx-ex));
						t=(px+r*(qx-px)-ex)/(fx-ex);
					}
					provv=NAP->vertex(2,NAP->edge[j][0])+r*(NAP->vertex(2,NAP->edge[j][1])-NAP->vertex(2,NAP->edge[j][0]))
											-(vertex(2,edge[i][0])+t*(vertex(2,edge[i][1])-vertex(2,edge[i][0])));
					if(provv<dist)dist=provv;
				}
			}
			/// fine sezione da rivedere ///
		}
	}
	return dist;
}

void pristine::rotate(double x, double y, double z, double alpha) {
	double invnorma=1.0/sqrt(x*x+y*y+z*z);
	x*=invnorma;
	y*=invnorma;
	z*=invnorma;
	Matrix3d R;
	double sina=sin(alpha);
	double cosa=cos(alpha);

	R << x*x+(1-x*x)*cosa,(1-cosa)*x*y-sina*z,(1-cosa)*x*z+sina*y,
		 (1-cosa)*y*x+sina*z,y*y+(1-y*y)*cosa,(1-cosa)*y*z-sina*x,
		 (1-cosa)*z*x-sina*y,(1-cosa)*z*y+sina*x,z*z+(1-z*z)*cosa;
	CM=R*CM;
	vertex=R*vertex;
}

void pristine::rotate(Matrix3d& R) {
	CM=R*CM;
	vertex=R*vertex;
}

void pristine::translate(double dx, double dy, double dz) {
	Vector3d delta;
	delta<<dx,dy,dz;
	translate(delta);
}

void pristine::translate(Vector3d& T) {
	CM+=T;
	for(int i=0;i<N_vertexes;i++) {
		vertex.col(i)+=T;
	}
}

void pristine::radius_of_gyration(double d) {
	Rgyr=1.;
}

void pristine::fill(MatrixX3i& point,double d) {
	double minX,maxX,minY,maxY,minZ,maxZ;
	minX=(vertex.row(0)).minCoeff();
	maxX=(vertex.row(0)).maxCoeff();
	minY=(vertex.row(1)).minCoeff();
	maxY=(vertex.row(1)).maxCoeff();
	minZ=(vertex.row(2)).minCoeff();
	maxZ=(vertex.row(2)).maxCoeff();
	int NminX,NmaxX,NminY,NmaxY,NminZ,NmaxZ;
	NminX=floor(minX/d)-1;
	NmaxX=ceil(maxX/d)+1;
	NminY=floor(minY/d)-1;
	NmaxY=ceil(maxY/d)+1;
	NminZ=floor(minZ/d)-1;
	NmaxZ=ceil(maxZ/d)+1;
	int N_point=(NmaxX-NminX)*(NmaxY-NminY)*(NmaxZ-NminZ);
	point.resize(N_point,3); // bounding box
	N_point=0;

	bool change=false;
	bool inside=false;
	double* nx;
	double* ny;
	double* nz;
	bool* pip;
	double* polyX;
	double* polyY;
	int N_vert_face;
	double* dcoeff;
	double Zplane;

	int count_crossed_faces=0; ///
	vector<int> crossed_faces_idx;
	vector<double> crossed_faces_Z;
	bool forced_inside=false;  /// perchè era true di default???
	bool forced_outside=false;

	vector<int> crossing_vertexes;
	bool alert_crossing_vertexes=false;
	vector<int> crossing_edges;
	vector<double> Z_crossing_edge;
	double vx,vy,vz,tx,ty;
	bool alert_crossing_edges=false;
	double dot_prod;
	double dot_prod1,dot_prod2;

	bool inside_projection=true;
	vector<int> nearest_faces;
	nearest_faces.push_back(999);
	nearest_faces.push_back(999);
	vector<double> dist_nearest_faces;
	double dist_face;
	double dist_provv;
	double dist_prov;
	double mx0,my0,mz0,mx1,my1,mz1;

	double t,xt,yt,zt;

	nx=new double[N_faces];
	ny=new double[N_faces];
	nz=new double[N_faces];
	pip=new bool[N_faces];
	dcoeff=new double[N_faces];

	for(int l=0;l<N_faces;l++) {
		plane_normal(vertex.col(face[l][0]),vertex.col(face[l][1]),vertex.col(face[l][2]),nx[l],ny[l],nz[l]);
		dcoeff[l]=-nx[l]*vertex(0,face[l][0])-ny[l]*vertex(1,face[l][0])-nz[l]*vertex(2,face[l][0]);
	}

	for(int i=NminX;i<NmaxX;i++) {
		for(int j=NminY;j<NmaxY;j++) {
			for(int nf=0;nf<N_faces;nf++) {
				N_vert_face=face[nf].size();
				polyX=new double[N_vert_face];
				polyY=new double[N_vert_face];
				for(int m=0;m<N_vert_face;m++) {
					polyX[m]=vertex(0,face[nf][m]);
					polyY[m]=vertex(1,face[nf][m]);
				}
				pip[nf]=point_in_polygon_convex(N_vert_face,(nz[nf]>0?1:-1),polyX,polyY,i*d,j*d);
				delete polyX;
				delete polyY;
			}
			// Find potential vertexes on the vertical
			for(int nv=0;nv<N_vertexes;nv++) {
				if((vertex(0,nv)==i*d)&&(vertex(1,nv)==j*d)) {
					crossing_vertexes.push_back(nv);
				}
				if(crossing_vertexes.size()) {
					alert_crossing_vertexes=true;
				}
			}
			// Find potential edges on the vertical NO VERTEXES!!!
			for(int ne=0;ne<N_edges;ne++) {
				vx=vertex(0,edge[ne][1])-vertex(0,edge[ne][0]);
				vy=vertex(1,edge[ne][1])-vertex(1,edge[ne][0]);
				if(vx&&vy) {
					tx=(i*d-vertex(0,edge[ne][0]))/vx;
					ty=(j*d-vertex(1,edge[ne][0]))/vy;
					if((tx)==(ty)) {
						if((0.<tx)&&(tx<1.)) {
							vz=vertex(2,edge[ne][1])-vertex(2,edge[ne][0]);
							crossing_edges.push_back(ne);
							Z_crossing_edge.push_back(vertex(2,edge[ne][0])+tx*vz);
						}
					}
				}
				else {
					if(vx) {
						if(j*d==vertex(1,edge[ne][0])) {
							tx=(i*d-vertex(0,edge[ne][0]))/vx;
							if((0.<tx)&&(tx<1.)) {
								vz=vertex(2,edge[ne][1])-vertex(2,edge[ne][0]);
								crossing_edges.push_back(ne);
								Z_crossing_edge.push_back(vertex(2,edge[ne][0])+tx*vz);
							}
						}
					}
					else if(vy) {
						if(i*d==vertex(0,edge[ne][0])) {
							ty=(j*d-vertex(1,edge[ne][0]))/vy;
							if((0.<ty)&&(ty<1.)) {
								vz=vertex(2,edge[ne][1])-vertex(2,edge[ne][0]);
								crossing_edges.push_back(ne);
								Z_crossing_edge.push_back(vertex(2,edge[ne][0])+ty*vz);
							}
						}
					}
					else {
						if((i*d==vertex(0,edge[ne][0]))&&(j*d==vertex(1,edge[ne][0]))) {
							cout<<"Trivial case. skipping"<<endl;
						}
					}
				}
				if(crossing_edges.size()) {
					alert_crossing_edges=true;
				}
			}
			for(int k=NminZ+1;k<NmaxZ;k++) {
				for(int l=0;l<N_faces;l++) {
					if(pip[l]) {
						Zplane=(-dcoeff[l]-nx[l]*i*d-ny[l]*j*d)/nz[l];
						if(((k*d)>Zplane)&&(((k-1)*d)<=Zplane)) {
							change=!change;
						}
					}
				}
				if(change) inside=!inside;
				if(alert_crossing_edges) {
					for(unsigned int nce=0; nce<crossing_edges.size(); nce++) {
						if((k*d)==Z_crossing_edge[nce]) forced_inside=true;
						else if((!forced_inside)&&(!forced_outside)&&((((k-1)*d)<=Z_crossing_edge[nce])&&(k*d>Z_crossing_edge[nce]))) {   
							dist_nearest_faces.push_back(999999.);
							dist_nearest_faces.push_back(9999999.);
							for(int nf=0;nf<N_faces;nf++) {
								inside_projection=true;
								t=-(dcoeff[nf]+d*(nx[nf]*i+ny[nf]*j+nz[nf]*k))/(nx[nf]*nx[nf]+ny[nf]*ny[nf]+nz[nf]*nz[nf]);
								xt=nx[nf]*t+i*d;
								yt=ny[nf]*t+j*d;
								zt=nz[nf]*t+k*d;
								
								Vector3d A=vertex.col(face[nf][face[nf].size()-2]);
								Vector3d B=vertex.col(face[nf][face[nf].size()-1]);
								Vector3d C=vertex.col(face[nf][0]);
								Vector3d P;
								P<<xt,yt,zt;
								Vector3d cross1;
								Vector3d cross2;
								Vector3d P_B=P-B;
								cross1=P_B.cross(C-B);
								cross2=P_B.cross(A-B);
								dot_prod=cross1.dot(cross2);
								if(dot_prod<0) inside_projection=false;
								else {
									A=vertex.col(face[nf][face[nf].size()-1]);
									B=vertex.col(face[nf][0]);
									C=vertex.col(face[nf][1]);
									P_B=P-B;
									cross1=P_B.cross(C-B);
									cross2=P_B.cross(A-B);
									dot_prod=cross1.dot(cross2);
									if(dot_prod<0) inside_projection=false;
									else
									{
										for(unsigned int nfv=1;nfv<(face[nf].size()-1);nfv++)
										{
											A=vertex.col(face[nf][nfv-1]);
											B=vertex.col(face[nf][nfv]);
											C=vertex.col(face[nf][nfv+1]);
									P_B=P-B;
									cross1=P_B.cross(C-B);
									cross2=P_B.cross(A-B);
									dot_prod=cross1.dot(cross2);
									if(dot_prod<0) inside_projection=false;
										}
									}
								}
								///////
								if(inside_projection) // distanza punto piano
								{
									//if((j==4)&&(i==0))cout<<"Sono dentro la proiezione della faccia "<<nf<<endl;
									dist_face=pow((nx[nf]*i+ny[nf]*j+nz[nf]*k)*d-dcoeff[nf],2)/(nx[nf]*nx[nf]+ny[nf]*ny[nf]+nz[nf]*nz[nf]);
									if(dist_face<dist_nearest_faces[0])
									{
										dist_nearest_faces[1]=dist_nearest_faces[0];
										nearest_faces[1]=nearest_faces[0];
										dist_nearest_faces[0]=dist_face;
										nearest_faces[0]=nf;
									}
									else if(dist_face<dist_nearest_faces[1])
									{
										dist_nearest_faces[1]=dist_face;
										nearest_faces[1]=nf;
									}
								}
								else // distanza punto segmento
								{
									//if((j==4)&&(i==0))cout<<"Sono fuori la proiezione della faccia "<<nf<<endl;
									dot_prod1=((vertex(0,face[nf][0])-vertex(0,face[nf][face[nf].size()-1]))*(i*d-vertex(0,face[nf][face[nf].size()-1]))
											 +(vertex(1,face[nf][0])-vertex(1,face[nf][face[nf].size()-1]))*(j*d-vertex(1,face[nf][face[nf].size()-1]))
											 +(vertex(2,face[nf][0])-vertex(2,face[nf][face[nf].size()-1]))*(k*d-vertex(2,face[nf][face[nf].size()-1])));
									dot_prod2=((vertex(0,face[nf][face[nf].size()-1])-vertex(0,face[nf][0]))*(i*d-vertex(0,face[nf][0]))
											 +(vertex(1,face[nf][face[nf].size()-1])-vertex(1,face[nf][0]))*(j*d-vertex(1,face[nf][0]))
											 +(vertex(2,face[nf][face[nf].size()-1])-vertex(2,face[nf][0]))*(k*d-vertex(2,face[nf][0])));
									dot_prod=dot_prod1*dot_prod2;
									if(dot_prod>0)
									{
										dist_face=(pow(i*d-vertex(0,face[nf][face[nf].size()-1]),2)+pow(j*d-vertex(1,face[nf][face[nf].size()-1]),2)+pow(k*d-vertex(2,face[nf][face[nf].size()-1]),2))
												 -(dot_prod1*dot_prod1)/(pow(vertex(0,face[nf][0])-vertex(0,face[nf][face[nf].size()-1]),2)+pow(vertex(1,face[nf][0])-vertex(1,face[nf][face[nf].size()-1]),2)+pow(vertex(2,face[nf][0])-vertex(2,face[nf][face[nf].size()-1]),2));
									}
									else if(dot_prod<0) // distanza dai vertici
									{
										dist_face=pow(i*d-vertex(0,face[nf][0]),2)+pow(j*d-vertex(1,face[nf][0]),2)+pow(k*d-vertex(2,face[nf][0]),2);
										dist_provv=pow(i*d-vertex(0,face[nf][face[nf].size()-1]),2)+pow(j*d-vertex(1,face[nf][face[nf].size()-1]),2)+pow(k*d-vertex(2,face[nf][face[nf].size()-1]),2);
										if(dist_provv<dist_face)
										{
											dist_face=dist_provv;
										}
									}
									else
									{
										if(dot_prod1==0)
										{
											dist_face=pow(vertex(0,face[nf][face[nf].size()-1])-i*d,2)+pow(vertex(1,face[nf][face[nf].size()-1])-j*d,2)+pow(vertex(2,face[nf][face[nf].size()-1])-k*d,2);
										}
										else if(dot_prod2==0)
										{
											dist_face=pow(vertex(0,face[nf][0])-i*d,2)+pow(vertex(1,face[nf][0])-j*d,2)+pow(vertex(2,face[nf][0])-k*d,2);
										}
										else
										{
											cout<<"PANIC"<<endl;
										}
									}
									for(unsigned int nfe=1;nfe<face[nf].size();nfe++)
									{
										dot_prod1=((vertex(0,face[nf][nfe])-vertex(0,face[nf][nfe-1]))*(i*d-vertex(0,face[nf][nfe-1]))
												 +(vertex(1,face[nf][nfe])-vertex(1,face[nf][nfe-1]))*(j*d-vertex(1,face[nf][nfe-1]))
												 +(vertex(2,face[nf][nfe])-vertex(2,face[nf][nfe-1]))*(k*d-vertex(2,face[nf][nfe-1])));
										dot_prod2=((vertex(0,face[nf][nfe-1])-vertex(0,face[nf][nfe]))*(i*d-vertex(0,face[nf][nfe]))
												 +(vertex(1,face[nf][nfe-1])-vertex(1,face[nf][nfe]))*(j*d-vertex(1,face[nf][nfe]))
												 +(vertex(2,face[nf][nfe-1])-vertex(2,face[nf][nfe]))*(k*d-vertex(2,face[nf][nfe])));
										dot_prod=dot_prod1*dot_prod2;
										if(dot_prod>0)
										{
											dist_provv=(pow(i*d-vertex(0,face[nf][nfe-1]),2)+pow(j*d-vertex(1,face[nf][nfe-1]),2)+pow(k*d-vertex(2,face[nf][nfe-1]),2))
												 -(dot_prod1*dot_prod1)/(pow(vertex(0,face[nf][nfe])-vertex(0,face[nf][nfe-1]),2)+pow(vertex(1,face[nf][nfe])-vertex(1,face[nf][nfe-1]),2)+pow(vertex(2,face[nf][nfe])-vertex(2,face[nf][nfe-1]),2));
										}
										else if(dot_prod<0) // distanza dai vertici
										{
											dist_prov=pow(i*d-vertex(0,face[nf][nfe]),2)+pow(j*d-vertex(1,face[nf][nfe]),2)+pow(k*d-vertex(2,face[nf][nfe]),2);
											dist_provv=pow(i*d-vertex(0,face[nf][nfe-1]),2)+pow(j*d-vertex(1,face[nf][nfe-1]),2)+pow(k*d-vertex(2,face[nf][nfe-1]),2);
											if(dist_prov<dist_provv)
											{
												dist_provv=dist_prov;
											}
										}                              
										else
										{
											if(dot_prod1==0)
										{
											dist_provv=pow(vertex(0,face[nf][nfe-1])-i*d,2)+pow(vertex(1,face[nf][nfe-1])-j*d,2)+pow(vertex(2,face[nf][nfe-1])-k*d,2);
										}
										else if(dot_prod2==0)
										{
											dist_provv=pow(vertex(0,face[nf][nfe])-i*d,2)+pow(vertex(1,face[nf][nfe])-j*d,2)+pow(vertex(2,face[nf][nfe])-k*d,2);
										}
										else
										{
											cout<<"PANIC"<<endl;
										}
										}
										if(dist_provv<dist_face)
										{
											dist_face=dist_provv;
										}
									}
								}
								if(dist_face<dist_nearest_faces[0])
								{
									dist_nearest_faces[1]=dist_nearest_faces[0];
									nearest_faces[1]=nearest_faces[0];
									dist_nearest_faces[0]=dist_face;
									nearest_faces[0]=nf;
								}
								else if(dist_face<dist_nearest_faces[1])
								{
									dist_nearest_faces[1]=dist_face;
									nearest_faces[1]=nf;
								}
							}
							dist_nearest_faces.clear();
							mx0=0.;
							my0=0.;
							mz0=0.;
							for(unsigned int nfv=0;nfv<face[nearest_faces[0]].size();nfv++)
							{
								mx0+=vertex(0,face[nearest_faces[0]][nfv]);
								my0+=vertex(1,face[nearest_faces[0]][nfv]);
								mz0+=vertex(2,face[nearest_faces[0]][nfv]);
							}
							mx0=mx0/((double)face[nearest_faces[0]].size());
							my0=my0/((double)face[nearest_faces[0]].size());
							mz0=mz0/((double)face[nearest_faces[0]].size());
							mx1=0.;
							my1=0.;
							mz1=0.;
							for(unsigned int nfv=0;nfv<face[nearest_faces[1]].size();nfv++)
							{
								mx1+=vertex(0,face[nearest_faces[1]][nfv]);
								my1+=vertex(1,face[nearest_faces[1]][nfv]);
								mz1+=vertex(2,face[nearest_faces[1]][nfv]);
							}
							mx1=mx1/((double)face[nearest_faces[1]].size());
							my1=my1/((double)face[nearest_faces[1]].size());
							mz1=mz1/((double)face[nearest_faces[1]].size());
							dot_prod1=(mx1-mx0)*nx[nearest_faces[0]]+(my1-my0)*ny[nearest_faces[0]]+(mz1-mz0)*nz[nearest_faces[0]];
							dot_prod2=(mx0-mx1)*nx[nearest_faces[1]]+(my0-my1)*ny[nearest_faces[1]]+(mz0-mz1)*nz[nearest_faces[1]];
							if((dot_prod1<0)&&(dot_prod2<0))
							{
								//if((j==4)&&(i==0))cout<<"sono facce divergenti ricordo che"<<nearest_faces[0]<<" e "<<nearest_faces[1]<<endl;
								dot_prod1=(i*d-mx0)*nx[nearest_faces[0]]+(j*d-my0)*ny[nearest_faces[0]]+(k*d-mz0)*nz[nearest_faces[0]];
								dot_prod2=(i*d-mx1)*nx[nearest_faces[1]]+(j*d-my1)*ny[nearest_faces[1]]+(k*d-mz1)*nz[nearest_faces[1]];
								dot_prod=dot_prod1*dot_prod2;
								if(dot_prod)
								{
									if((dot_prod1<0)&&(dot_prod2<0))
									{
										forced_inside=true;
									}
									else
									{
										forced_outside=true;
									}
								}
								else
								{
									if((dot_prod1<0)||(dot_prod2<0)||(dot_prod1==dot_prod2))
									{
										forced_inside=true;
										//if((j==4)&&(i==0))cout<<"< && <forzato ad entrare <||<"<<endl;
									}
									else
									{
										//if((j==4)&&(i==0))cout<<"< && <forzato ad uscire else"<<endl;
										forced_outside=true;
									}
								}
							}
							else if ((dot_prod1>0)&&(dot_prod2>0))
							{
								//if((j==4)&&(i==0))cout<<"Facce convergenti dot_prod1 "<<dot_prod1<<"dot_prod2"<<dot_prod2<<endl;
								dot_prod1=(i*d-mx0)*nx[nearest_faces[0]]+(j*d-my0)*ny[nearest_faces[0]]+(k*d-mz0)*nz[nearest_faces[0]];
								dot_prod2=(i*d-mx1)*nx[nearest_faces[1]]+(j*d-my1)*ny[nearest_faces[1]]+(k*d-mz1)*nz[nearest_faces[1]];
								dot_prod=dot_prod1*dot_prod2;
								if(dot_prod)
								{
									if((dot_prod1>0)&&(dot_prod2>0))
										{
											forced_outside=true;
										//if((j==4)&&(i==0))cout<<"> && > Forzato ad uscire >&&>"<<endl;
									}
									else
									{
											forced_inside=true;
											//if((j==4)&&(i==0))cout<<"> && > forzato ad entrare ##"<<endl;
										}
								}
								else
								{
									if((dot_prod1>0)||(dot_prod2>0)||(dot_prod1==dot_prod2))
									{
										forced_inside=true;
										//if((j==4)&&(i==0))cout<<"> && > Forzato ad entrare >||>"<<endl;
									}
									else
									{
										//if((j==4)&&(i==0))cout<<"> && > forzato ad uscire else"<<endl;
										forced_outside=true;
									}
								}
							}
							else cout<<"Damn edge  dot1 "<<dot_prod1<<"  dot2 "<<dot_prod2<<" coord"<<i<<" "<<j<<" "<<k<<endl;
						}  /// qui si chiude l'if(!forced_inside) che sarà cambiato (ora è stato identato)
					}///agginta questa parentesi per includere l'if di prima nel for su tutti i bordi
				}  /// qui si chiude l'alert crossing edges

				if(alert_crossing_vertexes) {
					for(unsigned int ncv=0;ncv<crossing_vertexes.size();ncv++) {
						if((k*d)==vertex(2,crossing_vertexes[ncv])) forced_inside=true; // I am on a vertex
						else if((!forced_inside)&&(!forced_outside)&&((((k-1)*d)<=vertex(2,crossing_vertexes[ncv]))&&(k*d>vertex(2,crossing_vertexes[ncv])))) {
							// I have crossed a vertex
							dist_nearest_faces.push_back(99999999.);
							dist_nearest_faces.push_back(999999999.);
							for(int nf=0;nf<N_faces;nf++) {                         
								inside_projection=true;
								t=-(dcoeff[nf]+d*(nx[nf]*i+ny[nf]*j+nz[nf]*k))/(nx[nf]*nx[nf]+ny[nf]*ny[nf]+nz[nf]*nz[nf]);
								xt=nx[nf]*t+i*d;
								yt=ny[nf]*t+j*d;
								zt=nz[nf]*t+k*d;
								
								Vector3d A=vertex.col(face[nf][face[nf].size()-2]);
								Vector3d B=vertex.col(face[nf][face[nf].size()-1]);
								Vector3d C=vertex.col(face[nf][0]);
								Vector3d P;
								P<<xt,yt,zt;
								Vector3d cross1;
								Vector3d cross2;
								Vector3d P_B=P-B;
								cross1=P_B.cross(C-B);
								cross2=P_B.cross(A-B);
								dot_prod=cross1.dot(cross2);
								if(dot_prod<0) inside_projection=false;
								else {
									A=vertex.col(face[nf][face[nf].size()-1]);
									B=vertex.col(face[nf][0]);
									C=vertex.col(face[nf][1]);
									P_B=P-B;
									cross1=P_B.cross(C-B);
									cross2=P_B.cross(A-B);
									dot_prod=cross1.dot(cross2);
									if(dot_prod<0) inside_projection=false;
									else {
										for(unsigned int nfv=1;nfv<(face[nf].size()-1);nfv++) {
											A=vertex.col(face[nf][nfv-1]);
											B=vertex.col(face[nf][nfv]);
											C=vertex.col(face[nf][nfv+1]);
											P_B=P-B;
											cross1=P_B.cross(C-B);
											cross2=P_B.cross(A-B);
											dot_prod=cross1.dot(cross2);
											if(dot_prod<0) inside_projection=false;
										}
									}
								}
								if(inside_projection) { // Inside the face projection compute distance point-plane
									dist_face=pow((nx[nf]*i+ny[nf]*j+nz[nf]*k)*d-dcoeff[nf],2)/(nx[nf]*nx[nf]+ny[nf]*ny[nf]+nz[nf]*nz[nf]);
									if(dist_face<dist_nearest_faces[0])	{
										dist_nearest_faces[1]=dist_nearest_faces[0];
										nearest_faces[1]=nearest_faces[0];
										dist_nearest_faces[0]=dist_face;
										nearest_faces[0]=nf;
									}
									else if(dist_face<dist_nearest_faces[1]) {
										dist_nearest_faces[1]=dist_face;
										nearest_faces[1]=nf;
									}
								}
								else { // point segment distance
									dot_prod1=((vertex(0,face[nf][0])-vertex(0,face[nf][face[nf].size()-1]))*(i*d-vertex(0,face[nf][face[nf].size()-1]))
											 +(vertex(1,face[nf][0])-vertex(1,face[nf][face[nf].size()-1]))*(j*d-vertex(1,face[nf][face[nf].size()-1]))
											 +(vertex(2,face[nf][0])-vertex(2,face[nf][face[nf].size()-1]))*(k*d-vertex(2,face[nf][face[nf].size()-1])));
									dot_prod2=((vertex(0,face[nf][face[nf].size()-1])-vertex(0,face[nf][0]))*(i*d-vertex(0,face[nf][0]))
											 +(vertex(1,face[nf][face[nf].size()-1])-vertex(1,face[nf][0]))*(j*d-vertex(1,face[nf][0]))
											 +(vertex(2,face[nf][face[nf].size()-1])-vertex(2,face[nf][0]))*(k*d-vertex(2,face[nf][0])));
									dot_prod=dot_prod1*dot_prod2;
									if(dot_prod>0) {
										dist_face=(pow(i*d-vertex(0,face[nf][face[nf].size()-1]),2)+pow(j*d-vertex(1,face[nf][face[nf].size()-1]),2)+pow(k*d-vertex(2,face[nf][face[nf].size()-1]),2))
												 -(dot_prod1*dot_prod1)/(pow(vertex(0,face[nf][0])-vertex(0,face[nf][face[nf].size()-1]),2)+pow(vertex(1,face[nf][0])-vertex(1,face[nf][face[nf].size()-1]),2)+pow(vertex(2,face[nf][0])-vertex(2,face[nf][face[nf].size()-1]),2));                                               
									}
									else if(dot_prod<0) { // distance from the vertexes
										dist_face=pow(i*d-vertex(0,face[nf][0]),2)+pow(j*d-vertex(1,face[nf][0]),2)+pow(k*d-vertex(2,face[nf][0]),2);
										dist_provv=pow(i*d-vertex(0,face[nf][face[nf].size()-1]),2)+pow(j*d-vertex(1,face[nf][face[nf].size()-1]),2)+pow(k*d-vertex(2,face[nf][face[nf].size()-1]),2);
										if(dist_provv<dist_face) dist_face=dist_provv;
									}
									else {
										if(dot_prod1==0) dist_face=pow(vertex(0,face[nf][face[nf].size()-1])-i*d,2)+pow(vertex(1,face[nf][face[nf].size()-1])-j*d,2)+pow(vertex(2,face[nf][face[nf].size()-1])-k*d,2);
										else if(dot_prod2==0) dist_face=pow(vertex(0,face[nf][0])-i*d,2)+pow(vertex(1,face[nf][0])-j*d,2)+pow(vertex(2,face[nf][0])-k*d,2);
										else cout<<"PANIC"<<endl;
									}                                    
									for(unsigned int nfe=1;nfe<face[nf].size();nfe++) {
										dot_prod1=((vertex(0,face[nf][nfe])-vertex(0,face[nf][nfe-1]))*(i*d-vertex(0,face[nf][nfe-1]))
												 +(vertex(1,face[nf][nfe])-vertex(1,face[nf][nfe-1]))*(j*d-vertex(1,face[nf][nfe-1]))
												 +(vertex(2,face[nf][nfe])-vertex(2,face[nf][nfe-1]))*(k*d-vertex(2,face[nf][nfe-1])));
										dot_prod2=((vertex(0,face[nf][nfe-1])-vertex(0,face[nf][nfe]))*(i*d-vertex(0,face[nf][nfe]))
												 +(vertex(1,face[nf][nfe-1])-vertex(1,face[nf][nfe]))*(j*d-vertex(1,face[nf][nfe]))
												 +(vertex(2,face[nf][nfe-1])-vertex(2,face[nf][nfe]))*(k*d-vertex(2,face[nf][nfe])));
										dot_prod=dot_prod1*dot_prod2;
										if(dot_prod>0) {
											dist_provv=(pow(i*d-vertex(0,face[nf][nfe-1]),2)+pow(j*d-vertex(1,face[nf][nfe-1]),2)+pow(k*d-vertex(2,face[nf][nfe-1]),2))
												 -(dot_prod1*dot_prod1)/(pow(vertex(0,face[nf][nfe])-vertex(0,face[nf][nfe-1]),2)+pow(vertex(1,face[nf][nfe])-vertex(1,face[nf][nfe-1]),2)+pow(vertex(2,face[nf][nfe])-vertex(2,face[nf][nfe-1]),2));
										}
										else if(dot_prod<0) { // distanza dai vertici
											dist_prov=pow(i*d-vertex(0,face[nf][nfe]),2)+pow(j*d-vertex(1,face[nf][nfe]),2)+pow(k*d-vertex(2,face[nf][nfe]),2);
											dist_provv=pow(i*d-vertex(0,face[nf][nfe-1]),2)+pow(j*d-vertex(1,face[nf][nfe-1]),2)+pow(k*d-vertex(2,face[nf][nfe-1]),2);
											if(dist_prov<dist_provv) dist_provv=dist_prov;
										}
										else {
											if(dot_prod1==0) dist_provv=pow(vertex(0,face[nf][nfe-1])-i*d,2)+pow(vertex(1,face[nf][nfe-1])-j*d,2)+pow(vertex(2,face[nf][nfe-1])-k*d,2);
											else if(dot_prod2==0) dist_provv=pow(vertex(0,face[nf][nfe])-i*d,2)+pow(vertex(1,face[nf][nfe])-j*d,2)+pow(vertex(2,face[nf][nfe])-k*d,2);
											else cout<<"PANIC"<<endl;
										}
										if(dist_provv<dist_face) dist_face=dist_provv;
									}                                    
								}
								if(dist_face<dist_nearest_faces[0]) {
									dist_nearest_faces[1]=dist_nearest_faces[0];
									nearest_faces[1]=nearest_faces[0];
									dist_nearest_faces[0]=dist_face;
									nearest_faces[0]=nf;
								}
								else if(dist_face<dist_nearest_faces[1]) {
									dist_nearest_faces[1]=dist_face;
									nearest_faces[1]=nf;
								}
							}
							dist_nearest_faces.clear();
							mx0=0.;
							my0=0.;
							mz0=0.;
							for(unsigned int nfv=0;nfv<face[nearest_faces[0]].size();nfv++) {
								mx0+=vertex(0,face[nearest_faces[0]][nfv]);
								my0+=vertex(1,face[nearest_faces[0]][nfv]);
								mz0+=vertex(2,face[nearest_faces[0]][nfv]);
							}
							mx0=mx0/((double)face[nearest_faces[0]].size());
							my0=my0/((double)face[nearest_faces[0]].size());
							mz0=mz0/((double)face[nearest_faces[0]].size());
							mx1=0.;
							my1=0.;
							mz1=0.;
							for(unsigned int nfv=0;nfv<face[nearest_faces[1]].size();nfv++) {
								mx1+=vertex(0,face[nearest_faces[1]][nfv]);
								my1+=vertex(1,face[nearest_faces[1]][nfv]);
								mz1+=vertex(2,face[nearest_faces[1]][nfv]);
							}
							mx1=mx1/((double)face[nearest_faces[1]].size());
							my1=my1/((double)face[nearest_faces[1]].size());
							mz1=mz1/((double)face[nearest_faces[1]].size());
							dot_prod1=(mx1-mx0)*nx[nearest_faces[0]]+(my1-my0)*ny[nearest_faces[0]]+(mz1-mz0)*nz[nearest_faces[0]];
							dot_prod2=(mx0-mx1)*nx[nearest_faces[1]]+(my0-my1)*ny[nearest_faces[1]]+(mz0-mz1)*nz[nearest_faces[1]];
							if((dot_prod1<0)&&(dot_prod2<0)) {
								dot_prod1=(i*d-mx0)*nx[nearest_faces[0]]+(j*d-my0)*ny[nearest_faces[0]]+(k*d-mz0)*nz[nearest_faces[0]];
								dot_prod2=(i*d-mx1)*nx[nearest_faces[1]]+(j*d-my1)*ny[nearest_faces[1]]+(k*d-mz1)*nz[nearest_faces[1]];
								dot_prod=dot_prod1*dot_prod2;
								if(dot_prod) {
									if((dot_prod1<0)&&(dot_prod2<0)) forced_inside=true;
									else forced_outside=true;
								}
								else {
									if((dot_prod1<0)||(dot_prod2<0)||(dot_prod1==dot_prod2)) forced_inside=true;
									else forced_outside=true;
								}
							}
							else if ((dot_prod1>0)&&(dot_prod2>0)) {
								dot_prod1=(i*d-mx0)*nx[nearest_faces[0]]+(j*d-my0)*ny[nearest_faces[0]]+(k*d-mz0)*nz[nearest_faces[0]];
								dot_prod2=(i*d-mx1)*nx[nearest_faces[1]]+(j*d-my1)*ny[nearest_faces[1]]+(k*d-mz1)*nz[nearest_faces[1]];
								dot_prod=dot_prod1*dot_prod2;
								if(dot_prod) {
									if((dot_prod1>0)&&(dot_prod2>0)) forced_outside=true;
									else forced_inside=true;
								}
								else {
									if((dot_prod1>0)||(dot_prod2>0)||(dot_prod1==dot_prod2)) forced_inside=true;
									else forced_outside=true;
								}
							}
							else cout<<"Damn vertex  dot1 "<<dot_prod1<<"  dot2 "<<dot_prod2<<endl;
						}
					}
				}
				if(forced_inside) inside=true;
				if(forced_outside) inside=false;
				if(inside) {
					point(N_point,0)=i;
					point(N_point,1)=j;
					point(N_point,2)=k;
					N_point++;
				}
				change=false;
//                degeneration=false;
				//forced_inside=true;
				forced_inside=false; /// nuovo metodo
				forced_outside=false;
				count_crossed_faces=0;
				crossed_faces_idx.clear();
				crossed_faces_Z.clear();
			}
			inside=false;
			change=false;
			crossing_vertexes.clear();
			crossing_edges.clear();
			Z_crossing_edge.clear();
			alert_crossing_edges=false;
			alert_crossing_vertexes=false;
		}
	}
	point.conservativeResize(N_point,3); // cut the crap not filled
	delete nx;
	delete ny;
	delete nz;
	delete pip;
	delete dcoeff;
}


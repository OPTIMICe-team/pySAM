#include "geom_lib.h"

double cross_prod_2D(double vx, double vy, double wx, double wy)
{
    return vx*wy-vy*wx;
}

void plane_normal(Vector3d A, Vector3d B, Vector3d C, double& nx, double &ny, double& nz)
{
    // it's basically a cross product (B-A)X(C-A)
    nx=(B(1)-A(1))*(C(2)-A(2))-(C(1)-A(1))*(B(2)-A(2));
    ny=(C(0)-A(0))*(B(2)-A(2))-(B(0)-A(0))*(C(2)-A(2));
    nz=(B(0)-A(0))*(C(1)-A(1))-(B(1)-A(1))*(C(0)-A(0));
}

bool point_in_polygon_convex(int N_poly, int sign, double* X, double* Y, double x, double y)
{
    // if sign is > 0 => polygon vertexes are defined counterclockwise
    double* Xt=new double[N_poly];
    double* Yt=new double[N_poly];
    int end=N_poly-1;
    //cout<<"la fine della faccia è"<<end<<"il punto da controllare è "<<x<<"  "<<y<<endl;

    for(int i=0;i<N_poly; i++)  // traslate polygon, (x,y) is the new origin
    {
        //cout<<"I punti della faccia passati al pip sono "<<X[i]<<"   "<<Y[i]<<endl;
        Xt[i]=X[i]-x;
        Yt[i]=Y[i]-y;
    }

    double tangent;
    bool pip=true;
    //cout<<"xt[end] "<<Xt[end]<<"  yt[end]"<<Yt[end]<<"  xt[0] "<<Xt[0]<<"  yt[0]"<<Yt[0]<<endl;
    tangent=Xt[end]*Yt[0]-Xt[0]*Yt[end];
    //cout<<"La prima tangente è"<<tangent<<endl;
    if(sign*tangent<0) pip=false;

    for(int i=0; i<(N_poly-1) && pip; i++)
    {
        tangent=Xt[i]*Yt[i+1]-Xt[i+1]*Yt[i];
        //cout<<"Le altre tangenti sono "<<tangent<<endl;
        if(sign*tangent<0) pip=false;
    }

    delete Xt;
    delete Yt;

    return pip;
}

void min_encl_ellipsoid(Matrix3Xd& points, double tol, Vector3d& center, Vector3d& axes, Matrix3d& rotation)
{
	//cout<<"Sono arrivato all'inizio di kachyan"<<endl;
	const unsigned int N=points.cols();
	MatrixXd Q(4,N);
	Q.block(0,0,3,N)=points;
	Q.row(3)=RowVectorXd::Constant(N,1.);
	MatrixXd Qt(N,4);
	Qt=Q.transpose();
	unsigned int count=0;
	double err=1000;
	VectorXd u(N);
	u=VectorXd::Constant(N,1./((double)N));
	VectorXd u_new(N);
	Matrix4d X;
	MatrixXd NN(N,N);
	VectorXd M(N);
	double maxM,step;
	unsigned int maxMidx;
	
	//cout<<"inizio l'algoritmo"<<endl;
	while((err>tol)&&(count<10000))
	{
		X=Q*u.asDiagonal()*Qt;
		NN=Qt*X.inverse()*Q;
		M=NN.diagonal();
		maxM=M.maxCoeff(&maxMidx);
		step=(maxM-4)/(4*(maxM-1));
		u_new=(1-step)*u;
		u_new(maxMidx)+=step;
		err=(u-u_new).norm();
		u=u_new;
		count++;		
	}
	//cout<<"finisco l'algoritmo"<<endl;
	NN=u.asDiagonal();
	//cout<<"assegno NN"<<endl;
	Matrix3d PUPt;
	//cout<<"calcolo PUPt"<<endl;
	center=points*u;
	PUPt=points*NN*points.transpose()-center*center.transpose();
	//cout<<"calcolo A"<<endl;
	Matrix3d A;
	A=(PUPt.inverse())/3.0;
	//cout<<"Faccio SVD"<<endl;
	JacobiSVD<MatrixXd> svd(A,ComputeThinV);
	axes=svd.singularValues().cwiseSqrt().cwiseInverse();
	rotation=svd.matrixV();
	//cout<<"Sono arrivato alla fine di kachyan"<<endl;
}

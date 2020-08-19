#ifndef GEOM_LIB_INCLUDED
#define GEOM_LIB_INCLUDED

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <iostream>

using namespace std;
using namespace Eigen;

double cross_prod_2D(double, double, double, double);

void plane_normal(Vector3d, Vector3d, Vector3d, double&, double&, double&);     // riesco a passarli solo per copia dopo vertex.col(i)

bool point_in_polygon_convex(int, int, double*, double*, double, double);

bool point_in_polyhedron();     // DA SCRIVERE, ma per ora non serve

void min_encl_ellipsoid(Matrix3Xd& points, double tol, Vector3d& center, Vector3d& axes, Matrix3d& rotation);

#endif // GEOM_LIB_INCLUDED

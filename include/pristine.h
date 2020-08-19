#ifndef PRISTINE_H
#define PRISTINE_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <math.h>
#include <iostream>

using namespace Eigen;
using namespace std;

class pristine
{
    public:
        pristine();
        virtual ~pristine();

        Vector3d CM;                // column vector of center of mass coordinates
        double vol;
        double Rgyr;

        int N_vertexes;             // number of vertexes
        Matrix3Xd vertex;         // contains vertexes cartesian coordinates as column vector X=N_vertexes

        int N_edges;                // number of edges
        vector<vector<int> > edge; // contains couples of vertex inedexes as row vector

        int N_faces;                // number of faces
        vector<vector<int> > face; // contains vector of vertexes indexes as row vector

        double distance(const pristine*);      //get distance

        void rotate(double, double, double, double);
        void rotate(Matrix3d& R);

        void translate(double, double, double);
        void translate(Vector3d&);

        virtual void radius_of_gyration(double);    // da scrivere, ora scrive 1
        void fill(MatrixX3i&,double);
};

#endif // PRISTINE_H

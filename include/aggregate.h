#ifndef AGGREGATE_H
#define AGGREGATE_H

#include "pristine.h"
#include "column.h"
#include "plate.h"
#include "dendrite.h"
#include "rosette.h"
#include "geom_lib.h"

#include <fstream>
#include <sstream>
#include <string>
#include <time.h>

#include <eigen3/Eigen/Core>

using namespace Eigen;

class aggregate {
  public:
    static constexpr double g=9.81;
    static constexpr double rho_air=0.7364;
    static constexpr double rho_ice=0.917e-15;   /// kg/micron^3
    static constexpr double ni=2.2344e-5;

    aggregate();
    aggregate(int,double);
    virtual ~aggregate();


    double mass;    // mass
    double v;       // terminal fall velocity
    double Rmax;    // maximum radius
    double Rgyr;
    double Maximum_dimension;

    Vector3d CM;    // center of mass

    vector<pristine*> crystal;
    int N;          // number of pristine

    void set_mass();    // calculate mass
    void set_v();       // calculate terminal fall velocity
    void set_Rmax();    // calculate maximum radius of the particle
    void set_Rgyr();
    void set_CM();      // calculate center of mass

    void rotate(double, double, double, double);
    void rotate(Matrix3d&);
    void translate(double, double, double);
    void translate(Vector3d&);
    bool collide(aggregate&);
    void update_mic();      // update all microphysical members
    void reset_CM();        // traslate particle in order to reset CM to (0,0,0)
    void update_maximum_dimension();
    void orient_horizontally();

    void save_vtk_mesh(string);
    void save_aggregate(double);

    protected:
    private:
};

#endif // AGGREGATE_H

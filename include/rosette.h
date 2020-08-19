#ifndef ROSETTE_H
#define ROSETTE_H

#include <pristine.h>


class rosette6b : public pristine
{
    public:
        /** Default constructor */
        rosette6b();
        /** Default destructor */
        virtual ~rosette6b();
        rosette6b(double);

        double f(double, double, double);
        double df(double, double);

        virtual void radius_of_gyration(double);
    protected:
    private:
};

#endif // ROSETTE_H

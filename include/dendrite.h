#ifndef DENDRITE_H
#define DENDRITE_H

#include <pristine.h>


class dendrite : public pristine
{
    public:
        /** Default constructor */
        dendrite();
        /** Default destructor */
        virtual ~dendrite();

        dendrite(double);

        virtual void radius_of_gyration(double);

    protected:
    private:
};

#endif // DENDRITE_H

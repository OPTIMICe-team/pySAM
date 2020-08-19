#ifndef PLATE_H
#define PLATE_H

#include "pristine.h"

class plate : public pristine
{
    public:
        /** Default constructor */
        plate();
        /** Default destructor */
        virtual ~plate();
        /** Copy constructor
         *  \param other Object to copy from
         */
        plate(const plate& other);
        /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */
        plate& operator=(const plate& other);

        plate(double);

        virtual void radius_of_gyration(double);

    protected:
    private:
};

#endif // PLATE_H

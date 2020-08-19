#ifndef COLUMN_H
#define COLUMN_H

#include "pristine.h"

class column : public pristine
{
    public:
        /** Default constructor */
        column();
        /** Default destructor */
        virtual ~column();
        /** Copy constructor
         *  \param other Object to copy from
         */
        column(const column& other);
        /** Assignment operator
         *  \param other Object to assign from
         *  \return A reference to this
         */
        column& operator=(const column& other);

        column(double);

        virtual void radius_of_gyration(double);

    protected:
    private:
};

#endif // COLUMN_H

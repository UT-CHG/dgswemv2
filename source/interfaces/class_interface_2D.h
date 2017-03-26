#ifndef CLASS_INTERFACE_2D_H
#define CLASS_INTERFACE_2D_H

#include "../class_interface.h"

class INTERFACE_2D : public INTERFACE {
public:
    INTERFACE_2D(int, double**, double*, double*, bool, bool boundary = false);
    ~INTERFACE_2D();
};

#endif
#ifndef EHDG_GN_PRE_DBATH_SERIAL_HPP
#define EHDG_GN_PRE_DBATH_SERIAL_HPP

#ifdef D_PROJECTION
#include "grad_methods/projection_dbath_serial.hpp"
#endif

#ifdef D_GREENGAUSS
#include "grad_methods/greengauss_dbath_serial.hpp"
#endif

#ifdef D_LEASTSQUARES
#include "grad_methods/leastsquares_dbath_serial.hpp"
#endif

#endif

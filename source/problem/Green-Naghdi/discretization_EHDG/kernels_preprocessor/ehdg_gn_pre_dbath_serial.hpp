#ifndef EHDG_GN_PRE_DBATH_SERIAL_HPP
#define EHDG_GN_PRE_DBATH_SERIAL_HPP

#ifdef B_RECON_INT
#include "grad_methods/recon_int_dbath.hpp"
#endif

#ifdef B_RECON_LS
#include "grad_methods/recon_ls_dbath.hpp"
#endif

#ifdef B_RECON_AVG
#include "grad_methods/recon_avg_dbath.hpp"
#endif

#ifdef B_PROJECTION
#include "grad_methods/projection_dbath_serial.hpp"
#endif

#ifdef B_GREENGAUSS
#include "grad_methods/greengauss_dbath_serial.hpp"
#endif

#ifdef B_LEASTSQUARES
#include "grad_methods/leastsquares_dbath_serial.hpp"
#endif

#endif

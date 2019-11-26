#ifndef EHDG_GN_PROC_DERIVATIVES_SERIAL_HPP
#define EHDG_GN_PROC_DERIVATIVES_SERIAL_HPP

#ifdef D_RECON_INT
#include "grad_methods/recon_int_derivatives.hpp"
#endif

#ifdef D_RECON_LS
#include "grad_methods/recon_ls_derivatives.hpp"
#endif

#ifdef D_RECON_AVG
#include "grad_methods/recon_avg_derivatives.hpp"
#endif

#ifdef D_PROJECTION
#include "grad_methods/projection_derivatives_serial.hpp"
#endif

#ifdef D_GREENGAUSS
#include "grad_methods/greengauss_derivatives_serial.hpp"
#endif

#ifdef D_LEASTSQUARES
#include "grad_methods/leastsquares_derivatives_serial.hpp"
#endif

#endif

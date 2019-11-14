#ifndef EHDG_GN_PRE_DBATH_OMPI_HPP
#define EHDG_GN_PRE_DBATH_OMPI_HPP

#ifdef D_RECON_INT
#include "grad_methods/recon_int_dbath.hpp"
#endif

#ifdef D_RECON_LS
#include "grad_methods/recon_ls_dbath.hpp"
#endif

#ifdef D_RECON_AVG
#include "grad_methods/recon_avg_dbath.hpp"
#endif

#ifdef D_PROJECTION
#include "grad_methods/projection_dbath_ompi.hpp"
#endif

#ifdef D_GREENGAUSS
#include "grad_methods/greengauss_dbath_ompi.hpp"
#endif

#ifdef D_LEASTSQUARES
#include "grad_methods/leastsquares_dbath_ompi.hpp"
#endif

#endif
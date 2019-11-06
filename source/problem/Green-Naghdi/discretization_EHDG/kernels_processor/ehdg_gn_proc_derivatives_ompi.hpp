#ifndef EHDG_GN_PROC_DERIVATIVES_OMPI_HPP
#define EHDG_GN_PROC_DERIVATIVES_OMPI_HPP

#ifdef D_PROJECTION
#include "grad_methods/projection_derivatives_ompi.hpp"
#endif

#ifdef D_GREENGAUSS
#include "grad_methods/greengauss_derivatives_ompi.hpp"
#endif

#ifdef D_LEASTSQUARES
#include "grad_methods/leastsquares_derivatives_ompi.hpp"
#endif

#endif
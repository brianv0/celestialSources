//
// File: GRBobsConstants.h

#ifndef GRB_OBS_CONSTANTS_H
#define GRB_OBS_CONSTANTS_H 

#include <cmath>


namespace grbcst
{
	const int     ndecs       = 3;
	const int     binsperdec  = 5;
	const int     nbins       = ndecs * binsperdec;
//	const int     ntotal      = 810;
	const int     ntotal      = 667;

	const double  degInPi     = 180.0 / M_PI;
	const double  dynrange    = 3;
	const double  ethres      = 0.03;
	const double  zenNorm     = cos(0) - cos(75. / degInPi);
	const double  minwid      = 0.01;
	const double  maxwid      = 10.0;
	const double  logfac0     = pow((maxwid/minwid), (1./nbins));
	const double  yr_frac     = 1.0;
	const double  alpha       = 1.0;
	const double  emax        = 100.;
	const double  logdyn      = log10(dynrange);
	const double  nu          = 1.5;
	const double  timres      = 0.001;
	const double  frackeep    = 0.01;

	const long    seed        = 83498275L;
	const long    nglast_yr   = long(ntotal * (zenNorm/2));
	const long    nbsim       = long(yr_frac * nglast_yr + 0.5);	// number of bursts

}


namespace grbobstypes
{
	 typedef std::vector<long>::iterator LongIter;
	 typedef std::vector<double>::iterator DoubleIter;
	 typedef std::vector<int>::iterator IntIter;
	 typedef std::vector<int>::size_type IntSize;
	 typedef std::vector<long>::const_iterator LongConstIter;
	 typedef std::vector<double>::const_iterator DoubleConstIter;
	 typedef std::vector<double>::size_type DoubleSize;
	 typedef std::vector<long>::size_type LongSize;
}

#endif



/**
 * @file test.cxx
 * @brief Simple program to exercise EblAtten class.  Real unit tests are
 * forthcoming.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */

#include <cmath>

#include <iostream>

#include "eblAtten/EblAtten.h"

int main() {

   IRB::EblAtten tau0(IRB::Kneiske), tau1(IRB::Primack05),
      tau2(IRB::Kneiske_HighUV), tau3(IRB::Salamon_Stecker);

   double emin(1e4), emax(5e5);
   int npts(50);

   double z [5] = {0., 0.5, 1., 2., 4.};

   double estep(std::log(emax/emin)/float(npts));

  for (int j = 0; j < 5; j++)
   for (int i = 0; i < npts; i++) {
      double energy(emin*std::exp(i*estep));
      std::cout << energy/1e3 << " GeV   "<< "  "
                << "redshift = "<< z[j] << "   Tau = "
                << tau0(energy, z[j]) << "  "
                << tau1(energy, z[j]) << "  "
                << tau2(energy, z[j]) << "  "
		<< tau3(energy, z[j]) << std::endl;
   }
}

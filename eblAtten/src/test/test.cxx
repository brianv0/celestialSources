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

   IRB::EblAtten tau0(IRB::Kneiske), tau1(IRB::Salamon_Stecker),
      tau2(IRB::Primack_Bullock99);

   double emin(1e4), emax(5e5);
   int npts(50);
   
   double z(0.1);

   double estep(std::log(emax/emin)/float(npts));
   for (int i = 0; i < npts; i++) {
      double energy(emin*std::exp(i*estep));
      std::cout << energy << "  " 
                << tau0(energy, z) << "  "
                << tau1(energy, z) << "  "
                << tau2(energy, z) << std::endl;
   }
}

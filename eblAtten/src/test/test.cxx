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
      tau2(IRB::Kneiske_HighUV), tau3(IRB::Stecker05), tau4(IRB::Finke), tau5(IRB::Franceschini), tau6(IRB::Gilmore),
      tau7(IRB::Stecker05_FE), tau8(IRB::SalamonStecker), tau9(IRB::Generic);

   double emin(1e4), emax(5e5);
   int npts(50);

   double z [5] = {0., 0.5, 1., 2., 4.};

   double estep(std::log(emax/emin)/float(npts));

  for (int j = 0; j < 5; j++)
   for (int i = 0; i < npts; i++) {
      double energy(emin*std::exp(i*estep));
      std::cout << energy/1e3 << " GeV   "<< "  "
                << "redshift = "<< z[j] << "   Tau = "
                << tau0(energy, z[j]) << " Kneiske BF |  "
                << tau1(energy, z[j]) << " Primack 05 | "
                << tau2(energy, z[j]) << " Kneiske UV |  "
                << tau3(energy, z[j]) << " Stecker 05 |  "
                << tau4(energy, z[j]) << " Finke  |  "
                << tau5(energy, z[j]) << " Franceschini  | "
		<< tau6(energy, z[j]) << " Gilmore  |" 
                << tau7(energy, z[j]) << " Stecker05 - FE |"
                << tau8(energy, z[j]) << " Salamon & Stecker 98  |"
                << tau9(energy, z[j]) << " Generic " << std::endl;
   }
}

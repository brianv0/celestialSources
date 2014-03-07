/**
 * @file test.cxx
 * @brief Simple program to exercise EblAtten class.  Real unit tests are
 * forthcoming.
 * @author J. Chiang <jchiang@slac.stanford.edu>
 *
 * $Header$
 */
#ifdef TRAP_FPE
#include <fenv.h>
#endif

#include <cmath>

#include <iostream>

#include "facilities/commonUtilities.h"

#include "eblAtten/EblAtten.h"

int main() {
#ifdef TRAP_FPE
   feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

   facilities::commonUtilities::setupEnvironment();

   IRB::EblAtten tau0(IRB::Kneiske), tau1(IRB::Primack05),
      tau2(IRB::Kneiske_HighUV), tau3(IRB::Stecker05), 
      tau4(IRB::Franceschini), tau5(IRB::Finke), tau6(IRB::Gilmore),
      tau7(IRB::Stecker05_FE), tau8(IRB::SalamonStecker), tau9(IRB::Generic),
      tau10(IRB::Dominguez11), tau11(IRB::KneiskeDole10), tau12(IRB::KneiskeDole10_CMB),
      tau13(IRB::Gilmore09), tau14(IRB::Gilmore12_fiducial), tau15(IRB::Gilmore12_fixed), 
      tau16(IRB::Scully14_lowOp), tau17(IRB::Scully14_highOp), tau18(IRB::Inoue13);

   double emin(1e4), emax(5e5);
   int npts(50);

   double z [5] = {0., 0.5, 1., 2., 4.};

   double estep(std::log(emax/emin)/float(npts));


   std::cout << "Energy (GeV)  Tau  KneiskeBF  Primack05  KneiskeUV  "
             << "Stecker05  Finke  Franceschini  Gilmore  Stecker05_FE  "
             << "Salamon&Stecker98  Generic  Gilmore12_fixed  "
             << "Gilmore12_fiducial Inoue13 Dominguez11 "
             << "Scully14_highOp  Scully14_lowOp KneiskeDole10 KneiskeDole10_CMB\n";
   for (int j = 0; j < 5; j++) {
      for (int i = 0; i < npts; i++) {
         double energy(emin*std::exp(i*estep));
         std::cout << energy/1e3 << "  "
                   << z[j] << "  "
                   << tau0(energy, z[j]) << "  "
                   << tau1(energy, z[j]) << "  "
                   << tau2(energy, z[j]) << "  "
                   << tau3(energy, z[j]) << "  "
                   << tau4(energy, z[j]) << "  "
                   << tau5(energy, z[j]) << "  "
                   << tau6(energy, z[j]) << "  " 
                   << tau7(energy, z[j]) << "  "
                   << tau8(energy, z[j]) << "  "
                   << tau9(energy, z[j]) << "  " 
                   << tau10(energy, z[j]) << "  " 
                   << tau11(energy, z[j]) << "  "
                   << tau12(energy, z[j]) << "  "
                   << tau13(energy, z[j]) << "  "
                   << tau14(energy, z[j]) << "  "
                   << tau15(energy, z[j]) << "  "
                   << tau16(energy, z[j]) << "  "
                   << tau17(energy, z[j]) << "  "
                   << tau18(energy, z[j]) << "  "
                   << std::endl;
      }
   }
}

/**
 * @file TestUtil.h
 * @brief Test code for Util static functions.
 * @author J. Chiang
 * 
 * $Header$
 */

#include "Util.h"

bool assert_equals(double x1, double x2, double tol=1e-6) {
   return (std::fabs(x1 - x2)/(std::fabs(x1) + std::fabs(x2)) < tol);
}

bool Util_tests() {
   std::vector<double> x;
   std::vector<double> y;
   x.push_back(1);
   x.push_back(2);
   x.push_back(3);
   y.push_back(1);
   y.push_back(4);
   y.push_back(13.5);

   double value(genericSources::Util::logInterpolate(x, y, 1.5));
   if (!assert_equals(value, 2.25)) {
      std::cout << "logInterpolate fails for x = " << 1.5 << "\n"
                << "value = " << value << std::endl;
      return false;
   }

   value = genericSources::Util::logInterpolate(x, y, 2.5);
   if (!assert_equals(value, 7.8125)) {
      std::cout << "logInterpolate fails for x = " << 2.5 << "\n"
                << "value = " << value << std::endl;
      return false;
   }

   value = genericSources::Util::powerLawIntegral(1, 2, 1, 0.5);
   if (!assert_equals(value, std::log(2.))) {
      std::cout << "powerLawIntegral fails for x1, x2, y1, y2 = 1, 2, 1, 0.5"
                << std::endl;
      return false;
   }

   value = genericSources::Util::powerLawIntegral(2, 3, 4, 9);
   if (!assert_equals(value, 19./3.)) {
      std::cout << "powerLawIntegral fails for x1, x2, y1, y2 = 1, 2, 1, 4"
                << std::endl;
      return false;
   }
   return true;

}

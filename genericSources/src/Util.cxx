/**
 * @file Util.cxx
 * @brief
 * @author J. Chiang
 *
 * $Header$
 */

#include <cmath>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoisson.h"

#include "facilities/Util.h"

#include "Util.h"

namespace genericSources {

   bool Util::fileExists(const std::string & filename) {
      std::ifstream file(filename.c_str());
      return file.is_open();
   }

   void Util::file_ok(std::string filename) {
      facilities::Util::expandEnvVar(&filename);
      if (fileExists(filename)) {
         return;
      } else {
         throw std::runtime_error("File not found: " + filename);
      }
   }

   void Util::readLines(std::string inputFile, 
                        std::vector<std::string> &lines,
                        const std::string &skip) {
      facilities::Util::expandEnvVar(&inputFile);
      std::ifstream file(inputFile.c_str());
      lines.clear();
      std::string line;
      while (std::getline(file, line, '\n')) {
         if (line != "" && line != " "             //skip (most) blank lines 
             && line.find_first_of(skip) != 0) {   //and commented lines
            lines.push_back(line);
         }
      }
   }

   void Util::resolve_fits_files(std::string filename, 
                                 std::vector<std::string> &files) {
      facilities::Util::expandEnvVar(&filename);
      files.clear();
// Read the first line of the file and see if the first 6 characters
// are "SIMPLE".  If so, then we assume it's a FITS file.
      std::ifstream file(filename.c_str());
      std::string firstLine;
      std::getline(file, firstLine, '\n');
      if (firstLine.find("SIMPLE") == 0) {
// This is a FITS file. Return that as the sole element in the files
// vector.
         files.push_back(filename);
         return;
      } else {
// filename contains a list of fits files.
         readLines(filename, files);
         return;
      }
   }

   bool Util::isXmlFile(std::string filename) {
      std::vector<std::string> tokens;
      facilities::Util::stringTokenize(filename, ".", tokens);
      if (*(tokens.end()-1) == "xml") {
         return true;
      }
      return false;
   }

   double Util::interpolate(const std::vector<double> &x,
                            const std::vector<double> &y,
                            double xx) {
      if (xx < x.front() || xx > x.back()) {
         std::ostringstream message;
         message << "Util::interpolate:\n"
                 << "abscissa value out-of-range, "
                 << xx << " is not in (" 
                 << x.front() << ", "
                 << x.back() << ")";
         throw std::range_error(message.str());
      }
      std::vector<double>::const_iterator it 
         = std::upper_bound(x.begin(), x.end(), xx) - 1;
      unsigned int indx = it - x.begin();
      double yy;
      if (*(it+1) != *it) {
         yy = (xx - *it)/(*(it+1) - *it)*(y[indx+1] - y[indx]) + y[indx];
      } else {
         yy = (y[indx+1] + y[indx])/2.;
      }
      return yy;
   }

   double Util::bilinear(const std::vector<double> &xx, double x, 
                         const std::vector<double> &yy, double y, 
                         const std::vector<double> &z) {

      std::vector<double>::const_iterator ix;
      if (x < *(xx.begin())) {
         ix = xx.begin() + 1;
      } else if (x >= *(xx.end()-1)) {
         ix = xx.end() - 1;
      } else {
         ix = std::upper_bound(xx.begin(), xx.end(), x);
      }
      int i = ix - xx.begin();
      
      std::vector<double>::const_iterator iy;
      if (y < *(yy.begin())) {
         iy = yy.begin() + 1;
      } else if (y >= *(yy.end()-1)) {
         iy = yy.end() - 1;
      } else {
         iy = std::upper_bound(yy.begin(), yy.end(), y);
      }
      int j = iy - yy.begin();
      
      double tt = (x - *(ix-1))/(*(ix) - *(ix-1));
      double uu = (y - *(iy-1))/(*(iy) - *(iy-1));
      
      double y1 = z[yy.size()*(i-1) + (j-1)];
      double y2 = z[yy.size()*(i) + (j-1)];
      double y3 = z[yy.size()*(i) + (j)];
      double y4 = z[yy.size()*(i-1) + (j)];

      double value = (1. - tt)*(1. - uu)*y1 + tt*(1. - uu)*y2 
         + tt*uu*y3 + (1. - tt)*uu*y4; 
      if (value < 0.) {
         std::ostringstream message;
         message << "irfUtil::Util::bilinear:\n"
                 << "value = " << value << " < 0\n";
         message << xx[i-1] << "  " << *(ix-1) << "  " 
                 << x << "  " << *ix << "\n";
         message << yy[j-1] << "  " << *(iy-1) << "  " 
                 << y << "  " << *iy << "\n";
         message << tt << "  " << uu << "  " 
                 << y1 << "  " << y2 << "  "
                 << y3 << "  " << y4;
         throw std::runtime_error(message.str());
      }
      return value;
   }

   bool Util::expectedException(const std::exception & eObj, 
                                const std::string & targetMessage) {
      std::string message(eObj.what());
      return message.find_first_of(targetMessage.c_str()) 
         != std::string::npos;
   }

   double Util::drawFromPowerLaw(double emin, double emax, double gamma) {
      double xi = RandFlat::shoot();
      double one_m_gamma = 1. - gamma;
      double arg = xi*(std::pow(emax, one_m_gamma) - 
                       std::pow(emin, one_m_gamma)) 
         + std::pow(emin, one_m_gamma);
      double energy = std::pow(arg, 1./one_m_gamma);
      return energy;
   }

} // namespace genericSources

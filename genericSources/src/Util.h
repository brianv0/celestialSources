/**
 * @file Util.h
 * @brief Some basic utility functions.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef genericSources_Util_h
#define genericSources_Util_h

#include <string>
#include <vector>

#include "facilities/Util.h"

namespace genericSources {

/**
 * @class Util
 * @brief Various static functions of general use for Science Tools 
 * applications.
 *
 * @author J. Chiang
 *
 * $Header$
 */

class Util {

public:

   /// @return true if the file exists.  Environment variables are 
   ///         *not* expanded.
   static bool fileExists(const std::string & filename);

   /// This expands any environment variables in filename and checks
   /// if the file exists.  If it doesn't, a runtime_error exception
   /// is thrown.  Otherwise, it does nothing.
   static void file_ok(std::string filename);

   /// @brief Read lines, separated by "\n", from a file.  Blank lines
   ///        are skipped.
   /// @param inputFile file to be read; environment variables are expanded
   /// @param lines On return, this vector is filled with each line read 
   ///        from the file.
   /// @param skip The comment string. Lines beginning with this string are
   ///        not put into lines.
   static void readLines(std::string inputFile, 
                         std::vector<std::string> &lines,
                         const std::string &skip = "#");

   /// @brief Determine if a file is a FITS file by looking for the "SIMPLE"
   ///        keyword as the first six characters of the file.  If it is
   ///        not a FITS file, then it is assumed to be a list if FITS files.
   /// @param filename The name of the candidate file; enviroment 
   ///        variables are expanded.
   /// @param files If filename is a FITS file, this is filled with
   ///        the name of that file; otherwise, it is filled with the 
   ///        "\n"-separated lines in the file.
   static void resolve_fits_files(std::string filename, 
                                  std::vector<std::string> &files);

   /// @return true if the filename ends in ".xml" extension
   static bool isXmlFile(std::string filename);

   /// Linear interpolation.
   static double interpolate(const std::vector<double> &x,
                             const std::vector<double> &y,
                             double xx);

   /// A zeroth order bilinear interpolater.
   static double bilinear(const std::vector<double> &xx, double x,
                          const std::vector<double> &yy, double y, 
                          const std::vector<double> &z);

   /// @return true if eObj.what() contains the targetMessage
   ///         as a substring.
   static bool expectedException(const std::exception & eObj,
                                 const std::string & targetMessage);

   static double drawFromPowerLaw(double emin, double emax, double gamma);

};

} // namespace genericSources

#endif // genericSources_Util_h

/** 
 * @file FitsImage.h
 * @brief Declaration of FitsImage class
 * @authors J. Chiang
 *
 * $Header$
 *
 */

#ifndef genericSources_FitsImage_h
#define genericSources_FitsImage_h

#include <sstream>
#include <string>
#include <vector>

#include "fitsio.h"

#include "astro/SkyDir.h"

namespace genericSources {

/** 
 * @class FitsImage
 *
 * @brief A class for accessing FITS image data.
 *
 * @author J. Chiang
 *    
 * $Header$
 *
 */

class FitsImage {
    
public:

   FitsImage() {}

   FitsImage(const std::string & filename);

   virtual ~FitsImage() {}

   /// @return The image value at the SkyDir location. Zero if outside
   ///         the domain of the map.  
   /// @param dir The sky location for which the map value is returned.
   /// @param iz For 3D images, this is the index of the image plane
   ///        to be sampled.  Ignored for 2D maps.
   double operator()(const astro::SkyDir & dir, size_t iz=0) const;

   /// A vector of the image axes dimensions
   virtual void getAxisDims(std::vector<int> & axisDims);

   /// The names (CTYPEs) of the image axes
   virtual void getAxisNames(std::vector<std::string> & axisNames);

   /// Get a vector filled with axis abscissa points for the naxis-th
   /// coordinate.
   virtual void getAxisVector(unsigned int naxis, 
                              std::vector<double> & axisVector) const;

   /// This method computes arrays of longitude and latitude obtained
   /// by traversing the image plane by column number then row.
   virtual void getCelestialArrays(std::vector<double> & lonArray,
                                   std::vector<double> & latArray);
   
   /// Get the pixel values.  They will be indexed by column, row,
   /// then plane, i.e., indx = i + j*NAXIS1 + k*NAXIS1*NAXIS2.  Note
   /// that each image plane is indexed starting at the lower left
   /// (South-East) corner.
   virtual void getImageData(std::vector<double> & imageData);

   /// This returns the pixel solid angles.  Use of this method assumes
   /// that m_axis[0] represents a longitudinal coordinate and that
   /// m_axis[1] represents a latitudinal coordinate.  The pixel values
   /// will be indexed by column then row, indx = i + j*NAXIS1.
   virtual void getSolidAngles(std::vector<double> & solidAngles) const;

   /// @return The integral over solid angle of the map.
   double mapIntegral() const;

   /// @return The coordinate system used (Galactic or Equatorial)
   astro::SkyDir::CoordSystem coordSys() const {
      return m_coordSys;
   }

   /// @return The HDU number of the specified extension.
   static int findHdu(const std::string & fitsFile,
                      const std::string & extension);

   template <typename Functor>
   static FitsImage sampledImage(const Functor & functor, 
                                 const std::vector<double> & longitudes,
                                 const std::vector<double> & latitudes,
                                 astro::SkyDir::CoordSystem coordSys);

#ifndef SWIG
   static void readColumn(fitsfile * fptr, const std::string & colname,
                          std::vector<double> & coldata);

   static void readRowVector(fitsfile * fptr, const std::string & colname,
                             int row, int nelements,
                             std::vector<double> & data);
#endif // SWIG

   static void fitsReportError(int status, std::string routine="");

protected:

/** 
 * @class AxisParams
 * @brief Nested class to represent FITS image axis information
 */
   class AxisParams {
   public:
      AxisParams() {}
      AxisParams(const std::vector<double> & coords,
                 const std::string & axis_type) : axisType(axis_type) {
         size = coords.size();
         refVal = coords.at(0);
         step = coords.at(1) - coords.at(0);
         refPixel = 1.;
         logScale = false;
      }
      ~AxisParams() {}
      int size;
      float refVal;
      float step;
      float refPixel;
      std::string axisType;
      std::string comment;
      bool logScale;

      /// Returns a vector of abscissa values based on the axis parameters.
      void computeAxisVector(std::vector<double> &axisVector);
   };

   /// Interface to cfitsio routines.
   void read_fits_image(std::string &filename, std::vector<AxisParams> &axes,
                        std::vector<double> &image);

   /// FITS file name.
   std::string m_filename;

   /// Descriptions for each image axis.
   std::vector<AxisParams> m_axes;

   /// Vectors of abscissa values for each axis.
   std::vector< std::vector<double> > m_axisVectors;

   /// The FITS image data
   std::vector<double> m_image;

   astro::SkyDir::CoordSystem m_coordSys;

   static std::string s_fitsRoutine;

};

template <typename Functor>
FitsImage FitsImage::sampledImage(const Functor & image, 
                                  const std::vector<double> & longitudes,
                                  const std::vector<double> & latitudes,
                                  astro::SkyDir::CoordSystem coordSys) {
   FitsImage my_image;

   my_image.m_filename = "from Functor";

   my_image.m_coordSys = coordSys;

   std::string lonType("GLON-CAR");
   std::string latType("GLAT-CAR");

   if (coordSys == astro::SkyDir::EQUATORIAL) {
      lonType = "RA---CAR";
      latType = "DEC--CAR";
   }

   my_image.m_axes.push_back(AxisParams(longitudes, lonType));
   my_image.m_axes.push_back(AxisParams(latitudes, latType));

   my_image.m_axisVectors.push_back(longitudes);
   my_image.m_axisVectors.push_back(latitudes);

   size_t nz(1);
/// @bug Using image.m_axes and indexing over k breaks polymorphism.
   if (image.m_axes.size() == 3) {
      my_image.m_axes.push_back(image.m_axes.at(2));
      my_image.m_axisVectors.push_back(image.m_axisVectors.at(2));
      nz = my_image.m_axes.at(2).size;
   }
   for (size_t k = 0; k < nz; k++) {
      for (size_t j = 0; j < latitudes.size(); j++) {
         for (size_t i = 0; i < longitudes.size(); i++) {
            astro::SkyDir dir(longitudes.at(i), latitudes.at(j), coordSys);
            my_image.m_image.push_back(image(dir, k));
         }
      }
   }
   return my_image;
}

} // namespace genericSources

#endif // genericSources_FitsImage.h

/** 
 * @file FitsImage.cxx
 * @brief Implementation of FitsImage member functions
 * @authors J. Chiang
 *
 * $Header$
 *
 */

#include <cassert>
#include <cmath>

#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "FitsImage.h"

namespace genericSources {

#include "fitsio.h"

FitsImage::FitsImage(const std::string &fitsfile) {
   m_filename = fitsfile;
   read_fits_image(m_filename, m_axes, m_image);
   for (unsigned int i = 0; i < m_axes.size(); i++) {
      std::vector<double> axisVector;
      m_axes[i].computeAxisVector(axisVector);
      m_axisVectors.push_back(axisVector);
   }
}

void FitsImage::getAxisDims(std::vector<int> &axisDims) {
   axisDims.clear();
   for (unsigned int i = 0; i < m_axes.size(); i++) {
      axisDims.push_back(m_axes[i].size);
   }
}

void FitsImage::getAxisNames(std::vector<std::string> &axisNames) {   
   axisNames.clear();
   for (unsigned int i = 0; i < m_axes.size(); i++) {
      axisNames.push_back(m_axes[i].axisType);
   }
}

void FitsImage::getAxisVector(unsigned int naxis,
                              std::vector<double> &axisVector) {
   if (naxis >= m_axes.size()) {
      std::ostringstream message;
      message << "FitsImage::getAxisVector: Invalid axis number " << naxis;
      throw std::invalid_argument(message.str());
   }
   axisVector = m_axisVectors[naxis];
}

void FitsImage::getCelestialArrays(std::vector<double> &lonArray,
                                   std::vector<double> &latArray) {
   unsigned int npixels = m_axes[0].size*m_axes[1].size;
   lonArray.resize(npixels);
   latArray.resize(npixels);
   for (int j = 0; j < m_axes[1].size; j++) {
      for (int i = 0; i < m_axes[0].size; i++) {
         int indx = i + j*m_axes[0].size;
         lonArray[indx] = m_axisVectors[0][i];
         latArray[indx] = m_axisVectors[1][j];
      }
   }
}
         
void FitsImage::getSolidAngles(std::vector<double> &solidAngles){
// This solid angle calculation *assumes* that m_axes[0] is a
// longitudinal coordinate and that m_axes[1] is a latitudinal one.
// Furthermore, the axis units are assumed to be degrees, while the
// solid angles are returned as steradians.

   solidAngles.resize(m_axes[0].size*m_axes[1].size);
   for (int i = 0; i < m_axes[0].size; i++) {
      for (int j = 0; j < m_axes[1].size; j++) {
         int indx = i + j*m_axes[0].size;
         double thetamin = (m_axisVectors[1][j] - m_axes[1].step/2.)*M_PI/180;
         double thetamax = (m_axisVectors[1][j] + m_axes[1].step/2.)*M_PI/180;
         solidAngles[indx] = m_axes[0].step*M_PI/180
            *(sin(thetamax) - sin(thetamin));
      }
   }
}

void FitsImage::getImageData(std::vector<double> &imageData) {
   imageData.resize(m_image.size());
   imageData = m_image;
}

void FitsImage::
AxisParams::computeAxisVector(std::vector<double> &axisVector) {
   axisVector.clear();
   axisVector.reserve(size);
   for (int i = 0; i < size; i++) {
      double value = step*(i - refPixel + 1) + refVal;
      if (logScale) value = exp(value);
      axisVector.push_back(value);
   }
}

void FitsImage::read_fits_image(std::string &filename, 
                                std::vector<AxisParams> &axes,
                                std::vector<double> &image) {
   fitsfile * fptr = 0;
   char *file = const_cast<char *>(filename.c_str());
   int status = 0;

   fits_open_file(&fptr, file, READONLY, &status);
   fits_report_error(stderr, status);
   if (status != 0) {
      throw std::runtime_error("FitsImage::read_fits_image: cfitsio error");
   }

// Get dimensions of the data cube
   long naxes;
   char comment[80];
   fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status);
   fits_report_error(stderr, status);
   if (status != 0) {
      throw std::runtime_error("FitsImage::read_fits_image: cfitsio error");
   }

// Assume at least 1 image plane, but at most 3 dimensions...
   if (naxes != 2 && naxes != 3) {
      std::ostringstream errorMessage;
      errorMessage << "FitsImage::read_fits_image: \n"
                   << "FITS file " << filename 
                   << " does not have the expected number of dimensions:"
                   << " naxes = " << naxes << "\n";
      throw std::runtime_error(errorMessage.str());
   }

// prepare the axes vector 
   axes.clear();
   axes.resize(naxes);

// keyword names
   char *naxis[] = {"NAXIS1", "NAXIS2", "NAXIS3"};
   char *crval[] = {"CRVAL1", "CRVAL2", "CRVAL3"};
   char *cdelt[] = {"CDELT1", "CDELT2", "CDELT3"};
   char *crpix[] = {"CRPIX1", "CRPIX2", "CRPIX3"};
   char *ctype[] = {"CTYPE1", "CTYPE2", "CTYPE3"};

   long ivalue;
   double value;
   char svalue[40];

   long npixels = 1;
   for (int i = 0; i < naxes; i++) {
// axis size
      fits_read_key_lng(fptr, naxis[i], &ivalue, comment, &status);
      fits_report_error(stderr, status);
      if (status != 0) {
         throw std::runtime_error
            ("FitsImage::read_fits_image: cfitsio error");
      }
      axes[i].size = ivalue;

// Compute the number of pixels in the image.
      npixels *= ivalue;
   }

// account for degenerate case of NAXIS3 = 1

   if (naxes == 3 && axes[2].size == 1) {
      naxes = 2;
      axes.resize(naxes);
   }

   for (int i = 0; i < naxes; i++) {
// reference values
      fits_read_key_dbl(fptr, crval[i], &value, comment, &status);
      fits_report_error(stderr, status);
      if (status != 0) {
         throw std::runtime_error
            ("FitsImage::read_fits_image: cfitsio error");
      }
      axes[i].refVal = value;

// step sizes
      fits_read_key_dbl(fptr, cdelt[i], &value, comment, &status);
      fits_report_error(stderr, status);
      if (status != 0) {
         throw std::runtime_error
            ("FitsImage::read_fits_image: cfitsio error");
      }
      axes[i].step = value;

// reference pixels
      fits_read_key_lng(fptr, crpix[i], &ivalue, comment, &status);
      fits_report_error(stderr, status);
      if (status != 0) {
         throw std::runtime_error
            ("FitsImage::read_fits_image: cfitsio error");
      }
      axes[i].refPixel = ivalue;

// axis types and commentary
      fits_read_key_str(fptr, ctype[i], svalue, comment, &status);
      fits_report_error(stderr, status);
      if (status != 0) {
         throw std::runtime_error
            ("FitsImage::read_fits_image: cfitsio error");
      }
      axes[i].axisType = svalue;
      axes[i].comment = comment;
      
// Check for logarithmic scaling.
      int offset = axes[i].axisType.substr(0).find_first_of("log_");
      if (offset == 0) {
         axes[i].logScale = true;
      } else {
         axes[i].logScale = false;
      }
   } // naxes

   status = 0;

// Read in the image pixels.
   long group = 0;
   long fpixel = 1;
   double nullval = 0.;
   int anynull;
   double *tmpImage;
   tmpImage = new double[npixels];
   fits_read_img_dbl(fptr, group, fpixel, npixels, nullval, 
                     tmpImage, &anynull, &status);
   fits_report_error(stderr, status);
   if (status != 0) {
      throw std::runtime_error("FitsImage::read_fits_image: cfitsio error");
   }

   image.resize(npixels);

   for (int i = 0; i < npixels; i++)
      image[i] = tmpImage[i];

   delete [] tmpImage;

   fits_close_file(fptr, &status);
   fits_report_error(stderr, status);
   if (status != 0) {
      throw std::runtime_error("FitsImage::read_fits_image: cfitsio error");
   }
}

} // namespace genericSources

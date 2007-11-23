/** 
 * @file FitsImage.cxx
 * @brief Implementation of FitsImage member functions
 * @authors J. Chiang
 *
 * $Header$
 *
 */

#include <cmath>

#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "FitsImage.h"

namespace genericSources {

#include "fitsio.h"

std::string FitsImage::s_fitsRoutine("");

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
         solidAngles[indx] = std::fabs(m_axes[0].step*M_PI/180
            *(sin(thetamax) - sin(thetamin)));
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
   s_fitsRoutine = "read_fits_image";
   fitsfile * fptr = 0;
   char *file = const_cast<char *>(filename.c_str());
   int status = 0;

   fits_open_file(&fptr, file, READONLY, &status);
   fitsReportError(status);

// Get dimensions of the data cube
   long naxes;
   char comment[80];
   fits_read_key_lng(fptr, "NAXIS", &naxes, comment, &status);
   fitsReportError(status);

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
      fitsReportError(status);
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
      fitsReportError(status);
      axes[i].refVal = value;

// step sizes
      fits_read_key_dbl(fptr, cdelt[i], &value, comment, &status);
      fitsReportError(status);
      axes[i].step = value;

// reference pixels
      fits_read_key_lng(fptr, crpix[i], &ivalue, comment, &status);
      fitsReportError(status);
      axes[i].refPixel = ivalue;

// axis types and commentary
      fits_read_key_str(fptr, ctype[i], svalue, comment, &status);
      fitsReportError(status);
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
   fitsReportError(status);

   image.resize(npixels);

   for (int i = 0; i < npixels; i++)
      image[i] = tmpImage[i];

   delete [] tmpImage;

   fits_close_file(fptr, &status);
   fitsReportError(status);
}

int FitsImage::findHdu(const std::string & fitsFile, 
                       const std::string & extension) {
   s_fitsRoutine = "findHdu";
   
   int status(0);
   fitsfile * fptr = 0;

   fits_open_file(&fptr, fitsFile.c_str(), READONLY, &status);
   fitsReportError(status);
   
   int nhdus;
   fits_get_num_hdus(fptr, &nhdus, &status);
   fitsReportError(status);

   int hdutype(0);
   char extname[20];
   char comment[72];
   for (int hdu = 1; hdu < nhdus+1; hdu++) {
      fits_movabs_hdu(fptr, hdu, &hdutype, &status);
      fitsReportError(status);
      
      fits_read_key_str(fptr, "EXTNAME", extname, comment, &status);
      if (status == 202) {
         status = 0;
         continue;
      } else {
         fitsReportError(status);
      }
      
      if (extension == extname) {
         fits_close_file(fptr, &status);
         fitsReportError(status);
         return hdu;
      }
   }
   fits_close_file(fptr, &status);
   fitsReportError(status);

   std::ostringstream message;
   message << "FitsImage::findHdu: HDU number not found for file "
           << fitsFile << " and extension " << extension;
   throw std::runtime_error(message.str());
   return -1;
}

void FitsImage::readColumn(fitsfile * fptr, const std::string & colname,
                           std::vector<double> & coldata) {
   std::string routineName("FitsImage::readColumn");
   int status(0);
   int colnum(0);
   fits_get_colnum(fptr, CASEINSEN, const_cast<char *>(colname.c_str()),
                   &colnum, &status);
   fitsReportError(status, routineName);

   long nrows(0);
   fits_get_num_rows(fptr, &nrows, &status);
   fitsReportError(status, routineName);

   int anynul(0), nulval(0);
   coldata.resize(nrows);
   fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, &nulval, &coldata[0],
                 &anynul, &status);
   fitsReportError(status, routineName);
}

void FitsImage::
readRowVector(fitsfile * fptr, const std::string & colname, 
              int row, int nelements, std::vector<double> & data) {
   std::string routineName("FitsImage::readRowVector");
   int status(0);
   int colnum(0);
   fits_get_colnum(fptr, CASEINSEN, const_cast<char *>(colname.c_str()),
                   &colnum, &status);
   fitsReportError(status, routineName);

   long nrows(0);
   fits_get_num_rows(fptr, &nrows, &status);
   fitsReportError(status, routineName);
   if (row >= nrows) {
      std::ostringstream message;
      message << routineName << "\n"
              << "Request for row " << row << " from a table that "
              << "has only " << nrows << "rows.";
      throw std::runtime_error(message.str());
   }

   int anynul(0), nulval(0);
   data.resize(nelements);
   fits_read_col(fptr, TDOUBLE, colnum, row+1, 1, nelements, &nulval,
                 &data[0], &anynul, &status);
   fitsReportError(status, routineName);
}

void FitsImage::fitsReportError(int status, std::string routine) {
   if (status == 0) {
      return;
   }
   if (routine == "") {
      routine = "FitsImage::" + s_fitsRoutine;
   }
   fits_report_error(stderr, status);
   std::ostringstream message;
   message << routine << ": CFITSIO error " << status;
   throw std::runtime_error(message.str());
}

} // namespace genericSources

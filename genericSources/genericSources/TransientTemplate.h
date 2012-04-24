/**
 * @file TransientTemplate.h
 * @brief A flaring source whose light curve is given by a template file.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef mySpectrum_TransientTemplate_h
#define mySpectrum_TransientTemplate_h

#include "SimpleTransient.h"

/**
 * @class TransientTemplate
 *
 * @brief A flaring source whose light curve shape is given by a
 * template file.  The duration, mean flux, and spectral index are
 * given as parameters.
 *
 */

class TransientTemplate : public SimpleTransient {

public:

   TransientTemplate(const std::string & params);
   
   virtual ~TransientTemplate() {}

   static double drawTime(const std::vector<double> & tt,
                          const std::vector<double> & integralDist);


   virtual std::pair<double, double> dir(double) {
      return std::make_pair(0, 0);
   }


private:

   void createEventTimes(std::string templateFile);

};

#endif // mySpectrum_TransientTemplate_h

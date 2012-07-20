/**
 * @file ConstParMap.h
 * @brief const interface to std::map<std::string, std::string>.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef celestialSources_ConstParMap_h
#define celestialSources_ConstParMap_h

#include <map>
#include <stdexcept>

#include "facilities/Util.h"

namespace celestialSources {

class ConstParMap {

public:

   ConstParMap(const std::map<std::string, std::string> & parmap) 
      : m_parmap(parmap) {
      check_flux();
   }

   ConstParMap(const std::string & params) {
      facilities::Util::keyValueTokenize(params, ", ", m_parmap);
      check_flux();
   }

   const std::string & operator[](const std::string & name) const {
      std::map<std::string, std::string>::const_iterator item 
         = m_parmap.find(name);
      if (item == m_parmap.end()) {
         throw std::runtime_error("Cannot find item named " + name);
      }
      return item->second;
   }

   double value(const std::string & name) const {
      return std::atof(operator[](name).c_str());
   }

   size_t size() const {
      return m_parmap.size();
   }

   std::map<std::string, std::string>::const_iterator 
   find(const std::string & parname) {
      return m_parmap.find(parname);
   }

   std::map<std::string, std::string>::const_iterator end() {
      return m_parmap.end();
   }

private:

   std::map<std::string, std::string> m_parmap;

   // Ensure that any parameter named "flux" is positive since
   // (amazingly enough) there are users who will build xml model defs
   // with negative flux values specified.
   void check_flux() const {
      std::map<std::string, std::string>::const_iterator item 
         = m_parmap.find("flux");
      if (item != m_parmap.end()) {
         double value(std::atof(item->second.c_str()));
         if (value < 0) {
            throw std::runtime_error("negative flux values in xml "
                                     "source model definitions are not "
                                     "allowed.");
         }
      }
   }

};

} // namespace celestialSources

#endif // celestialSources_ConstParMap_h

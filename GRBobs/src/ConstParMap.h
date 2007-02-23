/**
 * @file ConstParMap.h
 * @brief const interface to std::map<std::string, std::string>.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef genericSources_ConstParMapGRB_h
#define genericSources_ConstParMapGRB_h

#include <map>
#include <stdexcept>

#include "facilities/Util.h"

namespace GRBobs {

class ConstParMap {

 public:
  
  ConstParMap(const std::map<std::string, std::string> & parmap) 
    : m_parmap(parmap) {}
  
  ConstParMap(const std::string & params) {
    facilities::Util::keyValueTokenize(params, ", ", m_parmap);
  }
  
  const std::string & operator[](const std::string & name) const {
    std::map<std::string, std::string>::const_iterator item = m_parmap.find(name);
    if (item == m_parmap.end()) 
      {
	return name;
      }
    return item->second;
  }
  
  double value(const std::string & name) const {
    std::string value_name = operator[](name);
    if (operator[](name)==name)
      return -999;
    else
      return std::atof(value_name.c_str());
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
  
};
 
} // namespace genericSources

#endif // genericSources_ConstParMapGRB_h

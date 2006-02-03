/** 
* @file FileSpectrum.cxx
* @brief Implementation of FileSpectrum
*
*  $Header$
*/

#include "genericSources/FileSpectrum.h"

#include <cmath>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>

#include "flux/SpectrumFactory.h"
#include "facilities/Util.h"

ISpectrumFactory &FileSpectrumFactory() {
  static SpectrumFactory<FileSpectrum> factory;
  return factory;
}

FileSpectrum::FileSpectrum(const std::string& params) : 
  m_inLog(false)
{
  std::vector<efpair> temp_vector;
  //retrieve and parse the parameter string
  std::map<std::string,std::string> parmap;
  facilities::Util::keyValueTokenize(params,",",parmap);

  m_flux = std::atof(parmap["flux"].c_str());
  if(parmap["log"].find("yes")!=std::string::npos)
    {
      m_inLog = true; 
    }

  //open the spectrum file
  std::string fileName = parmap["specFile"];
  facilities::Util::expandEnvVar(&fileName);
  std::ifstream input_file;
  input_file.open(fileName.c_str());
  
  if(! input_file.is_open())
    {
      std::cerr << "ERROR:  Unable to open:  " << fileName.c_str() << std::endl;
      throw(std::string("ERROR:  Unable to open: "+ fileName ));
    }
  else 
    {
      temp_vector=readFile(input_file);
    }
        
  double total_flux = 0.0;
  std::pair<double,double> ef;
  std::vector<std::pair<double,double> >::iterator it; 
  for(it = temp_vector.begin(); it != temp_vector.end(); it++)
    {
      ef = (*it);
      double factor;
      
      if(it == temp_vector.begin() && (it+1) != temp_vector.end() )
	factor = ( (it+1)->first - ef.first );
      else if( it == temp_vector.begin() && (it+1) == temp_vector.end() )
	factor = 1;
      else if( (it+1) != temp_vector.end() )
	factor = (((it+1)->first - (it-1)->first)/2);
      else if( (it+1) == temp_vector.end() )
	factor = ( ef.first - (it-1)->first );
      
//      std::cout<<factor<<" "<<ef.second<<std::endl;
      ef.second *= factor;
      total_flux += ef.second;
      
      ef.second = total_flux;
      
      integ_flux.push_back(ef);
    }

  m_fileflux = total_flux;

  //set the flux to the file integral, if not set by XML
  if(m_flux==0){
    float scale =1.;
    if(parmap["cm2"].find("yes")!=std::string::npos){
      scale=10000.;}
    m_flux = scale*m_fileflux;
  }
}


double FileSpectrum::flux() const
{
    return m_flux;
}

double FileSpectrum::flux (double time ) const{
    return flux();
}


float FileSpectrum::operator() (float r)
{
    /// Purpose: sample a single particle energy from the spectrum
    double target_flux = r * m_fileflux;
    
    std::vector<efpair>::const_iterator i;
    
    std::pair<double,double> previous;
    
    i = integ_flux.begin();
    previous = (*i);
    
    for(i = integ_flux.begin(); i != integ_flux.end(); i++)
    {
      if((*i).second >= target_flux){
	if(i==integ_flux.begin()) i++;
	break;
      }
        previous = (*i);
    }
    
    // Use linear interpolation between bins
    double m = ( (*i).first - previous.first ) / ( (*i).second - previous.second );
    double b = (*i).first - m * (*i).second;
    float scale = 1.;
    double raw_e = m * target_flux + b;
    if(m_inLog)
      {
        return scale*pow(10., raw_e);
      }
    else
      {
        return scale*raw_e;
      }
}


std::string FileSpectrum::title() const
{
    return "FileSpectrum";
}

const char * FileSpectrum::particleName() const
{
    return m_particle_name.c_str();
}


std::vector<std::pair<double,double> > 
FileSpectrum::readFile(std::ifstream& input_file)
{
  char buffer[256];
  std::vector<std::pair<double,double> > temp_vector;

  while(input_file.getline(buffer,256,'\n'))
    {
      std::string line(buffer);
      std::vector<std::string> entries;
      if(line.find('%')!= std::string::npos)
      	{
	  //This is taken to be a comment line
	  continue;
      	}
      facilities::Util::stringTokenize(line," \t",entries);
      if(entries.back() == "") entries.pop_back(); //this should be removed after a fix in stringTokenize
      int size = entries.size();
      if(size == 2)
	{
	  temp_vector.push_back(std::make_pair<double,double>(atof(entries[0].c_str()),atof(entries[1].c_str())));
	} 
      else if(size == 3)
	{
	  std::cerr<<"Line contained more than 2 columns"<<std::endl;
        }
      
    }
  return temp_vector;
}

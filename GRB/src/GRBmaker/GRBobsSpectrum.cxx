// FILE: GRBobsSpectrum.cxx

#include "GRBmaker.h"
#include "GRBobsSpectrum.h"
#include "GRBobsConstants.h"
#include "src/GPS.h"
#include "src/SpectrumFactory.h"
#include "CLHEP/Random/RandFlat.h"


//static SpectrumFactory<GRBobsSpectrum> factory;
//const ISpectrumFactory& GRBobsSpectrumFactory = factory;

#include <fstream>
static std::ofstream os("gleam.lis");


// Constructors

GRBobsSpectrum::GRBobsSpectrum(const std::string &params)
	:ISpectrum(),
	 m_title("GRBobsSpectrum"),
	 m_particleName("gamma")
//	 m_grbMaker(new GRBmaker("GRB_000.lis"))
{
	std::cout << "GRBOBSSPECTRUM: params: " << params << std::endl;
	std::vector<std::string> paramVector;
	//paramVector.push_back("$(GRBROOT)/GRBmakerOutput/");
	//paramVector.push_back("GRB_000.lis");
	parseParamList(params, paramVector);

	m_grbMaker = new GRBmaker(paramVector);
}


GRBobsSpectrum::GRBobsSpectrum(double duration, int npuls, double flux, double fraction, double powerLawIndex, bool flag)
	:ISpectrum(),
	 m_title("GRBobsSpectrum"),
	 m_particleName("gamma"),
	 m_grbMaker(new GRBmaker(duration, npuls, flux, fraction, powerLawIndex, flag))
{
}


// Copy Constructor
GRBobsSpectrum::GRBobsSpectrum(const GRBobsSpectrum &right)
	:m_title(right.m_title),
	 m_particleName(right.m_particleName)
{
	if (right.m_grbMaker)
		m_grbMaker = right.m_grbMaker->clone();
}


// Destructor
GRBobsSpectrum::~GRBobsSpectrum()
{
	delete m_grbMaker;
	m_grbMaker = 0;
}


// Helper function for the assignment operator
void GRBobsSpectrum::swap(GRBobsSpectrum &other) throw()
{
	std::swap(m_title, other.m_title);
	std::swap(m_particleName, other.m_particleName);
	std::swap(m_grbMaker, other.m_grbMaker);
}


// Assignment Operator
GRBobsSpectrum &GRBobsSpectrum::operator=(const GRBobsSpectrum &right)
{
	GRBobsSpectrum temp(right);   // does all the work
	swap(temp);   // this can't throw
	return *this;
}


// Parse input parameter list obtained from the xml file
void GRBobsSpectrum::parseParamList(const std::string &input, std::vector<std::string>& output) const
{ 
	std::cout << "input: " << input << std::endl;
	int i = input.find_last_of("/");
	std::cout << "i: " << i << std::endl;

	if (i > 0)
	{
		std::cout << "input.size: " << input.size() << std::endl;
		if (i == input.size()-1)
			output.push_back(input);
		else
		{
			std::string temp = input.substr(0,i+1);
			std::cout << "temp: " << temp << std::endl;
			output.push_back(temp);
			temp = input.substr(i+1);
			std::cout << "temp: " << temp << std::endl;
			output.push_back(temp);
		}
	}

	else
	{
		std::string temp = input;
	    output.push_back(temp);
	}

#ifdef FOO
	output.push_back("$(GRBROOT)/GRBmakerOutput/");
	std::cout << "output: " << output.back() << std::endl;
	output.push_back("GRB_000.lis");
	std::cout << "output: " << output.back() << std::endl;

	
	int i = input.find_first_of(",");

	if (i > 0)
	{
		std::string temp = input.substr(0,i);
		output.push_back(temp);
		temp = input.substr(i+1);
		output.push_back(temp);
	}

	else
	{
		std::string temp = input;
	    output.push_back(temp);
	}
#endif
}


double GRBobsSpectrum::flux(double time) const
{ 
	std::cout << "m_grbmaker->flux: " << m_grbMaker->flux() << std::endl;
	return m_grbMaker->flux();
}


float GRBobsSpectrum::fraction(float energy) 
{
	std::cout<<"fraction"<< std::endl;
	return m_grbMaker->fraction();
}


double GRBobsSpectrum::interval(double time)
{
	static double currentPTime = 0.0;
	double nextPTime = 0.0;
	double intv = 0.0;


	nextPTime = nextTime();
	intv = nextPTime - currentPTime;
	std::cout << "time: " << time << " nextPTime: " << nextPTime << " currentPTime: " << currentPTime << 
		" intv: " << intv << std::endl;
	currentPTime = nextPTime;

	return intv;
}


    //JCT needs const to match pure virtual method
float GRBobsSpectrum::operator () (float randomNumber) const  
{
	return (float) nextEnergy();
}


// returns next available energy
// expects to be called once per time - so iterator should never become invalid before times run out
double GRBobsSpectrum::nextEnergy() const
{
	static std::vector<TimeEnergy>::iterator it = m_grbMaker->photonlist().begin();

	if (it != m_grbMaker->photonlist().end())
		return (*it++).energy();
	else
		return -1;
}


// returns next available energy to the simulation
double GRBobsSpectrum::energySrc(HepRandomEngine *engine, double time)
{
	return nextEnergy();
}


// returns next available time
double GRBobsSpectrum::nextTime() const
{
	static std::vector<TimeEnergy>::iterator it = m_grbMaker->photonlist().begin();

	if (it != m_grbMaker->photonlist().end())
		return (*it++).time();
	else
		return 1.0e+6;
}


std::pair<float,float> GRBobsSpectrum::dir(float energy) const
{
  return m_grbMaker->dir();
}
    

std::pair<double,double> GRBobsSpectrum::dir(double energy, HepRandomEngine *engine)
{
  return dir(energy);
}
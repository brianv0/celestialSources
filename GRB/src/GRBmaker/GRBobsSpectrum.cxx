// FILE: GRBobsSpectrum.cxx


#include <memory>  // for auto_ptr
#include <fstream>



#include "GRBmaker.h"
#include "GRBurst.h"
#include "GRBobsSpectrum.h"
#include "GRBobsConstants.h"
#include "flux/GPS.h"
#include "flux/SpectrumFactory.h"
#include "CLHEP/Random/RandFlat.h"


namespace {
        const double INVTIME=1.0e8;
        const double INVENERGY=-1;
        int indx=-1;
}



// Constructors


GRBobsSpectrum::GRBobsSpectrum(const std::string &params)
              : ISpectrum(),
                m_title("GRBobsSpectrum"),
                m_particleName("gamma")
{
    std::vector<std::string> paramVector;
    parseParamList(params, paramVector);
    
    std::auto_ptr<GRBmaker> grbMaker(new GRBmaker);
    std::auto_ptr<GRBurst> p(grbMaker->create(paramVector));
    m_grb = p.release();
}



GRBobsSpectrum::GRBobsSpectrum(const double duration, const int npuls, 
                               const double flux, const double fraction, 
                               const double alpha, const double beta, 
                               const double epeak, const double specnorm, 
                               const bool flag)
              :ISpectrum(),
               m_title("GRBobsSpectrum"),
               m_particleName("gamma")
{
    std::auto_ptr<GRBmaker> grbMaker(new GRBmaker);
    std::auto_ptr<GRBurst> p(grbMaker->create(duration, npuls, flux, fraction, alpha,
        beta, epeak, specnorm, flag));
    m_grb = p.release();
}



// Copy Constructor
GRBobsSpectrum::GRBobsSpectrum(const GRBobsSpectrum &right)
:ISpectrum(),
m_title(right.m_title),
m_particleName(right.m_particleName)
{
    if (right.m_grb)
        m_grb = right.m_grb->clone();
}



// Destructor
GRBobsSpectrum::~GRBobsSpectrum()
{
    delete m_grb;
    m_grb = 0;
}



// Helper function for the assignment operator
void GRBobsSpectrum::swap(GRBobsSpectrum &other) throw()
{
    std::swap(m_title, other.m_title);
    std::swap(m_particleName, other.m_particleName);
    std::swap(m_grb, other.m_grb);
}



// Assignment Operator
GRBobsSpectrum &GRBobsSpectrum::operator=(const GRBobsSpectrum &right)
{
    GRBobsSpectrum temp(right);   // does all the work
    swap(temp);   // this can't throw
    return *this;
}



// Parse input parameter list obtained from the xml file
void GRBobsSpectrum::parseParamList(const std::string &input, 
                                    std::vector<std::string>& output) const
{ 
    std::string::size_type i = input.find_last_of("/");
    
    if (i > 0)
    {
        if (i == input.size()-1)
            output.push_back(input);
        else
        {
            std::string temp = input.substr(0,i+1);
            output.push_back(temp);
            temp = input.substr(i+1);
            output.push_back(temp);
        }
    }
    
    else
    {
        std::string temp = input;
        output.push_back(temp);
    }
}



double GRBobsSpectrum::flux(double time) const
{ 
 //       std::cout << "Returning Flux: " << m_grb->globalData()->flux() << std::endl;
    return m_grb->globalData()->flux();
}



float GRBobsSpectrum::fraction(float energy) 
{
 //       std::cout << "Returning Fraction: " << m_grb->globalData()->fraction() << std::endl;
    return m_grb->globalData()->fraction();
}



double GRBobsSpectrum::interval(double time)
{
        static const std::vector<TimeEnergy> &photonList = m_grb->photonlist();
        static int sz = photonList.size();

// find the next photon in the list that comes after _time_ - recall that indx starts at -1
// if the next guy in the list is greater than time, just increment indx for use here and in
// energy and dir functions.

		if (time > photonList[sz-1].time()) return INVTIME;

		if (photonList[indx+1].time() < time){
			while (photonList[++indx].time() < time) {};
		}
		else indx++;

        if (indx < sz)
                return photonList[indx].time() - time;

        else
                return INVTIME;



#ifdef FOO
    static double currentPTime = 0.0;
    double nextPTime = 0.0;
    double intv = 0.0;
    
    
    nextPTime = nextTime();


        if (nextPTime != INVTIME)
        {
            intv = nextPTime - currentPTime;
                currentPTime = nextPTime;
        }


        else
                currentPTime = 0.0;


        std::cout << "Returning Interval: " << intv << std::endl;
   
    return intv;
#endif
}



//JCT needs const to match pure virtual method
float GRBobsSpectrum::operator () (float randomNumber) const  
{
        //float en = (float) nextEnergy();
        //return en;
    return (float) nextEnergy();
}



// returns next available energy
// expects to be called once per time - so iterator should never become invalid before times run out
double GRBobsSpectrum::nextEnergy() const
{
        static const std::vector<TimeEnergy> &photonList = m_grb->photonlist();
        static int sz = photonList.size();


        if (indx < sz)
                return photonList[indx].energy();
        else
                return INVENERGY;


#ifdef FOO
        static bool qNew=true;
    static std::vector<TimeEnergy>::iterator it = m_grb->photonlist().begin();


        if (qNew)
        {
                qNew = false;
            it = m_grb->photonlist().begin();
        }
    
    if (it != m_grb->photonlist().end())
        return (*it++).energy();
    else
        {
                qNew = true;
        return INVENERGY;
        }
#endif
}



// returns next available energy to the simulation
double GRBobsSpectrum::energy(double time)
{
        return nextEnergy();



#ifdef FOO
        float en = nextEnergy();
        std::cout << "Returning Energy: " << en << std::endl;
        return en;
    //return nextEnergy();
#endif
}



// returns next available time
double GRBobsSpectrum::nextTime() const
{
        static bool qNew=true;
        static std::vector<TimeEnergy>::iterator it;


        if (qNew)
        {
                qNew = false;
            it = m_grb->photonlist().begin();
        }
    
    if (it != m_grb->photonlist().end())
        return (*it++).time();
    else
        {
                qNew = true;
        return INVTIME;
        }
}



std::pair<float,float> GRBobsSpectrum::dir(float energy) const
{
        static const std::vector<TimeEnergy> &photonList = m_grb->photonlist();
        static int sz = photonList.size();


        if (indx < sz)
                return photonList[indx].photonDir();
        else
                return std::make_pair<float,float>(1E8,1E8);


#ifdef FOO
        // >>>>>>>>>>>> TEMPORARILY RETURN PHOTON DIRECTION -- 12/08/03 <<<<<<<<<<<<<<<<<<<<<<<<
    //return m_grb->dir();



        static bool qNew=true;
    static std::vector<TimeEnergy>::iterator it = m_grb->photonlist().begin();


        if (qNew)
        {
                qNew = false;
            it = m_grb->photonlist().begin();
        }
    
    if (it != m_grb->photonlist().end())
        {
                std::pair<float,float> d = (*it++).photonDir();
                std::cout << "Returning Direction: (" << d.first << "," << d.second << ")" << std::endl;
                return d;
        //return (*it++).photonDir();
        }


    else
        {
                qNew = true;
                std::cout << "Returning Direction: (-1,-1)" <<  std::endl;
                return std::make_pair<float,float>(-1,-1);
        }
#endif
}



std::pair<double,double> GRBobsSpectrum::dir(double energy)
{
        static const std::vector<TimeEnergy> &photonList = m_grb->photonlist();
        static int sz = photonList.size();


        if (indx < sz)
                return photonList[indx].photonDir();
        else
                return std::make_pair<double,double>(1E8,1E8);


#ifdef FOO
        // >>>>>>>>>>>> TEMPORARILY RETURN PHOTON DIRECTION -- 12/08/03 <<<<<<<<<<<<<<<<<<<<<<<<
        //return m_grb->dir();



        static bool qNew=true;
    static std::vector<TimeEnergy>::iterator it = m_grb->photonlist().begin();


        if (qNew)
        {
                qNew = false;
            it = m_grb->photonlist().begin();
        }
    
    if (it != m_grb->photonlist().end())
        {
                std::pair<float,float> d = (*it++).photonDir();
                std::cout << "Returning Direction: (" << d.first << "," << d.second << ")" << std::endl;
                return d;
        //return (*it++).photonDir();
        }


    else
        {
                qNew = true;
                std::cout << "Returning Direction: (-1,-1)" <<  std::endl;
                return std::make_pair<float,float>(-1,-1);
        }
#endif
}
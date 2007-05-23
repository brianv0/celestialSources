/**
* @file microQuasar.h
* @brief A source class for the flux package that 
* simulates the gamma emission of microQuasars
* @author R. Dubois
*
* $Header$
*/

#ifndef microQuasar_h
#define microQuasar_h

#include <vector>

#include "flux/Spectrum.h"
#include "CLHEP/Random/Random.h"

/**
* @class microQuasar
*
* @brief 
* A source class for the flux package that
* simulates the gamma emission of microQuasars
* @author R. Dubois
*
*/

class microQuasar : public Spectrum {

public:

	/// This constructor is required and used in FluxSource for
	/// "SpectrumClass" sources.
	microQuasar(const std::string &params);

	virtual ~microQuasar(){};


	/// @return Particle type, "gamma".
	virtual const char * particleName() const {return "gamma";}

	/// @return average differential flux integrated over energy (photons/m^2/sr).
	/// @param time Simulation time in seconds.
	virtual double flux(double) const {return 0;};

	/// @return particle energy in MeV
	/// @param xi uniform random deviate on the unir interval
	virtual float operator()(float xi) const;

	/// @return "Effective" solid angle (sr).
	virtual double solidAngle() const {return 0;}

	/// @return Title describing the spectrum.
	virtual std::string title() const {return "microQuasar";}

	/// @return Interval to the next event (seconds)
	virtual double interval(double time);  

	/// @return Photon energy (MeV).
	virtual double energy(double time);

	/// @return Photon direction in (zenith angle (degrees), azimuth (degrees,East=0,North=90deg)).
	virtual std::pair<double, double> dir(double energy){
		return std::make_pair(0,0);
	}

	/// @return energy and photon direction (zenith angle (degrees), azimuth (degrees,East=0,North=90deg)).
	std::pair<double, std::pair<double, double> > photon();

	class OrbitalRegion {
	public:
		OrbitalRegion() {};
		~OrbitalRegion() {};
		float getMinOrbitalPhase() {return m_minOrbitalPhase;}
		float getMaxOrbitalPhase() {return m_maxOrbitalPhase;}
		float getSpectralIndex(int region) {return m_spectralIndex[region];}
		void setSpectralIndex(float p1, float p2) {
			m_spectralIndex[0] = p1;
			m_spectralIndex[1] = p2;
		}
		void setOrbitalPhase(float p1, float p2) {
			m_minOrbitalPhase = p1;
			m_maxOrbitalPhase = p2;
		}
		int findRegion(double time, float period);

	private:
		/// minimum orbital phase
		float m_minOrbitalPhase;
		/// max orbital phase
		float m_maxOrbitalPhase;
		/// number of orbital phases with different spectral indices
		int m_numOrbitalRegions;
		/// spectral index for this region
		float m_spectralIndex[2];
	};

	class DiskCycleProperties {
	public:
			DiskCycleProperties() {}
			~DiskCycleProperties() {};
			
			float getCycleDuration() {return m_cycleDuration;}
			float getCycleDurationFluct() {return m_cycleDurationFluctuation;}
			void setCycleDuration(float p) {m_cycleDuration = p;}
			void setCycleDurationFluct(float p) {m_cycleDurationFluctuation = p;}

	private:
		/// disk cycle duration
		float m_cycleDuration;
		/// fractional fluctuation in cycle duration
		float m_cycleDurationFluctuation;

	};

	class JetProperties {
	public:
		JetProperties() {}
		~JetProperties() {};

		float getJetOnCycle() {return m_jetOnCycle;}
		float getJetOnCycleFluct() {return m_jetOnCycleFluctuation;}
		float getJetOnDuration() {return m_jetOnDuration;}
		float getJetOnDurationFluct() {return m_jetOnDurationFluctuation;}

		void setJetOnCycle(float p) { m_jetOnCycle = p;}
		void setJetOnCycleFluct(float p) { m_jetOnCycleFluctuation = p;}
		void setJetOnDuration(float p) { m_jetOnDuration = p;}
		void setJetOnDurationFluct(float p) { m_jetOnDurationFluctuation = p;}

	private:

		/// onset of jet during cycle (in fraction of the duration)
		float m_jetOnCycle;
		/// fluctuation of jet on onset
		float m_jetOnCycleFluctuation;
		/// duration of jet
		float m_jetOnDuration;
		/// fluctuation of jet duration
		float m_jetOnDurationFluctuation;
		/// current spectral index

	};
protected:

	void makeGrid(unsigned int n, double xmin, double xmax, 
		std::vector<double> &x, bool makeLog=false);
	double interpolate(const std::vector<double> &x, 
		const std::vector<double> &y,
		double xx);


private:
	/// function for generating the orbit modulation
	void modulation(const double x, double& funcValue, double& derivValue);

	/// function to tokenize the parameters and get rid of whitespaces
	std::vector<std::string> tokenize(std::string params, char token);

	/// fiddle of Numerical Recipe's rtsafe to avoid function passing
	double rtsafe(const double x1, const double x2, const double xacc);

	class burstPairs {
	public:
		burstPairs() {};
		burstPairs(double start, double end) {
		m_start = start;
		m_end =end;
		}
		~burstPairs() {};
		const bool operator () (const burstPairs* before, const burstPairs* after) const 
        {
		    return before->m_start == after->m_start;
		}
		std::pair<double,double> getBurstPairs() const { return std::make_pair(m_start,m_end);}
	private:
		double m_start;
		double m_end;
	};

    class InBurst
    {
    public:
        InBurst(double time = 0.) : m_time(time) {}
        const bool operator()(const microQuasar::burstPairs& pair) const
        {
            const std::pair<double,double> pairVals = pair.getBurstPairs();

            if (pairVals.first <= m_time && pairVals.second >= m_time) return true;
            return false;
        }
    private:
        double m_time;
    };

    class NextBurst
    {
    public:
        NextBurst(double time = 0.) : m_time(time) {}
        const bool operator()(const microQuasar::burstPairs& pair) const
        {
            const std::pair<double,double> pairVals = pair.getBurstPairs();

            if (pairVals.first >= m_time) return true;
            return false;
        }
    private:
        double m_time;
    };

	std::vector<burstPairs>::iterator getJetStart(double time);
	bool inJet(double time);

	burstPairs calculateJetStart(double time);

	/// flux (ph s^-1 cm^-2)
	float m_ftot;
	float m_alt;
	// minimum energy
	double m_eMin;
	/// max energy
	double m_eMax;
	/// period used 
	float m_orbitalPeriod;
	/// modulation over the period
	float m_phi0;
	float m_orbitalModulation;
	float m_nTurns;
	// current spectral index in use
	float m_currentSpectralIndex;
	/// randomly generated phase for generating interval
	float m_randPhase;	
	/// current time
	double m_currentTime;
	double m_jetStart;
	double m_jetEnd;
	double m_nJet;
	double m_cycleStart;
	// random see for burst generation
	double m_burstSeed;
	CLHEP::HepRandom m_randGenBurst;
	std::vector<burstPairs > m_bursts;

	/// disk-cycle properties
	DiskCycleProperties m_diskProperties;
	/// jet timing properties
	JetProperties m_jetProperties;
	/// allow max two orbital regions for now
	OrbitalRegion m_orbitalRegion;

};

#endif 

// LatGRBAlert.h
//
//
#ifndef LATGRBALERT_H
#define LATGRBALERT_H

#include <vector>
#include <deque>
#include <functional>   // for unary_function


class PhotonInfo;

class LatGRBAlert
{
 public:
	 // Constructor;  no destructor is needed
	 LatGRBAlert(const std::vector<PhotonInfo> &photonData, const long nbckoff, const long nGrb, const long nBck);

	 
	 bool falseTrig() const                             { return m_falseTrig; }
	 long nbckoff() const                               { return m_nbckoff; }
	 double firstTriggerTime() const                    { return m_firstTriggerTime; }
	 double highbck() const                             { return m_highbck; }
	 double maxResult() const                           { return m_maxResult; }
	 const std::vector<PhotonInfo> &photonData() const  { return m_photonData; }
	 const std::vector<double> &jointLike() const       { return m_jointLike; }
	 const std::deque<bool> &grbTruth() const           { return m_grbTruth; }
	 const std::deque<bool> &grb1phot() const           { return m_grb1phot; }


 private:
	 double firstTriggerTime(const long region, const std::vector<PhotonInfo> &photonData);
	 bool isFalseTrigger(const long region);
	 bool isSignal();
	 void recordFalseTriggers(const std::vector<long> &vindex);


	 // Data members
	 bool                     m_falseTrig;
	 long                     m_nbckoff;
	 double                   m_firstTriggerTime;
	 double                   m_highbck;
	 double                   m_maxResult;
	 std::vector<PhotonInfo>  m_photonData;
	 std::vector<double>      m_jointLike;
	 std::deque<bool>         m_grbTruth;
	 std::deque<bool>         m_grb1phot;




	 struct isSignalSet : public std::unary_function<PhotonInfo, PhotonInfo>
	 {
		 isSignalSet() {}
		 bool operator() (const PhotonInfo &x);
	 };
};


#endif
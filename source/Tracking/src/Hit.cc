/**
 * @file Hit.cpp
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "Hit.h"

#include <global_constants.h>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "DetectorModule.h"
#include "InactiveElement.h"
#include "MessageLogger.h"
#include "SimParms.h"
#include "Track.h"

//using namespace ROOT::Math;
using namespace std;

// bool Track::debugRemoval = false; // debug
//#ifdef HIT_DEBUG_RZ
//bool Track::debugRZCovarianceMatrix = false;  // debug
//bool Track::debugRZCorrelationMatrix = false;  // debug
//bool Track::debugRZErrorPropagation = false;  // debug
//#endif

/**
 * This is a comparator for two Hit objects based on smaller radius.
 * @param h1 A pointer to the first hit
 * @param h2 A pointer to the second hit
 * @return The result of the comparison: <i>true</i> if the distance from the z-axis of h1 is smaller than that of h2, false otherwise
 */
bool Hit::sortSmallerR(const HitPtr& h1, const HitPtr& h2) { return (h1->getRPos() < h2->getRPos()); }

/**
 * This is a comparator for two Hit objects based on higher radius
 * @param h1 A pointer to the first hit
 * @param h2 A pointer to the second hit
 * @return The result of the comparison: <i>true</i> if the distance from the z-axis of h1 is smaller than that of h2, false otherwise
 */
bool Hit::sortHigherR(const HitPtr& h1, const HitPtr& h2) { return (h1->getRPos() > h2->getRPos()); }

/**
 * Nothing to do for the destructor, as a hit never owns any objects it has pointers to...
 */
Hit::~Hit() {}

/**
 * The default constructor sets the internal parameters to default values.
 */
Hit::Hit() {
    m_detName          = "Undefined";
    m_distance         = 0;
    m_A                = 0;
    m_C                = 1;
    m_xPos             = 0;
    m_yPos             = 0;
    m_rPos             = 0;
    m_zPos             = 0;
    m_phi              = 0;
    m_rPhiLength       = 0;
    m_sLength          = 0;
    m_activity         = HitActivity::Undefined;
    m_hitModuleType    = HitModuleType::NONE;
    m_hitModule        = nullptr;
    m_track            = nullptr;
    m_isTrigger        = false;
    m_activeHitType    = HitActiveType::Undefined;
    m_passiveHitType   = HitPassiveType::Undefined;
    m_hitPassiveElem   = nullptr;
    m_isPixel          = false;
    m_layerID          = -1;
    m_discID           = -1;
    m_resolutionRPhi   = 0;
    m_resolutionZ      = 0;
}

/**
 * The copy constructor makes sure the new object doesn't point to the old track (the track pointer needs to
 * be set explicitly later). The pointer to the module, on the other hand, stays the same as that of the original.
 */
Hit::Hit(const Hit& h) {
    m_detName           = h.m_detName;
    m_distance          = h.m_distance;
    m_A                 = h.m_A;
    m_C                 = h.m_C;
    m_xPos              = h.m_xPos;
    m_yPos              = h.m_yPos;
    m_rPos              = h.m_rPos;
    m_zPos              = h.m_zPos;
    m_phi               = h.m_phi;
    m_rPhiLength        = h.m_rPhiLength;
    m_sLength           = h.m_sLength;
    m_activity          = h.m_activity;
    m_hitModuleType     = h.m_hitModuleType;
    m_hitModule         = h.m_hitModule;
    m_track             = nullptr;
    m_correctedMaterial = h.m_correctedMaterial;
    m_isTrigger         = h.m_isTrigger;
    m_activeHitType     = h.m_activeHitType;
    m_passiveHitType    = h.m_passiveHitType;
    m_hitPassiveElem    = h.m_hitPassiveElem;
    m_isPixel           = h.m_isPixel;
    m_layerID           = h.m_layerID;
    m_discID            = h.m_discID;
    m_resolutionRPhi    = h.m_resolutionRPhi;
    m_resolutionZ       = h.m_resolutionZ;
}

/*
 * Constructor for a hit on an inactive surface at [rPos, zPos] (cylindrical position) from the origin
 */
Hit::Hit(double rPos, double zPos, const InactiveElement* myPassiveElem, HitPassiveType passiveHitType) {
    m_detName           = "Undefined";
    m_distance          = sqrt(rPos*rPos + zPos*zPos);
    m_A                 = 0;    // High-momentum limit -> recalculated once known track
    m_C                 = 1;    // High-momentum limit -> recalculated once known track
    m_xPos              = rPos; // High-momentum limit -> recalculated once known track
    m_yPos              = 0;    // High-momentum limit -> recalculated once known track
    m_rPos              = rPos;
    m_zPos              = zPos;
    m_phi               = 0;    // High-momentum limit -> recalculated once known track
    m_rPhiLength        = rPos; // High-momentum limit -> recalculated once known track
    m_sLength           = 0;    // Needs to be recalculated once known track;
    m_activity          = HitActivity::Inactive;
    m_hitModuleType     = HitModuleType::NONE;
    m_hitModule         = nullptr;
    m_track             = nullptr;
    m_isTrigger         = false;
    m_activeHitType     = HitActiveType::Undefined;
    m_passiveHitType    = passiveHitType;
    setHitPassiveElement(myPassiveElem);
    m_isPixel           = false;
    m_layerID           = -1;
    m_discID            = -1;
    m_resolutionRPhi    = 0;
    m_resolutionZ       = 0;

    if (m_passiveHitType==HitPassiveType::BeamPipe) m_detName = "BeamPipe";
    if (m_passiveHitType==HitPassiveType::IP)       m_detName = "IP";
    if (m_passiveHitType==HitPassiveType::Support)  m_detName = "Support";
    if (m_passiveHitType==HitPassiveType::Service)  m_detName = "Service";
}


/**
 * //! Constructor for a hit on a given module at [rPos, zPos] (cylindrical position) from the origin
 * @param myModule pointer to the module with the hit 
 */
Hit::Hit(double rPos, double zPos, const DetectorModule* myModule, HitModuleType hitModuleType) {
    m_detName          = "Undefined";
    m_distance         = sqrt(rPos*rPos + zPos*zPos);
    m_A                = 0;    // High-momentum limit -> recalculated once known track
    m_C                = 1;    // High-momentum limit -> recalculated once known track
    m_xPos             = rPos; // High-momentum limit -> recalculated once known track
    m_yPos             = 0;    // High-momentum limit -> recalculated once known track
    m_rPos             = rPos;
    m_zPos             = zPos;
    m_phi              = 0;    // High-momentum limit -> recalculated once known track
    m_rPhiLength       = rPos; // High-momentum limit -> recalculated once known track
    m_sLength          = 0;    // Needs to be recalculated once known track;
    m_activity         = HitActivity::Active;
    m_hitModuleType    = hitModuleType;
    setHitModule(myModule);
    m_track            = nullptr;
    m_isTrigger        = false;
    if (myModule->measurementType()==MeasurementType::POSITION)   m_activeHitType    = HitActiveType::Position;
    if (myModule->measurementType()==MeasurementType::TIME)       m_activeHitType    = HitActiveType::Time;
    if (myModule->measurementType()==MeasurementType::POSANDTIME) m_activeHitType    = HitActiveType::PosAndTime;
    m_passiveHitType   = HitPassiveType::Undefined;
    m_hitPassiveElem   = nullptr;
    m_isPixel          = false;
    m_layerID          = -1;
    m_discID           = -1;
    m_resolutionRPhi   = 0;
    m_resolutionZ      = 0;
}


/*
 * Setter for the pointer to the active surface that caused the hit.
 * @param myModule A pointer to a barrel or endcap module; may be <i>NULL</i>
 */
void Hit::setHitModule(const DetectorModule* myModule) {

  if (myModule) m_hitModule = myModule;
  else logWARNING("Hit::setHitModule -> can't set module to given hit, pointer null!");
}

/*
 * Setter for the pointer to the inactive surface that caused the hit.
 */
void Hit::setHitPassiveElement(const InactiveElement* myPassiveElem) {

  if (myPassiveElem) m_hitPassiveElem = myPassiveElem;
  //else logWARNING("Hit::setHitPassiveElement -> can't set inactive element to given hit, pointer null!");
}

/*
 * Set track and recalculate hitXPos, hitYPos etc.
 */
void Hit::setTrack(const Track* newTrack) {

  m_track      = newTrack;
  m_A          = m_rPos/2./m_track->getRadius(m_zPos); // r_i/2R
  m_C          = sqrt(1-m_A*m_A);
  m_xPos       = m_rPos*m_C;
  m_yPos       = m_rPos*m_A;
  m_phi        = (m_A==0) ? 0. : asin(m_A)*2;
  m_rPhiLength = (m_A==0) ? m_rPos : m_track->getRadius(m_zPos)*m_phi;
  m_sLength    = m_rPhiLength/sin(m_track->getTheta());
};

/**
 * Get the track angle theta.
 * @return The angle from the z-axis of the entire track
 */
double Hit::getTrackTheta() {

  if (m_track==nullptr) {

    logWARNING("Hit::getTrackTheta -> no track assigned, will return zero!");
    return 0;
  }
  return (m_track->getTheta());
};

/**
 * Getter for the final, angle corrected pair of radiation and interaction lengths.
 * @return A copy of the pair containing the requested values; radiation length first, interaction length second
 */
RILength Hit::getCorrectedMaterial() {
    return m_correctedMaterial;
}

/**
 * Getter for the rPhi resolution (local rphi coordinate for a module)
 * If the hit is not active it returns -1
 * If the hit is connected to a module, then the module's resolution
 * is retured (if the hit is trigger-type, then then module's trigger resultion is requested)
 * if there is not any hit module, then the hit's resolution property is read and returned
 * @return the hit's local resolution
 */
double Hit::getRphiResolution(double trackRadius) {

  if (!(this->isActive())) {

    logERROR("Hit::getResolutionRphi called on a non-active hit");
    return -1;
  }
  else {

    // Module hit
    if (m_hitModule) {

      // R-Phi-resolution calculated as for a virtual barrel-type module -> transform local R-Phi res. to a true module orientation (rotation by theta angle, skew, tilt)
      // In detail, take into account a propagation of MS error on virtual barrel plane, on which all measurements are evaluated for consistency (global chi2 fit applied) ->
      // in limit R->inf. propagation along line used, otherwise a very small correction factor coming from the circular shape of particle track is required (similar
      // approach as for local resolutions)
      // TODO: Currently, correction mathematicaly derived only for use case of const magnetic field -> more complex mathematical expression expected in non-const B field
      // (hence correction not applied in such case)
      double A = 0;
      double C = 1;
      if (SimParms::getInstance().isMagFieldConst()) A = m_A; // r_i / 2R
      if (SimParms::getInstance().isMagFieldConst()) C = m_C; // sqrt[1-(r_i / 2R)^2)]
      double B            = A/sqrt(1-A*A);
      double tiltAngle    = m_hitModule->tiltAngle();
      double skewAngle    = m_hitModule->skewAngle();
      double resLocalRPhi = m_hitModule->resLocalRPhi();
      double resLocalZ    = m_hitModule->resLocalZ();

      // All modules & its resolution propagated to the resolution of a virtual barrel module (endcap is tilted by 90 degrees, barrel is tilted by 0 degrees)
      double resolutionRPhi = 0.;

      // For consistency reasons use parabolic approx or full approach
      if (SimParms::getInstance().useParabolicApprox()) {

        // TODO: Skew angle not used at the moment
        // resolutionRPhi = sqrt(pow((B*sin(skewAngle)*cos(tiltAngle) + cos(skewAngle)) * resLocalRPhi,2) + pow(B*sin(tiltAngle) * resLocalZ,2));
        resolutionRPhi = sqrt( pow(resLocalRPhi,2) + pow(B*sin(tiltAngle)*resLocalZ,2) );
      }
      else resolutionRPhi = sqrt( pow(C*resLocalRPhi,2) + pow(A*sin(tiltAngle)*resLocalZ,2) );

      return resolutionRPhi;
    }
    // IP or beam-constraint etc. hit
    else return m_resolutionRPhi;
  }
}

/**
 * Getter for the z resolution (local z coordinate for a module)
 * This corresponds to z coord for barrel modules and r coord for end-caps
 * If the hit is not active it returns -1
 * If the hit is connected to a module, then the module's resolution
 * is retured (if the hit is trigger-type, then then module's trigger resultion is requested)
 * if there is not any hit module, then the hit's resolution property is read and returned
 * @return the hit's local resolution
 */
double Hit::getZResolution(double trackRadius) {

  if (!(this->isActive())) {

    logERROR("Hit::getResolutionZ called on a non-active hit");
    return -1;
  }
  else {

    // Module hit
    if (m_hitModule) {

      if (m_track==nullptr) {

        logWARNING("Hit::getResolutionZ -> no track assigned, will return zero!");
        return 0;
      }
      else {

        // Z-resolution calculated as for a virtual barrel-type module -> transform local Z res. to a true module orientation (rotation by theta angle, skew, tilt)
        // In detail, take into account a propagation of MS error on virtual barrel plane, on which all measurements are evaluated for consistency (global chi2 fit applied) ->
        // in limit R->inf. propagation along line used, otherwise a very small correction factor coming from the circular shape of particle track is required (similar
        // approach as for local resolutions)
        // TODO: Currently, correction mathematicaly derived only for use case of const magnetic field -> more complex mathematical expression expected in non-const B field
        // (hence correction not applied in such case)
        double A = 0;
        double C = 1;
        if (SimParms::getInstance().isMagFieldConst()) A = m_A; // r_i / 2R
        if (SimParms::getInstance().isMagFieldConst()) C = m_C; // sqrt[1-(r_i / 2R)^2)]
        double cotgTheta    = m_track->getCotgTheta();
        double D            = cotgTheta/C;
        double tiltAngle    = m_hitModule->tiltAngle();
        double skewAngle    = m_hitModule->skewAngle();
        double resLocalRPhi = m_hitModule->resLocalRPhi();
        double resLocalZ    = m_hitModule->resLocalZ();

        // All modules & its resolution propagated to the resolution of a virtual barrel module (endcap is a tilted module by 90 degrees, barrel is tilted by 0 degrees)
        double resolutionZ = 0;

        // For consistency reasons use parabolic approx or full approach
        if (SimParms::getInstance().useParabolicApprox()) {

          // TODO: Skew angle not used at the moment
          //double resolutionZ = sqrt(pow(((D*cos(tiltAngle) + sin(tiltAngle))*sin(skewAngle)) * resLocalX,2) + pow((D*sin(tiltAngle) + cos(tiltAngle)) * resLocalY,2));
          resolutionZ = sqrt( pow((D*sin(tiltAngle) + cos(tiltAngle))*resLocalZ,2) );
        }
        else resolutionZ = sqrt( pow((C*cotgTheta*sin(tiltAngle) + cos(tiltAngle))*resLocalZ,2) + pow(A*cotgTheta*resLocalRPhi,2) );

        return resolutionZ;
      }
    }
    // IP or beam-constraint etc. hit
    else return m_resolutionZ;
  }
}

/*
 * Get R-Phi resulution in local module coordinates if module measures position
 */
double Hit::getLocalRPhiResolution() {

  if (!m_hitModule || !isPosMeasured()) {

    logERROR("Hit::getLocalRPhiResolution() - required for either passive hit or hit assigned to module, which doesn't measure position!!!");
    return -1;
  }
  else return m_hitModule->resLocalRPhi();
}

/*
 * Get Z resulution in local module coordinates if module measures position
 */
double Hit::getLocalZResolution() {

  if (!m_hitModule || !isPosMeasured()) {

    logERROR("Hit::getLocalZResolution() - required for either passive hit or hit assigned to module, which doesn't measure position!!!");
    return -1;
  }
  else return m_hitModule->resLocalZ();
}

/*
 * Get time-stamp resulution if module measures time
 */
double Hit::getTimeResolution() {

  if (!m_hitModule || !isTimeMeasured()) {

    logERROR("Hit::getTimeResolution() - required for either passive hit or hit assigned to module, which doesn't measure time!!!");
    return -1;
  }
  else return m_hitModule->resTime();
}

/*
 * Checks wether a module belongs to the outer endcap (no pixel allowed)
 * and the hit module is made of a square sensor
 * @return true if the module is in outer endcap and square
 */
bool Hit::isSquareEndcap() {

  //std::cout << "Hit::isSquareEndcap() "; //debug

  if (m_hitModule) {
    //std::cout << " hitModule_!= NULL "; //debug
    if (m_hitModule->subdet() == ENDCAP && m_hitModule->shape() == RECTANGULAR) {
      //std::cout << " getSubdetectorType()==Endcap "; //debug
       //std::cout << " getShape()==Rectangular "; //debug
       return true;
    }
  }
  //std::cout << std::endl; // debug
  return false;
}

bool Hit::isStub() const
{
  return m_hitModuleType == HitModuleType::STUB;
}

/*
 * Retrieves the module's half width
 * for hit related to endcap modules only
 * @return Modules half width
 */
double Hit::getD() {

  double result = 0;
  //std::cout << "Hit::getD() "; //debug
  if (m_hitModule) {
    //std::cout << " hitModule_!= NULL "; //debug
    try {

      const EndcapModule* myECModule = dynamic_cast<const EndcapModule*>(m_hitModule);//->as<EndcapModule>();
      if (myECModule) {
        //std::cout << " myECModule!= NULL "; //debug
        result = (myECModule->minWidth() + myECModule->maxWidth()) / 2. / 2.;
        //std::cout << " result = " << result; //debug
      }
    }
    catch (exception& e) {}
  }
  //std::cout << std::endl; // debug
  return result;
}



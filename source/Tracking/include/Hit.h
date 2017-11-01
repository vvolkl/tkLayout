/**
 * @file Hit.h
 * @brief This header file defines the hit and track classes used for internal analysis
 */
#ifndef INCLUDE_HIT_H_
#define INCLUDE_HIT_H_

#include <cmath>
#include <vector>

#include "DetectorModule.h"
#include "MaterialProperties.h"

// Forward declaration
class DetectorModule;
class Hit;
class InactiveElement;
class Track;

#undef HIT_DEBUG
#undef HIT_DEBUG_RZ

// Typedefs
typedef std::unique_ptr<Hit> HitPtr;
typedef std::vector<HitPtr>  HitCollection;

enum class HitActivity    : short { Undefined, Active, Inactive };               //!< Hit defined as pure material (inactive) or measurement point (active)
enum class HitActiveType  : short { Undefined, Position, Time, PosAndTime };     //!< Hit may be coming from position and/or time measurement
enum class HitPassiveType : short { Undefined, BeamPipe, IP, Service, Support }; //!< Hit defined as IP or pure material hit coming from: support, service, etc.

/**
 * @class Hit
 * @brief The Hit class is used when analysing a tracker layout to record information about a volume that was hit by a test track.
 *
 * It is used for both active and inactive surface hits. In case of an inactive surface, the pointer to the hit module will be <i>NULL</i>
 * since those volumes only matter with respect to the radiation and interaction lengths they add to the total in the error calculations.
 * All the other information is available to both categories. For convenience, the scaled radiation and interaction lengths are stored in
 * here as well to avoid additional computation and callbacks to the material property objects.
 */
// TODO: Remove pointer to track and find another method without any pointers involved!
class Hit {

public:

  //! Copy constructor
  Hit(const Hit& h);

  //! Constructor for a hit on an inactive surface at [rPos, zPos] (cylindrical position) from the origin
  Hit(double rPos, double zPos, const InactiveElement* myPassiveElem, HitPassiveType passiveHitType);

  //! Constructor for a hit on a given module at [rPos, zPos] (cylindrical position) from the origin
  Hit(double rPos, double zPos, const DetectorModule* myModule, HitModuleType hitModuleType);

  //! Destructor
  ~Hit();

  //! Given two hits, compare the distance to the z-axis based on smaller R
  static bool sortSmallerR(const HitPtr& h1, const HitPtr& h2);

  //! Given two hits, compare the distance to the z-axis based on higher R
  static bool sortHigherR(const HitPtr& h1, const HitPtr& h2);

  // Setter methods
  void setTrack(const Track* newTrack);

  void setAsActive()                                 { m_activity = HitActivity::Active;};
  void setAsPassive()                                { m_activity = HitActivity::Inactive;};
  void setAsTimeMeasurement()                        { m_activeHitType = HitActiveType::Time;};
  void setAsPosMeasurement()                         { m_activeHitType = HitActiveType::Position;};
  void setAsPosAndTimeMeasurement()                  { m_activeHitType = HitActiveType::PosAndTime;};
  void setAsPixel()                                  { m_isPixel  = true;}
  void setHitModuleType(HitModuleType hitModuleType) { m_hitModuleType = hitModuleType; }
  void setCorrectedMaterial(RILength newMaterial)    { m_correctedMaterial = newMaterial;};

  void setTrigger(bool isTrigger)                 { m_isTrigger = isTrigger;}
  void setResolutionRphi(double newRes)           { m_resolutionRPhi = newRes; } // Only used for virtual hits on non-modules
  void setResolutionZ(double newRes)              { m_resolutionZ = newRes; }    // Only used for virtual hits on non-modules
  void setResolutionY(double newRes)              { setResolutionZ(newRes); } // Used for compatibility only -> use setResolutionZ(double newRes) instead

  void setDetName(std::string detName)            { m_detName = detName; }
  void setLayerID(int layerID)                    { m_layerID = layerID; }
  void setDiscID(int discID)                      { m_discID  = discID; }

  // Getter methods
  const DetectorModule*  getHitModule() const         { return m_hitModule; };
  const InactiveElement* getHitPassiveElement() const { return m_hitPassiveElem; }

  double        getDistance() const         { return m_distance;};
  double        getXPos() const             { return m_xPos;};
  double        getYPos() const             { return m_yPos;};
  double        getRPos() const             { return m_rPos;};
  double        getZPos() const             { return m_zPos;};
  double        getRPhiLength() const       { return m_rPhiLength; }
  double        getCosBeta() const          { return m_C; }
  double        getSinBeta() const          { return m_A; }
  double        getPhiMinusPhi0() const     { return 2*m_beta; }
  double        getBeta() const             { return m_beta;}
  double        getSLength()  const         { return m_sLength; }
  double        getTilt() const             { if (this->isMeasurable()) return m_hitModule->tiltAngle(); else return 0; };
  bool          isActive() const            { if (m_activity==HitActivity::Active) return true; else return false;};
  bool          isPassive() const           { if (m_activity==HitActivity::Inactive) return true; else return false;};
  bool          isActivityUndefined() const { if (m_activity==HitActivity::Undefined) return true; else return false;};
  bool          isTimeMeasured() const      { if (m_activeHitType==HitActiveType::Time     || m_activeHitType==HitActiveType::PosAndTime) return true; else return false;};
  bool          isPosMeasured() const       { if (m_activeHitType==HitActiveType::Position || m_activeHitType==HitActiveType::PosAndTime) return true; else return false;};
  HitModuleType getHitModuleType() const    { return m_hitModuleType; } // NONE, INNER, OUTER, BOTH or STUB -- only meaningful for hits on active elements
  RILength getCorrectedMaterial();
  double   getLocalRPhiResolution();              //!< Get hit resolution in R-Phi in local module coordinates
  double   getRphiResolution(double trackRadius); //!< Get hit resolution in direction normal to tangent of track, i.e. in direction of curvature
  double   getLocalZResolution();                 //!< Get hit resolution in Z in local module coordinates
  double   getZResolution(double trackRadius);    //!< Get hit resolution in direction [n x t], i.e. perpendicular to tangent and normal (direction of curvature)
  double   getTimeResolution();                   //!< Get hit resolution in time
  double   getD();

  std::string getDetName()       const { return m_detName; };
  int         getLayerOrDiscID() const { if(this->isBarrel()) return m_layerID; else if(this->isEndcap()) return m_discID; else return -1;}; //!< Return positive number for layer (barrel hit) or disc (end-cap hit), -1 for beam-pipe or IP

  bool     isBeamPipe() const  { if (m_passiveHitType==HitPassiveType::BeamPipe) return true; else return false; };
  bool     isService() const   { if (m_passiveHitType==HitPassiveType::Service) return true; else return false; };
  bool     isSupport() const   { if (m_passiveHitType==HitPassiveType::Support) return true; else return false; };
  bool     isIP() const        { if (m_passiveHitType==HitPassiveType::IP) return true; else return false; };
  bool     isPixel() const     { return m_isPixel; };
  bool     isBarrel() const    { if (m_hitModule && (m_hitModule->subdet()==BARREL)) return true; else return false;};
  bool     isEndcap() const    { if (m_hitModule && (m_hitModule->subdet()==ENDCAP)) return true; else return false;};
  bool     isMeasurable() const{ if (m_hitModule!=nullptr) return true; else return false;};
  bool     isTrigger() const   { return m_isTrigger; };


  bool     isSquareEndcap();
  bool     isStub() const;

protected:
  
  //! Default constructor
  Hit();

  //! Set pointer to hit module in the constructor
  void setHitModule(const DetectorModule* myModule);

  //! Set pointer to inactive element in the constructor
  void setHitPassiveElement(const InactiveElement* myPassiveElem);

  double         m_distance;      //!< Distance of hit from origin in 3D = sqrt(rPos*rPos + zPos*zPos)
  double         m_xPos;          //!< Given track radius, distance of hit from origin in x = r.cos(beta-beta0) = r.cos(Phi/2-Phi0/2)
  double         m_yPos;          //!< Given track radius, distance of hit from origin in y = r.sin(beta-beta0) = r.sin(Phi/2-Phi0/2)
  double         m_rPos;          //!< Distance of hit from origin in the x/y plane (cylindrical coordinates -> r)
  double         m_zPos;          //!< Distance of hit from origin in z (cylindrical coordinates -> z)
  double         m_A;             //!< Given track radius: rPos/2R (sin[(Phi/2-Phi0/2)])
  double         m_C;             //!< Given track radius: cos[(Phi/2-Phi0/2)]
  double         m_beta;          //!< Given track radius: Phi/2-Phi0/2, where Phi & (rPos) represent circle polar coordinates with respect to PCR at Phi0
  double         m_sLength;       //!< Given track radius, R(Phi/2-Phi0/2)/sin(theta) -> track length with respect to PCR at Phi0
  double         m_rPhiLength;    //!< Given track radius, arc length with respect to PCR at phi0 : R(Phi/2-Phi0/2)
  HitActivity    m_activity;      //!< Hit defined as pure material (inactive) or measurement point (active)
  HitModuleType  m_hitModuleType; //!< Hit coming from inner, outer, stub, ... module
  HitActiveType  m_activeHitType; //!< Hit may be coming from position and/or time measurement
  HitPassiveType m_passiveHitType;//!< Hit coming from which passive part: beam-pipe, service, support etc.
  
  const DetectorModule*  m_hitModule;     //!< Const pointer to the hit module
  const InactiveElement* m_hitPassiveElem;//!< Const pointer to the hit inactive element
  const Track*           m_track;         //!< Const pointer to the track, into which the hit was assigned
  
  RILength m_correctedMaterial; //!< Material in the way of particle shot at m_track direction, i.e. theta, module tilt angles corrected
  
  bool m_isTrigger; //!< Hit comint from the trigger module?
  bool m_isPixel;   //!< Hit coming from the pixel module?

  std::string m_detName; //!< Detector name, in which the hit has been measured
  int         m_layerID; //!< Corresponding layer ID
  int         m_discID;  //!< Corresponding disc ID

private:
  
  double getTrackTheta();

  double m_resolutionRPhi; // Only used for virtual hits on non-modules
  double m_resolutionZ;    // Only used for virtual hits on non-modules

}; // Class

#endif /* INCLUDE_HIT_H_ */

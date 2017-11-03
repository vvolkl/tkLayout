/**
 * @file Track.cpp
 * @brief This file implements the hit and track classes used for internal analysis
 */

#include "Track.h"

#include <algorithm>
#include <cstdlib>

#include <global_constants.h>
#include "Hit.h"
#include "MessageLogger.h"
#include "MaterialProperties.h"
#include "SimParms.h"
#include "Units.h"


using namespace ROOT::Math;
using namespace std;

//
// Track constructor -> need to use setter methods to set: 2 of these [theta, phi, eta, cot(theta)] & 2 of these [mag. field, transv. momentum, radius]
//
Track::Track() :
  m_theta(0),
  m_phi(0),
  m_cotgTheta(0),
  m_eta(0),
  m_pt(0),
  m_p(0),
  m_reSortHits(true),
  m_covDone(false),
  m_refPointRPosCache(0),
  m_propagOutInCache(true)
{}

//
// Track copy-constructor -> creates deep copy of hit vector
//
Track::Track(const Track& track) {

  m_theta        = track.m_theta;
  m_phi          = track.m_phi;
  m_cotgTheta    = track.m_cotgTheta;
  m_eta          = track.m_eta;
  m_pt           = track.m_pt;
  m_p            = track.m_p;

  m_origin       = track.m_origin;
  m_direction    = track.m_direction;

  m_reSortHits   = track.m_reSortHits;
  m_covDone      = track.m_covDone;
  m_refPointRPosCache = track.m_refPointRPosCache;
  m_propagOutInCache  = track.m_propagOutInCache;

  m_covMatrixRPhi.ResizeTo(track.m_covMatrixRPhi);
  m_covMatrixRPhi = track.m_covMatrixRPhi;

  m_covMatrixRZ.ResizeTo(track.m_covMatrixRZ);
  m_covMatrixRZ = track.m_covMatrixRZ;

  m_covMatrixFull.ResizeTo(track.m_covMatrixFull);
  m_covMatrixFull = track.m_covMatrixFull;

  for (auto& iHit : track.m_hits) {
    HitPtr hit(new Hit(*iHit));
    addHit(std::move(hit));
  }
  m_tags = track.m_tags;
}

//
// Assign operator with deep copy of hit vector
//
Track& Track::operator= (const Track& track) {

  // check for self-assignment by comparing the address of the
  // implicit object and the parameter
  if (this == &track) return *this;
  
  // Do the copy
  m_theta        = track.m_theta;
  m_phi          = track.m_phi;
  m_cotgTheta    = track.m_cotgTheta;
  m_eta          = track.m_eta;
  m_pt           = track.m_pt;
  m_p            = track.m_p;

  m_origin       = track.m_origin;
  m_direction    = track.m_direction;

  m_reSortHits   = track.m_reSortHits;
  m_covDone      = track.m_covDone;
  m_refPointRPosCache = track.m_refPointRPosCache;
  m_propagOutInCache  = track.m_propagOutInCache;

  m_covMatrixRPhi.ResizeTo(track.m_covMatrixRPhi);
  m_covMatrixRPhi = track.m_covMatrixRPhi;

  m_covMatrixRZ.ResizeTo(track.m_covMatrixRZ);
  m_covMatrixRZ = track.m_covMatrixRZ;

  m_covMatrixFull.ResizeTo(track.m_covMatrixFull);
  m_covMatrixFull = track.m_covMatrixFull;

  for (auto& iHit : track.m_hits) {
    HitPtr hit(new Hit(*iHit));
    addHit(std::move(hit));
  }
  m_tags = track.m_tags;

  // Return the existing object
  return *this;
}

//
// Destructor
//
Track::~Track() {

  // Clear memory
  m_hits.clear();
}

//
// Calculate magnetic field at given z, assuming B = B(z).e_z + 0.e_x + 0 e_y
//
double Track::getMagField(double z) const {

  double magField = 0;

  // Option 1: Const mag. field across the detector: Bz = const
  if (SimParms::getInstance().isMagFieldConst()) {

    magField = SimParms::getInstance().magField[0];
  }
  // Option 2: Mag. field is a function in Z: B = B(z)
  else {

    for (unsigned int i=0; i<SimParms::getInstance().getNMagFieldRegions(); i++) {

      // Magnetic field regions are considered to be defined in metres
      if (z<SimParms::getInstance().magFieldZRegions[i]) {

        // Magnetic field is considered to be in Tesla
        magField = SimParms::getInstance().magField[i];
        break;
      }
    }
  }

  return magField;
}

//
// Get DeltaRho (error on 1/R) at refPoint [rPos, zPos] in "parabolic" approximation.
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
// Using 3x3 covariance propagator in case [r,z]!=[0,0]
//
double Track::getDeltaRho(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in R-Phi if something changed
  if (!m_covDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeCovarianceMatrix(refPointRPos, propagOutIn);

  double deltaRho = -1.;
  if (m_covDone && m_covMatrixRPhi(0,0)>=0) deltaRho = sqrt(m_covMatrixRPhi(0,0));

  // TODO: Not working correctly for B = B(z) & [r,z]!=[0,0], so print warning ...
  if (refPointRPos!=0. && !SimParms::getInstance().isMagFieldConst()) {

    logWARNING("Track::getDeltaRho(): Mathematical method to get deltaRho at [r,z]!=[0,0] in non const. B field not implemented, hence returned value at [r,z]=[0,0].");
  }

  return deltaRho;
}

//
// Get DeltaPtOvePt at refPoint [rPos, zPos] in "parabolic" approximation (utilize the calculated deltaRho quantity).
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
//
double Track::getDeltaPtOverPt(double refPointRPos, bool propagOutIn/*=true*/) {

  double deltaPtOverPt = -1.;

  // delta(1/R) & delta(pT) -> estimated at point [r,z] = [0,0] (important for use case, when B != const -> B = B(z))
  double deltaRho = getDeltaRho(refPointRPos,propagOutIn);
  double radius   = getRadius(refPointRPos*m_cotgTheta);       // Approximative transformation from rPos to zPos using tan(theta)
  if (deltaRho!=-1) deltaPtOverPt = deltaRho * radius; // dpT(z)/pT(z) = dRho(z) / Rho(z) = dRho(z) * R(z)

  return deltaPtOverPt;
}

// Get DeltaPtOverPt at [0,0] in full math. approach, i.e. using full 5x5 covariance matrix
double Track::getDeltaPtOverPtFull() {

  double deltaPtOverPt = -1.;

  // (Re)compute cov. matrix if something changed
  double refPointRPos = 0;
  bool   propagOutIn  = true;
  if (!m_covDone) computeCovarianceMatrix(refPointRPos, propagOutIn);

  double deltaRho = -1.;

  if (m_covDone && m_covMatrixFull(0,0)>=0) deltaRho = sqrt(m_covMatrixFull(0,0));
  double radius   = getRadius(0.);
  if (deltaRho!=-1) deltaPtOverPt = deltaRho * radius; // dpT(z)/pT(z) = dRho(z) / Rho(z) = dRho(z) * R(z)

  return deltaPtOverPt;
}

//
// Get DeltaPOverP at refPoint [rPos, zPos] in "parabolic" approximation (utilize deltaRho & deltaCotgTheta quantities)
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
//
double Track::getDeltaPOverP(double refPointRPos, bool propagOutIn/*=true*/) {

  double deltaPOverP = -1.;

  // Combining into p measurement
  // dp/p = dp_t/p_t + A / (1+A^2) * dA // with A = ctg(theta)
  // dp/p = dp_t/p_t + sin(theta)*cos(theta)*dcotg(theta)
  double deltaPtOverPt = getDeltaPtOverPt(refPointRPos,propagOutIn);
  double deltaCtgTheta = getDeltaCtgTheta(refPointRPos,propagOutIn);
  if (deltaPtOverPt!=-1 && deltaCtgTheta!=-1) deltaPOverP = sqrt(deltaPtOverPt*deltaPtOverPt + sin(m_theta)*sin(m_theta) * cos(m_theta)*cos(m_theta) * deltaCtgTheta*deltaCtgTheta);

  return deltaPOverP;
}

//
// Get DeltaPhi(Phi0) at refPoint [rPos, zPos] ([0,0]) in "parabolic" approximation
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
// Using 3x3 covariance propagator in case [r,z]!=[0,0]
//
double Track::getDeltaPhi(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in R-Phi if something changed
  if (!m_covDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeCovarianceMatrix(refPointRPos, propagOutIn);

  double deltaPhi0 = -1.;
  // No covariance propagation necessary at [0,0,0] point
  if (refPointRPos==0.) {

    if (m_covDone && m_covMatrixRPhi(1,1)>=0) deltaPhi0 = sqrt(m_covMatrixRPhi(1, 1));
  }
  else {

    if (m_covDone) {

      double covRhoRho    = m_covMatrixRPhi(0,0);
      double covRhoPhi0   = m_covMatrixRPhi(0,1);
      double covRhoD0     = m_covMatrixRPhi(0,2);
      double covPhi0Phi0  = m_covMatrixRPhi(1,1);
      double covPhi0D0    = m_covMatrixRPhi(1,2);
      double covD0D0      = m_covMatrixRPhi(2,2);

      double rho          = getRho(refPointRPos*m_cotgTheta);
      double rho2         = rho*rho;
      double rho4         = rho2*rho2;
      double refPointRPos2= refPointRPos*refPointRPos;

      double deltaPhi0Sq  = refPointRPos2*covRhoRho   + covPhi0Phi0                   + rho4*refPointRPos2*covD0D0;
             deltaPhi0Sq += 2*refPointRPos*covRhoPhi0 - 2*rho2*refPointRPos2*covRhoD0 - 2*rho2*refPointRPos*covPhi0D0;

      if (deltaPhi0Sq>=0) deltaPhi0 = sqrt(deltaPhi0Sq);
    }
  }

  // TODO: Not working correctly for B = B(z) & [r,z]!=[0,0], so print warning ...
  if (refPointRPos!=0. && !SimParms::getInstance().isMagFieldConst()) {

    logWARNING("Track::getDeltaPhi(): Mathematical method to get deltaPhi at [r,z]!=[0,0] in non const. B field not implemented, hence returned value at [r,z]=[0,0].");
  }

  return deltaPhi0;
}

//
// GetDeltaPhi at [0,0] in full math. approach, i.e. using full 5x5 covariance matrix
//
double Track::getDeltaPhi0Full() {

  double deltaPhi0 = -1.;

  // (Re)compute cov. matrix if something changed
  double refPointRPos = 0;
  bool   propagOutIn  = true;
  if (!m_covDone) computeCovarianceMatrix(refPointRPos, propagOutIn);

  if (m_covDone && m_covMatrixFull(1,1)>=0) deltaPhi0 = sqrt(m_covMatrixFull(1,1));

  return deltaPhi0;
}

//
// Get DeltaD (D0) at refPoint [rPos, zPos] ([0,0]) in "parabolic" approximation
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
// Using 3x3 covariance propagator in case [r,z]!=[0,0]
//
double Track::getDeltaD(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in R-Phi if something changed
  if (!m_covDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeCovarianceMatrix(refPointRPos, propagOutIn);

  double deltaD0 = -1.;
  // No covariance propagation necessary at [0,0,0] point
  if (refPointRPos==0.) {

    if (m_covDone && m_covMatrixRPhi(2,2)>=0) deltaD0 = sqrt(m_covMatrixRPhi(2,2));
  }
  else {

    if (m_covDone) {

      double covRhoRho     = m_covMatrixRPhi(0,0);
      double covRhoPhi0    = m_covMatrixRPhi(0,1);
      double covRhoD0      = m_covMatrixRPhi(0,2);
      double covPhi0Phi0   = m_covMatrixRPhi(1,1);
      double covPhi0D0     = m_covMatrixRPhi(1,2);
      double covD0D0       = m_covMatrixRPhi(2,2);

      double refPointRPos2 = refPointRPos*refPointRPos;
      double refPointRPos3 = refPointRPos2*refPointRPos;
      double refPointRPos4 = refPointRPos3*refPointRPos;

      double deltaD0Sq  = refPointRPos4/4.*covRhoRho + refPointRPos3*covRhoPhi0 + refPointRPos2*covRhoD0;
             deltaD0Sq += refPointRPos2*covPhi0Phi0  + 2*refPointRPos*covPhi0D0 + covD0D0;

      if (deltaD0Sq>=0) deltaD0 = sqrt(deltaD0Sq);
    }
  }

  // TODO: Not working correctly for B = B(z) & [r,z]!=[0,0], so print warning ...
  if (refPointRPos!=0. && !SimParms::getInstance().isMagFieldConst()) {

    logWARNING("Track::getDeltaD(): Mathematical method to get deltaD0 at [r,z]!=[0,0] in non const. B field not implemented, hence returned value at [r,z]=[0,0].");
  }

  return deltaD0;
}

//
// Get DeltaD0, i.e. at [0,0] in full math. approach, i.e. using full 5x5 covariance matrix
//
double Track::getDeltaD0Full() {

  double deltaD0 = -1.;

  // (Re)compute cov. matrix if something changed
  double refPointRPos = 0;
  bool   propagOutIn  = true;
  if (!m_covDone) computeCovarianceMatrix(refPointRPos, propagOutIn);

  if (m_covDone && m_covMatrixFull(2,2)>=0) deltaD0 = sqrt(m_covMatrixFull(2,2));

  return deltaD0;
}

//
// Get DeltaCtgTheta at refPoint [rPos, zPos] in "parabolic" approximation
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
// Using 2x2 covariance propagator in case [r,z]!=[0,0]
//
double Track::getDeltaCtgTheta(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in s-Z if something changed
  if (!m_covDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeCovarianceMatrix(refPointRPos, propagOutIn);

  double deltaCtgTheta = -1.;
  if (m_covDone && m_covMatrixRZ(0, 0)>=0) deltaCtgTheta = sqrt(m_covMatrixRZ(0, 0));

  return deltaCtgTheta;
}

//
// Get DeltaCtgTheta at [0,0] in full math. approach, i.e. using fulll 5x5 covariance matrix
//

double Track::getDeltaCtgThetaFull() {

  double deltaCtgTheta = -1.;

  // (Re)compute cov. matrix if something changed
  double refPointRPos = 0;
  bool   propagOutIn  = true;
  if (!m_covDone) computeCovarianceMatrix(refPointRPos, propagOutIn);

  if (m_covDone && m_covMatrixFull(3,3)>=0) deltaCtgTheta = sqrt(m_covMatrixFull(3,3));

  return deltaCtgTheta;
}

//
// Get DeltaZ (Z0) at refPoint [rPos, zPos] ([0,0]) in "parabolic" approximation
// Propagator direction defines, which part of tracker (at higher radii or lower radii from the ref. point) is going to be used.
// Using 2x2 covariance propagator in case [r,z]!=[0,0]
//
double Track::getDeltaZ(double refPointRPos, bool propagOutIn/*=true*/) {

  // (Re)compute cov. matrix in s-Z if something changed
  if (!m_covDone || refPointRPos!=m_refPointRPosCache || propagOutIn!=m_propagOutInCache) computeCovarianceMatrix(refPointRPos, propagOutIn);

  double deltaZ0 = -1.;

  // No covariance propagation necessary at [0,0,0] point
  if (refPointRPos==0.) {

    if (m_covDone && m_covMatrixRZ(1, 1)>=0) deltaZ0 = sqrt(m_covMatrixRZ(1, 1));
  }
  else {

    if (m_covDone) {

      double covZ0Z0       = m_covMatrixRZ(1,1);
      double covZ0CtgTheta = m_covMatrixRZ(0,1);
      double covCtgThCtgTh = m_covMatrixRZ(0,0);

      //std::cout << "<< " << rPos << " " << covZ0Z0 << " " << covZ0CtgTheta << "  " << covCtgThCtgTh << std::endl;

      double deltaZ0Sq = covZ0Z0 + 2*refPointRPos*covZ0CtgTheta + refPointRPos*refPointRPos*covCtgThCtgTh;
      if (deltaZ0Sq>=0) deltaZ0 = sqrt(deltaZ0Sq);
    }
  }

  return deltaZ0;
}

//
// Get DeltaZ0, i.e. at [0,0] in full math. approach, i.e. using full 5x5 covariance matrix
//

double Track::getDeltaZ0Full() {

  double deltaZ0 = -1.;

  // (Re)compute cov. matrix if something changed
  double refPointRPos = 0;
  bool   propagOutIn  = true;
  if (!m_covDone) computeCovarianceMatrix(refPointRPos, propagOutIn);

  if (m_covDone && m_covMatrixFull(4,4)>=0) deltaZ0 = sqrt(m_covMatrixFull(4, 4));

  return deltaZ0;
}

//
// Get DeltaT0 at refPoint [rPos, zPos] ([0,0]) combining all time-stamps along the track. Use track path length, defined as
// difference between time layer and ref. point to calculated the time of flight correction factor to individual time measurements.
//
double Track::getDeltaT(double refPointRPos) {

  double       deltaT = 0;

  for (auto& iHit : m_timeHits) {

    double radius    = getRadius(iHit->getZPos());
    //double deltaR    = iHit->getRPos() - refPointRPos;
    double deltaZ    = iHit->getZPos() - refPointRPos*getCotgTheta();
    //double sinTheta  = sin(this->getTheta());
    //double cotgTheta = m_cotgTheta;
    double cosTilt   = cos(iHit->getTilt());
    double sinTilt   = sin(iHit->getTilt());

    double sinPhiHalf    = iHit->getRPos()/2./radius;
    double sinPhiHalfRef = refPointRPos/2./radius;
    double cosPhiHalf    = sqrt(1-sinPhiHalf*sinPhiHalf);
    double deltaPhi      = asin(sinPhiHalf)*2 - asin(sinPhiHalfRef)*2;
    double pathLength    = radius*deltaPhi/sin(m_theta);

    double sigmaT         = iHit->getTimeResolution();
    double sigmaLocRPhiDet= iHit->getLocalRPhiResolution();
    double sigmaLocZDet   = iHit->getLocalZResolution();
    double sigmaPhiRef    = getDeltaPhi(refPointRPos);
    double sigmaZRef      = getDeltaZ(refPointRPos);

    double sigmaTRefSq = -1;
    if (sigmaT!=-1 && sigmaLocRPhiDet!=-1 && sigmaLocZDet!=-1 && sigmaPhiRef!=-1) {

      sigmaTRefSq = pow( (deltaZ*cosTilt - radius*deltaPhi*cosPhiHalf*sinTilt) * sigmaLocZDet,2);
      sigmaTRefSq+= pow( (deltaZ) * sigmaZRef,2);
      sigmaTRefSq+= pow( (radius*deltaPhi*sinPhiHalf) *sigmaLocRPhiDet,2);
      sigmaTRefSq+= pow( (radius*radius*deltaPhi) * sigmaPhiRef,2);

      sigmaTRefSq*= 1/Units::c/Units::c/pathLength/pathLength;
      sigmaTRefSq += sigmaT*sigmaT;

      // Do weighted average sigmaTot = sqrt[ 1/(Sum_i 1/sigma_i^2) ]
      deltaT   += 1/sigmaTRefSq;
    }
  }

  // Get final sigma from weighted average
  if (deltaT>=0) deltaT = sqrt(1./deltaT);
  else           deltaT = -1;

  return deltaT;
}

//
// Get DeltaCTau for secondary particles coming from the primary vertex at ~ [0,0] -> an important quantity to estimate the
// resolution of secondary vertices
//
double Track::getDeltaCTau() {

  double deltaCTau = -1;
  double deltaD0   = getDeltaD0();
  double deltaZ0   = getDeltaZ0();

  if (deltaD0>0 && deltaZ0>0) {
    deltaCTau = (cos(m_theta)*cos(m_theta)+1)/deltaD0/deltaD0 + sin(m_theta)*sin(m_theta)/deltaZ0/deltaZ0;
    deltaCTau = sqrt(2/deltaCTau);
  }

  return deltaCTau;
}

//
// Adds a new hit to the track (hit radius automatically updated)
//
void Track::addHit(HitPtr newHit) {

  // Add tracking tags
  if (newHit->getHitModule() != nullptr) {
    m_tags.insert(newHit->getHitModule()->trackingTags.begin(), newHit->getHitModule()->trackingTags.end());
  }
  newHit->setTrack(this);

  // Check that radial pos of new hit below track curling radius
  if (isHitRPosLowerThanCurlingRadius(newHit->getRPos(), newHit->getZPos())) {

    // Clone hit for timing information
    if (newHit->isTimeMeasured()) {

      HitPtr cloneHit(new Hit(*newHit));
      m_timeHits.push_back(std::move(cloneHit));
    }

    // Save position hit into hits
    m_hits.push_back(std::move(newHit));

  }
  else newHit.reset(nullptr);

  // Hits need to be re-sorted & cov. matrices recalculated
  m_reSortHits  = true;
  m_covDone     = false;
}

//
// Add IP constraint to the track, technically new hit is assigned: with no material and hit resolution in R-Phi as dr, in s-Z as dz
//
void Track::addIPConstraint(double dr, double dz) {

  // This modeling of the IP constraint was validated:
  // By placing dr = 0.5 mm and dz = 1 mm one obtains
  // sigma(d0) = 0.5 mm and sigma(z0) = 1 mm
  HitPtr newHit(new Hit(0,0,nullptr,HitPassiveType::IP));

  RILength emptyMaterial;
  emptyMaterial.radiation   = 0;
  emptyMaterial.interaction = 0;

  newHit->setCorrectedMaterial(emptyMaterial);
  newHit->setAsActive();
  newHit->setResolutionRphi(dr);
  newHit->setResolutionZ(dz);

  // Check that radial pos of new hit below track curling radius
  if (isHitRPosLowerThanCurlingRadius(newHit->getRPos(), newHit->getZPos())) m_hits.push_back(std::move(newHit));
  else newHit.reset(nullptr);

  // Hits need to be re-sorted & cov. matrices recalculated
  m_reSortHits  = true;
  m_covDone     = false;
}

//
// Set track polar angle - theta, azimuthal angle - phi, particle transverse momentum - pt
// (magnetic field obtained automatically from SimParms singleton class)Setter for the track azimuthal angle.
//
const Polar3DVector& Track::setThetaPhiPt(const double& newTheta, const double& newPhi, const double& newPt) {

  m_theta     = newTheta;
  m_cotgTheta = 1/tan(newTheta);
  m_eta       = -log(tan(m_theta/2));
  m_phi       = newPhi;
  m_pt        = newPt;
  m_p         = m_pt/sin(m_theta);

  if (m_pt>=0) m_direction.SetCoordinates(+1, m_theta, m_phi); // Particle inside-out
  else         m_direction.SetCoordinates(-1, m_theta, m_phi); // Particle outside-in

  // Clear all previously assigned hits -> hits need to be recalculated
  m_hits.clear();

  // Hits need to be re-sorted & cov. matrices recalculated
  m_reSortHits  = true;
  m_covDone     = false;

  return m_direction;
}

//
// Re-set transverse momentum (recalculate hit quantities dependent on track radius) + resort hits (if changing direction) +
// initiate recalc of cov matrices + prune hits (otherwise they may not lie on the new track, originally found at high pT limit)
void Track::resetPt(double newPt) {

  if (newPt*m_pt<0) m_reSortHits = true;
  m_covDone     = false;

  // Set pt, recalculate track radius related hit quantities & prune hits
  m_pt = newPt;
  m_p  = m_pt/sin(m_theta);
  for (auto& iHit : m_hits) { iHit->setTrack(this); }
  pruneHits();
}

//
// Sort internally all hits assigned to this track -> sorting algorithm based on hit radius - by smaller radius sooner or vice-versa (inner-2-outer approach or vice-versa)
//
void Track::sortHits(bool bySmallerR) { bySmallerR ? std::stable_sort(m_hits.begin(), m_hits.end(), Hit::sortSmallerR) : std::stable_sort(m_hits.begin(), m_hits.end(), Hit::sortHigherR); }

//
// Remove hits that don't follow the parabolic approximation used in tracking - TODO: still needs to be updated (not all approximations taken into account here)
//
bool Track::pruneHits() {

  bool isPruned = false;

  HitCollection newHits;
  for (auto& iHit : m_hits) {

    // Check that radial pos of new hit below track curling radius
    if (isHitRPosLowerThanCurlingRadius(iHit->getRPos(),iHit->getZPos())) newHits.push_back(std::move(iHit));
    else {

      // Clear memory
      iHit.reset();

      isPruned = true;
    }
  }

  m_hits.clear();
  for (auto& iHit : newHits) m_hits.push_back(std::move(iHit));

  return isPruned;
}

//
// Set active only trigger hits, so all other hits are made as inactive
//
void Track::keepTriggerHitsOnly() {

  for (auto& iHit : m_hits) {

    // Hit needs to be measurable, i.e. is linked to module
    if (iHit->isMeasurable()) {

      if (iHit->isActive()) {
        if      (iHit->isPixel()) iHit->setAsPassive();
        else if (iHit->getHitModule()->sensorLayout()!=PT) iHit->setAsPassive();
      }
    }
  } // For
}

//
// Set active only hits with the given tag
//
void Track::keepTaggedHitsOnly(const string& tag, bool useIP /*=true*/) {

  for (auto& iHit : m_hits) {

    // IP constraint hit
    if (tag=="all" && iHit->isIP() && useIP) iHit->setAsActive();

    // Measurement hit
    if (iHit->isMeasurable()) {
      if (tag=="all") iHit->setAsActive();
      else {
        if (std::count_if(iHit->getHitModule()->trackingTags.begin(), iHit->getHitModule()->trackingTags.end(), [&tag](const string& s){ return s == tag; })) iHit->setAsActive();
        else iHit->setAsPassive();
      }
    }
  }

  // Cov. matrices need to be recalculated
  m_covDone = false;
}

//
// Remove material from all assigned hits -> modify all hits such as they are without any material
//
void Track::removeMaterial() {

  // Material object with no material assigned
  RILength nullMaterial;

  // Reset all material assigned to hits
  for (auto& iHit : m_hits) iHit->setCorrectedMaterial(nullMaterial);

  // Cov. matrices need to be recalculated
  m_covDone = false;
}

//
// Helper method printing track covariance matrices in R-Phi
//
void Track::printErrors() {

  std::cout << "Covariance matrix: " << std::endl;
  m_covMatrixRPhi.Print();

  // Print errors @ [r,z]=[0,0]
  double rPos = 0.0;

  std::cout << "Rho errors by momentum: " << getDeltaRho(rPos) << std::endl;
  std::cout << "Phi0 errors by momentum: "<< getDeltaPhi0()    << std::endl;
  std::cout << "D0 errors by momentum: "  << getDeltaD0()      << std::endl;
}

//
// Helper method printing symmetric matrix
//
void Track::printSymMatrix(const TMatrixTSym<double>& matrix) const {

  std::cout << std::endl;

  int nCols = matrix.GetNcols();
  int nRows = matrix.GetNrows();

  for (int i = 0; i<nRows; i++) {
    std::cout << "(";
    for (int j=0; j<nCols;j++) {

      std::cout << " " << std::showpos << std::scientific << std::setprecision(4) << matrix(i,j);
    }
    std::cout << ")" << std::endl;
  }
  std::cout << std::fixed << std::noshowpos << std::endl;
}

//
// Helper method printing matrix
//
void Track::printMatrix(const TMatrixT<double>& matrix) const {

  std::cout << std::endl;

  int nCols = matrix.GetNcols();
  int nRows = matrix.GetNrows();

  for (int i = 0; i<nRows; i++) {
    std::cout << "(";
    for (int j=0; j<nCols;j++) {

      std::cout << " " << std::showpos << std::scientific << std::setprecision(5) << matrix(i,j);
    }
    std::cout << ")" << std::endl;
  }
  std::cout << std::fixed << std::noshowpos << std::endl;
}

//
// Helper method printing track hits
//
void Track::printHits() const {

  std::cout << "******************" << std::endl;
  std::cout << "Track eta=" << m_eta << std::endl;

  for (const auto& it : m_hits) {
    std::cout << "    Hit";
    if (it->isActive())   std::cout << " r="  << it->getRPos() << " +- " << it->getRphiResolution(getRadius(it->getZPos()));
    else                  std::cout << " r="  << it->getRPos();
    if (it->isActive())   std::cout << " z="  << it->getZPos() << " +- " << it->getZResolution(getRadius(it->getZPos()));
    else                  std::cout << " z="  << it->getZPos();
    std::cout << " d="  << it->getDistance()
              << " rl=" << it->getCorrectedMaterial().radiation
              << " il=" << it->getCorrectedMaterial().interaction;
    if (it->isActive())   std::cout << " active";
    else                  std::cout << " inactive";
    if (it->isBarrel())   std::cout << " barrel";
    if (it->isEndcap())   std::cout << " endcap";
    if (it->isBeamPipe()) std::cout << " beam-pipe";
    if (it->isIP())       std::cout << " ip";
    if (it->getLayerOrDiscID()!=-1) std::cout << " " << it->getDetName() << " L/D_id= " << it->getLayerOrDiscID();

    if (it->isActive()) {
      std::cout << " hitModuleType_=" << static_cast<short>(it->getHitModuleType());
    }
    std::cout << std::endl;
  }
}

//
// Helper method printing track hits
//
void Track::printActiveHits() const {

  std::cout << "******************" << std::endl;
  std::cout << "Track eta=" << m_eta << std::endl;

  for (const auto& it : m_hits) {
    if (it->isActive()) {

      std::cout << "    Hit";
      std::cout << " r="  << it->getRPos() << " +- " << it->getRphiResolution(getRadius(it->getZPos()));
      std::cout << " z="  << it->getZPos() << " +- " << it->getZResolution(getRadius(it->getZPos()));
      std::cout << " d="  << it->getDistance()
                << " rl=" << it->getCorrectedMaterial().radiation
                << " il=" << it->getCorrectedMaterial().interaction;
      if (it->isActive())   std::cout << " active";
      else                  std::cout << " inactive";
      if (it->isBarrel())   std::cout << " barrel";
      if (it->isEndcap())   std::cout << " endcap";
      if (it->isBeamPipe()) std::cout << " beam-pipe";
      if (it->isIP())       std::cout << " ip";
      if (it->getLayerOrDiscID()!=-1) std::cout << " " << it->getDetName() << " L/D_id= " << it->getLayerOrDiscID();
      std::cout << " hitModuleType_=" << static_cast<short>(it->getHitModuleType());
      std::cout << std::endl;
    }
  }
}

//
// Get number of active hits assigned to track for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file). If tag specified as "all" no extra tag required
//
int Track::getNActiveHits (std::string tag, bool useIP /* = true */ ) const {

  // Result variable
  int nHits=0;

  for (auto& iHit : m_hits) {
    if (iHit && iHit->isActive()){
      if (iHit->isIP() && useIP) {
        nHits++;
      }
      else if (!iHit->isIP()) {

        // Check tag for non-IP assigned hits
        bool tagOK = false;
        for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
          if (tag==*it || tag=="all") tagOK = true;
        }

        if (tagOK) nHits++;
      }
    }
  } // For

  return nHits;
}

//
// Get number of active hits coming from measurement planes or IP constraint assigned to track for given tag. If tag specified as "all", all module & IP hits assigned.
//
int Track::getNMeasuredHits(std::string tag, bool useIP /*=true*/) const {

  // Result variable
  int nHits=0;

  for (auto& iHit : m_hits) {
    if (iHit && iHit->isActive()) {
      if (iHit->isIP() && useIP) {
        nHits++;
      }
      else if (iHit->isMeasurable()) {

        // Check tag for non-IP assigned hits
        bool tagOK = false;
        for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
          if (tag==*it || tag=="all") tagOK = true;
        }
        if (tagOK) nHits++;
      }
    }
  } // For

  return nHits;

}

//
// Get reference to a hit, which can be measured, i.e. coming from measurement plane (active or inactive) or IP constraint
//
const Hit* Track::getMeasurableOrIPHit(int iHit) {

  int   hitCounter = 0;
  const Hit* pHit  = nullptr;

  // Sort hits based on particle direction: in-out or out-in
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  for (auto& hit : m_hits) {
    if (hit && (hit->isIP() || hit->isMeasurable())) {

      // Hit we're looking for!
      if (hitCounter==iHit) {

        pHit = hit.get();
        break;
      }
      hitCounter++;
    }
  }

  return pHit;
}

//
// Reverse search - Get reference to a hit, which can be measured, i.e. coming from measurement plane (active or inactive) or IP constraint
//
const Hit* Track::getRMeasurableOrIPHit(int iHit) {

  int   hitCounter = 0;
  const Hit* pHit  = nullptr;

  // Sort hits based on particle direction: in-out or out-in
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  for (auto hit = m_hits.rbegin(); hit != m_hits.rend(); ++hit) {
    if (*hit && ((*hit)->isIP() || (*hit)->isMeasurable())) {

      // Hit we're looking for!
      if (hitCounter==iHit) {

        pHit = (*hit).get();
        break;
      }
      hitCounter++;
    }
  }

  return pHit;
}

//
// Get the probabilty of having "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
// If tag specified as "all" no extra tag required
//
std::vector<double> Track::getHadronActiveHitsProbability(std::string tag) {

  // Result variable
  std::vector<double> probabilities;
  double probability = 1;

  // Sort hits based on particle direction: in-out or out-in
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  for (auto& iHit : m_hits) {
    if (iHit) {
      if (iHit->isActive()){

        // Check tag
        bool tagOK = false;
        for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
          if (tag==*it || tag=="all") tagOK = true;
        }

        if (tagOK) probabilities.push_back(probability);
      }

      // Decrease the probability that the next hit is a clean one
      RILength myMaterial = iHit->getCorrectedMaterial();
      probability /= exp(myMaterial.interaction);
    }
  } // For

  return probabilities;
}

//
// Get the probabilty of having a given number of "clean" hits for nuclear-interacting particles for given tag: pixel, strip, tracker, etc. (as defined in the geometry config file)
// If tag specified as "all" no extra tag required
//
double Track::getHadronActiveHitsProbability(std::string tag, int nHits) {

  // Probability
  double probability = 1;

  // Number of clean hits
  int goodHits = 0;

  // Sort hits based on particle direction: in-out or out-in
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  for (auto& iHit : m_hits) {

    if (iHit) {
      if (iHit->isActive()) {

        // Check tag
        bool tagOK = false;
        for (auto it=iHit->getHitModule()->trackingTags.begin(); it!=iHit->getHitModule()->trackingTags.end(); it++) {
          if (tag==*it || tag=="all") tagOK = true;
        }

        if (tagOK) goodHits++;
      }

      // If I reached the requested number of hits
      if (goodHits==nHits) return probability;

      // Decrease the probability that the
      // next hit is a clean one
      RILength myMaterial = iHit->getCorrectedMaterial();
      probability /= exp(myMaterial.interaction);
    }
  }

  // If I did not reach the requested number of active hits
  // The probability is zero
  return 0;
}

//
// Get track material
//
RILength Track::getMaterial() const {

  RILength totalMaterial;
  totalMaterial.radiation   = 0;
  totalMaterial.interaction = 0;

  for (auto& iHit : m_hits) totalMaterial += iHit->getCorrectedMaterial();

  return totalMaterial;
}

//
// Get a vector of pairs: Detector module & hit type for Trigger hits
//
std::vector<std::pair<const DetectorModule*, HitModuleType>> Track::getHitModules() const {

  std::vector<std::pair<const DetectorModule*, HitModuleType>> result;

  for (auto& iHit : m_hits) {

    if ((iHit) && (iHit->isTrigger()) && (!iHit->isIP()) && (iHit->isActive())) {

      // We've got a possible trigger here
      // Let's find the corresponding module
      const DetectorModule* myModule = iHit->getHitModule();
      if (myModule) result.push_back(std::make_pair(myModule, iHit->getHitModuleType()));
      else {
        // Whoops: problem here: an active hit is not linked to any module
        logERROR("Track::getHitModules: This SHOULD NOT happen. In expectedTriggerPoints() an active hit does not correspond to any module!");
      }
    }
  }
  return result;
}

//
// Main method calculating track parameters in both r-phi & s-z planes, using linear fit with parameters: 1/R, d0, phi0, cotg(theta) & z0. Within the calculation an NxN variance
// matrix is used (N hits = K+L: K active hits on detectors + L passive (artificial) hits due to material). Ref. point dictates, whether Multiple scattering effects need to be
// calculated inside-out or outside-in. MS effect is symmetric as regards track fitting. PropagOutIn variable defines whether detectors at higher R than the ref. point (true)
// should be used for error calculation/propagation (e.g. d0,z0) or whether detectors at lower R (false). E.g. for standard estimation of D0,Z0 parameters one calculates MS
// inside->out from rPos=0 (zPos can be calculated from rPos using theta).
// Return true if covariance matrix correctly calculated
//
bool Track::computeCovarianceMatrix(double refPointRPos, bool propagOutIn) {

  // Sort hits based on particle direction: in-out or out-in (if needed)
  if (m_reSortHits) {

    bool bySmallerRadius = true;
    if (m_pt>=0) sortHits(bySmallerRadius);
    else         sortHits(!bySmallerRadius);
    m_reSortHits = false;
  }

  m_refPointRPosCache = refPointRPos;
  m_propagOutInCache  = propagOutIn;

  // Variance matrices size
  int nHits = m_hits.size();

  TMatrixTSym<double> varMatrixRR;
  TMatrixTSym<double> varMatrixRZ;
  TMatrixTSym<double> varMatrixZR;
  TMatrixTSym<double> varMatrixZZ;
  TMatrixTSym<double> varMatrixFull;

  varMatrixRR.ResizeTo(nHits,nHits);
  varMatrixRZ.ResizeTo(nHits,nHits);
  varMatrixZR.ResizeTo(nHits,nHits);
  varMatrixZZ.ResizeTo(nHits,nHits);

  //
  // Find hit index ranges relevant for MS effects calculation based on requirements on refPoint & propagation direction
  int iStart          = nHits;
  int iEnd            = nHits-1;
  int nHitsUsed       = 0;
  int nActiveHitsUsed = 0;

  bool bySmallerR = true;
  //bool useIP      = false;

  // Particle traverses inside-out with error propagation outside-in -> hits already sorted correctlly
  if (m_pt>=0 && propagOutIn) {

    for (auto it=m_hits.rbegin(); it!=m_hits.rend(); it++) {

      if (refPointRPos<(*it)->getRPos()) {
        if ((*it)->isMeasurable()) nActiveHitsUsed++;
        iStart--;
      }
      else break;
    }
  }
  // Particle traverses inside-out with error propagation inside-out -> sort hits by higher radius & resort backwards after matrix calculated
  else if (m_pt>=0 && !propagOutIn) {

    sortHits(!bySmallerR);
    for (auto it=m_hits.rbegin(); it!=m_hits.rend(); it++) {

      if (refPointRPos>(*it)->getRPos()) {
        if ((*it)->isMeasurable()) nActiveHitsUsed++;
        iStart--;
      }
      else break;
    }
  }
  // TODO: Implement
  else {

    logWARNING("Variance matrix V(NxN) -> calculations of particle traversing outside-in not yet implemented.");
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); // Important -> resort back
    return false;
  }

  // Check that enough hits for track calculation, i.e. >=3
  nHitsUsed = nHits - iStart;
  if ((nActiveHitsUsed)<3) {

    std::string message = "Variance matrix V(NxN) -> refPointRPos[mm]="+any2str(refPointRPos/Units::mm,1);
                message+= " in combination with propagator direction doesn't provide sufficient number of hits for tracking!";
    logWARNING(message);
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); // Important -> resort back
    return false;
  }

  // Pre-compute the squares of the scattering angles projected to the virtual measurement plane
  std::vector<double> msThetaSqProj;

  for (int i=iStart; i<iEnd; i++) {

    // MS theta
    double msThetaSq = 0.0;

    // Material in terms of rad. lengths -> corrected by 1/cos(phi/2-phi0), i.e. by true material crossed, unless parabolic approximation required
    // Correction to 1/sin(theta) for BRL or 1/cos(theta) for ECAP already applied by getCorrectedMaterial() function
    double XtoX0 = m_hits[i]->getCorrectedMaterial().radiation;

    // Correct crossed material accounting circular shape (in non-parabolic approximation)
    if (!(SimParms::getInstance().useParabolicApprox())) XtoX0 *= 1./m_hits[i]->getCosBeta();

    if (XtoX0>0) {

      // Comment (see further in code): MS error depends on total path, i.e. projected path in XY plane scales like 1/sin^2(theta)
      msThetaSq = (13.6*Units::MeV * 13.6*Units::MeV) / (m_p/Units::MeV * m_p/Units::MeV) * XtoX0 * (1 + 0.038 * log(XtoX0)) * (1 + 0.038 * log(XtoX0));
    }
    else {
      msThetaSq = 0;
    }
    msThetaSqProj.push_back(msThetaSq);
  }

  //
  // Calculate the variance matrix with all correlation terms: c is column, r is row (hits are assumed to be sorted)
  for (int c=iStart ; c<=iEnd; c++) {

    double beta_c       = m_hits[c]->getBeta();
    double sinBeta_c    = m_hits[c]->getSinBeta();
    //double cosBeta_c    = m_hits[c]->getCosBeta();
    double rPos_c       = m_hits[c]->getRPos();
    double rPhiLength_c = m_hits[c]->getRPhiLength();
    //double rPosSq_c     = rPos_c*rPos_c;

    // Dummy value for correlations involving inactive surfaces
    if (m_hits[c]->isPassive()) {
      for (int r = 0; r <= c; r++) {

        varMatrixRR(r, c) = 0.0;
        varMatrixRZ(r, c) = 0.0;
        varMatrixZR(r, c) = 0.0;
        varMatrixZZ(r, c) = 0.0;
      }
    }
    // One of the correlation factors refers to an active surface
    else {

      for (int r=iStart; r <= c; r++) {

        // Dummy value for correlation involving an inactive surface
        if (m_hits[r]->isPassive()) {

          varMatrixRR(r, c) = 0.0;
          varMatrixRZ(r, c) = 0.0;
          varMatrixZR(r, c) = 0.0;
          varMatrixZZ(r, c) = 0.0;
        }

        // Correlations between two active surfaces
        else {

          double sumRR        = 0.0;
          double sumRZ        = 0.0;
          double sumZR        = 0.0;
          double sumZZ        = 0.0;

          double beta_r       = m_hits[r]->getBeta();
          double sinBeta_r    = m_hits[r]->getSinBeta();
          //double cosBeta_r    = m_hits[r]->getCosBeta();
          double rPos_r       = m_hits[r]->getRPos();
          double rPhiLength_r = m_hits[r]->getRPhiLength();
          //double rPosSq_r     = rPos_r*rPos_r;

          for (int i=iStart; i<r; i++) {

            double beta_i       = m_hits[i]->getBeta();
            //double sinBeta_i    = m_hits[i]->getSinBeta();
            //double cosBeta_i    = m_hits[i]->getCosBeta();
            double rPos_i       = m_hits[i]->getRPos();
            double rPhiLength_i = m_hits[i]->getRPhiLength();
            //double rPosSq_i     = rPos_i*rPos_i;

            if (SimParms::getInstance().useParabolicApprox()) {

              sumRR += msThetaSqProj[i-(nHits-nHitsUsed)] * (rPos_c-rPos_i) * (rPos_r-rPos_i) /sin(m_theta)/sin(m_theta);                           // Correction to path length
              sumZZ += msThetaSqProj[i-(nHits-nHitsUsed)] * (rPos_c-rPos_i) * (rPos_r-rPos_i) /sin(m_theta)/sin(m_theta)/sin(m_theta)/sin(m_theta); // Correction to path length + error propagation to barrel virtual plane
            }
            else {

              //double trackR = getRadius(0.);

              double beta_ci = beta_c - beta_i;
              double beta_ri = beta_r - beta_i;
              double rPos_ci = rPos_c/sinBeta_c*sin(beta_ci);
              double rPos_ri = rPos_r/sinBeta_r*sin(beta_ri);

              sumRR += msThetaSqProj[i-(nHits-nHitsUsed)] * rPos_ci * cos(beta_ci)     // Project MS in XY plane in direction of curvature @ rPos_c (all with respect to rPos_i)
                                                          * rPos_ri * cos(beta_ri)     // Project MS in XY plane in direction of curvature @ rPos_r (all with respect to rPos_i)
                                                          / sin(m_theta)/sin(m_theta); // MS increases with path length -> projection to XY plane ~1/sin^2(theta)

              sumZZ += msThetaSqProj[i-(nHits-nHitsUsed)] * ( (rPhiLength_c-rPhiLength_i) )
                                                          * ( (rPhiLength_r-rPhiLength_i) )
                                                          / sin(m_theta)/sin(m_theta)/sin(m_theta)/sin(m_theta); // MS increases with path length ~1/sin^2(theta) x error propagation to barrel virtual plane ~1/sin^2(theta);

            }
          }
          if (r == c) {

            double resRPhi = 0;
            double resZ    = 0;
            resRPhi = m_hits[r]->getRphiResolution(getRadius(m_hits[r]->getZPos()));
            resZ    = m_hits[r]->getZResolution(getRadius(m_hits[r]->getZPos()));

            sumRR = sumRR + resRPhi*resRPhi;
            sumZZ = sumZZ + resZ*resZ;

          }
          varMatrixRR(r, c) = sumRR;
          if (r != c) varMatrixRR(c, r) = sumRR;
          varMatrixRZ(r, c) = sumRZ;
          if (r != c) varMatrixRZ(c, r) = sumZR;
          varMatrixZR(r, c) = sumZR;
          if (r != c) varMatrixZR(c, r) = sumRZ;
          varMatrixZZ(r, c) = sumZZ;
          if (r != c) varMatrixZZ(c, r) = sumZZ;
        }
      }
    }
  } // Correlation terms: c is column, r is row

  // Print variance matrix
  //std::cout << "Variance matrix in R-Phi (with zero cols/rows): " << std::endl;
  //printSymMatrix(varMatrixRR);
  //std::cout << "Variance matrix in S-Z (with zero cols/rows): " << std::endl;
  //printSymMatrix(varMatrixZZ);

  //
  // Remove zero rows and columns in covariance matrix
  int  rActual = -1;          // Row, at which to move the active row due to a sequence of zero rows or inactive hits inbetween
  bool lookForActive = false; // Start looking for shift of active rows, after first passive row found

  for (int r=0; r<nHits; r++) {

    // Keep actual row @ zero for first N passive (zero) layers
    if ((m_hits[r]->isPassive() || r<iStart) && (!lookForActive)) {

      // Next hit has to be active (set as active and considered in track fitting (see iStart))
      if ((r+1)<nHits && m_hits[r+1]->isActive() && (r+1)>=iStart) lookForActive = true;

      // Previous hit has to be passive or not being considered in track fitting (see iStart))
      if (!((r-1)>=0 && (m_hits[r-1]->isPassive() || (r-1)<iStart)) ) rActual = r;
    }
    // Shift active layer to zero-th row + i active layers, which have already been shifted by number of zero layers
    else if ((m_hits[r]->isActive()) && (lookForActive)) {
      for (int c=0; c<nHits; c++) {

        varMatrixRR(rActual, c) = varMatrixRR(r, c);
        varMatrixRR(c, rActual) = varMatrixRR(c, r);

        varMatrixRZ(rActual, c) = varMatrixRZ(r, c);
        varMatrixRZ(c, rActual) = varMatrixRZ(c, r);

        varMatrixZR(rActual, c) = varMatrixZR(r, c);
        varMatrixZR(c, rActual) = varMatrixZR(c, r);

        varMatrixZZ(rActual, c) = varMatrixZZ(r, c);
        varMatrixZZ(c, rActual) = varMatrixZZ(c, r);
      }

      varMatrixRR(rActual, rActual) = varMatrixRR(r, r);
      varMatrixRZ(rActual, rActual) = varMatrixRZ(r, r);
      varMatrixZR(rActual, rActual) = varMatrixZR(r, r);
      varMatrixZZ(rActual, rActual) = varMatrixZZ(r, r);

      rActual++;
    }
  }
  // If some rows/colums were zero -> matrix rank needs to be adjusted
  int nDim = rActual;
  if (nDim!=0) {

    varMatrixRR.ResizeTo(nDim, nDim);
    varMatrixRZ.ResizeTo(nDim, nDim);
    varMatrixZR.ResizeTo(nDim, nDim);
    varMatrixZZ.ResizeTo(nDim, nDim);
  }

  varMatrixFull.ResizeTo(2*nDim, 2*nDim);
  for (int iRow=0; iRow<nDim; iRow++) {
    for (int iCol=0; iCol<nDim; iCol++) {

      varMatrixFull(iRow     , iCol)      = varMatrixRR(iRow, iCol);
      varMatrixFull(iRow+nDim, iCol)      = varMatrixRZ(iRow, iCol);
      varMatrixFull(iRow     , iCol+nDim) = varMatrixZR(iRow, iCol);
      varMatrixFull(iRow+nDim, iCol+nDim) = varMatrixZZ(iRow, iCol);

    }
  }

  // Print variance matrix
  //std::cout << "Variance matrix in R-Phi: " << std::endl;
  //printSymMatrix(varMatrixRR);
  //std::cout << "Variance matrix in S-Z: " << std::endl;
  //printSymMatrix(varMatrixZZ);
  //std::cout << "Variance matrix in R-Z: " << std::endl;
  //printSymMatrix(varMatrixRZ);
  //std::cout << "Variance matrix in Z-R: " << std::endl;
  //printSymMatrix(varMatrixZR);
  //std::cout << "Full matrix: " << std::endl;
  //printSymMatrix(varMatrixFull);
  //std::cout << "Det: " << std::scientific << varMatrixFull.Determinant() << std::fixed << std::endl;

  // Check if matrix is sane and worth keeping
  if (!((varMatrixRR.GetNoElements() > 0) && (varMatrixRR.Determinant() != 0.0))) {
    logWARNING("Variance matrix Vrr(NxN) -> zero determinat or zero number of elements");
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); //  Important -> resort back
    return false;
  }
  if (!((varMatrixZZ.GetNoElements() > 0) && (varMatrixZZ.Determinant() != 0.0))) {
    logWARNING("Variance matrix Vzz(NxN) -> zero determinat or zero number of elements");
    if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR); //  Important -> resort back
    return false;
  }

  //
  // Compute covariance matrix of the track parameters
  unsigned int offset = iStart;

  // Get inverse of variance (measurement) matrices (r-phi versus s-z are assumed not to be correlated)
  TMatrixT<double> VRPhiInv(varMatrixRR);   // Local copy of R-Phi variance matrix to be inverted
  VRPhiInv.Invert();
  TMatrixT<double> VSZInv(varMatrixZZ);     // Local copy of S-Z variance matrix to be inverted
  VSZInv.Invert();
  TMatrixT<double> VFullInv(varMatrixFull); // Local copy of full variance matrix to be inverted
  VFullInv.Invert();

  // Derivative matrices 3x3 (R-Phi only), 2x2 (S-Z only), 5x5 (full matrix)
  TMatrixT<double> diffAT(3, nDim);      // Derivatives of track parameters transposed (in R-Phi -> 3 track parameters)
  TMatrixT<double> diffA(nDim, 3);       // Derivatives of track parameters (in R-Phi -> 3 track parameters)

  TMatrixT<double> diffBT(2, nDim);      // Derivatives of track parameters transposed (in R-Phi -> 3 track parameters)
  TMatrixT<double> diffB(nDim, 2);       // Derivatives of track parameters (in R-Phi -> 3 track parameters)

  TMatrixT<double> diffFullT(5, 2*nDim); // Derivatives of track parameters (transposed) in full representation (R-Phi equation)
  TMatrixT<double> diffFull(2*nDim, 5);  // Derivatives of track parameters in full representation (R-Phi equation)

  for (auto i=iStart; i<=iEnd; i++) {

    if (m_hits[i]->isActive()) {

      double hitRPos    = m_hits[i]->getRPos();
      double hitZPos    = m_hits[i]->getZPos();

      if (SimParms::getInstance().useParabolicApprox()) {

        diffA(i - offset, 0) = computeDfOverDRho(hitRPos, hitZPos);
        diffA(i - offset, 1) = m_hits[i]->getRPos();
        diffA(i - offset, 2) = 1;

        diffB(i - offset, 0) = m_hits[i]->getRPos();
        diffB(i - offset, 1) = 1;
      }
      else {

        diffA(i - offset, 0) = computeDfOverDRho(hitRPos, hitZPos);
        diffA(i - offset, 1) = m_hits[i]->getXPos();
        diffA(i - offset, 2) = 1;

        diffFull(i - offset, 0) = computeDfOverDRho(hitRPos, hitZPos);
        diffFull(i - offset, 1) = m_hits[i]->getXPos();
        diffFull(i - offset, 2) = 1;
        diffFull(i - offset, 3) = 0;
        diffFull(i - offset, 4) = 0;

        diffB(i - offset, 0) = m_hits[i]->getRPhiLength();
        diffB(i - offset, 1) = 1;


        diffFull(i - offset+nDim, 0) = 0;
        diffFull(i - offset+nDim, 1) = 0;
        diffFull(i - offset+nDim, 2) = 0;
        diffFull(i - offset+nDim, 3) = m_hits[i]->getRPhiLength();
        diffFull(i - offset+nDim, 4) = 1;
      }
    }
    else offset++;
  }

  // Transpose
  diffAT.Transpose(diffA);
  diffBT.Transpose(diffB);
  diffFullT.Transpose(diffFull);

  // Print A, B matrices
  //std::cout << "DiffA matrix in R-Phi: " << std::endl;
  //printMatrix(diffA);
  //std::cout << "DiffB matrix in S-Z: " << std::endl;
  //printMatrix(diffB);
  //std::cout << "DiffFull matrix: " << std::endl;
  //printMatrix(diffFull);

  // Get covariance matrix using global chi2 fit: C = cov(i,j) = (D^T * V^-1 * D)^-1
  m_covMatrixRPhi.ResizeTo(3,3);
  m_covMatrixRPhi = diffAT * VRPhiInv * diffA;
  m_covMatrixRPhi.Invert();

  m_covMatrixRZ.ResizeTo(2,2);
  m_covMatrixRZ = diffBT * VSZInv * diffB;
  m_covMatrixRZ.Invert();

  if (!SimParms::getInstance().useParabolicApprox()) {

    m_covMatrixFull.ResizeTo(5,5);
    m_covMatrixFull = diffFullT * VFullInv * diffFull;
    m_covMatrixFull.Invert();
  }

  //std::cout << "Full covariance matrix: " << std::endl;
  //printMatrix(m_covMatrixFull);

  // Sort-back hits based on particle direction if they were resorted
  if (m_pt>=0 && !propagOutIn) sortHits(bySmallerR);

  m_covDone = true;
  return m_covDone;
}

//
// Helper fce returning derivative: df(rho, d0, phi0)/drho, where f approximates
// a helix by set of parabolas. In general, N connected parabolas used, for const B
// field only one parabola applied.
//
double Track::computeDfOverDRho(double rPos, double zPos) {

  double DfOverDRho = 0;

  // Option 1: Const mag. field across the detector: Bz = const
  if (SimParms::getInstance().isMagFieldConst()) {

    DfOverDRho = 0.5 * rPos*rPos;
  }
  // Option 2: Mag. field is a function in Z: B = B(z)
  else {

    int nRegions = SimParms::getInstance().getNMagFieldRegions();

    // Find i-th region corresponding to the current zPos
    int iRegion = 0;
    for (iRegion=0; iRegion < nRegions; iRegion++) {

      if (zPos<(SimParms::getInstance().magFieldZRegions[iRegion])) break;
    }

    // Check that zPos not beyond Z-range, in which B field has been defined
    if (iRegion==nRegions) {

      std::ostringstream message;
      message << "Track::computeDfOverDRho(): Hit z-position: " << zPos/Units::mm << " beyond defined B field Z-range: [0," << SimParms::getInstance().magFieldZRegions[nRegions-1]/Units::mm << "]!";
      logERROR(message.str());
      exit(1);
    }

    // Z pos. in the first region or only 1 region defined (const mag. field)
    if (iRegion==0) {
      DfOverDRho = 0.5 * rPos*rPos;
    }
    // Z pos in i-th region (generally N regions defined)
    else {

      // Get reference magnetic field B0 (at [r,z] = [0,0]
      double B0   = SimParms::getInstance().magField[0];

      double Bi   = 0.; // B-field in ith z-region
      double Bi_1 = 0.; // B-field in (i-1)th z-region
      double xi   = 0.; // x-position corresponding to the ith z-region

      // Sum-up all contributions across the regions: 0th - ith
      for (int i=1; i<=iRegion; i++) {

        // Get current value of magnetic field B_i & B_i-1
        Bi   = SimParms::getInstance().magField[i];
        Bi_1 = SimParms::getInstance().magField[i-1];
        xi   = SimParms::getInstance().magFieldZRegions[i-1] * tan(m_theta); // (z0,z1,z2...) -> intervals defined as 0-z0, z0-z1, z1-z2

        // Add dB/dz terms
        DfOverDRho+= -1.0*(Bi-Bi_1)/B0 * xi*rPos;
        DfOverDRho+= +0.5*(Bi-Bi_1)/B0 * xi*xi;

        // Add Bi/B0 term
        if (i==iRegion) DfOverDRho += +0.5*Bi/B0 * rPos*rPos;
      }
    }
  }

  return DfOverDRho;
}

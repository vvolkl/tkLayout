/*
 * VisitorMatTrack.cc
 *
 *  Created on: 17. 1. 2017
 *      Author: drasal
 */
#include "../include/VisitorMatTrack.h"

#include "Barrel.h"
#include "BeamPipe.h"
#include "DetectorModule.h"
#include "Disk.h"
#include "Endcap.h"
#include "InactiveElement.h"
#include "Layer.h"
#include "MaterialProperties.h"
#include "ModuleCap.h"
#include "SupportStructure.h"
#include "Track.h"
#include "Tracker.h"

//
// Material track visitor - constructor
//
VisitorMatTrack::VisitorMatTrack(Track& matTrack) :
    m_matTrack(matTrack)
{
  m_detName  = "Undefined";
  m_brlName  = "Undefined";
  m_ecapName = "Undefined";
  m_layerID  = -1;
  m_discID   = -1;
}

//
// Destructor
//
VisitorMatTrack::~VisitorMatTrack()
{
}

//
// Visit BeamPipe -> update track with beam pipe hit
//
void VisitorMatTrack::visit(const BeamPipe& bp)
{
  //
  // Add hit corresponding with beam-pipe
  double theta = m_matTrack.getTheta();
  //double eta   = m_matTrack.getEta();
  double z0    = m_matTrack.getOrigin().Z();

  // Eta assumed to be positive so, theta between 0-90deg, hence material is always given as firstUpper-firstLower hit, secondUpper-secondLower etc. -> sort from lowest to highest Z
  std::vector<std::pair<double,double>> bpZRLower;
  std::vector<std::pair<double,double>> bpZRUpper;
  std::vector<std::pair<double,double>> bpRadIntLength;

  // Check all segments from negative to positive part of the beam-pipe & create crossing points with beam-pipe contours (lower & upper)
  // Upper contour is asssumed to close-up the beam-pipe segments containing a certain material, e.g. beam-pipe constructed from central Be & aluminium flanges
  // If beam-pipe consists of 1 material only, the upper contour will close-up just the end of beam-pipe, so that particles at very shallow angle are forced to
  // create two cross-points, i.e. 1 hit point!
  for (signed int i=-bp.getNSegments(); i<bp.getNSegments(); i++) {

    float sgni = (i>=0) ? +1 : -1;

    double zi         = sgni*bp.zPos[abs(i)];
    double zipls1     = sgni*bp.zPos[abs(i+1)];
    double riLower    = bp.radiusLower[abs(i)];
    double ripls1Lower= bp.radiusLower[abs(i+1)];
    double riUpper    = bp.radiusUpper[abs(i)];
    double ripls1Upper= bp.radiusUpper[abs(i+1)];
    double tiltLower  = (sgni<0) ? -bp.getTiltLower(abs(i+1)) : +bp.getTiltLower(abs(i));
    double tiltUpper  = (sgni<0) ? -bp.getTiltUpper(abs(i+1)) : +bp.getTiltUpper(abs(i));

    double rPosLower = -1;
    double rPosUpper = -1;
    double zPosLower = -1;
    double zPosUpper = -1;

    // Tilt differen from +-pi/2 (@ pi/2 formulae based on zPos diverge)
    if (zi!=zipls1) {

      zPosUpper = (tan(theta)*z0 - tan(tiltUpper)*zi + riUpper)/(tan(theta) -tan(tiltUpper));
      zPosLower = (tan(theta)*z0 - tan(tiltLower)*zi + riLower)/(tan(theta) -tan(tiltLower));

      if ( ((zPosLower>=0) && (zPosLower>=zi) && (zPosLower<zipls1)) ||
         ((zPosLower<0)  && (zPosLower<zi)  && (zPosLower>=zipls1)) ) rPosLower = tan(theta)*(zPosLower-z0);

      if ( ((zPosUpper>=0) && (zPosUpper>=zi) && (zPosUpper<zipls1)) ||
         ((zPosUpper<0)  && (zPosUpper<zi)  && (zPosUpper>=zipls1)) ) rPosUpper = tan(theta)*(zPosUpper-z0);
    }
    // Tilt = +-pi/2, formulae based on zPos diverge -> use different approach
    else {

      double rPos = tan(theta)*(zi-z0);

      if ( (rPos>=MIN(riLower,ripls1Lower)) && (rPos<MAX(riLower,ripls1Lower)) ) {

        zPosLower = zi;
        rPosLower = rPos;
      }
      if ( (rPos>=MIN(riUpper,ripls1Upper)) && (rPos<MAX(riUpper,ripls1Upper)) ) {

        zPosUpper = zi;
        rPosUpper = rPos;
      }
    }

    // Save individual points of crossings with beam-pipe contours
    if (rPosLower!=-1 || rPosUpper!=-1) {

      bpZRUpper.push_back(std::pair<double,double>(zPosUpper,rPosUpper));
      bpZRLower.push_back(std::pair<double,double>(zPosLower,rPosLower));

      // Material
      //
      // Positive z -> use directly index i to find material of corresponding segment
      if (i>=0) bpRadIntLength.push_back(std::pair<double,double>(bp.radLength[i],bp.intLength[i]));
      // Negative z -> use abs(i) -1 (indicis start at -(N-segments+1) and end-up at (-1) to find material of corresponding segment
      else      bpRadIntLength.push_back(std::pair<double,double>(bp.radLength[abs(i)-1],bp.intLength[abs(i)-1]));
    }
  }

  //
  // Form hits from crossed points -> hit defined by average of two consecutive crossed points (maybe cross-points from lower & upper
  // contours or lower only/upper only contours)
  //
  // Check that cross-points algorithm created points correctly first
  if (bpZRLower.size()!=bpZRUpper.size()) {

    std::ostringstream message;
    message << "Beam-pipe algorithm assigning material not working correctly, check the code!";
    throw PathfulException(message.str(), "VisitorMatBudget");
  }
  // Form hits
  else {

    for (unsigned int i=0; i<bpZRLower.size(); i++) {

      double pathLength  = 0;
      double rPosAvg     = 0;
      double zPosAvg     = 0;
      double radMaterial = -1;
      double intMaterial = -1;

      double zPosUpper = bpZRUpper[i].first;
      double zPosLower = bpZRLower[i].first;
      double rPosUpper = bpZRUpper[i].second;
      double rPosLower = bpZRLower[i].second;
      double radLength = bpRadIntLength[i].first;
      double intLength = bpRadIntLength[i].second;

      // Upper & Lower crossed
      if (rPosLower!=-1 && rPosUpper!=-1) {

        pathLength  = sqrt(pow(rPosUpper-rPosLower,2)+pow(zPosUpper-zPosLower,2));
        rPosAvg     = (rPosUpper+rPosLower)/2.;
        zPosAvg     = (zPosUpper+zPosLower)/2.;
        radMaterial = pathLength/radLength;
        intMaterial = pathLength/intLength;

      }
      // Upper or lower & next(next) upper or lower crossed -> update current point index by the found next(next) crossing point
      else if ( (rPosLower==-1 && rPosUpper!=-1) || (rPosLower!=-1 && rPosUpper==-1)) {

        // Look for next, next-next, next-next-next etc., set once found!
        bool found = false;
        while (!found) {

          i++; // Get next point

          if (i>=bpZRLower.size()) {
            logERROR("BeamPipe not defined as a closed object, can't find all hit part!!!");
            break;
          }
          else {

            // Only one of these must be different from -1
            double zPosUpperNext = bpZRUpper[i].first;
            double zPosLowerNext = bpZRLower[i].first;
            double rPosUpperNext = bpZRUpper[i].second;
            double rPosLowerNext = bpZRLower[i].second;
            double radLengthNext = bpRadIntLength[i].first;
            double intLengthNext = bpRadIntLength[i].second;

            // Different combinations may appear, but at least one of these: upper with upper or lower with upper or lower with lower
            if      (rPosUpperNext!=-1 && rPosLowerNext!=-1) logERROR("Algorithmic problem when looking for material hits, check the code!!!");
            else if (rPosUpperNext!=-1) {

              if (rPosLower!=-1) {
                pathLength = sqrt(pow(rPosUpperNext-rPosLower,2)+pow(zPosUpperNext-zPosLower,2));
                rPosAvg    = (rPosUpperNext+rPosLower)/2.;
                zPosAvg    = (zPosUpperNext+zPosLower)/2.;
                found      = true;
              }
              else {
                pathLength = sqrt(pow(rPosUpperNext-rPosUpper,2)+pow(zPosUpperNext-zPosUpper,2));
                rPosAvg    = (rPosUpperNext+rPosUpper)/2.;
                zPosAvg    = (zPosUpperNext+zPosUpper)/2.;
                found      = true;
              }
            }
            else if (rPosLowerNext!=-1) {

              if (rPosLower!=-1) {
                pathLength = sqrt(pow(rPosLowerNext-rPosLower,2)+pow(zPosLowerNext-zPosLower,2));
                rPosAvg    = (rPosLowerNext+rPosLower)/2.;
                zPosAvg    = (zPosLowerNext+zPosLower)/2.;
                found      = true;
              }
              else {
                pathLength = sqrt(pow(rPosLowerNext-rPosUpper,2)+pow(zPosLowerNext-zPosUpper,2));
                rPosAvg    = (rPosLowerNext+rPosUpper)/2.;
                zPosAvg    = (zPosLowerNext+zPosUpper)/2.;
                found      = true;
              }
            }

            // Check that material consistently defined: crossed-points forming the potential hit must be within the same material
            if (radLength!=radLengthNext || intLength!=intLengthNext) {
              logERROR("Algorithmic problem when looking for material hits: 2 crossed-points forming a potential hit not within the same material -> Check that Beam-pipe upper contour closing object & material definition!!!");
            }
            else {
              radMaterial = pathLength/radLength;
              intMaterial = pathLength/intLength;
            }
          }
        } // While
      }
      else {

        continue;
      }

      HitPtr hit(new Hit(rPosAvg, zPosAvg, nullptr, HitPassiveType::BeamPipe));

      Material material;
      material.radiation   = radMaterial;
      material.interaction = intMaterial;

      hit->setCorrectedMaterial(material);
      m_matTrack.addHit(std::move(hit));

    } // For
  } // Form hits algorithm
}

//
// Visit Tracker
//
void VisitorMatTrack::visit(const Tracker& t)
{
  m_detName = t.myid();
  m_layerID    = -1;
  m_discID     = -1;
}

//
// Visit Barrel
//
void VisitorMatTrack::visit(const Barrel& b)
{
  m_brlName  = b.myid();
  m_layerID  = -1;
  m_discID   = -1;

  for (auto& element : b.services()) {
    analyzeInactiveElement(element);
  }
}

//
// Visit Endcap
//
void VisitorMatTrack::visit(const Endcap& e)
{
  m_ecapName = e.myid();
  m_layerID  = -1;
  m_discID   = -1;
}

//
// Visit Layer
//
void VisitorMatTrack::visit(const Layer& l)
{
  m_layerID = l.myid();
  m_discID  = -1;
}

//
// Visit Disc
//
void VisitorMatTrack::visit(const Disk& d)
{
  m_layerID = -1;
  m_discID  = d.myid();
}

//
// Visit BarrelModule (no limits on Rods, Layers or Barrels)
//
void VisitorMatTrack::visit(const BarrelModule& m)
{
  analyzeModuleMB(m);
}

//
// Visit EndcapModule (no limits on Rings or Endcaps)
//
void VisitorMatTrack::visit(const EndcapModule& m)
{
  analyzeModuleMB(m);
}

//
// Visit Support strucutre
//
void VisitorMatTrack::visit(const SupportStructure& s)
{
  for (auto& elem : s.inactiveElements()) {
    analyzeInactiveElement(elem);
  }
}

//
// Analyze if module crossed by given track & how much material is in the way
//
void VisitorMatTrack::analyzeModuleMB(const DetectorModule& m)
{
  // Collision detection: material tracks being shot in z+ only, so consider only modules that lie on +Z side after correction on primary vertex origin
  if ((m.maxZ()-m_matTrack.getOrigin().Z())>0) {

    XYZVector     direction(m_matTrack.getDirection());
    Material      hitMaterial;
    XYZVector     hitPosition;
    HitModuleType hitModuleType;

    if (m.checkTrackHits(m_matTrack.getOrigin(), direction, hitMaterial, hitModuleType, hitPosition)) {

      auto hitRPos = hitPosition.rho();
      auto hitZPos = hitPosition.z();

      // Create Hit object with appropriate parameters, add to Track t
      HitPtr hit(new Hit(hitRPos, hitZPos, &m, hitModuleType));
      hit->setCorrectedMaterial(hitMaterial);
      if (m.subdet() == BARREL) {
        hit->setDetName(m_detName+"_"+m_brlName);
        hit->setLayerID(m_layerID);
      }
      else {
        hit->setDetName(m_detName+"_"+m_ecapName);
        hit->setDiscID(m_discID);
      }
      m_matTrack.addHit(std::move(hit));
    }
  } // (module_Z-orig_Z)>0
}

//
// Helper method - analyse inactive element & estimate how much material is in the way
//
void VisitorMatTrack::analyzeInactiveElement(const InactiveElement& e)
{

  XYZVector direction(m_matTrack.getDirection());
  Material  hitMaterial;
  XYZVector hitPosition;

  // Hit found -> create
  if (e.checkTrackHits(m_matTrack.getOrigin(),direction, hitMaterial, hitPosition)) {

    auto hitRPos = hitPosition.rho();
    auto hitZPos = hitPosition.z();

    // Create Hit object with appropriate parameters, add to Track t
    if (e.getCategory() == MaterialProperties::b_sup ||
        e.getCategory() == MaterialProperties::e_sup ||
        e.getCategory() == MaterialProperties::u_sup ||
        e.getCategory() == MaterialProperties::t_sup ) {

      HitPtr hit(new Hit(hitRPos, hitZPos, &e, HitPassiveType::Support));
      hit->setCorrectedMaterial(hitMaterial);
      m_matTrack.addHit(std::move(hit));
    }
    else if (e.getCategory() == MaterialProperties::b_ser ||
             e.getCategory() == MaterialProperties::e_ser ) {

      HitPtr hit(new Hit(hitRPos, hitZPos, &e, HitPassiveType::Service));
      hit->setCorrectedMaterial(hitMaterial);
      m_matTrack.addHit(std::move(hit));
    }
  } // Hit found
}




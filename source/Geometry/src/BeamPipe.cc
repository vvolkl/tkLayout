/*
 *  BeamPipe.cc
 *
 *  Created on: 11. 5. 2016
 *      Author: drasal
 */

#include "BeamPipe.h"

#include "math_functions.h"

#include "Math/DisplacementVector3D.h"
#include "Math/Cylindrical3D.h"
#include "InactiveTube.h"
#include "MessageLogger.h"
#include "SimParms.h"
#include "Visitor.h"

//
// Constructor -> object set through build method
//
BeamPipe::BeamPipe(const PropertyTree& treeProperty) :
 radiusUpper(     "radiusUpper", parsedOnly()),
 radiusLower(     "radiusLower", parsedOnly()),
 radius(          "radius"     , parsedOnly()),
 thickness(       "thickness"  , parsedOnly()),
 minRadius(string("minRadius")),
 maxRadius(string("maxRAdius")),
 zPos(            "zPos"       , parsedOnly()),
 radLength(       "radLength"  , parsedOnly()),
 intLength(       "intLength"  , parsedOnly())
{
  // Set the geometry config parameters
  this->myid(treeProperty.data());
  this->store(treeProperty);
}

BeamPipe::~BeamPipe()
{
  //m_tube.reset();
}

//
// Buidl method - setting all parameters
//
void BeamPipe::build()
{
  try {
    check();

    // Build beam pipe as inactive material
//    double zLength = 2*maxZ[maxZ.size()-1];
//    double zOffset = 0.0;
//
//    m_tube.reset(new InactiveTube(-zLength/2.,zLength,radius(),thickness()));
//    m_tube->setRadiationLength(radLength());
//    m_tube->setInteractionLength(intLength());
  }
  catch (PathfulException& pe) {

    throw;
  }

  cleanup();
  builtok(true);
}

//
// Cross-check parameters provided from geometry configuration file
//
void BeamPipe::check() {

  PropertyObject::check();

  // Check for beam-pipe defined in 3D version
  if (SimParms::getInstance().use3DBeamPipe()) {

    if (radiusLower.size()<2) {

      std::ostringstream message;
      message << "Number of defined beam-pipe points for lower contour " << radiusLower.size() << " must be at least 2 in order to define at least 1 segment.";
      message << " Check whether not defining radius only for 2DBeamPipe, but assuming 3DBeamPipe in SimParms config file!";

      throw PathfulException(message.str(),"BeamPipe3D");
    }
    if (radiusLower.size()!=zPos.size()) {
      throw PathfulException("Number of defined beam-pipe segments not the same in lower radius & z position!","BeamPipe3D");
    }
    if (radiusUpper.size()!=zPos.size()) {
      throw PathfulException("Number of defined beam-pipe segments not the same in  radius & z position!","BeamPipe3D");
    }
    if (radLength.size()!=(zPos.size()-1)) {
      throw PathfulException("Number of defined rad. length parameters must equal to number of segments, i.e. number of radii/zpos points minus one!","BeamPipe3D");
    }
    if (intLength.size()!=(zPos.size()-1)) {
      throw PathfulException("Number of defined int. length parameters must equal to number of segments, i.e. number of radii/zpos points minus one!","BeamPipe3D");
    }
  }
  // otherwise 2D version expected
  else {

    if (radius.size()<2) {

      std::ostringstream message;
      message << "Number of defined beam-pipe points " << radius.size() << " must be at least 2 in order to define at least 1 segment.";
      message << " Check whether not defining lower & upper radius for 3DBeamPipe, but assuming 2DBeamPipe in SimParms config file!";

      throw PathfulException(message.str(),"BeamPipe2D");
    }
    if (radius.size()!=zPos.size()) {
      throw PathfulException("Number of defined beam-pipe segments not the same in radius & z position!","BeamPipe2D");
    }
    if (thickness.size()!=(zPos.size()-1)) {
      throw PathfulException("Number of thickness parameters must equal to number of segments, i.e. number of radii/zpos points minus one!","BeamPipe2D");
    }
    if (radLength.size()!=(zPos.size()-1)) {
      throw PathfulException("Number of defined rad. length parameters must equal to number of segments, i.e. number of radii/zpos points minus one!","BeamPipe2D");
    }
    if (intLength.size()!=(zPos.size()-1)) {
      throw PathfulException("Number of defined int. length parameters must equal to number of segments, i.e. number of radii/zpos points minus one!","BeamPipe2D");
    }
  }
}

//
// Setup: link lambda functions to various beampipe related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
//
void BeamPipe::setup() {

  maxRadius.setup([&]() { double max = 0;

    if (SimParms::getInstance().use3DBeamPipe()) { for (unsigned int i=0; i<radiusUpper.size(); i++) { max = MAX(max, radiusUpper[i]); };}
    else                                         { for (unsigned int i=0; i<radius.size(); i++)      { max = MAX(max, radius[i]); }; }
    return max; });
  minRadius.setup([&]() { double min = std::numeric_limits<double>::max();

    if (!SimParms::getInstance().use3DBeamPipe()){ for (unsigned int i=0; i<radiusLower.size(); i++) { min = MIN(min, radiusLower[i]); };}
    else                                         { for (unsigned int i=0; i<radius.size(); i++)      { min = MIN(min, radius[i]); }; }
    return min; });
}

//
// Check if track hit the beam-pipe -> if yes, return true with passed material & hit position vector
//
bool BeamPipe::checkTrackHits(const ROOT::Math::XYZVector& trackOrig, const ROOT::Math::XYZVector& trackDir, std::vector<RILength*>& hitMaterialVec, std::vector<ROOT::Math::RhoZPhiVector*>& hitPosVec) const {

  // Hit containers
  Material*                  hitMaterial = nullptr;
  ROOT::Math::RhoZPhiVector* hitPos      = nullptr;

  // Found at least one hit?
  bool   hitFound = false;

  // Track parameters
  double theta    = trackDir.Theta();
  double z0       = trackOrig.Z();

  //
  // Beam-pipe defined approximatively as a 2D object only (phi-symmetric object with thickness  projected)
  if (!SimParms::getInstance().use3DBeamPipe()) {

    // Check all segments in positive part of the beam-pipe
    for (unsigned int i=0; i<getNSegments(); i++) {

      // Initialize
      hitMaterial = nullptr;
      hitPos      = nullptr;

      double zi       = zPos[i];
      double zipls1   = zPos[i+1];
      double ri       = radius[i];
      double tilt     = getAvgTilt(i);
      double thickness= getAvgThickness(i);

      double zPos = (tan(theta)*z0 - tan(tilt)*zi + ri)/(tan(theta) -tan(tilt));

      if ( ((zPos>=0) && (zPos>=zi) && (zPos<zipls1)) ||
           ((zPos<0)  && (zPos<zi)  && (zPos>=zipls1)) ) {

        double rPos = tan(theta)*(zPos-z0);

        hitMaterial = new Material();
        hitMaterial->radiation   = thickness/radLength(i)/fabs(sin(theta-tilt));
        hitMaterial->interaction = thickness/intLength(i)/fabs(sin(theta-tilt));

        hitPos = new ROOT::Math::RhoZPhiVector();
        hitPos->SetZ(zPos);
        hitPos->SetRho(rPos);
        hitPos->SetPhi(0); // Assumed phi symmetry, why not e.g. phi=0 for completness

        hitMaterialVec.push_back(hitMaterial);
        hitPosVec.push_back(hitPos);

        hitFound = true;
      }
    }

    // Check all segments in negative part of the beam-pipe
    for (unsigned int i=0; i<getNSegments(); i++) {

      // Initialize
      hitMaterial = nullptr;
      hitPos      = nullptr;

      double zi       = -zPos[i];
      double zipls1   = -zPos[i+1];
      double ri       = +radius[i];
      double tilt     = -getAvgTilt(i);
      double thickness= +getAvgThickness(i);

      double zPos = (tan(theta)*z0 - tan(tilt)*zi + ri)/(tan(theta) -tan(tilt));

      if ( ((zPos>=0) && (zPos>=zi) && (zPos<zipls1)) ||
           ((zPos<0)  && (zPos<zi)  && (zPos>=zipls1)) ) {

        double rPos = tan(theta)*(zPos-z0);

        hitMaterial = new Material();
        hitMaterial->radiation   = thickness/radLength(i)/fabs(sin(theta-tilt));
        hitMaterial->interaction = thickness/intLength(i)/fabs(sin(theta-tilt));

        hitPos = new ROOT::Math::RhoZPhiVector();
        hitPos->SetZ(zPos);
        hitPos->SetRho(rPos);
        hitPos->SetPhi(0); // Assumed phi symmetry, why not e.g. phi=0 for completness

        hitMaterialVec.push_back(hitMaterial);
        hitPosVec.push_back(hitPos);

        hitFound = true;
      }
    }
  } // If beam-pipe defined approximately like 2D object

  //
  // Beam-pipe defined as a 3D object: phi-symmetric object defined by 2 contours upper & lower, 2 contours being defined by a sequence of points
  else {

    // Eta assumed to be positive so, theta between 0-90deg, hence material is always given as firstUpper-firstLower hit, secondUpper-secondLower etc. -> sort from lowest to highest Z
    std::vector<std::pair<double,double>> bpZRLower;
    std::vector<std::pair<double,double>> bpZRUpper;
    std::vector<std::pair<double,double>> bpRadIntLength;

    // Check all segments from negative to positive part of the beam-pipe & create crossing points with beam-pipe contours (lower & upper)
    // Upper contour is asssumed to close-up the beam-pipe segments containing a certain material, e.g. beam-pipe constructed from central Be & aluminium flanges
    // If beam-pipe consists of 1 material only, the upper contour will close-up just the end of beam-pipe, so that particles at very shallow angle are forced to
    // create two cross-points, i.e. 1 hit point!
    for (signed int i=-getNSegments(); i<getNSegments(); i++) {

      // Initialize
      hitMaterial = nullptr;
      hitPos      = nullptr;

      float sgni = (i>=0) ? +1 : -1;

      double zi         = sgni*zPos[abs(i)];
      double zipls1     = sgni*zPos[abs(i+1)];
      double riLower    = radiusLower[abs(i)];
      double ripls1Lower= radiusLower[abs(i+1)];
      double riUpper    = radiusUpper[abs(i)];
      double ripls1Upper= radiusUpper[abs(i+1)];
      double tiltLower  = (sgni<0) ? -getTiltLower(abs(i+1)) : +getTiltLower(abs(i));
      double tiltUpper  = (sgni<0) ? -getTiltUpper(abs(i+1)) : +getTiltUpper(abs(i));

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
        if (i>=0) bpRadIntLength.push_back(std::pair<double,double>(radLength[i],intLength[i]));
        // Negative z -> use abs(i) -1 (indicis start at -(N-segments+1) and end-up at (-1) to find material of corresponding segment
        else      bpRadIntLength.push_back(std::pair<double,double>(radLength[abs(i)-1],intLength[abs(i)-1]));
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

        hitMaterial = new Material();
        hitMaterial->radiation   = radMaterial;
        hitMaterial->interaction = intMaterial;

        hitPos = new ROOT::Math::RhoZPhiVector();
        hitPos->SetZ(zPosAvg);
        hitPos->SetRho(rPosAvg);
        hitPos->SetPhi(0); // Assumed phi symmetry, why not e.g. phi=0 for completness

        hitMaterialVec.push_back(hitMaterial);
        hitPosVec.push_back(hitPos);

        hitFound = true;
      } // For
    } // Form hits algorithm

  } // If/Else 2D/3D beam-pipe geometry

  return hitFound;
}

//
// Get radius of upper wall (contour) of i-th segment
//
double BeamPipe::getRadiusUpper(unsigned int i) const {

  double radiusValue = -1;

  if (i>getNSegments()) {

    std::ostringstream message;
    message << "Upper radius not defined for beam-pipe countour point i=" << i;
    throw PathfulException(message.str(), "BeamPipe");
  }
  else {

    // 3D version of beam-pipe geometry
    if (SimParms::getInstance().use3DBeamPipe()) radiusValue = radiusUpper[i];
    else                                         radiusValue = radius[i];
  }

  return radiusValue;
}

//
// Get radius of lower wall (contour) of i-th segment
//
double BeamPipe::getRadiusLower(unsigned int i) const {

  double radiusValue = -1;

  if (i>getNSegments()) {

    std::ostringstream message;
    message << "Lower radius not defined for beam-pipe countour point i=" << i;
    throw PathfulException(message.str(), "BeamPipe");
  }
  else {

    // 3D version of beam-pipe geometry
    if (SimParms::getInstance().use3DBeamPipe()) radiusValue = radiusLower[i];
    else                                         radiusValue = radius[i];
  }

  return radiusValue;
}

//
// Get tilt of upper wall of i-th segment
//
double BeamPipe::getTiltUpper(unsigned int i) const {

  double tiltAngle = -1;

  if (i>=getNSegments()) {

    std::ostringstream message;
    message << "Tilt not defined for segment indexed by required i=" << i;
    throw PathfulException(message.str(), "BeamPipe");
  }
  else {

    // 3D version of beam-pipe geometry
    if (SimParms::getInstance().use3DBeamPipe()) {
      double radiusIPls1 = radiusUpper[i+1];
      double radiusI     = radiusUpper[i]  ;

      if (radiusIPls1==radiusI) tiltAngle = 0;
      else                      tiltAngle = atan((radiusIPls1-radiusI)/(zPos[i+1]-zPos[i]));
    }
    // Simplified 2D version of beam-pipe
    else {

      double radiusIPls1 = radius[i+1];
      double radiusI     = radius[i]  ;

      if (radiusIPls1==radiusI) tiltAngle = 0;
      else                      tiltAngle = atan((radiusIPls1-radiusI)/(zPos[i+1]-zPos[i]));
    }
  }

  return tiltAngle;
}

//
// Get tilt of lower wall of i-th segment
//
double BeamPipe::getTiltLower(unsigned int i) const {

  double tiltAngle = -1;

  if (i>=getNSegments()) {

    std::ostringstream message;
    message << "Tilt not defined for segment indexed by required i=" << i;
    throw PathfulException(message.str(), "BeamPipe");
  }
  else {

    // 3D version of beam-pipe geometry
    if (SimParms::getInstance().use3DBeamPipe()) {
      double radiusIPls1 = radiusLower[i+1];
      double radiusI     = radiusLower[i]  ;

      if (radiusIPls1==radiusI) tiltAngle = 0;
      else                      tiltAngle = atan((radiusIPls1-radiusI)/(zPos[i+1]-zPos[i]));
    }
    // Simplified 2D version of beam-pipe
    else {

      double radiusIPls1 = radius[i+1];
      double radiusI     = radius[i]  ;

      if (radiusIPls1==radiusI) tiltAngle = 0;
      else                      tiltAngle = atan((radiusIPls1-radiusI)/(zPos[i+1]-zPos[i]));
    }
  }

  return tiltAngle;
}

//
// Get average thickness of i-th segment
//
double BeamPipe::getAvgThickness(unsigned int i) const {

  double avgThickness = -1;

  if (i>=getNSegments()) {

    std::ostringstream message;
    message << "Thickness not defined for segment indexed by required i=" << i;
    throw PathfulException(message.str(), "BeamPipe");
  }
  else {

    if (SimParms::getInstance().use3DBeamPipe()) avgThickness = (radiusUpper[i]+radiusUpper[i+1]-radiusLower[i]-radiusLower[i+1])/2;
    else                                         avgThickness = thickness[i];
  }

  return avgThickness;
}


//
// GeometryVisitor pattern -> beam pipe visitable
//
void BeamPipe::accept(GeometryVisitor& v)
{
  v.visit(*this);
}

//
// GeometryVisitor pattern -> beam pipe visitable (const. option)
//
void BeamPipe::accept(ConstGeometryVisitor& v) const
{
  v.visit(*this);
}

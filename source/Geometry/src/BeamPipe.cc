/*
 * BeamPipe.cc
 *
 *  Created on: 11. 5. 2016
 *      Author: drasal
 */

#include "BeamPipe.h"

#include "math_functions.h"
#include "InactiveTube.h"
#include "Visitor.h"

//
// Constructor -> object set through build method
//
BeamPipe::BeamPipe(const PropertyTree& treeProperty) :
 radiusUpper(     "radiusUpper", parsedOnly()),
 radiusLower(     "radiusLower", parsedOnly()),
 minRadius(string("minRadius")),
 maxRadius(string("maxRAdius")),
 zPos(            "zPos"     , parsedOnly()),
 radLength(       "radLength", parsedOnly()),
 intLength(       "intLength", parsedOnly())
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

  // Check that number of defined beam-pipe segments the same in radius & Z-length
  if (radiusLower.size()!=zPos.size()) {
    throw PathfulException("Number of defined beam-pipe segments not the same in lower radius & z position!","BeamPipe");
  }
  if (radiusUpper.size()!=zPos.size()) {
    throw PathfulException("Number of defined beam-pipe segments not the same in  radius & z position!","BeamPipe");
  }
  if (radiusLower.size()<2) {
    throw PathfulException("Number of defined beam-pipe points must be at least 2 in order to define at least 1 segment!","BeamPipe");
  }
  if (radLength.size()!=(zPos.size()-1)) {
    throw PathfulException("Number of defined rad. length parameters must equal to number of segments, i.e. number of radii/zpos points minus one!","BeamPipe");
  }
  if (intLength.size()!=(zPos.size()-1)) {
    throw PathfulException("Number of defined int. length parameters must equal to number of segments, i.e. number of radii/zpos points minus one!","BeamPipe");
  }
}

//
// Setup: link lambda functions to various beampipe related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
//
void BeamPipe::setup() {

  maxRadius.setup([&]() { double max = 0;                                  for (auto i=0; i<radiusUpper.size(); i++) { max = MAX(max, radiusUpper[i]); } return max; });
  minRadius.setup([&]() { double min = std::numeric_limits<double>::max(); for (auto i=0; i<radiusLower.size(); i++) { min = MIN(min, radiusLower[i]); } return min; });;
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

    double radiusIPls1 = radiusUpper[i+1];
    double radiusI     = radiusUpper[i]  ;

    if (radiusIPls1==radiusI) tiltAngle = 0;
    else                      tiltAngle = atan((radiusIPls1-radiusI)/(zPos[i+1]-zPos[i]));
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

    double radiusIPls1 = radiusLower[i+1];
    double radiusI     = radiusLower[i]  ;

    if (radiusIPls1==radiusI) tiltAngle = 0;
    else                      tiltAngle = atan((radiusIPls1-radiusI)/(zPos[i+1]-zPos[i]));
  }

  return tiltAngle;
}

//
// Get average thickness of i-th segment
//
double BeamPipe::getAvgThickness(unsigned int i) const {

  double thickness = -1;

  if (i>=getNSegments()) {

    std::ostringstream message;
    message << "Thickness not defined for segment indexed by required i=" << i;
    throw PathfulException(message.str(), "BeamPipe");
  }
  else {

    thickness = (radiusUpper[i]+radiusUpper[i+1]-radiusLower[i]-radiusLower[i+1])/2;
  }

  return thickness;
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

/*
 * Detector.cc
 *
 *  Created on: 25. 7. 2016
 *      Author: Z.Drasal (CERN)
 */
#include "Detector.h"

#include "BeamPipe.h"
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include "InactiveSurfaces.h"
#include "Materialway.h"
#include "MessageLogger.h"
#include "StopWatch.h"
#include "Tracker.h"

//
// Default constructor
//
Detector::Detector(std::string name) : m_name(name)
{}

//
// Default destructor
//
Detector::~Detector()
{
  m_activeTrackers.clear();
  m_passiveTrackers.clear();
  m_beamPipe.reset();
}

//
// Build active part of the tracker (a bare-bones geometry of active modules). The resulting tracker object consists of individual sub-trackers
// (pixel, strip, central, forward, ...). This procedure replaces the previously registered tracker (applying correct memory managment), if such an object
// existed.
// Return True if there were no errors during processing, false otherwise
//
bool Detector::buildActiveTracker(boost::property_tree::ptree& geomTree)
{
  // Clear out content of trackers container
  m_activeTrackers.clear();

  // Build active tracker -> look for tag "Tracker"
  startTaskClock("Building active trackers");
  try {

    auto childRange = getChildRange(geomTree, "Tracker");
    for (auto it=childRange.first; it!=childRange.second; ++it) {

      // Message
      if (it==childRange.first) std::cout << std::endl;
      std::cout << " " << it->second.data() << " tracker" << std::endl;

      // Create tracker
      std::unique_ptr<Tracker> trk(GeometryFactory::make<Tracker>(it->second));
      trk->build();

      m_activeTrackers.push_back(std::move(trk));

    } // Trackers

    // Declare unmatched properties
    std::set<string> unmatchedProperties = PropertyObject::reportUnmatchedProperties();
    if (!unmatchedProperties.empty()) {

      std::ostringstream message;
      message << "Detector::buildActiveTracker(): The following unknown properties were ignored:" << std::endl;

      for (const string& property  : unmatchedProperties) {
        message << "  " << property << std::endl;
      }
      logERROR(message);
    }
  } // Try
  catch (PathfulException& e) {

    logERROR(e.path()+ " : " +e.what());
    stopTaskClock();
    return false;
  }

  stopTaskClock();
  return true;
}

//
// Build all passive components related to individual active sub-trackers. This procedure replaces the previously registered support (applying correct
// memory managment), if such an object existed.
//
bool Detector::buildPassiveTracker()
{
  if (m_activeTrackers.size()>0) {

    startTaskClock("Building passive materials of active trackers");

    for (auto& iTracker : m_activeTrackers) {

      // Build new inactive surface
      std:unique_ptr<insur::InactiveSurfaces> iPassive(new insur::InactiveSurfaces());

      // Get all materials in the way for given tracker
      material::Materialway materialway;
      materialway.build(*iTracker, *iPassive);

      m_passiveTrackers.push_back(std::move(iPassive));
    }

    stopTaskClock();
    return true;
  }
  else return false;
}

//
// Build beam pipe. This procedure replaces the previously registered beam pipe
//
bool Detector::buildBeamPipe(boost::property_tree::ptree& geomTree)
{
  // Clear out memory if beam pipe built before
  m_beamPipe.reset();

  // Build beam pipe -> look for tag "BeamPipe"
  startTaskClock("Building beam pipe");
  try {

    auto childRange = getChildRange(geomTree, "BeamPipe");
    for (auto it=childRange.first; it!=childRange.second; ++it) {

      // Declare if more than one beam pipes defined
      if (m_beamPipe.get()!=nullptr) {

        std::ostringstream message;
        message << "Detector::buildBeamPipe(): More than one beam pipe defined -> only first one found used!" << std::endl;
        logERROR(message);
      }
      else {

        std::cout << "Building new beam pipe" << std::endl;
        m_beamPipe = std::unique_ptr<BeamPipe>(GeometryFactory::make<BeamPipe>(it->second));
        m_beamPipe->build();
      }

    }; // Beam pipes

    // Declare unmatched properties
    std::set<string> unmatchedProperties = PropertyObject::reportUnmatchedProperties();
    if (!unmatchedProperties.empty()) {

      std::ostringstream message;
      message << "Detector::buildBeamPipe(): The following unknown properties were ignored:" << std::endl;

      for (const string& property  : unmatchedProperties) {
        message << "  " << property << std::endl;
      }
      logERROR(message);
    }
  } // Try
  catch (PathfulException& e) {

    logERROR(e.path()+ " : " +e.what());
    stopTaskClock();
    return false;
  }
  stopTaskClock();
  return true;
}

//
// GeometryVisitor pattern -> detector visitable
//
void Detector::accept(GeometryVisitor& v)
{
  v.visit(*this);
  for (auto& iTrk : m_activeTrackers) iTrk->accept(v);
  m_beamPipe->accept(v);
}

//
// GeometryVisitor pattern -> detector visitable (const. option)
//
void Detector::accept(ConstGeometryVisitor& v) const
{
  v.visit(*this);
  for (auto& iTrk : m_activeTrackers) iTrk->accept(v);
  m_beamPipe->accept(v);
}

//
// Get active sub-trackers as references
//
std::vector<const Tracker*> Detector::getActiveTrackers() const
{
  std::vector<const Tracker*> activeTrackers;

  for (auto& iTracker : m_activeTrackers) {

    const Tracker* trk = iTracker.get();
    activeTrackers.push_back(trk);
  }

  return activeTrackers;
}

//
// Get passive components related to active sub-trackers as references
//
std::vector<const insur::InactiveSurfaces*> Detector::getPassiveTrackers() const
{
  std::vector<const insur::InactiveSurfaces*> passiveTrackers;

  for (auto& iPassive : m_passiveTrackers) {

    const insur::InactiveSurfaces* passive = iPassive.get();
    passiveTrackers.push_back(passive);
  }

  return passiveTrackers;
}

//
// Get reference to beam pipe
//
const BeamPipe* Detector::getBeamPipe() const
{
  const BeamPipe* beamPipe = m_beamPipe.get();

  if (beamPipe==nullptr) throw std::invalid_argument( "Detector::getBeamPipe() - beam pipe not defined (null pointer), check!!!" );

  return beamPipe;
}



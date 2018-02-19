#ifndef CABLE_HH
#define CABLE_HH

#include <vector>
#include <string>

#include "Property.hh"
#include "Bundle.hh"


namespace insur { class DTC; }
using insur::DTC;


class Cable : public PropertyObject, public Buildable, public Identifiable<int> {
  typedef PtrVector<Bundle> Container;
public:
  Cable(const int id, const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);
  ~Cable();

  // BUNDLES CONNECTED TO THE CABLE.
  const Container& bundles() const { return bundles_; }

  void addBundle(Bundle* b) { bundles_.push_back(b); }
  int numBundles() const { return bundles_.size(); }

  // DTC THE CABLE IS CONNECTED TO.
  const DTC* getDTC() const {
    if (!myDTC_) throw PathfulException("myDTC_ is nullptr");
    return myDTC_; 
  }

  const Category& type() const { return type_; }
  const double phiSectorWidth() const { return phiSectorWidth_; }
  const int phiSectorRef() const { return phiSectorRef_; }
  const int slot() const { return slot_; }
  const bool isPositiveCablingSide() const { return isPositiveCablingSide_; }

  // SERVICES CHANNELS INFORMATION
  // Optical
  const ChannelSection* opticalChannelSection() const {
    if (!opticalChannelSection_) throw PathfulException("opticalChannelSection_ is nullptr");
    return opticalChannelSection_; 
  }
  // Power
  void assignPowerChannelSections();

private:
  void buildDTC(const double phiSectorWidth, const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide);
  const std::string computeDTCName(const int phiSectorRef, const Category& type, const int slot, const bool isPositiveCablingSide) const;
  
  Container bundles_;

  DTC* myDTC_ = nullptr;

  double phiSectorWidth_;
  int phiSectorRef_;
  Category type_;
  int slot_;
  bool isPositiveCablingSide_;

  ChannelSection* opticalChannelSection_ = nullptr;
};



#endif
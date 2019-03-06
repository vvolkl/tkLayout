#include "InnerCabling/inner_cabling_functions.hh"


namespace inner_cabling_functions {

  /* Compute the int n for which we have: phi ~= (n * phiSubUnitWidth + phiSubUnitStart).
   * n = 0 for the module with the lowest Phi > -Pi/2.
   * n = 0 for the module with the lowest Phi > Pi/2.
   */
  const int computePhiUnitRef(const double phi, const int numPhiUnits, const bool isPositiveZEnd) {

    const double projectedPhi = computePhiFromMinY(phi, isPositiveZEnd);

    const double phiUnitWidth = computePhiUnitWidth(numPhiUnits);
    const double phiUnitStart = computePhiUnitStart(projectedPhi, phiUnitWidth);

    const double phiRelative = femod(projectedPhi - phiUnitStart, M_PI);

    int phiUnitRef = 0;
    const double phiUnitRefExact = phiRelative / phiUnitWidth;

    // In case phiUnitRefExact is an integer, round it to an int!
    if (fabs((phiUnitRefExact - round(phiUnitRefExact))) < inner_cabling_roundingTolerance) phiUnitRef = fabs(round(phiUnitRefExact));
    else phiUnitRef = std::floor(phiUnitRefExact);

    return phiUnitRef;
  }


  /*
   * Starts counting Phi from -Pi/2.
   * This is useful as the module should be localized in Phi by IT half.
   */
  const double computePhiFromMinY(const double phi, const bool isPositiveZEnd) {
    const double stereoPhi = computeStereoPhi(phi, isPositiveZEnd);
    const double phiFromMinY = femod(stereoPhi + M_PI / 2., M_PI);
    return phiFromMinY;  
  }


  /*
   * If a module is on the negative cabling side, it is seen here as from a rotation of 180 deg around CMS_Y.
   */
  const double computeStereoPhi(const double phi, const bool isPositiveZEnd) {
    const double stereoPhi = (isPositiveZEnd ? phi : femod(M_PI - phi, 2.*M_PI) );
    return stereoPhi;
  }


  /*
   * Computes deltaPhi, depending on the number of phi units.
   */
  const double computePhiUnitWidth(const int numPhiUnits) {
    const double phiUnitWidth = (2.*M_PI) / numPhiUnits;
    return phiUnitWidth;
  }


  /* Compute the offset in Phi with respect to PhiUnitWidth.
   */
  const double computePhiUnitStart(const double phi, const double phiUnitWidth) {
    double phiUnitStart = femod(phi, phiUnitWidth);
    return phiUnitStart;
  }


  /*
   * Is the subdetector of barrel type?
   */
  const bool isBarrel(const std::string subDetectorName) {
    if (subDetectorName == inner_cabling_tbpx) return true;
    else if (subDetectorName == inner_cabling_tfpx || subDetectorName == inner_cabling_tepx) return false;
    else { 
      logERROR(any2str("Unknown subDetector name : ")
	       + any2str(subDetectorName)
	       );
      return false;
    }
  }


  /*
   * Split IT per (Z) end and (X) side: hence per quarter.
   */
  const int computeInnerTrackerQuarterIndex(const bool isPositiveZEnd, const bool isPositiveXSide) {
    int innerTrackerQuarterIndex = 0;
    if (isPositiveZEnd) {
      innerTrackerQuarterIndex  = (isPositiveXSide ? 1 : 2);
    }
    else {
      innerTrackerQuarterIndex  = (isPositiveXSide ? 3 : 4);
    }

    return innerTrackerQuarterIndex;
  }


  /*
   * Each subdetector is identified by a unique index.
   */
  const int computeSubDetectorIndex(const std::string subDetectorName) {
    int subDetectorIndex = 0;
    if (subDetectorName == inner_cabling_tbpx) subDetectorIndex = 1;
    else if (subDetectorName == inner_cabling_tfpx) subDetectorIndex = 2;
    else if (subDetectorName == inner_cabling_tepx) subDetectorIndex = 3;
    else {  logINFO(any2str("Unknown subDetector name : ")
		    + any2str(subDetectorName)
		    );
    }

    return subDetectorIndex;
  }


  /*
   * A given IT ring is divided by (X) side AND per dee side: hence notion of ring quarter.
   * Each ring quarter is identified by a unique index.
   */
  const int computeRingQuarterIndex(const int ringNumber, const bool isRingInnerEnd) {
    const int isRingInnerEndIndex = (!isRingInnerEnd);
    const int ringQuarterIndex = (ringNumber < 1 ? 0 : (ringNumber - 1) * 2 + isRingInnerEndIndex);
    return ringQuarterIndex;
  }


  /*
   * Retrieve the ring index from the ring quarter idex.
   * This is possible because a ring quarter is a part of a ring.
   */
  const int computeRingNumber(const int ringQuarterIndex) {
    const int ringNumber = 1 + ringQuarterIndex / 2;
    return ringNumber;
  }


  /*
   * Retrieve, from the index, whether one is on one dee side or the other.
   */
  const bool isRingInnerEnd(const int ringQuarterIndex) {
    const bool isRingInnerEnd = (ringQuarterIndex % 2 ? true : false);
    return isRingInnerEnd;
  }


  /*
   * Compute the number of ELinks per module, depending on the module location.
   */
  const int computeNumELinksPerModule(const std::string subDetectorName, const int layerOrRingNumber) {   
    int numELinksPerModule = 0;

    // tbpx
    if (subDetectorName == inner_cabling_tbpx) {
      if (layerOrRingNumber == 1) numELinksPerModule = inner_cabling_numELinksPerModuleBarrelLayer1;
      else if (layerOrRingNumber == 2) numELinksPerModule = inner_cabling_numELinksPerModuleBarrelLayer2;
      else if (layerOrRingNumber == 3) numELinksPerModule = inner_cabling_numELinksPerModuleBarrelLayer3;
      else if (layerOrRingNumber == 4) numELinksPerModule = inner_cabling_numELinksPerModuleBarrelLayer4;
      else { 
	logERROR(any2str("Found layer number ") + any2str(layerOrRingNumber)
		 + any2str(" in ") + any2str(inner_cabling_tbpx)
		 + any2str(". This is not supported.")
		 );
      }
    }
    // tfpx
    else if (subDetectorName == inner_cabling_tfpx) {
      if (layerOrRingNumber == 1) numELinksPerModule = inner_cabling_numELinksPerModuleForwardRing1;
      else if (layerOrRingNumber == 2) numELinksPerModule = inner_cabling_numELinksPerModuleForwardRing2;
      else if (layerOrRingNumber == 3) numELinksPerModule = inner_cabling_numELinksPerModuleForwardRing3;
      else if (layerOrRingNumber == 4) numELinksPerModule = inner_cabling_numELinksPerModuleForwardRing4;
      else { 
	logERROR(any2str("Found ring number ") + any2str(layerOrRingNumber)
		 + any2str(" in ") + any2str(inner_cabling_tfpx)
		 + any2str(". This is not supported.")
		 );
      }
    }
    // tepx
    else if (subDetectorName == inner_cabling_tepx) {
      numELinksPerModule = inner_cabling_numELinksPerModuleEndcap;
    }
    // other
    else { 
      logERROR(any2str("Unknown subDetector ") + any2str(subDetectorName)
	       );
    }

    return numELinksPerModule;
  }

} // namespace

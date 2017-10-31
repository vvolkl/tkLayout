/*
 * BeamPipe.h
 *
 *  Created on: 11. 5. 2016
 *      Author: Z.Drasal (CERN)
 */

#ifndef BEAMPIPE_H_
#define BEAMPIPE_H_

#include <memory>

#include "Property.h"
#include "Visitable.h"

// Forward declaration
class ConstGeometryVisitor;
class GeometryVisitor;
class InactiveTube;

/*
 * @class BeamPipe
 * @details Geometry & material object holding information about the beam pipe. Beam-pipe is defined by 2 contours: upper (@ higher R) & lower (@ lower R).
 * Both contours are defined by a sequence of [z,r] points, R-Phi symmetry is being assumed. The beam-pipe shape is then defined by 2 lines connecting these
 * points. In order to enclose the object and form thus a 3D object for material & tracking purposes, the upper line is assumed to end-up at the same [z,r]
 * point as the lower line. If the beam-pipe consists of several parts being constructed from different materials, e.g. Berrylium and Aluminium, the upper
 * contour needs to be enclosed at each border between these components. The final result is a 3D object consisting of several neighbouring 3D objects, each
 * of which once crossed by particle forms one or more material hits.
 *
 * Illustration:         --------              /
 *                       |      |             / tilted part
 * ----------------------|      |------------/
 *         Be           ||  Al  ||     Be     /
 * ------------------------------------------/
 *      Central BP        Flange          Conical BP
 *
 *
 * Finally, material is defined for each line segment, i.e. number of defined segments is (N-1), where N is number of [z,r] points. Tilt of each segment is
 * automatically calculated from [z,r] positions of individual upper/lower contour points, so no input from that point of view is needed.
 *
 */
class BeamPipe : public PropertyObject, public Identifiable<string>, public Buildable, Visitable {

 public:

  //! Constructor -> object set through build method
  BeamPipe(const PropertyTree& treeProperty);

  //! Destructor
  ~BeamPipe();

  //! Build method - setting all parameters
  void build();

  //! Setup: link lambda functions to various beampipe related properties (use setup functions for ReadOnly Computable properties -> use UncachedComputable if everytime needs to be recalculated)
  void setup();

  //! Get number of defined beam-pipe segments
  double getNSegments() const { return zPos.size()-1;};

  //! Get tilt of upper wall (contour) of i-th segment
  double getTiltUpper(unsigned int i) const;

  //! Get tilt of lower wall (contour) of i-th segment
  double getTiltLower(unsigned int i) const;

  //! Get average thickness of i-th segment
  double getAvgThickness(unsigned int i) const;

  //! GeometryVisitor pattern -> beam pipe visitable
  void accept(GeometryVisitor& v);

  //! GeometryVisitor pattern -> beam pipe visitable (const. option)
  void accept(ConstGeometryVisitor& v) const;

  // Beam pipe radius, thickness, thickness in rad. length, in int. length
  ReadonlyPropertyVector<double>      radiusLower;      //!< Beam pipe radii: lower wall/contour (beam-pipe shape described by sequence of segments defined by: lower & upper radii & zPos )
  ReadonlyPropertyVector<double>      radiusUpper;      //!< Beam pipe radii: upper wall/contour (beam-pipe shape described by sequence of segments defined by: lower & upper radii & zPos )
  ReadonlyPropertyVector<double>      zPos;             //!< Beam pipe zPos  (beam-pipe shape described by sequence of segments defined by: low & high radii & zPos )
  ReadonlyProperty<double, Computable>minRadius;        //!< Beam pipe minimal radius -> used e.g. for drawing purposes
  ReadonlyProperty<double, Computable>maxRadius;        //!< Beam pipe maximum radius -> used e.g. for drawing purposes
  ReadonlyPropertyVector<double>      radLength;        //!< Beam pipe material rad. length for each beam-pipe segment -> material calculated as: pathLength in given segment/radLength
  ReadonlyPropertyVector<double>      intLength;        //!< Beam pipe material int. lenght for each beam-pipe segment -> material calculated as: pathLength in given segment/radLength

 private:

  //! Cross-check parameters provided from geometry configuration file
  void check() override;

}; // Class

#endif /* BEAMPIPE_H_ */

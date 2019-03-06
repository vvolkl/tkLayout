/**
 * @file ConversionStation.h
 *
 * @date Jul 21, 2014
 * @author Stefano Martina
 */

#ifndef CONVERSIONSTATION_H_
#define CONVERSIONSTATION_H_

#include "Property.hh"
#include "MaterialObject.hh"

namespace insur {
  class InactiveElement;
}

using insur::InactiveElement;

namespace material {
  //class MaterialObject;
  
  class ConversionStation :public MaterialObject {
  public:
    enum Type {ERROR, FLANGE, SECOND};
    
    ConversionStation() :
      MaterialObject(MaterialObject::Type::STATION),
      stationName_ ("stationName", parsedAndChecked()),
      type_ ("type", parsedAndChecked()),
      minZ_ ("minZ", parsedOnly()),
      maxZ_ ("maxZ", parsedOnly()),
      stationType_ (ERROR),
      conversionsNode_ ("Conversion", parsedOnly())
    {};
    virtual ~ConversionStation() {};

    void build();
    void routeConvertedElements(MaterialObject& localOutput, MaterialObject& serviceOutput, InactiveElement& inactiveElement);
    //void routeConvertedServicesTo(MaterialObject& outputObject) const;
    //void routeConvertedLocalsTo(MaterialObject& outputObject) const;
    Type stationType() const;

    ReadonlyProperty<std::string, NoDefault> stationName_;
    ReadonlyProperty<std::string, NoDefault> type_;
    ReadonlyProperty<double, NoDefault> minZ_;
    ReadonlyProperty<double, NoDefault> maxZ_;
    const double meanZ() const { return (minZ_() + maxZ_()) / 2.; }

  private:
    static const std::map<std::string, Type> typeString;
    Type stationType_;
    bool valid_;
    PropertyNodeUnique<std::string> conversionsNode_;

    void buildConversions();

    /*
    class Element : public PropertyObject {
    public:
      ReadonlyProperty<std::string, NoDefault> elementName;
      ReadonlyProperty<long, NoDefault> quantity;
      ReadonlyProperty<std::string, NoDefault> unit;
      ReadonlyProperty<bool, Default> service;

      Element() :
        elementName ("elementName", parsedAndChecked()),
        quantity ("quantity", parsedAndChecked()),
        unit ("unit", parsedAndChecked()),
        service ("service", parsedOnly(), false) {};
      virtual ~Element() {};

      //void build();
    };
    */

    class Inoutput : public PropertyObject {
    public:
      PropertyNodeUnique<std::string> elementsNode_;

      Inoutput() :
        elementsNode_ ("Element", parsedOnly()),
        elementMaterialType(MaterialObject::Type::STATION) {};

        virtual ~Inoutput() {};

      void build();

      std::vector<MaterialObject::Element*> elements;
      MaterialObject::Type elementMaterialType;
    };

    class Conversion : public PropertyObject {
    public:
      PropertyNode<std::string> inputNode_;
      PropertyNode<std::string> outputNode_;

      Conversion() :
        inputNode_ ("Input", parsedAndChecked()),
        outputNode_ ("Output", parsedAndChecked()) {};
      virtual ~Conversion() {};

      void build();

      Inoutput* input;
      Inoutput* outputs;
    };

    std::vector<Conversion *> conversions;
  };

} /* namespace material */

#endif /* CONVERSIONSTATION_H_ */

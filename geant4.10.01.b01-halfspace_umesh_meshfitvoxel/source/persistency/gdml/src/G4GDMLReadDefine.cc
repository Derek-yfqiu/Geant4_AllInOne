//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4GDMLReadDefine.cc 68053 2013-03-13 14:39:51Z gcosmo $
//
// class G4GDMLReadDefine Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLReadDefine.hh"
#include "G4PolyhedronArbitrary.hh" //qiu
#include "G4VisAttributes.hh"

G4GDMLMatrix::G4GDMLMatrix()
  : m(0), rows(0), cols(0)
{
}

G4GDMLMatrix::G4GDMLMatrix(size_t rows0, size_t cols0)
{   
   if ((rows0==0) || (cols0==0))
   {
     G4Exception("G4GDMLMatrix::G4GDMLMatrix(r,c)", "InvalidSetup",
                 FatalException, "Zero indeces as arguments!?");
   }
   rows = rows0;
   cols = cols0;
   m = new G4double[rows*cols];
}

G4GDMLMatrix::G4GDMLMatrix(const G4GDMLMatrix& rhs)
  : m(0), rows(0), cols(0)
{
   if (rhs.m)
   {
     rows = rhs.rows;
     cols = rhs.cols;
     m = new G4double[rows*cols];
     for (size_t i=0; i<rows*cols; i++)  { m[i] = rhs.m[i]; }
   }
}

G4GDMLMatrix& G4GDMLMatrix::operator=(const G4GDMLMatrix& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy data
   //
   rows = rhs.rows;
   cols = rhs.cols;
   if (rhs.m)
   {
     m = new G4double[rows*cols];
     for (size_t i=0; i<rows*cols; i++)  { m[i] = rhs.m[i]; }
   }
   else
   {
     m = 0;
   }

   return *this;
}

G4GDMLMatrix::~G4GDMLMatrix()
{
   delete [] m;
}

void G4GDMLMatrix::Set(size_t r,size_t c,G4double a)
{   
   if (r>=rows || c>=cols)
   {
     G4Exception("G4GDMLMatrix::set()", "InvalidSetup",
                 FatalException, "Index out of range!");
   }
   m[cols*r+c] = a;
}

G4double G4GDMLMatrix::Get(size_t r,size_t c) const
{   
   if (r>=rows || c>=cols)
   {
     G4Exception("G4GDMLMatrix::get()", "InvalidSetup",
                 FatalException, "Index out of range!");
   }
   return m[cols*r+c];
}

size_t G4GDMLMatrix::GetRows() const
{
   return rows;
}

size_t G4GDMLMatrix::GetCols() const
{
   return cols;
}

G4GDMLReadDefine::G4GDMLReadDefine() : G4GDMLRead()
{
}

G4GDMLReadDefine::~G4GDMLReadDefine()
{
}

G4RotationMatrix
G4GDMLReadDefine::GetRotationMatrix(const G4ThreeVector& angles)
{
   G4RotationMatrix rot;

   rot.rotateX(angles.x());
   rot.rotateY(angles.y());
   rot.rotateZ(angles.z());

   return rot;
}

void
G4GDMLReadDefine::ConstantRead(const xercesc::DOMElement* const constantElement)
{
   G4String name  = "";
   G4double value = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = constantElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::ConstantRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name")  { name = attValue; }  else
      if (attName=="value") { value = eval.Evaluate(attValue); }
   }

   eval.DefineConstant(name,value);
}

void
G4GDMLReadDefine::ExpressionRead(const xercesc::DOMElement* const expElement)
{
   G4String name  = "";
   G4double value = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = expElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::ExpressionRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name")  { name = attValue; }
   }

   const G4String expValue = Transcode(expElement->getTextContent());
   value = eval.Evaluate(expValue);
   eval.DefineConstant(name,value);
}

void
G4GDMLReadDefine::MatrixRead(const xercesc::DOMElement* const matrixElement) 
{
   G4String name = "";
   G4int coldim  = 0;
   G4String values = "";

   const xercesc::DOMNamedNodeMap* const attributes
         = matrixElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::MatrixRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name")   { name  = GenerateName(attValue); } else
      if (attName=="coldim") { coldim = eval.EvaluateInteger(attValue); } else
      if (attName=="values") { values = attValue; }
   }

   std::stringstream MatrixValueStream(values);
   std::vector<G4double> valueList;

   while (!MatrixValueStream.eof())
   {
      G4String MatrixValue;
      MatrixValueStream >> MatrixValue;
      valueList.push_back(eval.Evaluate(MatrixValue));
   }

   eval.DefineMatrix(name,coldim,valueList);

   G4GDMLMatrix matrix(valueList.size()/coldim,coldim);

   for (size_t i=0;i<valueList.size();i++)
   {
      matrix.Set(i/coldim,i%coldim,valueList[i]);
   }

   matrixMap[name] = matrix;
}

void
G4GDMLReadDefine::PositionRead(const xercesc::DOMElement* const positionElement)
{
   G4String name = "";
   G4double unit = 1.0;
   G4ThreeVector position(0.,0.,0.);

   const xercesc::DOMNamedNodeMap* const attributes
         = positionElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::PositionRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") { name = GenerateName(attValue); }  else
      if (attName=="unit") { unit = eval.Evaluate(attValue); } else
      if (attName=="x") { position.setX(eval.Evaluate(attValue)); } else
      if (attName=="y") { position.setY(eval.Evaluate(attValue)); } else
      if (attName=="z") { position.setZ(eval.Evaluate(attValue)); }
   }

   positionMap[name] = position*unit;
}

void
G4GDMLReadDefine::RotationRead(const xercesc::DOMElement* const rotationElement)
{
   G4String name = "";
   G4double unit = 1.0;
   G4ThreeVector rotation(0.,0.,0.);

   const xercesc::DOMNamedNodeMap* const attributes
         = rotationElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::RotationRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") { name = GenerateName(attValue); }  else
      if (attName=="unit") { unit = eval.Evaluate(attValue); } else
      if (attName=="x") { rotation.setX(eval.Evaluate(attValue)); } else
      if (attName=="y") { rotation.setY(eval.Evaluate(attValue)); } else
      if (attName=="z") { rotation.setZ(eval.Evaluate(attValue)); }
   }

   rotationMap[name] = rotation*unit;
}

void G4GDMLReadDefine::ScaleRead(const xercesc::DOMElement* const scaleElement)
{
   G4String name = "";
   G4ThreeVector scale(1.0,1.0,1.0);

   const xercesc::DOMNamedNodeMap* const attributes
         = scaleElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::ScaleRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") { name = GenerateName(attValue); }    else
      if (attName=="x") { scale.setX(eval.Evaluate(attValue)); } else
      if (attName=="y") { scale.setY(eval.Evaluate(attValue)); } else
      if (attName=="z") { scale.setZ(eval.Evaluate(attValue)); }
   }

   scaleMap[name] = scale;
}

void
G4GDMLReadDefine::VariableRead(const xercesc::DOMElement* const variableElement)
{
   G4String name  = "";
   G4double value = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = variableElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::VariableRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name")  { name = attValue; } else
      if (attName=="value") { value = eval.Evaluate(attValue); }
   }

   eval.DefineVariable(name,value);
}

void G4GDMLReadDefine::QuantityRead(const xercesc::DOMElement* const element)
{
   G4String name = "";
   G4double unit = 1.0;
   G4double value = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::QuantityRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") { name = attValue; } else
      if (attName=="value") { value = eval.Evaluate(attValue); } else
      if (attName=="unit") { unit = eval.Evaluate(attValue); }
   }

   quantityMap[name] = value*unit;
   eval.DefineConstant(name,value*unit);
}

void
G4GDMLReadDefine::DefineRead(const xercesc::DOMElement* const defineElement)
{
   G4cout << "G4GDML: Reading definitions..." << G4endl;

   for (xercesc::DOMNode* iter = defineElement->getFirstChild();
        iter != 0;iter = iter->getNextSibling())
   {
      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

      const xercesc::DOMElement* const child
            = dynamic_cast<xercesc::DOMElement*>(iter);
      if (!child)
      {
        G4Exception("G4GDMLRead::DefineRead()", "InvalidRead",
                    FatalException, "No child found!");
        return;
      }
      const G4String tag = Transcode(child->getTagName());

      if (tag=="constant") { ConstantRead(child); } else
      if (tag=="matrix")   { MatrixRead(child); }   else
      if (tag=="position") { PositionRead(child); } else
      if (tag=="rotation") { RotationRead(child); } else
      if (tag=="scale")    { ScaleRead(child); }    else
      if (tag=="variable") { VariableRead(child); } else
      if (tag=="quantity") { QuantityRead(child); } else
      if (tag=="expression") { ExpressionRead(child); }
      //qiu
      else if (tag=="Polyhedron") { PolyhedronRead(child); }
      else
      {
        G4String error_msg = "Unknown tag in define: "+tag;
        G4Exception("G4GDMLReadDefine::defineRead()", "ReadError",
                    FatalException, error_msg);
      }
   }
}

void
G4GDMLReadDefine::VectorRead(const xercesc::DOMElement* const vectorElement,
                             G4ThreeVector& vec)
{
   G4double unit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes
         = vectorElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::VectorRead()", "InvalidRead",
                    FatalException, "No attribute found!");
        return;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="unit") { unit = eval.Evaluate(attValue); } else
      if (attName=="x") { vec.setX(eval.Evaluate(attValue)); } else
      if (attName=="y") { vec.setY(eval.Evaluate(attValue)); } else
      if (attName=="z") { vec.setZ(eval.Evaluate(attValue)); }
   }

   vec *= unit;
}

G4String G4GDMLReadDefine::RefRead(const xercesc::DOMElement* const element)
{
   G4String ref;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;
        attribute_index<attributeCount; attribute_index++)
   {
      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
      { continue; }

      const xercesc::DOMAttr* const attribute
            = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      if (!attribute)
      {
        G4Exception("G4GDMLRead::Read()", "InvalidRead",
                    FatalException, "No attribute found!");
        return ref;
      }
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="ref") { ref = attValue; }
   }

   return ref;
}

G4bool G4GDMLReadDefine::IsValidID(const G4String& ref) const
{
   return eval.IsVariable(ref);
}

G4double G4GDMLReadDefine::GetConstant(const G4String& ref)
{
   return eval.GetConstant(ref);
}

G4double G4GDMLReadDefine::GetVariable(const G4String& ref)
{
   return eval.GetVariable(ref);
}

G4double G4GDMLReadDefine::GetQuantity(const G4String& ref)
{
   if (quantityMap.find(ref) == quantityMap.end())
   {
     G4String error_msg = "Quantity '"+ref+"' was not found!";
     G4Exception("G4GDMLReadDefine::getQuantity()", "ReadError",
                 FatalException, error_msg);
   }
   return quantityMap[ref];
}

G4ThreeVector G4GDMLReadDefine::GetPosition(const G4String& ref)
{
   if (positionMap.find(ref) == positionMap.end())
   {
     G4String error_msg = "Position '"+ref+"' was not found!";
     G4Exception("G4GDMLReadDefine::getPosition()", "ReadError",
                 FatalException, error_msg);
   }
   return positionMap[ref];
}

G4ThreeVector G4GDMLReadDefine::GetRotation(const G4String& ref)
{
   if (rotationMap.find(ref) == rotationMap.end())
   {
     G4String error_msg = "Rotation '"+ref+"' was not found!";
     G4Exception("G4GDMLReadDefine::getRotation()", "ReadError",
                 FatalException, error_msg);
   }
   return rotationMap[ref];
}

G4ThreeVector G4GDMLReadDefine::GetScale(const G4String& ref)
{
   if (scaleMap.find(ref) == scaleMap.end())
   {
     G4String error_msg = "Scale '"+ref+"' was not found!";
     G4Exception("G4GDMLReadDefine::getScale()", "ReadError",
                 FatalException, error_msg);
   }
   return scaleMap[ref];
}

G4GDMLMatrix G4GDMLReadDefine::GetMatrix(const G4String& ref)
{
   if (matrixMap.find(ref) == matrixMap.end())
   {
     G4String error_msg = "Matrix '"+ref+"' was not found!";
     G4Exception("G4GDMLReadDefine::getMatrix()", "ReadError",
                 FatalException, error_msg);
   }
   return matrixMap[ref];
}

//qiu start here to add some codes
//method to read ONE node
G4ThreeVector  G4GDMLReadDefine::NodeRead(const xercesc::DOMElement* const NodeElement)
{
    G4double x=0., y=0., z=0.;
    const xercesc::DOMNamedNodeMap* const attributes
          = NodeElement->getAttributes();
    XMLSize_t attributeCount = attributes->getLength();
    for (XMLSize_t attribute_index=0;
         attribute_index<attributeCount; attribute_index++)
    {
       xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

       if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
         { continue; }

       const xercesc::DOMAttr* const attribute
             = dynamic_cast<xercesc::DOMAttr*>(attribute_node);
       if (!attribute)
       {
         G4Exception("G4GDMLReadDefine::NodeRead()",
                     "InvalidRead", FatalException, "No attribute found!");
         return G4ThreeVector() ;
       }
       const G4String attName = Transcode(attribute->getName());
       const G4String attValue = Transcode(attribute->getValue());

       if (attName=="x") { x = eval.Evaluate(attValue); } else
       if (attName=="y") { y = eval.Evaluate(attValue); } else
       if (attName=="z") { z = eval.Evaluate(attValue); }
    }
    return G4ThreeVector(x,y,z);
}



//method to read the nodes, defined by three
std::vector < G4ThreeVector >  G4GDMLReadDefine::NodesRead(const xercesc::DOMElement* const NodesElement)
{
    std::vector <G4ThreeVector > aNodeList;

    //first read the attribute, only one attribute possible: lunit
    G4double lunit = 1.0;
    const xercesc::DOMNamedNodeMap* const attributes
          = NodesElement->getAttributes();
    XMLSize_t attributeCount = attributes->getLength();
    for (XMLSize_t attribute_index=0;
         attribute_index<attributeCount; attribute_index++)
    {
       xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

       if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
         { continue; }

       const xercesc::DOMAttr* const attribute
             = dynamic_cast<xercesc::DOMAttr*>(attribute_node);
       if (!attribute)
       {
           //use the default unit: mm
           break;
       }
       const G4String attName = Transcode(attribute->getName());
       const G4String attValue = Transcode(attribute->getValue());

       if (attName=="lunit")
           lunit = eval.Evaluate(attValue);
    }

    //then read the nodes


    for (xercesc::DOMNode* iter = NodesElement->getFirstChild();
         iter != 0;iter = iter->getNextSibling())
    {
       if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

       const xercesc::DOMElement* const child
             = dynamic_cast<xercesc::DOMElement*>(iter);
       if (!child)
       {
         G4Exception("G4GDMLRead::NodesRead()", "InvalidRead",
                     FatalException, "No child found!");
         return aNodeList;
       }
       const G4String tag = Transcode(child->getTagName());

       if (tag=="Node") {
           G4ThreeVector aPoint = NodeRead(child);
           aPoint *= lunit; //use the unit
           aNodeList.push_back(aPoint);
       }
       else
       {
         G4String error_msg = "Unknown tag in define: "+tag;
         G4Exception("G4GDMLReadDefine::NodesRead()", "ReadError",
                     FatalException, error_msg);
       }
    }
    return aNodeList;
}

std::vector < G4int> G4GDMLReadDefine::FacetRead(const xercesc::DOMElement* const FacetElement)
{
    G4int n1=0, n2=0, n3=0, n4=0;
    std::vector <G4int> aFacet;


    const xercesc::DOMNamedNodeMap* const attributes
          = FacetElement->getAttributes();
    XMLSize_t attributeCount = attributes->getLength();
    for (XMLSize_t attribute_index=0;
         attribute_index<attributeCount; attribute_index++)
    {
       xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

       if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
         { continue; }

       const xercesc::DOMAttr* const attribute
             = dynamic_cast<xercesc::DOMAttr*>(attribute_node);
       if (!attribute)
       {
         G4Exception("G4GDMLReadDefine::FacetRead()",
                     "InvalidRead", FatalException, "No attribute found!");
         return aFacet;
       }
       const G4String attName = Transcode(attribute->getName());
       const G4String attValue = Transcode(attribute->getValue());
       G4int aTmpInt = eval.EvaluateInteger(attValue);

       if (attName=="n1") { n1 = aTmpInt; } else
       if (attName=="n2") { n2 = aTmpInt; } else
       if (attName=="n3") { n3 = aTmpInt; } else
       if (attName=="n4") { n4 = aTmpInt; }
    }

    //for return the result
    aFacet.push_back(n1);
    aFacet.push_back(n2);
    aFacet.push_back(n3);
    aFacet.push_back(n4);
    return aFacet;
}


//method to read the facets, defined by three or four integer which refer the node id
std::vector < std::vector <G4int> > G4GDMLReadDefine::FacetsRead(const xercesc::DOMElement* const FacetsElement)
{

    std::vector < std::vector <G4int>  > aFacetList;

    for (xercesc::DOMNode* iter = FacetsElement->getFirstChild();
         iter != 0;iter = iter->getNextSibling())
    {
       if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

       const xercesc::DOMElement* const child
             = dynamic_cast<xercesc::DOMElement*>(iter);
       if (!child)
       {
         G4Exception("G4GDMLRead::FacetsRead()", "InvalidRead",
                     FatalException, "No child found!");
         return aFacetList;
       }
       const G4String tag = Transcode(child->getTagName());

       if (tag=="Facet") {
           std::vector < G4int> aFacet = FacetRead(child);
           aFacetList.push_back(aFacet);
       }
       else
       {
         G4String error_msg = "Unknown tag in define: "+tag;
         G4Exception("G4GDMLReadDefine::FacetsRead()", "ReadError",
                     FatalException, error_msg);
       }
    }
    return aFacetList;
}

//for reading the color of this polyhedron
G4Color G4GDMLReadDefine::ColorRead (const xercesc::DOMElement* const ColorElement)
{
    G4int r = 175, g=175, b=175;
    const xercesc::DOMNamedNodeMap* const attributes
          = ColorElement->getAttributes();
    XMLSize_t attributeCount = attributes->getLength();
    for (XMLSize_t attribute_index=0;
         attribute_index<attributeCount; attribute_index++)
    {
       xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

       if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
         { continue; }

       const xercesc::DOMAttr* const attribute
             = dynamic_cast<xercesc::DOMAttr*>(attribute_node);
       if (!attribute)
       {
         G4Exception("G4GDMLReadDefine::ColorRead()",
                     "InvalidRead", FatalException, "No attribute found!");
         return G4Color() ;
       }
       const G4String attName = Transcode(attribute->getName());
       const G4String attValue = Transcode(attribute->getValue());

       if (attName=="r") { r = eval.Evaluate(attValue); } else
       if (attName=="g") { g = eval.Evaluate(attValue); } else
       if (attName=="b") { b = eval.Evaluate(attValue); }
    }
    return G4Color(r,g,b);
}

//for read the polyhedron,used in visualization of the HalfSpaceSolid
void G4GDMLReadDefine::PolyhedronRead (const xercesc::DOMElement* const PolyhedronElement)
{
    //first read the attribute, only one attribute possible: name
    G4String aName;
    const xercesc::DOMNamedNodeMap* const attributes
          = PolyhedronElement->getAttributes();
    XMLSize_t attributeCount = attributes->getLength();
    for (XMLSize_t attribute_index=0;
         attribute_index<attributeCount; attribute_index++)
    {
       xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

       if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE)
         { continue; }

       const xercesc::DOMAttr* const attribute
             = dynamic_cast<xercesc::DOMAttr*>(attribute_node);
       if (!attribute)
       {
         G4Exception("G4GDMLReadDefine::PolyhedronRead()",
                     "InvalidRead", FatalException, "No attribute found!");
         return;
       }
       const G4String attName = Transcode(attribute->getName());
       const G4String attValue = Transcode(attribute->getValue());

       if (attName=="name")
           aName = attValue;
    }

    if (aName.empty())        {
        G4Exception("G4GDMLReadDefine::PolyhedronRead()",
                    "InvalidRead", FatalException, "name should not be empty!");
        return;
    }

    //then read the nodes and facets
    std::vector <G4ThreeVector > aNodeList;
    std::vector < std::vector <G4int>  > aFacetList;
    G4Color aColor;

    for (xercesc::DOMNode* iter = PolyhedronElement->getFirstChild();
         iter != 0;iter = iter->getNextSibling())
    {
       if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) { continue; }

       const xercesc::DOMElement* const child
             = dynamic_cast<xercesc::DOMElement*>(iter);
       if (!child)
       {
         G4Exception("G4GDMLRead::PolyhedronRead()", "InvalidRead",
                     FatalException, "No child found!");
         return;
       }
       const G4String tag = Transcode(child->getTagName());

       if (tag=="Nodes") { aNodeList = NodesRead(child);} else
       if (tag=="Facets") { aFacetList = FacetsRead(child);} else
       if (tag=="Color"){ aColor = ColorRead(child);}
       else
       {
         G4String error_msg = "Unknown tag in define: "+tag;
         G4Exception("G4GDMLReadDefine::PolyhedronRead()", "ReadError",
                     FatalException, error_msg);
       }
    }

    //At the end, we make the polyhedron
    if (aNodeList.empty() || aFacetList.empty())
        G4Exception("G4GDMLRead::PolyhedronRead()", "InvalidRead",
                    FatalException, "Node or facet list empty!");
/*
    G4Polyhedron *ph=new G4Polyhedron;
    typedef G4int G4int4[4];
    typedef G4double G4double3[3];

    G4double3 * aPointArray = new G4double3 [aNodeList.size()];
    for (unsigned int i=0; i<aNodeList.size(); i++) {
        aPointArray[i][0] = aNodeList[i].x();
        aPointArray[i][1] = aNodeList[i].y();
        aPointArray[i][2] = aNodeList[i].z();
    }
    G4int4 * aFacetArray = new G4int4 [aFacetList.size()];
    for (unsigned int i=0; i<aFacetList.size(); i++) {
        aFacetArray[i][0] = aFacetList[i][0];
        aFacetArray[i][1] = aFacetList[i][1];
        aFacetArray[i][2] = aFacetList[i][2];
        aFacetArray[i][3] = aFacetList[i][3];
    }
    ph->createPolyhedron(aNodeList.size(), aFacetList.size(),
                         aPointArray, aFacetArray);
    delete [] aPointArray;
    delete [] aFacetArray;
*/


    G4PolyhedronArbitrary *ph =
      new G4PolyhedronArbitrary (aNodeList.size(), aFacetList.size());
    for (unsigned int i=0; i<aNodeList.size(); i++) {
        ph->AddVertex(aNodeList[i]);
    }
    for (unsigned int i=0; i<aFacetList.size(); i++) {
        ph->AddFacet(aFacetList[i][0],aFacetList[i][1], aFacetList[i][2], aFacetList[i][3] );
    }

    G4VisAttributes aVis;
    aVis.SetColor(aColor);
    ph->SetVisAttributes(aVis);

    polyhedronMap[aName] = ph;

}

G4Polyhedron * G4GDMLReadDefine::GetPolyhedron(const G4String& ref)
{
   if (polyhedronMap.find(ref) == polyhedronMap.end())
   {
     G4String error_msg = "polyhedron '"+ref+"' was not found!";
     G4Exception("G4GDMLReadDefine::getpolyhedron()", "ReadError",
                 FatalException, error_msg);
   }
   return polyhedronMap[ref];
}

#include "G4ScoringUMesh.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4UMeshParser.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4VScoreColorMap.hh"


G4ScoringUMesh::G4ScoringUMesh(G4String wName)
    :G4VScoringMesh(wName)
{
    fShape = unstrMesh;
    fUMeshParser = 0;

}

G4ScoringUMesh::~G4ScoringUMesh()
{
    if (fUMeshParser) delete fUMeshParser;
}

/*!
 * \brief G4ScoringUMesh::Construct
 * \param fWorldPhys
 */
void  G4ScoringUMesh::Construct(G4VPhysicalVolume* fWorldPhys)
{
    if(fConstructed) {

        if(verboseLevel > 0)
            G4cout << fWorldPhys->GetName() << " --- All quantities are reset."
                   << G4endl;
        ResetScore();

    } else {
        fConstructed = true;

        if(verboseLevel > 9) G4cout << " G4ScoringUMesh::Construct mesh..." << G4endl;

        // World
        G4VPhysicalVolume * scoringWorld = fWorldPhys;
        G4LogicalVolume * worldLogical = scoringWorld->GetLogicalVolume();

        //ask the mesh parser to read the mesh
        fUMeshParser->setLogicalEnvelop(worldLogical);
        fUMeshParser->setMaterial(0);
        fUMeshParser->setDimFactor(fDimensionFactor); //in case factor set after reading mesh
        fUMeshParser->setTranslation(fCenterPosition);
        fUMeshParser->setRotationMatrix(fRotationMatrix);
        IsUMeshRead = fUMeshParser->buildMesh();  //generate the mesh by creating element solids
        if (!IsUMeshRead){ //if read failed
            G4cerr << "ERROR : G4VScoringMesh::Construct() : "
                   << "Build mesh failed!" << G4endl;
            return;
        }
        //to assign a MultifunctionalDetector for all the mesh elements
        // and also vis. attributes
        G4VisAttributes * visatt = new G4VisAttributes(G4Colour(.5,.5,.5));
        visatt->SetVisibility(true);
        fElementPVList = fUMeshParser->getElementPVList();
        for (unsigned int i=0; i<fElementPVList.size(); i++) {
            fElementPVList[i]->GetLogicalVolume()->SetSensitiveDetector(fMFD);
            fElementPVList[i]->GetLogicalVolume()->SetVisAttributes(visatt);
        }
    }
}

void  G4ScoringUMesh::List() const
{
    G4cout << "G4ScoringUMesh : " << fWorldName << " --- Shape: unstructured mesh" << G4endl;
}

/*!
 * \brief G4ScoringUMesh::Draw
 *  draw the result of the mesh
 *  here is an important issue about how to associate the key of the G4THitsMap with
 *  mesh elements. we use a trick, assigning a sequential copyNo to elements, for example,
 *  0,1,2,..., then use the default GetIndex() method in G4PrimitiveScores to generate the
 *  key for G4THitsMap, in this way the key of the G4THitsMap is matching with Element id
 * \param map G4THitsMap
 * \param colorMap color legend
 */
void G4ScoringUMesh:: Draw(std::map<G4int, G4double*> * map, G4VScoreColorMap* colorMap, G4int )
{
    G4VVisManager * pVisManager = G4VVisManager::GetConcreteInstance();

    std::vector <G4double> results;
    results.assign(fElementPVList.size(), 0.0); //allocate and initiate
    if(pVisManager)
    {
        //get the result from the map
        std::map<G4int, G4double*>::iterator itr = map->begin();
        for(; itr != map->end(); itr++) {

//            G4cout << "Map key:"<<itr->first<<"\t\t value:"<<*(itr->second)/fDrawUnitValue <<G4endl; //for test
            results.at(itr->first) = *(itr->second)/fDrawUnitValue;
        }
        //get the min and max value
        G4double theMin = DBL_MAX, theMax = DBL_MIN;
        for (unsigned int i=0; i<results.size(); i++) {
            if (results[i] < theMin) theMin = results[i];
            if (results[i] > theMax) theMax = results[i];
        }
        //set the visualization attibute
        G4VisAttributes att;
        att.SetForceSolid(true);
        att.SetForceAuxEdgeVisible(true);
        //set the min and max of the color map
        if(colorMap->IfFloatMinMax()) { colorMap->SetMinMax(theMin ,theMax); }
        //for each element, get the solid
        for (unsigned int i=0; i<fElementPVList.size(); i++) {
            G4VSolid * aSolid = fElementPVList[i]->GetLogicalVolume()->GetSolid();
            G4Transform3D trans;
            if(fRotationMatrix) {
              trans = G4Rotate3D(*fRotationMatrix).inverse(); //?? invert?
              trans = G4Translate3D(fCenterPosition)*trans;
            } else {
              trans = G4Translate3D(fCenterPosition);
            }
            //set the color
            G4double c[4];
            colorMap->GetMapColor(results[i], c);
            att.SetColour(c[0], c[1], c[2]);//, c[3]);
            pVisManager->Draw(*aSolid, att, trans);
        }
    }
    colorMap->SetPSUnit(fDrawUnit);
    colorMap->SetPSName(fDrawPSName);
    colorMap->DrawColorChart();
}

void  G4ScoringUMesh::DrawColumn(std::map<G4int, G4double*> * , G4VScoreColorMap* ,
                         G4int , G4int )
{
    G4cerr <<" Not possible for unstructured mesh!"<<G4endl;
}

/*!
 * \brief G4ScoringUMesh::openUMeshFile
 * open the unstructured mesh file,
 */
void G4ScoringUMesh::openUMeshFile()
{
    if (fMeshFileName.isNull())        {
        G4cerr << "ERROR: file name should be set first!" << G4endl;
        return;
    }
    if (fMeshFormat == G4UM::Mesh_UnKnown) {
        G4cerr << "ERROR: mesh file format should be set !" << G4endl;
        return;
    }

    //create a new parser
    if (fUMeshParser == 0)
        fUMeshParser = new G4UMeshParser;
    fUMeshParser->setVerboseLevel(verboseLevel);
    fUMeshParser->setDimFactor(fDimensionFactor);  //in case set factor before reading mesh
    IsUMeshRead = fUMeshParser->readMesh(fMeshFileName, fMeshFormat);
    if (!IsUMeshRead){ //if read failed
        G4cerr << "ERROR : G4VScoringMesh::SetPrimitiveScorer() : "
               << "Read mesh failed!" << G4endl;
        return;
    }
}



















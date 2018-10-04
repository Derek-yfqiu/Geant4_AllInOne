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
//

#include "G4WorkerRunManager.hh"
#include "G4WorkerRunManagerKernel.hh"
#include "G4UImanager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4MTRunManager.hh"
#include "G4ScoringManager.hh"
#include "G4TransportationManager.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4WorkerThread.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserWorkerInitialization.hh"
#include "G4UserWorkerThreadInitialization.hh"
#include "G4UserRunAction.hh"
#include "G4RNGHelper.hh"
#include "G4Run.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4VScoringMesh.hh"
#include <sstream>

G4WorkerRunManager* G4WorkerRunManager::GetWorkerRunManager()
{ return static_cast<G4WorkerRunManager*>(G4RunManager::GetRunManager()); }

G4WorkerRunManagerKernel* G4WorkerRunManager::GetWorkerRunManagerKernel()
{ return static_cast<G4WorkerRunManagerKernel*>(GetWorkerRunManager()->kernel); }

G4WorkerRunManager::G4WorkerRunManager() : G4RunManager(workerRM) {
    //This constructor should never be called in non-multithreaded mode
#ifndef G4MULTITHREADED
    G4ExceptionDescription msg;
    msg<<"Geant4 code is compiled without multi-threading support (-DG4MULTITHREADED is set to off).";
    msg<<" This type of RunManager can only be used in mult-threaded applications.";
    G4Exception("G4WorkerRunManager::G4WorkerRunManager()","Run0035",FatalException,msg);
#endif
    G4ParticleTable::GetParticleTable()->WorkerG4ParticleTable();
    G4ScoringManager* masterScM = G4MTRunManager::GetMasterScoringManager();
    if(masterScM) G4ScoringManager::GetScoringManager(); //TLS instance for a worker

    eventLoopOnGoing = false;
    nevModulo = -1;
    currEvID = -1;
    workerContext = 0;

    G4UImanager::GetUIpointer()->SetIgnoreCmdNotFound(true);

#ifdef G4MULTITHREADED
    G4VVisManager* pVVis = G4VVisManager::GetConcreteInstance();
    if(pVVis)
    {
      pVVis->SetUpForAThread();
      visIsSetUp = true;
    }
    else
    { visIsSetUp = false; }
#endif
}

#include "G4MTRunManager.hh"

G4WorkerRunManager::~G4WorkerRunManager() {
    // Delete thread-local process manager objects
    physicsList->RemoveProcessManager();

    //Put these pointers to zero: owned by master thread
    //If not to zero, the base class destructor will attempt to
    //delete them
    userDetector = 0;
    userWorkerInitialization = 0;
    userWorkerThreadInitialization = 0;
    userActionInitialization = 0;
    physicsList = 0;
    if(verboseLevel>0) G4cout<<"Destroying WorkerRunManager ("<<this<<")"<<G4endl;
}


void G4WorkerRunManager::InitializeGeometry() {
    if(!userDetector)
    {
        G4Exception("G4RunManager::InitializeGeometry", "Run0033",
                    FatalException, "G4VUserDetectorConstruction is not defined!");
        return;
    }
    //Step1: Get pointer to the physiWorld (note: needs to get the "super pointer, i.e. the one shared by all threads"
    G4RunManagerKernel* masterKernel = G4MTRunManager::GetMasterRunManagerKernel();
    G4VPhysicalVolume* worldVol = masterKernel->GetCurrentWorld();
    //Step2:, Call a new "WorkerDefineWorldVolume( pointer from 2-, false); 
    kernel->WorkerDefineWorldVolume(worldVol,false);
    kernel->SetNumberOfParallelWorld(masterKernel->GetNumberOfParallelWorld());
    //Step3: Call user's ConstructSDandField()
    userDetector->ConstructSDandField();
    userDetector->ConstructParallelSD();
    geometryInitialized = true;
}

void G4WorkerRunManager::RunInitialization()
{
#ifdef G4MULTITHREADED
  if(!visIsSetUp)
  {
    G4VVisManager* pVVis = G4VVisManager::GetConcreteInstance();
    if(pVVis)
    {
      pVVis->SetUpForAThread();
      visIsSetUp = true;
    }
  }
#endif

  if(!(kernel->RunInitialization(fakeRun))) return;

  //Signal this thread can start event loop.
  //Note this will return only when all threads reach this point
  G4MTRunManager::GetMasterRunManager()->ThisWorkerReady();
  if(fakeRun) return;

  const G4UserWorkerInitialization* uwi
       = G4MTRunManager::GetMasterRunManager()->GetUserWorkerInitialization();
  if(currentRun) delete currentRun;
  currentRun = 0;

  //Call a user hook: this is guaranteed all threads are "synchronized"
  if(uwi) uwi->WorkerRunStart();

  if(userRunAction) currentRun = userRunAction->GenerateRun();
  if(!currentRun) currentRun = new G4Run();

  currentRun->SetRunID(runIDCounter);
  currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed);

  currentRun->SetDCtable(DCtable);
  G4SDManager* fSDM = G4SDManager::GetSDMpointerIfExist();
  if(fSDM)
  { currentRun->SetHCtable(fSDM->GetHCtable()); }

  std::ostringstream oss;
    G4Random::saveFullState(oss);
  randomNumberStatusForThisRun = oss.str();
  currentRun->SetRandomNumberStatus(randomNumberStatusForThisRun);

  previousEvents->clear();
  for(G4int i_prev=0;i_prev<n_perviousEventsToBeStored;i_prev++)
  { previousEvents->push_back((G4Event*)0); }

  if(printModulo>0 || verboseLevel>0)
  {
    G4cout << "### Run " << currentRun->GetRunID() << " starts on worker thread "
           << G4Threading::G4GetThreadId() << "." << G4endl;
  }
  if(userRunAction) userRunAction->BeginOfRunAction(currentRun);

  if(storeRandomNumberStatus) {
      G4String fileN = "currentRun";
      if ( rngStatusEventsFlag ) {
          std::ostringstream os;
          os << "run" << currentRun->GetRunID();
          fileN = os.str();
      }
      StoreRNGStatus(fileN);
  }

  runAborted = false;
  numberOfEventProcessed = 0;
}

void G4WorkerRunManager::DoEventLoop(G4int n_event, const char* macroFile , G4int n_select)
{
    if(!userPrimaryGeneratorAction)
    {
      G4Exception("G4RunManager::GenerateEvent()", "Run0032", FatalException,
                "G4VUserPrimaryGeneratorAction is not defined!");
    }

    //This is the same as in the sequential case, just the for-loop indexes are
    //different
    InitializeEventLoop(n_event,macroFile,n_select);

    // Reset random number seeds queue
    while(seedsQueue.size()>0)
    { seedsQueue.pop(); }

    // Event loop
    eventLoopOnGoing = true;
///////    G4int i_event = workerContext->GetThreadId();
    G4int i_event = -1;
    nevModulo = -1;
    currEvID = -1;

    while(eventLoopOnGoing)
    {
      ProcessOneEvent(i_event);
      if(eventLoopOnGoing)
      {
        TerminateOneEvent();
        if(runAborted)
        { eventLoopOnGoing = false; }
//////        else
//////        {
//////          i_event += workerContext->GetNumberThreads();
//////          eventLoopOnGoing = i_event<n_event;
//////        }
      }
    }
     
    TerminateEventLoop();
}

void G4WorkerRunManager::ProcessOneEvent(G4int i_event)
{
  currentEvent = GenerateEvent(i_event);
  if(eventLoopOnGoing)
  {  
    eventManager->ProcessOneEvent(currentEvent);
    AnalyzeEvent(currentEvent);
    UpdateScoring();
    if(currentEvent->GetEventID()<n_select_msg) G4UImanager::GetUIpointer()->ApplyCommand(msgText);
  }
}

G4Event* G4WorkerRunManager::GenerateEvent(G4int i_event)
{
  G4Event* anEvent = new G4Event(i_event);
  long s1, s2;
  long s3 = 0;

  if(i_event<0)
  {
    G4int nevM = G4MTRunManager::GetMasterRunManager()->GetEventModulo();
    if(nevM==1)
    {
      eventLoopOnGoing = G4MTRunManager::GetMasterRunManager()
                       ->SetUpAnEvent(anEvent,s1,s2,s3);
    }
    else
    {
      if(nevModulo<=0)
      {
        G4int nevToDo = G4MTRunManager::GetMasterRunManager()
                         ->SetUpNEvents(anEvent,&seedsQueue);
        if(nevToDo==0)
        { eventLoopOnGoing = false; }
        else
        {
          currEvID = anEvent->GetEventID();
          nevModulo = nevToDo - 1;
        }
      }
      else
      {
        anEvent->SetEventID(++currEvID);
        nevModulo--;
      }
      if(eventLoopOnGoing)
      {
        s1 = seedsQueue.front(); seedsQueue.pop();
        s2 = seedsQueue.front(); seedsQueue.pop();
      }
    }

    if(!eventLoopOnGoing)
    {
      delete anEvent;
      return 0;
    }
  }
  else
  {
    //Need to reseed random number generator
    G4RNGHelper* helper = G4RNGHelper::GetInstance();
    s1 = helper->GetSeed(i_event*2);
    s2 = helper->GetSeed(i_event*2+1);
  }

  long seeds[3] = { s1, s2, 0 };
  G4Random::setTheSeeds(seeds,-1);

  if(storeRandomNumberStatusToG4Event==1 || storeRandomNumberStatusToG4Event==3)
  {
    std::ostringstream oss;
    G4Random::saveFullState(oss);
    randomNumberStatusForThisEvent = oss.str();
    anEvent->SetRandomNumberStatus(randomNumberStatusForThisEvent);
  }

  if(storeRandomNumberStatus) {
      G4String fileN = "currentEvent";
      if ( rngStatusEventsFlag ) {
          std::ostringstream os;
          os << "run" << currentRun->GetRunID() << "evt" << anEvent->GetEventID();
          fileN = os.str();
      }
      StoreRNGStatus(fileN);
  }

  if(printModulo > 0 && anEvent->GetEventID()%printModulo == 0 )
  {
    G4cout << "--> Event " << anEvent->GetEventID() << " starts with initial seeds ("
           << s1 << "," << s2 << ")." << G4endl;
  }
  userPrimaryGeneratorAction->GeneratePrimaries(anEvent);
  return anEvent;
}

void G4WorkerRunManager::MergePartialResults()
{
    //Merge partial results into global run
    G4MTRunManager* mtRM = G4MTRunManager::GetMasterRunManager();
    G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
    if(ScM) mtRM->MergeScores(ScM);
    mtRM->MergeRun(currentRun);
}

void G4WorkerRunManager::RunTermination()
{
  if(!fakeRun)
  {
    MergePartialResults();
    
    //Call a user hook: note this is before the next barrier
    //so threads execute this method asyncrhonouzly
    //(TerminateRun allows for synch via G4RunAction::EndOfRun)
    const G4UserWorkerInitialization* uwi
       = G4MTRunManager::GetMasterRunManager()->GetUserWorkerInitialization();
    if(uwi) uwi->WorkerRunEnd();
  }

  G4RunManager::RunTermination();
  //Signal this thread has finished envent-loop.
  //Note this will return only whan all threads reach this point
  G4MTRunManager::GetMasterRunManager()->ThisWorkerEndEventLoop();

}

/****************************
void G4WorkerRunManager::BeamOn(G4int n_event,const char* macroFile,G4int n_select)
{ 
  if(n_event>0) 
  { G4RunManager::BeamOn(n_event,macroFile,n_select); }
  else
  {
    // fake BeamOn.
    G4MTRunManager::GetMasterRunManager()->ThisWorkerReady();
    G4MTRunManager::GetMasterRunManager()->ThisWorkerEndEventLoop();
  }
}
******************************/

#include "G4AutoLock.hh"
namespace { G4Mutex ConstructScoringWorldsMutex = G4MUTEX_INITIALIZER; }
void G4WorkerRunManager::ConstructScoringWorlds()
{
    // Return if unnecessary
    G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
    if(!ScM) return;
    G4int nPar = ScM->GetNumberOfMesh();
    if(nPar<1) return;

    // Update thread-local G4TransportationManager of all the world volumes
    kernel->WorkerUpdateWorldVolume();

    G4ScoringManager* masterScM = G4MTRunManager::GetMasterScoringManager();
    assert( masterScM != NULL );
    
    G4ParticleTable::G4PTblDicIterator* particleIterator
     = G4ParticleTable::GetParticleTable()->GetIterator();
    
    for(G4int iw=0;iw<nPar;iw++)
    {
      G4VScoringMesh* mesh = ScM->GetMesh(iw);
      G4VPhysicalVolume* pWorld
       = G4TransportationManager::GetTransportationManager()
         ->IsWorldExisting(ScM->GetWorldName(iw));
      if(!pWorld)
      {
        G4ExceptionDescription ed;
        ed<<"Mesh name <"<<ScM->GetWorldName(iw)<<"> is not found in the masther thread.";
        G4Exception("G4WorkerRunManager::ConstructScoringWorlds()","RUN79001",
                      FatalException,ed);
      }
      if(!(mesh->GetMeshElementLogical()) && mesh->GetMeshElementPVList()->empty())
      {
        G4AutoLock l(&ConstructScoringWorldsMutex);
        G4VScoringMesh* masterMesh = masterScM->GetMesh(iw);
        mesh->SetMeshElementLogical(masterMesh->GetMeshElementLogical());
        mesh->SetMeshElementPVList(*masterMesh->GetMeshElementPVList());
        l.unlock();
        
        G4ParallelWorldProcess* theParallelWorldProcess
          = new G4ParallelWorldProcess(ScM->GetWorldName(iw));
        theParallelWorldProcess->SetParallelWorld(ScM->GetWorldName(iw));

        particleIterator->reset();
        while( (*particleIterator)() ){
          G4ParticleDefinition* particle = particleIterator->value();
          G4ProcessManager* pmanager = particle->GetProcessManager();
          if(pmanager)
          {
            pmanager->AddProcess(theParallelWorldProcess);
            if(theParallelWorldProcess->IsAtRestRequired(particle))
            { pmanager->SetProcessOrdering(theParallelWorldProcess, idxAtRest, 9999); }
            pmanager->SetProcessOrderingToSecond(theParallelWorldProcess, idxAlongStep);
            pmanager->SetProcessOrdering(theParallelWorldProcess, idxPostStep, 9999);
          } //if(pmanager)
        }//while
      }
      mesh->WorkerConstruct(pWorld);
    }
}

void G4WorkerRunManager::SetUserInitialization(G4UserWorkerInitialization*)
{
    G4Exception("G4RunManager::SetUserInitialization(G4UserWorkerInitialization*)", "Run3021",
                FatalException, "This method should be used only with an instance of G4MTRunManager");
}

void G4WorkerRunManager::SetUserInitialization(G4UserWorkerThreadInitialization*)
{
    G4Exception("G4RunManager::SetUserInitialization(G4UserWorkerThreadInitialization*)", "Run3021",
                FatalException, "This method should be used only with an instance of G4MTRunManager");
}

void G4WorkerRunManager::SetUserInitialization(G4VUserActionInitialization*)
{
    G4Exception("G4RunManager::SetUserInitialization(G4VUserActionInitialization*)", "Run3021",
                FatalException, "This method should be used only with an instance of G4MTRunManager");
}

void G4WorkerRunManager::SetUserInitialization(G4VUserDetectorConstruction*)
{
    G4Exception("G4RunManager::SetUserInitialization(G4VUserDetectorConstruction*)", "Run3021",
                FatalException, "This method should be used only with an instance of G4MTRunManager");
}

void G4WorkerRunManager::SetUserInitialization(G4VUserPhysicsList* pl)
{
    pl->InitializeWorker();
    G4RunManager::SetUserInitialization(pl);
}

void G4WorkerRunManager::SetUserAction(G4UserRunAction* userAction)
{
    G4RunManager::SetUserAction(userAction);
    userAction->SetMaster(false);
}

void G4WorkerRunManager::SetupDefaultRNGEngine()
{
    const CLHEP::HepRandomEngine* mrnge = G4MTRunManager::GetMasterRunManager()->getMasterRandomEngine();
    assert(mrnge);//Master has created RNG
    const G4UserWorkerThreadInitialization* uwti
      =G4MTRunManager::GetMasterRunManager()->GetUserWorkerThreadInitialization();
    uwti->SetupRNGEngine(mrnge);
}


//Forward calls (avoid GCC compilation warnings)
void G4WorkerRunManager::SetUserAction(G4UserEventAction* ua)
{
    G4RunManager::SetUserAction(ua);
}

void G4WorkerRunManager::SetUserAction(G4VUserPrimaryGeneratorAction* ua)
{
    G4RunManager::SetUserAction(ua);
}

void G4WorkerRunManager::SetUserAction(G4UserStackingAction* ua)
{
    G4RunManager::SetUserAction(ua);
}

void G4WorkerRunManager::SetUserAction(G4UserTrackingAction* ua)
{
    G4RunManager::SetUserAction(ua);
}

void G4WorkerRunManager::SetUserAction(G4UserSteppingAction* ua)
{
    G4RunManager::SetUserAction(ua);
}

void G4WorkerRunManager::StoreRNGStatus(const G4String& fn )
{
    std::ostringstream os;
    os << randomNumberStatusDir << "G4Worker"<<workerContext->GetThreadId()<<"_"<<fn <<".rndm";
    G4Random::saveEngineStatus(os.str().c_str());    
}

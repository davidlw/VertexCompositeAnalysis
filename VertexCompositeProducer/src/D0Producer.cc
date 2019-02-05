// -*- C++ -*-
//
// Package:    VertexCompositeProducer
//
// Class:      D0Producer
// 
/**\class D0Producer D0Producer.cc VertexCompositeAnalysis/VertexCompositeProducer/src/D0Producer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li
//
//


// system include files
#include <memory>

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/D0Producer.h"

// Constructor
D0Producer::D0Producer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
{
  useAnyMVA_ = false;
  if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
 
  produces< reco::VertexCompositeCandidateCollection >("D0");
  produces< TrkIndexCollection >("negTrkIndex"); // mtd
  produces< TrkIndexCollection >("posTrkIndex"); // mtd
  if(useAnyMVA_) produces<MVACollection>("MVAValuesD0");
}

// (Empty) Destructor
D0Producer::~D0Producer() {
}


//
// Methods
//

// Producer Method
void D0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create D0Fitter object which reconstructs the vertices and creates
//   D0Fitter theVees(theParams, iEvent, iSetup);

   theVees.fitAll(iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
//   std::auto_ptr< reco::VertexCompositeCandidateCollection >
//     d0Candidates( new reco::VertexCompositeCandidateCollection );
//
   auto d0Candidates = std::make_unique<reco::VertexCompositeCandidateCollection>();
   d0Candidates->reserve( theVees.getD0().size() );

   std::copy( theVees.getD0().begin(),
              theVees.getD0().end(),
              std::back_inserter(*d0Candidates) );

   // create auto ptr for indices of tracks
   
   auto negTrkIndexPtr = std::make_unique<TrkIndexCollection>( // mtd
                 theVees.getNegTrkIndex().begin(),theVees.getNegTrkIndex().end());
   auto posTrkIndexPtr = std::make_unique<TrkIndexCollection>( // mtd
                 theVees.getPosTrkIndex().begin(),theVees.getPosTrkIndex().end());
   

   // Write the collections to the Event
   iEvent.put( std::move(d0Candidates), std::string("D0") );
   iEvent.put( std::move(negTrkIndexPtr), std::string("negTrkIndex") ); // mtd
   iEvent.put( std::move(posTrkIndexPtr), std::string("posTrkIndex") );
    
   if(useAnyMVA_) 
   {
     auto mvas = std::make_unique<MVACollection>(theVees.getMVAVals().begin(),theVees.getMVAVals().end());
     iEvent.put(std::move(mvas), std::string("MVAValuesD0"));
   }

   theVees.resetAll();
}


//void D0Producer::beginJob() {
void D0Producer::beginJob() {
}


void D0Producer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(D0Producer);

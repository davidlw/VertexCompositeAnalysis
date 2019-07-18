// -*- C++ -*-
//
// Package:    VertexCompositeProducer
//
// Class:      D0ProducerNew
// 
/**\class D0ProducerNew D0ProducerNew.cc VertexCompositeAnalysis/VertexCompositeProducer/src/D0ProducerNew.cc

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

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/D0ProducerNew.h"

// Constructor
D0ProducerNew::D0ProducerNew(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
{
  useAnyMVA_ = false;
  if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
 
  produces< reco::VertexCompositeCandidateCollection >("D0");
  if(useAnyMVA_) produces<MVACollection>("MVAValuesD0");
}

// (Empty) Destructor
D0ProducerNew::~D0ProducerNew() {
}


//
// Methods
//

// Producer Method
void D0ProducerNew::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create D0FitterNew object which reconstructs the vertices and creates
//   D0FitterNew theVees(theParams, iEvent, iSetup);

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

   // Write the collections to the Event
   iEvent.put( std::move(d0Candidates), std::string("D0") );
    
   if(useAnyMVA_) 
   {
     auto mvas = std::make_unique<MVACollection>(theVees.getMVAVals().begin(),theVees.getMVAVals().end());
     iEvent.put(std::move(mvas), std::string("MVAValuesD0"));
   }

   theVees.resetAll();
}


//void D0ProducerNew::beginJob() {
void D0ProducerNew::beginJob() {
}


void D0ProducerNew::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(D0ProducerNew);

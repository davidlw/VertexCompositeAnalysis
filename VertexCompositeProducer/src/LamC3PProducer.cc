// -*- C++ -*-
//
// Package:    VertexCompositeProducer
//
// Class:      LamC3PProducer
// 
/**\class LamC3PProducer LamC3PProducer.cc VertexCompositeAnalysis/VertexCompositeProducer/src/LamC3PProducer.cc

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

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/LamC3PProducer.h"

// Constructor
LamC3PProducer::LamC3PProducer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
{
  useAnyMVA_ = false;
  if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
 
  produces< pat::CompositeCandidateCollection >("LamC3P");
  if(useAnyMVA_) produces<MVACollection>("MVAValuesLamC3P");
}

// (Empty) Destructor
LamC3PProducer::~LamC3PProducer() {
}


//
// Methods
//

// Producer Method
void LamC3PProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create LamC3PFitter object which reconstructs the vertices and creates
//   LamC3PFitter theVees(theParams, iEvent, iSetup);

   theVees.fitAll(iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
//   std::auto_ptr< reco::VertexCompositeCandidateCollection >
//     lamCCandidates( new reco::VertexCompositeCandidateCollection );
//
   auto lamCCandidates = std::make_unique<pat::CompositeCandidateCollection>();
   lamCCandidates->reserve( theVees.getLamC3P().size() );

   std::copy( theVees.getLamC3P().begin(),
              theVees.getLamC3P().end(),
              std::back_inserter(*lamCCandidates) );

   // Write the collections to the Event
   iEvent.put( std::move(lamCCandidates), std::string("LamC3P") );
    
   if(useAnyMVA_) 
   {
     auto mvas = std::make_unique<MVACollection>(theVees.getMVAVals().begin(),theVees.getMVAVals().end());
     iEvent.put(std::move(mvas), std::string("MVAValuesLamC3P"));
   }

   theVees.resetAll();
}


//void LamC3PProducer::beginJob() {
void LamC3PProducer::beginJob() {
}


void LamC3PProducer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(LamC3PProducer);

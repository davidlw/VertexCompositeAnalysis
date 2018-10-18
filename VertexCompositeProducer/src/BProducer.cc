// -*- C++ -*-
//
// Package:    VertexCompositeProducer
//
// Class:      BProducer
// 
/**\class BProducer BProducer.cc VertexCompositeAnalysis/VertexCompositeProducer/src/BProducer.cc

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

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/BProducer.h"

// Constructor
BProducer::BProducer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
{
//  useAnyMVA_ = false;
//  if(iConfig.exists("useAnyMVA")) useAnyMVA_ = iConfig.getParameter<bool>("useAnyMVA");
 
  produces< reco::VertexCompositeCandidateCollection >("B");
//  if(useAnyMVA_) produces<MVACollection>("MVAValuesB");
}

// (Empty) Destructor
BProducer::~BProducer() {
}


//
// Methods
//

// Producer Method
void BProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create BFitter object which reconstructs the vertices and creates
//   BFitter theVees(theParams, iEvent, iSetup);

   theVees.fitAll(iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
//   std::auto_ptr< reco::VertexCompositeCandidateCollection >
//     bCandidates( new reco::VertexCompositeCandidateCollection );
//
   auto bCandidates = std::make_unique<reco::VertexCompositeCandidateCollection>();
   bCandidates->reserve( theVees.getB().size() );

   std::copy( theVees.getB().begin(),
              theVees.getB().end(),
              std::back_inserter(*bCandidates) );

   // Write the collections to the Event
   iEvent.put( std::move(bCandidates), std::string("B") );
/*    
   if(useAnyMVA_) 
   {
     auto mvas = std::make_unique<MVACollection>(theVees.getMVAVals().begin(),theVees.getMVAVals().end());
     iEvent.put(std::move(mvas), std::string("MVAValuesB"));
   }
*/
   theVees.resetAll();
}


//void BProducer::beginJob() {
void BProducer::beginJob() {
}


void BProducer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(BProducer);

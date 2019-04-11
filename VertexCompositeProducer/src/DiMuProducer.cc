// -*- C++ -*-
//
// Package:    DiMuProducer
//
// Class:      DiMuProducer
// 
/**\class DiMuProducer DiMuProducer.cc VertexCompositeAnalysis/VertexCompositeProducer/src/DiMuProducer.cc

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

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/DiMuProducer.h"

// Constructor
DiMuProducer::DiMuProducer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
//  theParams(iConfig) {
{
  // Trying this with Candidates instead of the simple reco::Vertex
  produces< pat::CompositeCandidateCollection >("DiMu");
}

// (Empty) Destructor
DiMuProducer::~DiMuProducer() {
}


//
// Methods
//

// Producer Method
void DiMuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create DiMuFitter object which reconstructs the vertices and creates
   theVees.fitAll(iEvent, iSetup);

   // Create auto_ptr for each collection to be stored in the Event
   auto dimuCandidates = std::make_unique<pat::CompositeCandidateCollection>();

   dimuCandidates->reserve( theVees.getDiMu().size() );

   std::copy( theVees.getDiMu().begin(),
              theVees.getDiMu().end(),
              std::back_inserter(*dimuCandidates) );

   // Write the collections to the Event
   iEvent.put( std::move(dimuCandidates), std::string("DiMu") ); 

   theVees.resetAll();
}


//void DiMuProducer::beginJob() {
void DiMuProducer::beginJob() {
}


void DiMuProducer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(DiMuProducer);

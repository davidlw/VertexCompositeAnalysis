// -*- C++ -*-
//
// Package:    DiMuProducer
//
// Class:      DiMuProducer
// 
/**\class DiMuProducer DiMuProducer.cc VertexCompositeAnalysis/HadronCompositeProducer/src/DiMuProducer.cc

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

#include "VertexCompositeAnalysis/HadronCompositeProducer/interface/DiMuProducer.h"

// Constructor
DiMuProducer::DiMuProducer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
//  theParams(iConfig) {
{
  // Trying this with Candidates instead of the simple reco::Vertex
  produces< reco::VertexCompositeCandidateCollection >("DiMu");
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
   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     dimuCandidates( new reco::VertexCompositeCandidateCollection );
   dimuCandidates->reserve( theVees.getDiMu().size() );

   std::copy( theVees.getDiMu().begin(),
              theVees.getDiMu().end(),
              std::back_inserter(*dimuCandidates) );

   // Write the collections to the Event
   iEvent.put( dimuCandidates, std::string("DiMu") ); 

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

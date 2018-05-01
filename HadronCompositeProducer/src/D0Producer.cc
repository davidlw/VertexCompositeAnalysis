// -*- C++ -*-
//
// Package:    HadronCompositeProducer
//
// Class:      D0Producer
// 
/**\class D0Producer D0Producer.cc VertexCompositeAnalysis/HadronCompositeProducer/src/D0Producer.cc

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

#include "VertexCompositeAnalysis/HadronCompositeProducer/interface/D0Producer.h"

// Constructor
D0Producer::D0Producer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
//  theParams(iConfig) {
{
  // Trying this with Candidates instead of the simple reco::Vertex
  produces< reco::VertexCompositeCandidateCollection >("D0");

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
   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     d0Candidates( new reco::VertexCompositeCandidateCollection );
   d0Candidates->reserve( theVees.getD0().size() );

   std::copy( theVees.getD0().begin(),
              theVees.getD0().end(),
              std::back_inserter(*d0Candidates) );

   // Write the collections to the Event
   iEvent.put( d0Candidates, std::string("D0") );
    
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

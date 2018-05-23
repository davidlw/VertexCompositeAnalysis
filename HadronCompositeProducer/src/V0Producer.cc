// -*- C++ -*-
//
// Package:    HadronCompositeProducer
//
// Class:      V0Producer
// 
/**\class V0Producer V0Producer.cc  VertexCompositeAnalysis/HadronCompositeProducer/src/V0Producer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//


// system include files
#include <memory>

#include "VertexCompositeAnalysis/HadronCompositeProducer/interface/V0Producer.h"

// Constructor
V0Producer::V0Producer(const edm::ParameterSet& iConfig) :
 theVees(iConfig, consumesCollector())
//  theParams(iConfig) {
{
   // Registering V0 Collections
  //produces<reco::VertexCollection>("Kshort");
  //produces<reco::VertexCollection>("Lambda");
  //produces<reco::VertexCollection>("LambdaBar");

  // Trying this with Candidates instead of the simple reco::Vertex
  produces< reco::VertexCompositeCandidateCollection >("Kshort");
  produces< reco::VertexCompositeCandidateCollection >("Phi");
  produces< reco::VertexCompositeCandidateCollection >("Lambda");
  produces< reco::VertexCompositeCandidateCollection >("Xi");
  produces< reco::VertexCompositeCandidateCollection >("Omega");
  produces< reco::VertexCompositeCandidateCollection >("D0");
  produces< reco::VertexCompositeCandidateCollection >("DSToKsK");
  produces< reco::VertexCompositeCandidateCollection >("DSToPhiPi");
  produces< reco::VertexCompositeCandidateCollection >("DPM");
  produces< reco::VertexCompositeCandidateCollection >("LambdaCToLamPi");
  produces< reco::VertexCompositeCandidateCollection >("LambdaCToKsP");
  //produces< reco::VertexCompositeCandidateCollection >("LambdaBar");

}

// (Empty) Destructor
V0Producer::~V0Producer() {
}


//
// Methods
//

// Producer Method
void V0Producer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   using namespace edm;

   // Create V0Fitter object which reconstructs the vertices and creates
   //  (and contains) collections of Kshorts, Lambda0s
//   V0Fitter theVees(theParams, iEvent, iSetup);

   theVees.fitAll(iEvent, iSetup);

   // Create unique_ptr for each collection to be stored in the Event
   std::unique_ptr< reco::VertexCompositeCandidateCollection > 
     kShortCandidates( new reco::VertexCompositeCandidateCollection );
   kShortCandidates->reserve( theVees.getKshorts().size() ); 

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     phiCandidates( new reco::VertexCompositeCandidateCollection );
   phiCandidates->reserve( theVees.getPhis().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     lambdaCandidates( new reco::VertexCompositeCandidateCollection );
   lambdaCandidates->reserve( theVees.getLambdas().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     xiCandidates( new reco::VertexCompositeCandidateCollection );
   xiCandidates->reserve( theVees.getXis().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     omegaCandidates( new reco::VertexCompositeCandidateCollection );
   omegaCandidates->reserve( theVees.getOmegas().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     d0Candidates( new reco::VertexCompositeCandidateCollection );
   d0Candidates->reserve( theVees.getD0().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     dsCandidates1( new reco::VertexCompositeCandidateCollection );
   dsCandidates1->reserve( theVees.getDSToKsK().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     dsCandidates2( new reco::VertexCompositeCandidateCollection );
   dsCandidates2->reserve( theVees.getDSToPhiPi().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     dpmCandidates( new reco::VertexCompositeCandidateCollection );
   dpmCandidates->reserve( theVees.getDPM().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     lambdaCCandidates1( new reco::VertexCompositeCandidateCollection );
   lambdaCCandidates1->reserve( theVees.getLambdaCToLamPi().size() );

   std::unique_ptr< reco::VertexCompositeCandidateCollection >
     lambdaCCandidates2( new reco::VertexCompositeCandidateCollection );
   lambdaCCandidates2->reserve( theVees.getLambdaCToKsP().size() );

   std::copy( theVees.getKshorts().begin(),
	      theVees.getKshorts().end(),
	      std::back_inserter(*kShortCandidates) );
   std::copy( theVees.getPhis().begin(),
              theVees.getPhis().end(),
              std::back_inserter(*phiCandidates) );
   std::copy( theVees.getLambdas().begin(),
	      theVees.getLambdas().end(),
	      std::back_inserter(*lambdaCandidates) );
   std::copy( theVees.getXis().begin(),
              theVees.getXis().end(),
              std::back_inserter(*xiCandidates) );
   std::copy( theVees.getOmegas().begin(),
              theVees.getOmegas().end(),
              std::back_inserter(*omegaCandidates) );
   std::copy( theVees.getD0().begin(),
              theVees.getD0().end(),
              std::back_inserter(*d0Candidates) );
   std::copy( theVees.getDSToKsK().begin(),
              theVees.getDSToKsK().end(),
              std::back_inserter(*dsCandidates1) );
   std::copy( theVees.getDSToPhiPi().begin(),
              theVees.getDSToPhiPi().end(),
              std::back_inserter(*dsCandidates2) );
   std::copy( theVees.getDPM().begin(),
              theVees.getDPM().end(),
              std::back_inserter(*dpmCandidates) );
   std::copy( theVees.getLambdaCToLamPi().begin(),
              theVees.getLambdaCToLamPi().end(),
              std::back_inserter(*lambdaCCandidates1) );
   std::copy( theVees.getLambdaCToKsP().begin(),
              theVees.getLambdaCToKsP().end(),
              std::back_inserter(*lambdaCCandidates2) );

   // Write the collections to the Event
   iEvent.put( std::move(kShortCandidates), std::string("Kshort") );
   iEvent.put( std::move(phiCandidates), std::string("Phi") );
   iEvent.put( std::move(lambdaCandidates), std::string("Lambda") );
   iEvent.put( std::move(xiCandidates), std::string("Xi") );
   iEvent.put( std::move(omegaCandidates), std::string("Omega") );
   iEvent.put( std::move(d0Candidates), std::string("D0") );
   iEvent.put( std::move(dsCandidates1), std::string("DSToKsK") );
   iEvent.put( std::move(dsCandidates2), std::string("DSToPhiPi") );
   iEvent.put( std::move(dpmCandidates), std::string("DPM") );
   iEvent.put( std::move(lambdaCCandidates1), std::string("LambdaCToLamPi") );
   iEvent.put( std::move(lambdaCCandidates2), std::string("LambdaCToKsP") );

   theVees.resetAll();
}


//void V0Producer::beginJob() {
void V0Producer::beginJob() {
}


void V0Producer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(V0Producer);
//DEFINE_FWK_MODULE(V0finder);

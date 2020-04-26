// -*- C++ -*-
//
// Package:    ParticleProducer
//
// Class:      ParticleProducer
// 
/**\class ParticleProducer ParticleProducer.cc VertexCompositeAnalysis/VertexCompositeProducer/src/ParticleProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andre Stahl
//
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ParticleProducer.h"


// constructor
ParticleProducer::ParticleProducer(const edm::ParameterSet& iConfig) :
  fitter_(iConfig, consumesCollector())
{
  produces<pat::GenericParticleCollection>();
}

// dDestructor
ParticleProducer::~ParticleProducer() {
}


//
// Methods
//

// producer method
void ParticleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
   // fit particles
   fitter_.fitAll(iEvent, iSetup);
   // extract particles
   const auto& particles = fitter_.getParticles();
   // store particles
   auto output = std::make_unique<pat::GenericParticleCollection>(particles);
   iEvent.put(std::move(output));
   // reset fitter
   fitter_.resetAll();
}


void ParticleProducer::beginJob() {
}


void ParticleProducer::endJob() {
}

//define this as a plug-in
#include "FWCore/PluginManager/interface/ModuleDef.h"

DEFINE_FWK_MODULE(ParticleProducer);

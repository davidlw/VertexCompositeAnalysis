// -*- C++ -*-
//
// Package:    ParticleProducer
// Class:      ParticleProducer
// 
/**\class ParticleProducer ParticleProducer.h VertexCompositeAnalysis/VertexCompositeProducer/interface/ParticleProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Andre Stahl 
//
//

#ifndef VertexCompositeAnalysis__Particle_PRODUCER_H
#define VertexCompositeAnalysis__Particle_PRODUCER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/ParticleFitter.h"

class ParticleProducer : public edm::EDProducer {
public:
  explicit ParticleProducer(const edm::ParameterSet&);
  ~ParticleProducer();

private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  ParticleFitter fitter_;
  ParticleDaughter daughter_;
};

#endif

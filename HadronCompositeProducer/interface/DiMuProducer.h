// -*- C++ -*-
//
// Package:    DiMuProducer
// Class:      DiMuProducer
// 
/**\class DiMuProducer DiMuProducer.h VertexCompositeAnalysis/HadronCompositeProducer/interface/DiMuProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li 
//
//

#ifndef VertexCompositeAnalysis__DiMu_PRODUCER_H
#define VertexCompositeAnalysis__DiMu_PRODUCER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "VertexCompositeAnalysis/HadronCompositeProducer/interface/DiMuFitter.h"

class DiMuProducer : public edm::EDProducer {
public:
  explicit DiMuProducer(const edm::ParameterSet&);
  ~DiMuProducer();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  DiMuFitter theVees; 
};

#endif

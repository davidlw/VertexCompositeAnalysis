// -*- C++ -*-
//
// Package:    HadronCompositeProducer
// Class:      V0Producer
// 
/**\class V0Producer V0Producer.h VertexCompositeAnalysis/HadronCompositeProducer/interface/V0Producer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//

#ifndef VertexCompositeAnalysis__V0_PRODUCER_H
#define VertexCompositeAnalysis__V0_PRODUCER_H

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
//#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "VertexCompositeAnalysis/HadronCompositeProducer/interface/V0Fitter.h"

class V0Producer : public edm::EDProducer {
public:
  explicit V0Producer(const edm::ParameterSet&);
  ~V0Producer();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  V0Fitter theVees; 
//  edm::ParameterSet theParams;
};

#endif

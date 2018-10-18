// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      BProducer
// 
/**\class BProducer BProducer.h VertexCompositeAnalysis/VertexCompositeProducer/interface/BProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li 
//
//

#ifndef VertexCompositeAnalysis__B_PRODUCER_H
#define VertexCompositeAnalysis__B_PRODUCER_H

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

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/BFitter.h"

class BProducer : public edm::EDProducer {
public:
//  using MVACollection = std::vector<float>;

  explicit BProducer(const edm::ParameterSet&);
  ~BProducer();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

//  bool useAnyMVA_;

  BFitter theVees; 
//  edm::ParameterSet theParams;
};

#endif

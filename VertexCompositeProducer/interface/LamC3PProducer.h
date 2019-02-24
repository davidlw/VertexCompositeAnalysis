// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      LamC3PProducer
// 
/**\class LamC3PProducer LamC3PProducer.h VertexCompositeAnalysis/VertexCompositeProducer/interface/LamC3PProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Wei Li 
//
//

#ifndef VertexCompositeAnalysis__LAMC3P_PRODUCER_H
#define VertexCompositeAnalysis__LAMC3P_PRODUCER_H

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
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/LamC3PFitter.h"

class LamC3PProducer : public edm::EDProducer {
public:
  using MVACollection = std::vector<float>;

  explicit LamC3PProducer(const edm::ParameterSet&);
  ~LamC3PProducer();

private:
  //virtual void beginJob() ;
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  bool useAnyMVA_;

  LamC3PFitter theVees; 
//  edm::ParameterSet theParams;
};

#endif

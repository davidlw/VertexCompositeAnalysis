// -*- C++ -*-
//
// Package:    VertexCompositeAnalysis/VertexCompositeProducer
// Class:      PackedPtrToObjConverter
// 
/**\class PackedPtrToObjConverter PackedPtrToObjConverter.cc VertexCompositeAnalysis/VertexCompositeProducer/plugins/PackedPtrToObjConverter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yousen Zhang
//         Created:  Wed, 28 Apr 2021 23:25:52 GMT
//
//


// system include files
#include <memory>
#include <type_traits>
#include <typeinfo>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/FwdPtr.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

//
// class declaration
//

template <class T, bool isFwd>
class PackedPtrToObjConverter : public edm::stream::EDProducer<> {
public:
  using MyPtr = typename std::conditional<isFwd, edm::FwdPtr<T>, edm::Ptr<T>>::type;
  using Coll = std::vector<MyPtr>;
  //
  // constructors and destructor
  //
  explicit  PackedPtrToObjConverter(const edm::ParameterSet& iConfig)
    : srcToken_(consumes<Coll>(iConfig.getParameter<edm::InputTag>("src")))
  {
    //register your products
    produces<std::vector<T>>("constituents");
    //now do what ever other initialization is needed
  }
  ~PackedPtrToObjConverter()
  {
    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
  }

  //
  // member functions
  //

  // ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  static void
  fillDescriptions(edm::ConfigurationDescriptions& descriptions)
  {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
  }


  private:

  // ------------ method called to produce the data  ------------

  void
  produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
  {
    // using namespace edm;
    auto candsOut = std::make_unique<std::vector<T>>();

    edm::Handle<Coll> srcHandle;
    iEvent.getByToken(srcToken_, srcHandle);

    candsOut->reserve(srcHandle->size());

    for (const auto e : *srcHandle.product()) {
      candsOut->push_back(T(*e));
    }
    iEvent.put(std::move(candsOut), "constituents");

  }

  // ----------member data ---------------------------
  edm::EDGetTokenT<Coll> const srcToken_;

};

// using GenJetConstituentPackedFwdToObj = PackedPtrToObjConverter<pat::PackedGenParticle, true>;
using GenJetConstituentPackedPtrToObj = PackedPtrToObjConverter<pat::PackedGenParticle, false>;
using PatJetConstituentPackedFwdToObj = PackedPtrToObjConverter<pat::PackedCandidate, true>;
using PatJetConstituentPackedPtrToObj = PackedPtrToObjConverter<pat::PackedCandidate, false>;
using PFJetConstituentPackedFwdToObj = PackedPtrToObjConverter<pat::PackedCandidate, true>;
//using PFJetConstituentPackedPtrToObj = PackedPtrToObjConverter<pat::PackedCandidate, false>;


//define this as a plug-in

DEFINE_FWK_MODULE(GenJetConstituentPackedPtrToObj);
DEFINE_FWK_MODULE(PatJetConstituentPackedFwdToObj);
DEFINE_FWK_MODULE(PatJetConstituentPackedPtrToObj);
DEFINE_FWK_MODULE(PFJetConstituentPackedFwdToObj);

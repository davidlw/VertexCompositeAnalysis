// system include files
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "TTree.h"
#include "TRegexp.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Luminosity/interface/LumiInfo.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/OnlineMetaData/interface/OnlineLuminosityRecord.h"
#include "DataFormats/EgammaCandidates/interface/HIPhotonIsolation.h"
#include "DataFormats/Math/interface/angle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//
// constants, enums and typedefs
//

//
// class decleration
//

class ParticleAnalyzer : public edm::EDAnalyzer {
public:
  explicit ParticleAnalyzer(const edm::ParameterSet&);
  ~ParticleAnalyzer();

private:
  virtual void beginJob();
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void getEventData(const edm::Event&, const edm::EventSetup&);
  virtual void getTriggerData(const edm::Event&, const edm::EventSetup&);
  virtual void fillEventInfo(const edm::Event&);
  virtual void fillTriggerInfo(const edm::Event&, const edm::EventSetup&);
  virtual void fillLumiInfo(const edm::Event&);
  virtual void fillRecoParticleInfo(const edm::Event&);
  virtual void fillGenParticleInfo(const edm::Event&);
  virtual void endJob();
  virtual void initTree();

  UShort_t fillTriggerObjectInfo(const pat::TriggerObject&, const UShort_t&, const UShort_t& candIdx=USHRT_MAX);
  UShort_t fillRecoParticleInfo(const pat::GenericParticle&, const UShort_t& momIdx=USHRT_MAX);
  UShort_t fillTrackInfo(const pat::GenericParticle&, const UShort_t&, const bool& force=false);
  UShort_t fillSourceInfo(const pat::GenericParticle&, const int&, const UShort_t&, const bool& force=false);
  UShort_t fillMuonInfo(const pat::GenericParticle& cand, const UShort_t&, const bool& force=false);
  UShort_t fillElectronInfo(const pat::GenericParticle& cand, const UShort_t&, const bool& force=false);
  UShort_t fillPhotonInfo(const pat::GenericParticle& cand, const UShort_t&, const bool& force=false);
  UShort_t fillJetInfo(const pat::GenericParticle& cand, const UShort_t&, const bool& force=false);
  UShort_t fillTauInfo(const pat::GenericParticle& cand, const UShort_t&, const bool& force=false);
  UShort_t fillPFCandidateInfo(const pat::GenericParticle& cand, const UShort_t&, const bool& force=false);
  UShort_t fillGenParticleInfo(const reco::GenParticleRef&, const UShort_t& candIdx=USHRT_MAX, const bool& force=false);

  void initParticleInfo(const std::string&, const int& pId=0);
  void addTriggerObject(pat::GenericParticle&);
  bool addTriggerObject(pat::GenericParticle&, const pat::TriggerObjectCollection&, const std::string&, const std::string&, const int&);
  bool addGenParticle(pat::GenericParticle&, const reco::GenParticleRefVector&);
  void findGenDaughters(reco::GenParticleRefVector&, const reco::GenParticleRef&, const pat::GenericParticle&, const short& iter=0);
  reco::GenParticleRef findGenDaughter(const reco::GenParticleRef&, const size_t&);
  reco::GenParticleRef findGenMother(const reco::GenParticleRef&);

  int   getInt  (const pat::GenericParticle& c, const std::string& n, const int&   d=-99  ) const { return (c.hasUserInt(n)   ? c.userInt(n)   : d); }
  float getFloat(const pat::GenericParticle& c, const std::string& n, const float& d=-99.9) const { return (c.hasUserFloat(n) ? c.userFloat(n) : d); }
  template <class T> Char_t   getChar  (const T& n) const { if (!(n >= CHAR_MIN && n < CHAR_MAX )) { throw(cms::Exception("Overflow")<<"Invalid char: "  <<n); }; return Char_t(n);   }
  template <class T> UChar_t  getUChar (const T& n) const { if (!(n >= 0        && n < UCHAR_MAX)) { throw(cms::Exception("Overflow")<<"Invalid uchar: " <<n); }; return UChar_t(n);  }
  template <class T> Short_t  getShort (const T& n) const { if (!(n >= SHRT_MIN && n < SHRT_MAX )) { throw(cms::Exception("Overflow")<<"Invalid short: " <<n); }; return Short_t(n);  }
  template <class T> UShort_t getUShort(const T& n) const { if (!(n >= 0        && n < USHRT_MAX)) { throw(cms::Exception("Overflow")<<"Invalid ushort: "<<n); }; return UShort_t(n); }
  template <class T> Int_t    getInt   (const T& n) const { if (!(n >= INT_MIN  && n < INT_MAX  )) { throw(cms::Exception("Overflow")<<"Invalid int: "   <<n); }; return Int_t(n);    }
  template <class T> UInt_t   getUInt  (const T& n) const { if (!(n >= 0        && n < UINT_MAX )) { throw(cms::Exception("Overflow")<<"Invalid uint: "  <<n); }; return UInt_t(n);   }

  bool isCompatible(const reco::Candidate& p1, const reco::Candidate& p2) const
  {
    return ((p1.status()==1)==(p2.status()==1) && p1.charge()==p2.charge() && std::abs(p1.pdgId())!=std::abs(p2.pdgId()));
  }
  double deltaPt(const double& pT1, const double& pT2) const
  {
    return std::abs(pT1 - pT2)/std::sqrt(pT1*pT2);
  }
  bool isMatched(const reco::Candidate::LorentzVector& c1, const reco::Candidate::LorentzVector& c2, const double& maxDeltaR, const double& maxDeltaPtRel)
  {
    const auto deltaR = reco::deltaR(c1.Eta(), c1.Phi(), c2.Eta(), c2.Phi());
    const auto dPtRel = deltaPt(c1.Pt(), c2.Pt());
    return (deltaR < maxDeltaR && dPtRel < maxDeltaPtRel);
  }

  // ----------member data ---------------------------

  // input tokens
  const edm::EDGetTokenT<reco::BeamSpot> tok_offlineBS_;
  const edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
  const edm::EDGetTokenT<pat::GenericParticleCollection> tok_recParticle_;
  const edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
  const edm::EDGetTokenT<GenEventInfoProduct> tok_genInfo_;
  const edm::EDGetTokenT<int> tok_centBin_;
  const edm::EDGetTokenT<reco::Centrality> tok_centSrc_;
  const edm::EDGetTokenT<reco::EvtPlaneCollection> tok_eventPlaneSrc_;
  const edm::EDGetTokenT<edm::TriggerResults> tok_triggerResults_;
  const edm::EDGetTokenT<trigger::TriggerEvent> tok_triggerEvent_;
  const edm::EDGetTokenT<edm::TriggerResults> tok_filterResults_;
  const edm::EDGetTokenT<LumiInfo> tok_lumiInfo_;
  const edm::EDGetTokenT<LumiScalersCollection> tok_lumiScalers_;
  const edm::EDGetTokenT<OnlineLuminosityRecord> tok_lumiRecord_;

  // input data
  const std::vector<edm::ParameterSet> triggerInfo_;
  const std::vector<std::string> eventFilters_;
  const std::string selectEvents_;
  const bool saveTree_, addTrack_, addSource_;
  const int maxGenIter_;
  const double maxGenDeltaR_, maxGenDeltaPtRel_, maxTrgDeltaR_, maxTrgDeltaPtRel_;
  const std::vector<int> dauIDs_;

  HLTPrescaleProvider hltPrescaleProvider_;

  // attributes
  edm::Service<TFileService> fileService_;
  TTree* tree_ = 0;
  TTree* ntuple_ = 0;

  bool isMC_;
  reco::Vertex vertex_;
  reco::VertexCollection vertices_;
  reco::GenParticleRefVector genParticles_;
  const MagneticField* magField_;

  // containers
  class Container
  {
   public:
    Container() {};
    ~Container() {};

    // getters
    std::map<std::string, bool     > boolM()   const { return boolM_;   };
    std::map<std::string, Char_t   > charM()   const { return charM_;   };
    std::map<std::string, short    > shortM()  const { return shortM_;  };
    std::map<std::string, int      > intM()    const { return intM_;    };
    std::map<std::string, UChar_t  > ucharM()  const { return ucharM_;  };
    std::map<std::string, UShort_t > ushortM() const { return ushortM_; };
    std::map<std::string, UInt_t   > uintM()   const { return uintM_;   };
    std::map<std::string, float    > floatM()  const { return floatM_;  };
    std::map<std::string, std::vector<bool>     > boolVM()   const { return boolVM_;   };
    std::map<std::string, std::vector<Char_t>   > charVM()   const { return charVM_;   };
    std::map<std::string, std::vector<short>    > shortVM()  const { return shortVM_;  };
    std::map<std::string, std::vector<int>      > intVM()    const { return intVM_;    };
    std::map<std::string, std::vector<UChar_t>  > ucharVM()  const { return ucharVM_;  };
    std::map<std::string, std::vector<UShort_t> > ushortVM() const { return ushortVM_; };
    std::map<std::string, std::vector<UInt_t>   > uintVM()   const { return uintVM_;   };
    std::map<std::string, std::vector<float>    > floatVM()  const { return floatVM_;  };

    // setters
    template <class T>
    void add(const std::string& n, const T&        v) = delete; // avoid implicit conversion
    void add(const std::string& n, const bool&     v) { boolM_[n]   = v; };
    void add(const std::string& n, const Char_t&   v) { charM_[n]   = v; };
    void add(const std::string& n, const short&    v) { shortM_[n]  = v; };
    void add(const std::string& n, const int&      v) { intM_[n]    = v; };
    void add(const std::string& n, const UChar_t&  v) { ucharM_[n]  = v; };
    void add(const std::string& n, const UShort_t& v) { ushortM_[n] = v; };
    void add(const std::string& n, const UInt_t&   v) { uintM_[n]   = v; };
    void add(const std::string& n, const float&    v) { floatM_[n]  = v; };
    void add(const std::string& n, const double&   v) { floatM_[n]  = v; };
    void add(const std::string& n, const std::vector<bool>&     v) { boolVM_[n]   = v; };
    void add(const std::string& n, const std::vector<Char_t>&   v) { charVM_[n]   = v; };
    void add(const std::string& n, const std::vector<short>&    v) { shortVM_[n]  = v; };
    void add(const std::string& n, const std::vector<int>&      v) { intVM_[n]    = v; };
    void add(const std::string& n, const std::vector<UChar_t>&  v) { ucharVM_[n]  = v; };
    void add(const std::string& n, const std::vector<UShort_t>& v) { ushortVM_[n] = v; };
    void add(const std::string& n, const std::vector<UInt_t>&   v) { uintVM_[n]   = v; };
    void add(const std::string& n, const std::vector<float>&    v) { floatVM_[n]  = v; };

    template <class T>
    void push(const std::string& n, const T&        v, const bool& c=0) = delete; // avoid implicit conversion
    void push(const std::string& n, const bool&     v, const bool& c=0) { boolVM_[n].push_back(v);  };
    void push(const std::string& n, const Char_t&   v, const bool& c=0) { if(!c || v!=CHAR_MAX ) charVM_[n].push_back(v);   };
    void push(const std::string& n, const short&    v, const bool& c=0) { if(!c || v!=SHRT_MAX ) shortVM_[n].push_back(v);  };
    void push(const std::string& n, const int&      v, const bool& c=0) { if(!c || v!=INT_MAX  ) intVM_[n].push_back(v);    };
    void push(const std::string& n, const UChar_t&  v, const bool& c=0) { if(!c || v!=UCHAR_MAX) ucharVM_[n].push_back(v);  };
    void push(const std::string& n, const UShort_t& v, const bool& c=0) { if(!c || v!=USHRT_MAX) ushortVM_[n].push_back(v); };
    void push(const std::string& n, const UInt_t&   v, const bool& c=0) { if(!c || v!=UINT_MAX ) uintVM_[n].push_back(v);   };
    void push(const std::string& n, const float&    v, const bool& c=0) { floatVM_[n].push_back(v); };

    // clear
    void clear()
    {
      for (auto& d : boolM_   ) { d.second = false; }
      for (auto& d : charM_   ) { d.second = -99;   }
      for (auto& d : shortM_  ) { d.second = -99;   }
      for (auto& d : intM_    ) { d.second = -99;   }
      for (auto& d : ucharM_  ) { d.second = 0;     }
      for (auto& d : ushortM_ ) { d.second = 0;     }
      for (auto& d : uintM_   ) { d.second = 0;     }
      for (auto& d : floatM_  ) { d.second = -99.9; }
      for (auto& d : boolVM_  ) { d.second.clear(); }
      for (auto& d : charVM_  ) { d.second.clear(); }
      for (auto& d : shortVM_ ) { d.second.clear(); }
      for (auto& d : intVM_   ) { d.second.clear(); }
      for (auto& d : ucharVM_ ) { d.second.clear(); }
      for (auto& d : ushortVM_) { d.second.clear(); }
      for (auto& d : uintVM_  ) { d.second.clear(); }
      for (auto& d : floatVM_ ) { d.second.clear(); }
    };

    // tree initializer
    void initTree(TTree& tree, const std::string& n="")
    {
      for (auto& d : boolM_   ) { tree.Branch((n+d.first).c_str(), &d.second, (d.first+"/O").c_str()); }
      for (auto& d : ucharM_  ) { tree.Branch((n+d.first).c_str(), &d.second, (d.first+"/b").c_str()); }
      for (auto& d : charM_   ) { tree.Branch((n+d.first).c_str(), &d.second, (d.first+"/B").c_str()); }
      for (auto& d : ushortM_ ) { tree.Branch((n+d.first).c_str(), &d.second, (d.first+"/s").c_str()); }
      for (auto& d : shortM_  ) { tree.Branch((n+d.first).c_str(), &d.second, (d.first+"/S").c_str()); }
      for (auto& d : uintM_   ) { tree.Branch((n+d.first).c_str(), &d.second, (d.first+"/i").c_str()); }
      for (auto& d : intM_    ) { tree.Branch((n+d.first).c_str(), &d.second, (d.first+"/I").c_str()); }
      for (auto& d : floatM_  ) { tree.Branch((n+d.first).c_str(), &d.second, (d.first+"/F").c_str()); }
      for (auto& d : boolVM_  ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : charVM_  ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : shortVM_ ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : intVM_   ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : ucharVM_ ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : ushortVM_) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : uintVM_  ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : floatVM_ ) { tree.Branch((n+d.first).c_str(), &d.second); }
    };

   private:
    std::map<std::string, bool     > boolM_;
    std::map<std::string, Char_t   > charM_;
    std::map<std::string, short    > shortM_;
    std::map<std::string, int      > intM_;
    std::map<std::string, UChar_t  > ucharM_;
    std::map<std::string, UShort_t > ushortM_;
    std::map<std::string, UInt_t   > uintM_;
    std::map<std::string, float    > floatM_;
    std::map<std::string, std::vector<bool>     > boolVM_;
    std::map<std::string, std::vector<Char_t>   > charVM_;
    std::map<std::string, std::vector<short>    > shortVM_;
    std::map<std::string, std::vector<int>      > intVM_;
    std::map<std::string, std::vector<UChar_t>  > ucharVM_;
    std::map<std::string, std::vector<UShort_t> > ushortVM_;
    std::map<std::string, std::vector<UInt_t>   > uintVM_;
    std::map<std::string, std::vector<float>    > floatVM_;
  };

  class ParticleContainer
  {
   public:
    ParticleContainer() : size_(0) {};
    ~ParticleContainer() {};

    // getters
    bool getIndex(size_t& index, const reco::Candidate& cand) const
    {
      index = size_;
      const auto& parIt = parM_.find(getPar(cand));
      if (parIt==parM_.end()) return false;
      index = parIt->second;
      return true;
    };

    // setters
    template <class T>
    void add(const std::string& n, const T& v) { data_.add(n, v); };

    template <class T>
    void push(const std::string& n, const T& v, const bool& c=0) { data_.push(n, v, c); };

    template <class T>
    void pushV(const size_t& i, const std::string& n, const T&        v, const bool& c=0) = delete; // avoid implicit conversion
    void pushV(const size_t& i, const std::string& n, const UChar_t&  v, const bool& c=0) { if(!c || v!=UCHAR_MAX) ucharVVM_[n][i].push_back(v);  };
    void pushV(const size_t& i, const std::string& n, const UShort_t& v, const bool& c=0) { if(!c || v!=USHRT_MAX) ushortVVM_[n][i].push_back(v); };
    void pushV(const size_t& i, const std::string& n, const UInt_t&   v, const bool& c=0) { if(!c || v!=UINT_MAX ) uintVVM_[n][i].push_back(v);   };
    template <class T>
    void push(const size_t& i, const std::string& n, const T& v, const bool& c=0) { (i < size_ ? pushV(i, n, v, c) : push(n, v, c)); }

    void pushData(const reco::Candidate& cand)
    { 
      const auto par = getPar(cand);
      if (parM_.find(par)!=parM_.end()) { throw std::logic_error("[ERROR] Particle already present"); }
      for (const auto& d : data_.boolM()   ) { boolVM_[d.first].push_back(d.second);    }
      for (const auto& d : data_.charM()   ) { charVM_[d.first].push_back(d.second);    }
      for (const auto& d : data_.shortM()  ) { shortVM_[d.first].push_back(d.second);   }
      for (const auto& d : data_.intM()    ) { intVM_[d.first].push_back(d.second);     }
      for (const auto& d : data_.ucharM()  ) { ucharVM_[d.first].push_back(d.second);   }
      for (const auto& d : data_.ushortM() ) { ushortVM_[d.first].push_back(d.second);  }
      for (const auto& d : data_.uintM()   ) { uintVM_[d.first].push_back(d.second);    }
      for (const auto& d : data_.floatM()  ) { floatVM_[d.first].push_back(d.second);   }
      for (const auto& d : data_.boolVM()  ) { boolVVM_[d.first].push_back(d.second);   }
      for (const auto& d : data_.charVM()  ) { charVVM_[d.first].push_back(d.second);   }
      for (const auto& d : data_.shortVM() ) { shortVVM_[d.first].push_back(d.second);  }
      for (const auto& d : data_.intVM()   ) { intVVM_[d.first].push_back(d.second);    }
      for (const auto& d : data_.ucharVM() ) { ucharVVM_[d.first].push_back(d.second);  }
      for (const auto& d : data_.ushortVM()) { ushortVVM_[d.first].push_back(d.second); }
      for (const auto& d : data_.uintVM()  ) { uintVVM_[d.first].push_back(d.second);   }
      parM_[par] = size_++;
    };

    // clear
    void clear()
    {
      size_ = 0;
      data_.clear();
      parM_.clear();
      for (auto& d : boolVM_   ) { d.second.clear(); }
      for (auto& d : charVM_   ) { d.second.clear(); }
      for (auto& d : shortVM_  ) { d.second.clear(); }
      for (auto& d : intVM_    ) { d.second.clear(); }
      for (auto& d : ucharVM_  ) { d.second.clear(); }
      for (auto& d : ushortVM_ ) { d.second.clear(); }
      for (auto& d : uintVM_   ) { d.second.clear(); }
      for (auto& d : floatVM_  ) { d.second.clear(); }
      for (auto& d : boolVVM_  ) { d.second.clear(); }
      for (auto& d : charVVM_  ) { d.second.clear(); }
      for (auto& d : shortVVM_ ) { d.second.clear(); }
      for (auto& d : intVVM_   ) { d.second.clear(); }
      for (auto& d : ucharVVM_ ) { d.second.clear(); }
      for (auto& d : ushortVVM_) { d.second.clear(); }
      for (auto& d : uintVVM_  ) { d.second.clear(); }
    };

    // tree initializer
    void initTree(TTree& tree, const std::string& n="")
    {
      for (auto& d : boolVM_   ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : charVM_   ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : shortVM_  ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : intVM_    ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : ucharVM_  ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : ushortVM_ ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : uintVM_   ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : floatVM_  ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : boolVVM_  ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : charVVM_  ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : shortVVM_ ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : intVVM_   ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : ucharVVM_ ) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : ushortVVM_) { tree.Branch((n+d.first).c_str(), &d.second); }
      for (auto& d : uintVVM_  ) { tree.Branch((n+d.first).c_str(), &d.second); }
    };

   private:
    typedef std::tuple<float, float, float, float, char, int> Particle;

    Particle getPar(const reco::Candidate& cand) const
    {
      return Particle(cand.pt(), cand.eta(), cand.phi(), cand.mass(), cand.charge(), cand.pdgId());
    };

    size_t size_;
    Container data_;
    std::map<Particle, size_t> parM_;
    std::map<std::string, std::vector<bool>     > boolVM_;
    std::map<std::string, std::vector<Char_t>   > charVM_;
    std::map<std::string, std::vector<short>    > shortVM_;
    std::map<std::string, std::vector<int>      > intVM_;
    std::map<std::string, std::vector<UChar_t>  > ucharVM_;
    std::map<std::string, std::vector<UShort_t> > ushortVM_;
    std::map<std::string, std::vector<UInt_t>   > uintVM_;
    std::map<std::string, std::vector<float>    > floatVM_;
    std::map<std::string, std::vector<std::vector<bool> >     > boolVVM_;
    std::map<std::string, std::vector<std::vector<Char_t> >   > charVVM_;
    std::map<std::string, std::vector<std::vector<short> >    > shortVVM_;
    std::map<std::string, std::vector<std::vector<int> >      > intVVM_;
    std::map<std::string, std::vector<std::vector<UChar_t> >  > ucharVVM_;
    std::map<std::string, std::vector<std::vector<UShort_t> > > ushortVVM_;
    std::map<std::string, std::vector<std::vector<UInt_t> >   > uintVVM_;
  };
  typedef std::map<std::string, ParticleContainer> ParticleContainerMap;

  class TriggerInfo
  {
   public:
    TriggerInfo()
    {
      triggerIndex_ = -1;
      filterIndex_  = -1;
    };
    ~TriggerInfo() {};

    // getters
    int               triggerIndex()  const { return triggerIndex_;  }
    int               filterIndex()   const { return filterIndex_;   }
    int               minN()          const { return minN_;          }
    bool              validPrescale() const { return validPrescale_; }
    UShort_t          hltPrescale()   const { return hltPrescale_;   }
    UShort_t          l1Prescale()    const { return l1Prescale_;    }
    std::string       triggerName()   const { return triggerName_;   }
    std::string       filterName()    const { return filterName_;    }
    std::vector<bool> triggerBit()    const { return triggerBit_;    }
    bool              triggerBit(const size_t& i) const { return (i<triggerBit_.size() ? triggerBit_[i] : false); }
    pat::TriggerObjectCollection filterObjects() const { return filterObjects_; }

    // setters
    void setInfo(const int& triggerIndex, const int& filterIndex,
                 const std::string& triggerName, const std::string& filterName, const int& minN,
                 const bool& validPrescale, const UShort_t& hltPrescale, const UShort_t& l1Prescale,
                 const std::vector<bool>& bit, const pat::TriggerObjectCollection& filterObjects)
    {
      triggerIndex_  = triggerIndex;
      filterIndex_   = filterIndex;
      triggerName_   = triggerName;
      filterName_    = filterName;
      minN_          = minN;
      validPrescale_ = validPrescale;
      hltPrescale_   = hltPrescale;
      l1Prescale_    = l1Prescale;
      triggerBit_    = bit;
      filterObjects_ = filterObjects;
    }

   private:
    int triggerIndex_, filterIndex_, minN_;
    bool validPrescale_;
    UShort_t hltPrescale_, l1Prescale_;
    std::string triggerName_, filterName_;
    std::vector<bool> triggerBit_;
    pat::TriggerObjectCollection filterObjects_;
  };
  typedef std::vector<TriggerInfo> TriggerContainer;

  // class container attributes
  Container eventInfo_;
  TriggerContainer triggerData_;
  ParticleContainerMap particleInfo_;
};

//
// static data member definitions
//

//
// constructors and destructor
//

ParticleAnalyzer::ParticleAnalyzer(const edm::ParameterSet& iConfig) :
  tok_offlineBS_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  tok_offlinePV_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  tok_recParticle_(consumes<pat::GenericParticleCollection>(iConfig.getParameter<edm::InputTag>("recoParticles"))),
  tok_genParticle_(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genParticles", edm::InputTag("genParticles")))),
  tok_genInfo_(consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("genInfo", edm::InputTag("generator")))),
  tok_centBin_(consumes<int>(iConfig.getUntrackedParameter<edm::InputTag>("centralityBin"))),
  tok_centSrc_(consumes<reco::Centrality>(iConfig.getUntrackedParameter<edm::InputTag>("centrality"))),
  tok_eventPlaneSrc_(consumes<reco::EvtPlaneCollection>(iConfig.getUntrackedParameter<edm::InputTag>("eventPlane"))),
  tok_triggerResults_(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("triggerResults", edm::InputTag("TriggerResults::HLT")))),
  tok_triggerEvent_(consumes<trigger::TriggerEvent>(iConfig.getUntrackedParameter<edm::InputTag>("triggerEvent", edm::InputTag("hltTriggerSummaryAOD::HLT")))),
  tok_filterResults_(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("eventFilterResults", edm::InputTag("TriggerResults")))),
  tok_lumiInfo_(consumes<LumiInfo>(iConfig.getUntrackedParameter<edm::InputTag>("lumiInfo"))),
  tok_lumiScalers_(consumes<LumiScalersCollection>(iConfig.getUntrackedParameter<edm::InputTag>("lumiScalers", edm::InputTag("scalersRawToDigi")))),
  tok_lumiRecord_(consumes<OnlineLuminosityRecord>(iConfig.getUntrackedParameter<edm::InputTag>("lumiRecord", edm::InputTag("onlineMetaDataDigis")))),
  triggerInfo_(iConfig.getUntrackedParameter<std::vector<edm::ParameterSet> >("triggerInfo")),
  eventFilters_(iConfig.getUntrackedParameter<std::vector<std::string> >("eventFilterNames")),
  selectEvents_(iConfig.getParameter<std::string>("selectEvents")),
  saveTree_(iConfig.getUntrackedParameter<bool>("saveTree", true)),
  addTrack_(iConfig.getUntrackedParameter<bool>("addTrack", true)),
  addSource_(iConfig.getUntrackedParameter<bool>("addSource", false)),
  maxGenIter_(iConfig.getUntrackedParameter<int>("maxGenIter", 0)),
  maxGenDeltaR_(iConfig.getUntrackedParameter<double>("maxGenDeltaR", 0.03)),
  maxGenDeltaPtRel_(iConfig.getUntrackedParameter<double>("maxGenDeltaPtRel", 0.5)),
  maxTrgDeltaR_(iConfig.getUntrackedParameter<double>("maxTrgDeltaR", 0.3)),
  maxTrgDeltaPtRel_(iConfig.getUntrackedParameter<double>("maxTrgDeltaPtRel", 10.)),
  dauIDs_(iConfig.getUntrackedParameter<std::vector<int> >("dauIDs", {})),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
}


ParticleAnalyzer::~ParticleAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //check event
  if (selectEvents_!="")
  {
    edm::Handle<edm::TriggerResults> filterResults;
    iEvent.getByToken(tok_filterResults_, filterResults);
    const auto& filterNames = iEvent.triggerNames(*filterResults);
    const auto& index = filterNames.triggerIndex(selectEvents_);
    if (index < filterNames.size() && filterResults->wasrun(index) && !filterResults->accept(index)) return;
  }

  // clear information
  vertex_ = reco::Vertex();
  vertices_.clear();
  genParticles_.clear();
  eventInfo_.clear();
  triggerData_.clear();
  for (auto& p : particleInfo_) { p.second.clear(); }

  // get event data
  getEventData(iEvent, iSetup);
  getTriggerData(iEvent, iSetup);

  // fill information
  fillEventInfo(iEvent);
  fillTriggerInfo(iEvent, iSetup);
  fillLumiInfo(iEvent);
  fillRecoParticleInfo(iEvent);
  fillGenParticleInfo(iEvent);

  // fill tree
  if (saveTree_)
  {
    if (!tree_) initTree();
    tree_->Fill();
  }
}


void
ParticleAnalyzer::getEventData(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // general information
  isMC_ = !iEvent.isRealData();

  // magnetic field information
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
  magField_ = bFieldHandle.product();

  // primary vertex information
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_offlinePV_, vertices);
  for (const auto& pv : *vertices)
  {
    if (!pv.isFake() && pv.tracksSize() >= 2)
    {
      vertices_.push_back(pv);
    }
  }
  auto byTracksSize = [] (const reco::Vertex& v1, const reco::Vertex& v2) -> bool { return v1.tracksSize() > v2.tracksSize(); };
  std::sort(vertices_.begin(), vertices_.end(), byTracksSize);
  if (vertices_.empty())
  {
    edm::Handle<reco::BeamSpot> beamspot;
    iEvent.getByToken(tok_offlineBS_, beamspot);
    vertex_ = reco::Vertex(beamspot->position(), beamspot->rotatedCovariance3D());
  }
  else { vertex_ = vertices_[0]; }

  // generated particle information
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(tok_genParticle_, genParticles);
  if (isMC_ && genParticles.isValid())
  {
    // initialize generated particle container
    initParticleInfo("gen");
    // extract generated particles
    for (size_t i=0; i<genParticles->size(); i++)
    {
      const auto& p = reco::GenParticleRef(genParticles, i); 
      if (p->isLastCopy() && (p->isPromptFinalState() || p->isPromptDecayed() || p->statusFlags().fromHardProcess()))
      {
        genParticles_.push_back(p);
      }
    }
  }
}


void
ParticleAnalyzer::getTriggerData(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (triggerInfo_.empty()) return;

  // get trigger data
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  iEvent.getByToken(tok_triggerEvent_, triggerEvent);
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(tok_triggerResults_, triggerResults);
  if (triggerEvent.isValid() && triggerResults.isValid())
  {
    // initialize trigger container
    triggerData_.resize(triggerInfo_.size());
    // retrieve trigger information
    auto& l1tGlobalUtil = *const_cast<l1t::L1TGlobalUtil*>(&hltPrescaleProvider_.l1tGlobalUtil());
    const auto& hltConfig = hltPrescaleProvider_.hltConfigProvider();
    l1tGlobalUtil.retrieveL1Event(iEvent, iSetup);
    const bool validPrescale = (hltConfig.inited() && hltPrescaleProvider_.prescaleSet(iEvent, iSetup) >= 0);
    const auto& hltPaths = hltConfig.triggerNames();
    const auto& triggerObjects = triggerEvent->getObjects();
    // extract trigger information
    for (size_t iTrg=0; iTrg<triggerInfo_.size(); iTrg++)
    {
      const auto& pSet = triggerInfo_[iTrg];
      const auto& minN = (pSet.existsAs<std::string>("minN") ? pSet.getParameter<int>("minN") : 0);
      const auto& isL1OR = (pSet.existsAs<bool>("isL1OR") ? pSet.getParameter<int>("isL1OR") : false);
      const auto& pathLabel = (pSet.existsAs<std::string>("path") ? pSet.getParameter<std::string>("path") : std::string());
      const auto& filterLabel = (pSet.existsAs<std::string>("filter") ? pSet.getParameter<std::string>("filter") : std::string());
      // initialize trigger data
      int triggerIndex=-1, filterIndex=-1;
      // extract the trigger index
      std::vector<size_t> trgIdxFound;
      if (pathLabel!="")
      {
        for (size_t trgIdx=0; trgIdx<hltPaths.size(); trgIdx++)
        {
          const auto& hltPath = hltPaths.at(trgIdx);
          if (TString(hltPath).Contains(TRegexp(TString(pathLabel))) && triggerResults->wasrun(trgIdx)) { trgIdxFound.push_back(trgIdx); }
        }
        if (trgIdxFound.empty()) continue;
      }
      // extract the filter index
      std::vector<size_t> filterIdxFound;
      if (filterLabel!="")
      {
        for (size_t filterIdx=0; filterIdx<triggerEvent->sizeFilters(); filterIdx++)
        {
          const auto& filterName = triggerEvent->filterLabel(filterIdx);
          if (TString(filterName).Contains(TRegexp(TString(filterLabel)))) { filterIdxFound.push_back(filterIdx); }
        }
      }
      // match the trigger and filter
      if (pathLabel!="" && !filterIdxFound.empty())
      {
        std::vector<std::pair<size_t, size_t> > trgFilterIdxFound;
        for (const auto& trgIdx : trgIdxFound)
        {
          for (const auto& filterIdx : filterIdxFound)
          {
            const auto& filterName = triggerEvent->filterLabel(filterIdx);
            if (hltConfig.moduleIndex(trgIdx, filterName)!=hltConfig.size(trgIdx)) { trgFilterIdxFound.push_back({trgIdx, filterIdx}); }
          }
        }
        if (trgFilterIdxFound.empty()) continue;
        // determine the trigger and filter index
        auto triggerFilterIndex = trgFilterIdxFound[0];
        for (const auto& trgFilterIdx : trgFilterIdxFound) { if (triggerResults->accept(trgFilterIdx.first)) { triggerFilterIndex = trgFilterIdx; break; } }
        triggerIndex = triggerFilterIndex.first;
        filterIndex = triggerFilterIndex.second;
      }
      // if no filter provided, select the last filter module ran
      else if (pathLabel!="" && filterIdxFound.empty())
      {
        // determine the trigger index
        triggerIndex = trgIdxFound[0];
        for (const auto& trgIdx : trgIdxFound) { if (triggerResults->accept(trgIdx)) { triggerIndex = trgIdx; break; } }
        // find last filter module
        if (filterLabel=="" && triggerResults->accept(triggerIndex))
        {
          for (int j = hltConfig.moduleLabels(triggerIndex).size()-1; j >= 0; j--)
          {
            const auto& filter = hltConfig.moduleLabels(triggerIndex).at(j);
            const auto& filterIdx = triggerEvent->filterIndex(edm::InputTag(filter, "", hltConfig.processName()));
            if (filterIdx < triggerEvent->sizeFilters()) { filterIndex = filterIdx; break; }
          }
          if (filterIndex<0) { throw(cms::Exception("Trigger")<<"Invalid filter index for: "<<pathLabel<<" , "<<filterLabel); }
        }
      }
      // if no path provided, select the first filter module found
      else if (pathLabel=="" && !filterIdxFound.empty())
      {
        filterIndex = filterIdxFound[0];
      }
      else if (pathLabel=="" && filterIdxFound.empty()) continue;
      // determine the trigger/filter name
      const auto triggerName = (triggerIndex>=0 ? hltPaths.at(triggerIndex) : std::string());
      const auto filterName  = (filterIndex>=0  ? triggerEvent->filterLabel(filterIndex) : std::string());
      // determine the trigger decision
      std::vector<bool> bit(4);
      bit[0] = (triggerIndex>=0 ? triggerResults->accept(triggerIndex) : false);
      // extract prescale information
      UShort_t hltPrescale=0, l1Prescale=0;
      if (validPrescale && triggerIndex>=0)
      {
        // HLT information
        const auto& presInfo = hltPrescaleProvider_.prescaleValuesInDetail(iEvent, iSetup, triggerName);
        hltPrescale = getUShort(presInfo.second);
        bit[2] = !hltPrescaleProvider_.rejectedByHLTPrescaler(*triggerResults, triggerIndex);
        // L1 information
        if (!presInfo.first.empty())
        {
          // L1 prescale
          int l1Pres = (isL1OR ? 1E9 : -1E9);
          for (const auto p : presInfo.first)
          {
            int pres; if (!l1tGlobalUtil.getPrescaleByName(p.first, pres)) continue;
            if (pres==0) { pres = 1E9; }
            l1Pres = (isL1OR ? std::min(l1Pres, pres) : std::max(l1Pres, pres));
          }
          if (l1Pres==1E9) { l1Pres = 0; }
          l1Prescale = getUShort(l1Pres);
          // L1 decision
          for (const auto p : presInfo.first)
          {
            int l1Bit; if (!l1tGlobalUtil.getAlgBitFromName(p.first, l1Bit)) continue;
            bit[1] = l1tGlobalUtil.decisionsInitial()[l1Bit].second;
            bit[3] = (l1tGlobalUtil.decisionsInitial()[l1Bit].second == l1tGlobalUtil.decisionsFinal()[l1Bit].second);
            if (isL1OR==bit[1]) break;
          }
        }
      }
      // extract filter objects
      pat::TriggerObjectCollection filterObjects;
      const auto& filterKeys = (filterIndex>=0 ? triggerEvent->filterKeys(filterIndex) : std::vector<UShort_t>());
      const auto& filterIds = (filterIndex>=0 ? triggerEvent->filterIds(filterIndex) : std::vector<int>());
      for (size_t iKey=0; iKey<filterKeys.size(); iKey++)
      {
        pat::TriggerObject obj(triggerObjects[filterKeys[iKey]]);
        obj.addFilterId(filterIds[iKey]);
        filterObjects.push_back(obj);
      }
      // store trigger information
      triggerData_[iTrg].setInfo(triggerIndex, filterIndex, triggerName, filterName, minN, validPrescale, hltPrescale, l1Prescale, bit, filterObjects);
    }
  }
}


void
ParticleAnalyzer::fillEventInfo(const edm::Event& iEvent)
{
  // fill general information
  eventInfo_.add("RunNb", iEvent.id().run());
  eventInfo_.add("EventNb", getUInt(iEvent.id().event()));
  eventInfo_.add("LSNb", getUShort(iEvent.luminosityBlock()));
  eventInfo_.add("BXNb", getUShort(iEvent.bunchCrossing()));

  // fill vertex information
  eventInfo_.add("nPV", getUChar(vertices_.size()));
  eventInfo_.add("bestvtxX", vertex_.x());
  eventInfo_.add("bestvtxY", vertex_.y());
  eventInfo_.add("bestvtxZ", vertex_.z());

  // fill event selection information
  edm::Handle<edm::TriggerResults> filterResults;
  iEvent.getByToken(tok_filterResults_, filterResults);
  if (filterResults.isValid() && !eventFilters_.empty())
  {
    std::vector<bool> evtSel;
    const auto& filterNames = iEvent.triggerNames(*filterResults);
    for (const auto& eventFilter : eventFilters_)
    {
      const auto& index = filterNames.triggerIndex(eventFilter);
      evtSel.push_back(index < filterNames.size() && filterResults->wasrun(index) && filterResults->accept(index));
    }
    eventInfo_.add("evtSel", evtSel);
  }

  // fill centrality information
  edm::Handle<reco::Centrality> cent;
  iEvent.getByToken(tok_centSrc_, cent);
  if (cent.isValid())
  {
    eventInfo_.add("HFsumETPlus", cent->EtHFtowerSumPlus());
    eventInfo_.add("HFsumETMinus", cent->EtHFtowerSumMinus());
    eventInfo_.add("Npixel", cent->multiplicityPixel());
    eventInfo_.add("ZDCPlus", cent->zdcSumPlus());
    eventInfo_.add("ZDCMinus", cent->zdcSumMinus());
    eventInfo_.add("Ntrkoffline", cent->Ntracks());
  }
  edm::Handle<int> centBin;
  iEvent.getByToken(tok_centBin_, centBin);
  if (centBin.isValid())
  {
    eventInfo_.add("centrality", getUChar(*centBin));
  }

  // fill event plane information
  edm::Handle<reco::EvtPlaneCollection> eventPlanes;
  iEvent.getByToken(tok_eventPlaneSrc_, eventPlanes);
  if (eventPlanes.isValid())
  {
    for (const auto& iEP : {0, 6, 13})
    {
      eventInfo_.push("ephfmAngle",      eventPlanes->at(iEP  ).angle(2));
      eventInfo_.push("ephfpAngle",      eventPlanes->at(iEP+1).angle(2));
      eventInfo_.push("eptrackmidAngle", eventPlanes->at(iEP+3).angle(2));
      eventInfo_.push("ephfmQ",          eventPlanes->at(iEP  ).q(2));
      eventInfo_.push("ephfpQ",          eventPlanes->at(iEP+1).q(2));
      eventInfo_.push("eptrackmidQ",     eventPlanes->at(iEP+3).q(2));
    }
    eventInfo_.add("ephfmSumW",      eventPlanes->at(6).sumw());
    eventInfo_.add("ephfpSumW",      eventPlanes->at(7).sumw());
    eventInfo_.add("eptrackmidSumW", eventPlanes->at(9).sumw());
  }
}


void
ParticleAnalyzer::fillLumiInfo(const edm::Event& iEvent)
{
  // offline luminosity information
  edm::Handle<LumiInfo> lumiInfo;
  iEvent.getByToken(tok_lumiInfo_, lumiInfo);
  if (lumiInfo.isValid())
  {
    eventInfo_.add("totalLumi", lumiInfo->integLuminosity());
    eventInfo_.add("recordLumi", lumiInfo->recordedLuminosity());
  }

  // online luminosity information
  edm::Handle<LumiScalersCollection> lumiScalers;
  iEvent.getByToken(tok_lumiScalers_, lumiScalers);
  edm::Handle<OnlineLuminosityRecord> lumiRecord;
  iEvent.getByToken(tok_lumiRecord_, lumiRecord);
  // use SCAL data
  if (lumiScalers.isValid() && !lumiScalers->empty())
  {
    eventInfo_.add("pileup", lumiScalers->begin()->pileup());
    eventInfo_.add("rawInstLumi", lumiScalers->begin()->instantLumi());
  }
  // otherwise, use online record
  else if (lumiRecord.isValid())
  {
    eventInfo_.add("pileup", lumiRecord->avgPileUp());
    eventInfo_.add("rawInstLumi", lumiRecord->instLumi());
  }
}


void
ParticleAnalyzer::fillTriggerInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // initialize trigger object container
  initParticleInfo("trig");

  // fill trigger information
  for (UShort_t idx=0; idx<triggerData_.size(); idx++)
  {
    const auto& data = triggerData_[idx];
    eventInfo_.push("passHLT", data.triggerBit(0));
    eventInfo_.push("passL1", data.triggerBit(1));
    eventInfo_.push("validPrescale", data.validPrescale());
    eventInfo_.push("passHLTPrescaler", data.triggerBit(2));
    eventInfo_.push("passL1Prescaler", data.triggerBit(3));
    eventInfo_.push("hltPrescale", data.hltPrescale());
    eventInfo_.push("l1Prescale", data.l1Prescale());
    for (const auto& obj: data.filterObjects()) { fillTriggerObjectInfo(obj, idx); }
  }
}


void
ParticleAnalyzer::initParticleInfo(const std::string& type, const int& pId)
{
  // return if already initialized
  if (particleInfo_.find(type)!=particleInfo_.end()) return;
  // check type
  if      (type=="trig" && triggerData_.empty()) return;
  else if (type=="trk"  && !addTrack_) return;
  else if (type=="gen"  && !isMC_) return;
  else if (type=="src"  && !addSource_) return;
  // proceed to initialize with dummy value
  if      (type=="cand") fillRecoParticleInfo(pat::GenericParticle(), 0);
  else if (type=="trig") fillTriggerObjectInfo(pat::TriggerObject(), 0, 0);
  else if (type=="gen" ) fillGenParticleInfo(reco::GenParticleRef(), 0, true);
  else if (type=="trk" ) fillTrackInfo(pat::GenericParticle(), 0, true);
  else if (type=="src" ) fillSourceInfo(pat::GenericParticle(), pId, 0, true);
  // clear the initialized info
  if (type!="src") { particleInfo_[type].clear(); }
  else
  {
    for (auto& p : particleInfo_)
    {
      if (p.first!="cand" && p.first!="trig" && p.first!="gen" && p.first!="trk") { p.second.clear(); }
    }
  }
}


UShort_t
ParticleAnalyzer::fillTriggerObjectInfo(const pat::TriggerObject& obj, const UShort_t& trigIdx, const UShort_t& candIdx)
{
  if (triggerData_.empty()) return USHRT_MAX;

  // fill trigger information
  auto& info = particleInfo_["trig"];

  // add input information
  size_t index;
  const bool found = info.getIndex(index, obj);
  info.push(index, "trigIdx", trigIdx, true);
  info.push(index, "candIdx", candIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;

  // basic information
  info.add("pT", obj.pt());
  info.add("eta", obj.eta());
  info.add("phi", obj.phi());
  info.add("mass", obj.mass());
  info.add("pdgId", obj.pdgId());
  const auto filterId = (obj.filterIds().empty() ? 0 :obj.filterIds()[0]);
  info.add("filterId", getChar(filterId));

  // push data and return index
  info.pushData(obj);
  return idx;
}


void
ParticleAnalyzer::fillRecoParticleInfo(const edm::Event& iEvent)
{
  // fill reconstructed particle information
  edm::Handle<pat::GenericParticleCollection> particles;
  iEvent.getByToken(tok_recParticle_, particles);
  if (particles.isValid())
  {
    // initialize reconstructed particle containers
    initParticleInfo("cand");
    initParticleInfo("trk");
    for (const auto& pId : dauIDs_) { initParticleInfo("src", pId); }
    // loop over reconstructed particles
    for (const auto& cand : *particles) { fillRecoParticleInfo(cand); }
  }
}


UShort_t
ParticleAnalyzer::fillRecoParticleInfo(const pat::GenericParticle& cand, const UShort_t& momIdx)
{
  // fill reconstructed particle information
  auto& info = particleInfo_["cand"];

  // add input information
  size_t index;
  const bool found = info.getIndex(index, cand);
  info.push(index, "momIdx", momIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;

  // basic information
  info.add("p", cand.p());
  info.add("y", cand.rapidity());
  info.add("pT", cand.pt());
  info.add("eta", cand.eta());
  info.add("phi", cand.phi());
  info.add("mass", cand.mass());
  info.add("charge", getChar(cand.charge()));
  info.add("pdgId", cand.pdgId());
  info.add("status", getUChar(cand.status()));

  // decay vertex information
  info.add("vtxChi2", getFloat(cand, "normChi2", -1.));
  info.add("vtxProb", getFloat(cand, "vertexProb", -1.));

  // decay angle information
  info.add("angle3D", getFloat(cand, "angle2D", -10.));
  info.add("angle2D", getFloat(cand, "angle3D", -10.));

  // decay length information
  const auto& lVtxMag = getFloat(cand, "lVtxMag", -99.9);
  const auto& lVtxSig = getFloat(cand, "lVtxSig", -99.9);
  info.add("decayLength3D", lVtxMag);
  info.add("decayLengthError3D", ((lVtxMag==-99.9 || lVtxSig==-99.9) ? -1. : lVtxMag/lVtxSig));
  info.add("pseudoDecayLengthError3D", getFloat(cand, "pdlErr3D"));
  const auto& rVtxMag = getFloat(cand, "rVtxMag", -99.9);
  const auto& rVtxSig = getFloat(cand, "rVtxSig", -99.9);
  info.add("decayLength2D", rVtxMag);
  info.add("decayLengthError2D", ((rVtxMag==-99.9 || rVtxSig==-99.9) ? -1. : rVtxMag/rVtxSig));
  info.add("pseudoDecayLength2D", getFloat(cand, "pdlErr2D"));

  // trigger information
  if (!triggerData_.empty())
  {
    std::vector<UShort_t> trigIdx;
    addTriggerObject(*const_cast<pat::GenericParticle*>(&cand));
    for (UShort_t i=0; i<triggerData_.size(); i++)
    {
      std::vector<UShort_t> trigObjIdx;
      const auto& triggerObjects = cand.triggerObjectMatchesByFilter(triggerData_[i].filterName());
      info.add(Form("matchTRG%d",i), !triggerObjects.empty());
      if (cand.status()==1)
      {
        for (const auto& obj : triggerObjects)
        {
          trigIdx.push_back(fillTriggerObjectInfo(obj, i, idx));
        }
      }
    }
    info.add("trigIdx", trigIdx);
  }

  // track information
  if (addTrack_) info.add("trkIdx", fillTrackInfo(cand, idx));

  // source information
  if (addSource_) info.add("srcIdx", fillSourceInfo(cand, cand.pdgId(), idx));

  // generated particle information
  if (isMC_)
  {
    if (!cand.hasUserInt("isGenMatched")) { addGenParticle(*const_cast<pat::GenericParticle*>(&cand), genParticles_); } 
    const auto& genPar = cand.genParticleRef();
    info.add("genIdx", fillGenParticleInfo(genPar, idx));
    info.add("matchGEN", genPar.isNonnull());
    info.add("isSwap", (cand.hasUserInt("isSwap") ? cand.userInt("isSwap") : false));
    info.add("idmom_reco", (genPar.isNonnull() ? genPar->pdgId() : 0));
  }

  // initialize daughter information
  info.add("dauIdx", std::vector<UShort_t>());

  // push data
  info.pushData(cand);

  // add daughter information
  if (cand.hasUserData("daughters"))
  {
    const auto& dauColl = *cand.userData<pat::GenericParticleCollection>("daughters");
    for(const auto& dau : dauColl)
    {
      info.push(idx, "dauIdx", fillRecoParticleInfo(dau, idx));
    }
  }

  // return index
  return idx;
}


void
ParticleAnalyzer::addTriggerObject(pat::GenericParticle& cand)
{
  if (triggerData_.empty()) return;

  // check if already matched
  if (cand.hasUserInt("isTrigMatched") || !cand.triggerObjectMatches().empty()) return;

  // initialize information
  cand.addUserInt("isTrigMatched", 1);

  // loop over trigger data
  for (const auto& data : triggerData_)
  {
    addTriggerObject(cand, data.filterObjects(), data.triggerName(), data.filterName(), data.minN());
  }
}


bool
ParticleAnalyzer::addTriggerObject(pat::GenericParticle& cand, const pat::TriggerObjectCollection& filterObjects,
                                   const std::string& triggerName, const std::string& filterName, const int& minN)
{
  if (filterObjects.empty() || minN<=0) return false;

  // initialize trigger objects
  pat::TriggerObjectStandAloneCollection triggerObjects;

  // if already matched, copy them
  if (cand.hasUserData("src") && cand.userData<pat::Muon>("src"))
  {
    for (const auto& obj : cand.userData<pat::Muon>("src")->triggerObjectMatches())
    {
      if (obj.hasFilterLabel(filterName)) { triggerObjects.push_back(obj); }
    }
  }

  // do trigger-reco matching
  if (triggerObjects.empty())
  {
    // case: final state particle (direct matching)
    if (cand.status()==1)
    {
      const bool isL3Mu = (triggerName.find("L3Mu")!=std::string::npos);
      const auto& maxDeltaR = (isL3Mu ? 0.1 : maxTrgDeltaR_);
      for (const auto& obj : filterObjects)
      {
        if (isMatched(obj.p4(), cand.p4(), maxDeltaR, maxTrgDeltaPtRel_))
        {
          triggerObjects.push_back(pat::TriggerObjectStandAlone(obj));
        }
      }
    }
    // case: decayed particle (daughter matching)
    else if (cand.status()>1)
    {
      int n=0;
      const auto& dauColl = *cand.userData<pat::GenericParticleCollection>("daughters");
      for (const auto& d: dauColl)
      {
        auto& dau = *const_cast<pat::GenericParticle*>(&d);
        if (addTriggerObject(dau, filterObjects, triggerName, filterName, minN)) { n+=1; }
        if (n==minN) { triggerObjects.push_back(pat::TriggerObjectStandAlone()); break; }
      }
    }
    for (auto& trigObj : triggerObjects)
    {
      trigObj.addPathName(triggerName);
      trigObj.addFilterLabel(filterName);
    }
  }
  if (triggerObjects.empty()) return false;

  // add trigger objects 
  for (const auto& trigObj : triggerObjects)
  {
    cand.addTriggerObjectMatch(trigObj);
  }
  return true;
}


UShort_t
ParticleAnalyzer::fillTrackInfo(const pat::GenericParticle& cand, const UShort_t& candIdx, const bool& force)
{
  const bool hasTrack = cand.track().isNonnull();
  if (!force && !hasTrack) return USHRT_MAX;

  // fill track information
  auto& info = particleInfo_["trk"];

  // add input information
  size_t index;
  const bool found = info.getIndex(index, cand);
  info.push(index, "candIdx", candIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;

  const auto& track = (hasTrack ? *cand.track() : reco::Track());

  // general information
  info.add("isHP", track.quality(reco::TrackBase::highPurity));
  info.add("nChi2", track.normalizedChi2());
  info.add("pTErr", track.ptError());
  info.add("nHit", track.numberOfValidHits());

  // dca information
  info.add("zDCASignificance", getFloat(cand, "dzSig"));
  info.add("xyDCASignificance", getFloat(cand, "dxySig"));

  // dEdx information
  info.add("dEdxHarmonic2", getFloat(cand, "dEdx"));

  // push data and return index
  info.pushData(cand);
  return idx;
}


UShort_t
ParticleAnalyzer::fillSourceInfo(const pat::GenericParticle& cand, const int& pId, const UShort_t& candIdx, const bool& force)
{
  // fill source information
  const auto pdgId = std::abs(pId);
  if      (pdgId<=6 ) { return fillJetInfo(cand, candIdx, force);      }
  else if (pdgId==11) { return fillElectronInfo(cand, candIdx, force); }
  else if (pdgId==13) { return fillMuonInfo(cand, candIdx, force);     }
  else if (pdgId==15) { return fillTauInfo(cand, candIdx, force);      }
  else if (pdgId==22) { return fillPhotonInfo(cand, candIdx, force);   }
  // if pdgId not matched, use PF candidates
  return fillPFCandidateInfo(cand, candIdx, force);
}


UShort_t
ParticleAnalyzer::fillMuonInfo(const pat::GenericParticle& cand, const UShort_t& candIdx, const bool& force)
{
  const bool hasSrc = (std::abs(cand.pdgId())==13 && cand.hasUserData("src") && cand.userData<pat::Muon>("src"));
  if (!force && !hasSrc) return USHRT_MAX;

  // fill muon information
  auto& info = particleInfo_["muon"];

  // add input information
  size_t index;
  const bool found = info.getIndex(index, cand);
  info.push(index, "candIdx", candIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;

  const auto& vertex = (cand.hasUserData("primaryVertex") ? *cand.userData<reco::Vertex>("primaryVertex") : vertex_).position();
  const auto& muon   = (hasSrc ? *cand.userData<pat::Muon>("src") : pat::Muon());
  const auto& gTrack = muon.globalTrack();
  const auto& iTrack = muon.innerTrack();
  const auto& mTrack = muon.muonBestTrack();

  // general information
  UInt_t selector = 0;
  for (size_t i=0; i<27; i++) { selector += (muon.passed(i) ? 1U<<i : 0); }
  info.add("selector", selector);

  UInt_t selectionType = 0;
  for (size_t i=0; i<30; i++) { selectionType += (muon::isGoodMuon(muon, muon::SelectionType(i)) ? 1U<<i : 0); }
  info.add("selectionType", selectionType);

  // tight ID muon POG Run 2
  info.add("isGlobal", muon.isGlobalMuon());
  info.add("isPF", muon.isPFMuon());
  info.add("glbTrkNChi2", (gTrack.isNonnull() ? gTrack->normalizedChi2() : 99.9));
  info.add("nMuonHit", getChar(gTrack.isNonnull() ? gTrack->hitPattern().numberOfValidMuonHits() : -1));
  info.add("nMatchedStation", getUChar(muon.numberOfMatchedStations()));
  info.add("nPixelHit", getChar(iTrack.isNonnull() ? iTrack->hitPattern().numberOfValidPixelHits() : -1));
  info.add("nTrackerLayer", getChar(iTrack.isNonnull() ? iTrack->hitPattern().trackerLayersWithMeasurement() : -1));
  info.add("bestTrkdXY", (mTrack.isNonnull() ? mTrack->dxy(vertex) : -99.9));
  info.add("bestTrkdZ", (mTrack.isNonnull() ? mTrack->dz(vertex) : -99.9));
  const auto tightmuon = (
                          muon.isGlobalMuon() &&
                          muon.isPFMuon() &&
                          gTrack.isNonnull() &&
                          gTrack->normalizedChi2() < 10. &&
                          gTrack->hitPattern().numberOfValidMuonHits() > 0 &&
                          muon.numberOfMatchedStations() > 1 &&
                          iTrack.isNonnull() &&
                          iTrack->hitPattern().numberOfValidPixelHits() > 0 &&
                          iTrack->hitPattern().trackerLayersWithMeasurement() > 5 &&
                          mTrack.isNonnull() &&
                          std::abs(mTrack->dxy(vertex)) < 0.2 &&
                          std::abs(mTrack->dz(vertex)) < 0.5
                         );
  info.add("tightID", tightmuon);

  // soft ID muon POG Run 2
  info.add("isOneStTight", muon::isGoodMuon(muon, muon::SelectionType::TMOneStationTight));
  info.add("nPixelLayer", getChar(iTrack.isNonnull() ? iTrack->hitPattern().pixelLayersWithMeasurement() : -1));
  info.add("isHP", (iTrack.isNonnull() ? iTrack->quality(reco::TrackBase::highPurity) : false));
  info.add("dXY", (iTrack.isNonnull() ? iTrack->dxy(vertex) : -99.9));
  info.add("dZ", (iTrack.isNonnull() ? iTrack->dz(vertex) : -99.9));
  const auto softmuon = (
                          muon::isGoodMuon(muon, muon::SelectionType::TMOneStationTight) &&
                          iTrack.isNonnull() &&
                          iTrack->hitPattern().trackerLayersWithMeasurement() > 5 &&
                          iTrack->hitPattern().pixelLayersWithMeasurement() > 0 &&
                          iTrack->quality(reco::TrackBase::highPurity) &&
                          std::abs(iTrack->dxy(vertex)) < 0.3 &&
                          std::abs(iTrack->dz(vertex)) < 20.
                        );
  info.add("softID", softmuon);

  // hybrid soft ID HIN PAG Run 2 PbPb
  info.add("isTracker", muon.isTrackerMuon());
  const auto hybridmuon = (
                           muon.isGlobalMuon() &&
                           iTrack.isNonnull() &&
                           iTrack->hitPattern().trackerLayersWithMeasurement() > 5 &&
                           iTrack->hitPattern().pixelLayersWithMeasurement() > 0 &&
                           std::abs(iTrack->dxy(vertex)) < 0.3 &&
                           std::abs(iTrack->dz(vertex)) < 20.
                          );
  info.add("hybridSoftID", hybridmuon);    
                                          
  // muon extra information
  info.add("nMatches", getUChar(muon.numberOfMatches()));
  info.add("hadMax", muon.calEnergy().hadMax);

  float d_seg = 1.E9; 
  float dx_seg = -99.9, dy_seg = -99.9, dxSig_seg = -99.9, dySig_seg = -99.9;
  float ddxdz_seg = -99.9, ddydz_seg = -99.9, ddxdzSig_seg = -99.9, ddydzSig_seg = -99.9;
  for (const auto& exp : muon.matches())
  {
    for (const auto& seg : exp.segmentMatches)
    {
      const auto dseg = std::sqrt((seg.x-exp.x)*(seg.x-exp.x) + (seg.y-exp.y)*(seg.y-exp.y));
      if(dseg < d_seg)
      {
        d_seg = dseg;
        dx_seg = seg.x - exp.x;
        dy_seg = seg.y - exp.y;
        dxSig_seg = dx_seg / std::sqrt(seg.xErr*seg.xErr + exp.xErr*exp.xErr);
        dySig_seg = dy_seg / std::sqrt(seg.yErr*seg.yErr + exp.yErr*exp.yErr);
        ddxdz_seg = seg.dXdZ - exp.dXdZ;
        ddydz_seg = seg.dYdZ - exp.dYdZ;
        ddxdzSig_seg = ddxdz_seg / std::sqrt(seg.dXdZErr*seg.dXdZErr + exp.dXdZErr*exp.dXdZErr);
        ddydzSig_seg = ddydz_seg / std::sqrt(seg.dYdZErr*seg.dYdZErr + exp.dYdZErr*exp.dYdZErr);
      }
    }
  }
  info.add("dX_seg", dx_seg);
  info.add("dY_seg", dy_seg);
  info.add("dXSig_seg", dxSig_seg);
  info.add("dYSig_seg", dySig_seg);
  info.add("ddXdZ_seg", ddxdz_seg);
  info.add("ddYdZ_seg", ddydz_seg);
  info.add("ddXdZSig_seg", ddxdzSig_seg);
  info.add("ddYdZSig_seg", ddydzSig_seg);

  // push data and return index
  info.pushData(cand);
  return idx;
};


UShort_t
ParticleAnalyzer::fillElectronInfo(const pat::GenericParticle& cand, const UShort_t& candIdx, const bool& force)
{
  const bool hasSrc = (std::abs(cand.pdgId())==13 && cand.hasUserData("src") && cand.userData<pat::Electron>("src"));
  if (!force && !hasSrc) return USHRT_MAX;

  // fill electron information
  auto& info = particleInfo_["elec"];

  // add input information
  size_t index;
  const bool found = info.getIndex(index, cand);
  info.push(index, "candIdx", candIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;

  const auto& vertex = (cand.hasUserData("primaryVertex") ? *cand.userData<reco::Vertex>("primaryVertex") : vertex_);
  const auto& electron = (hasSrc ? *cand.userData<pat::Electron>("src") : pat::Electron());
  const auto& sCluster = (electron.core().isNonnull() ? electron.superCluster() : reco::SuperClusterRef());
  const auto& gsfTrack = (electron.core().isNonnull() ? electron.gsfTrack() : reco::GsfTrackRef());

  // electron information
  info.add("sigmaIEtaIEta", electron.full5x5_sigmaIetaIeta());
  info.add("sigmaIPhiIPhi", electron.full5x5_sigmaIphiIphi());
  info.add("dEtaSeedAtVtx", (sCluster.isNonnull() ? electron.deltaEtaSeedClusterTrackAtVtx() : -99.9));
  info.add("dEtaAtVtx", electron.deltaEtaSuperClusterTrackAtVtx());
  info.add("dPhiAtVtx", electron.deltaPhiSuperClusterTrackAtVtx());
  info.add("eOverPInv", (1./electron.ecalEnergy() - 1./electron.trackMomentumAtVtx().R()));
  info.add("hOverEBc", electron.hcalOverEcalBc());

  // super cluster information
  info.add("sCEta", (sCluster.isNonnull() ? sCluster->eta() : -99.9));
  info.add("sCPhi", (sCluster.isNonnull() ? sCluster->phi() : -99.9));

  // gsf track information
  const auto ip3D = (gsfTrack.isNonnull() ? IPTools::absoluteImpactParameter3D(reco::TransientTrack(*gsfTrack, magField_), vertex).second.value() : -99.9);
  info.add("ip3D", ip3D);
  info.add("nLostHits", getChar(gsfTrack.isNonnull() ? gsfTrack->numberOfLostHits() : -1));

  // push data and return index
  info.pushData(cand);
  return idx;
};


UShort_t
ParticleAnalyzer::fillPhotonInfo(const pat::GenericParticle& cand, const UShort_t& candIdx, const bool& force)
{ 
  const bool hasSrc = (std::abs(cand.pdgId())==22 && cand.hasUserData("src") && cand.userData<pat::Photon>("src"));
  if (!force && !hasSrc) return USHRT_MAX;
  
  // fill photon information
  auto& info = particleInfo_["pho"];

  // add input information
  size_t index;
  const bool found = info.getIndex(index, cand);
  info.push(index, "candIdx", candIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;
  
  const auto& photon   = (hasSrc ? *cand.userData<pat::Photon>("src") : pat::Photon());
  const auto& phoIso   = cand.userData<reco::HIPhotonIsolation>("iso");
  const auto& sCluster = (photon.photonCore().isNonnull() ? photon.superCluster() : reco::SuperClusterRef());

  // photon information
  info.add("sigmaIEtaIEta", photon.full5x5_sigmaIetaIeta());
  info.add("hoverE", photon.hadronicOverEm());

  // super cluster information
  info.add("sCEta", (sCluster.isNonnull() ? sCluster->eta() : -99.9));
  info.add("sCPhi", (sCluster.isNonnull() ? sCluster->phi() : -99.9));

  // photon isolation information
  info.add("swissCrx", (phoIso ? phoIso->swissCrx() : -99.9));
  info.add("seedTime", (phoIso ? phoIso->seedTime() : -99.9));
  info.add("trackIsoR3PtCut20", (phoIso ? phoIso->trackIsoR3PtCut20() : 99.9));
  info.add("ecalClusterIsoR3", (phoIso ? phoIso->ecalClusterIsoR3() : 99.9));
  info.add("hcalRechitIsoR3", (phoIso ? phoIso->hcalRechitIsoR3() : 99.9));

  // push data and return index
  info.pushData(cand);
  return idx;
};


UShort_t
ParticleAnalyzer::fillJetInfo(const pat::GenericParticle& cand, const UShort_t& candIdx, const bool& force)
{
  const bool hasSrc = (std::abs(cand.pdgId())<=6 && cand.hasUserData("src") && cand.userData<pat::Jet>("src"));
  if (!force && !hasSrc) return USHRT_MAX;

  // fill jet information
  auto& info = particleInfo_["jet"];

  // add input information
  size_t index;
  const bool found = info.getIndex(index, cand);
  info.push(index, "candIdx", candIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;

  const auto& jet = (hasSrc ? *cand.userData<pat::Jet>("src") : pat::Jet());
  const bool hasJEC = jet.jecSetsAvailable();
  const auto& rawP4 = (hasJEC ? jet.correctedJet("Uncorrected") : pat::Jet());
  const auto& l2P4 = (hasJEC ? jet.correctedJet("L2Relative") : pat::Jet());
  const auto& l3P4 = (hasJEC ? jet.correctedJet("L3Absolute") : pat::Jet());

  // jet information
  info.add("rawPt", rawP4.pt());
  info.add("l2", l2P4.pt()/rawP4.pt());
  info.add("l3", l3P4.pt()/l2P4.pt());
  info.add("area", jet.jetArea());
  info.add("pileup", jet.pileup());
  info.add("hadronFlavour", getChar(jet.hadronFlavour()));
  info.add("partonFlavour", getShort(jet.partonFlavour()));

  // PF jet information
  const auto& isPF = jet.isPFJet();
  info.add("isPF", isPF);
  info.add("CHF", (isPF ? jet.chargedHadronEnergyFraction() : -1.));
  info.add("NHF", (isPF ? jet.neutralHadronEnergyFraction() : -1.));
  info.add("CEF", (isPF ? jet.chargedEmEnergyFraction() : -1.));
  info.add("NEF", (isPF ? jet.neutralEmEnergyFraction() : -1.));
  info.add("MUF", (isPF ? jet.muonEnergyFraction() : -1.));
  info.add("CHM", getUShort(isPF ? jet.chargedHadronMultiplicity() : 0));
  info.add("NHM", getUShort(isPF ? jet.neutralHadronMultiplicity() : 0));
  info.add("CEM", getUShort(isPF ? jet.electronMultiplicity() : 0));
  info.add("NEM", getUShort(isPF ? jet.photonMultiplicity() : 0));
  info.add("MUM", getUShort(isPF ? jet.muonMultiplicity() : 0));

  // b-tagging information
  UInt_t bTag = 0;
  const auto& v = jet.getPairDiscri();
  const auto n = std::min(size_t(32), v.size());
  for (size_t i=0; i<n; i++) { bTag += (std::next(v.begin(), i)->second ? 1U<<i : 0); }
  info.add("bTagging", bTag);

  // push data and return index
  info.pushData(cand);
  return idx;
};


UShort_t
ParticleAnalyzer::fillTauInfo(const pat::GenericParticle& cand, const UShort_t& candIdx, const bool& force)
{
  const bool hasSrc = (std::abs(cand.pdgId())==15 && cand.hasUserData("src") && cand.userData<pat::Tau>("src"));
  if (!force && !hasSrc) return USHRT_MAX; 

  // fill tau information
  auto& info = particleInfo_["tau"];

  // add input information
  size_t index;
  const bool found = info.getIndex(index, cand);
  info.push(index, "candIdx", candIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;

  const auto& tau = (hasSrc ? *cand.userData<pat::Tau>("src") : pat::Tau());

  // tau ID information
  UInt_t selector = 0;
  const auto& v = tau.tauIDs();
  const size_t n = std::min(size_t(32), v.size());
  for (size_t i=0; i<n; i++) { selector += (std::next(v.begin(), i)->second>0.5 ? 1U<<i : 0); }
  info.add("selector", selector);
  info.add("isPF", tau.isPFTau());

  // push data and return index
  info.pushData(cand);
  return idx;
};


UShort_t
ParticleAnalyzer::fillPFCandidateInfo(const pat::GenericParticle& cand, const UShort_t& candIdx, const bool& force)
{
  const bool hasSrc = (cand.hasUserData("src") && cand.userData<reco::PFCandidate>("src"));
  if (!force && !hasSrc) return USHRT_MAX;

  // fill jet information
  auto& info = particleInfo_["pf"];

  // add input information
  size_t index;
  const bool found = info.getIndex(index, cand);
  info.push(index, "candIdx", candIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;

  const auto& pfCand = (hasSrc ? *cand.userData<reco::PFCandidate>("src") : reco::PFCandidate());

  // general information
  info.add("id", getUChar(pfCand.particleId()));

  // calo information
  info.add("ecalE", pfCand.ecalEnergy());
  info.add("rawEcalE", pfCand.rawEcalEnergy());
  info.add("hcalE", pfCand.hcalEnergy());
  info.add("rawHcalE", pfCand.rawHcalEnergy());

  // push data and return index
  info.pushData(cand);
  return idx;
};


void
ParticleAnalyzer::fillGenParticleInfo(const edm::Event& iEvent)
{
  if (!isMC_) return;

  // fill generated particle information
  if (!genParticles_.empty())
  {
    for (const auto& cand : genParticles_) { fillGenParticleInfo(cand); }
  }

  // fill generator weight
  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(tok_genInfo_, genInfo);
  if (genInfo.isValid())
  {
    eventInfo_.add("genWeight", genInfo->weight());
  }
}


UShort_t
ParticleAnalyzer::fillGenParticleInfo(const reco::GenParticleRef& candR, const UShort_t& candIdx, const bool& force)
{
  const bool hasGen = candR.isNonnull();
  if (!force && !hasGen) return USHRT_MAX;

  // fill generated particle information
  auto& info = particleInfo_["gen"];

  const auto& cand = (hasGen ? *candR : reco::GenParticle());

  // add input information
  size_t index;
  const bool found = info.getIndex(index, cand);
  info.push(index, "candIdx", candIdx, true);
  const auto& idx = getUShort(index);
  // return if already added
  if (found) return idx;

  // basic information
  info.add("p", cand.p());
  info.add("y", cand.rapidity());
  info.add("pT", cand.pt());
  info.add("eta", cand.eta());
  info.add("phi", cand.phi());
  info.add("mass", cand.mass());
  info.add("charge", getChar(cand.charge()));
  info.add("pdgId", cand.pdgId());
  info.add("status", getUChar(cand.status()));

  UShort_t statusBit = 0;
  const auto& flags = cand.statusFlags().flags_;
  for (size_t i=0; i<flags.size(); i++) { i += (flags[i] ? 1U<<i : 0); }
  info.add("statusBit", statusBit);

  // decay vertex information
  float dl3D = -99.9, angle3D = -10., dl2D = -99.9, angle2D = -10.;
  if (cand.numberOfDaughters()>0 && cand.daughter(0))
  {
    auto mom = candR;
    while (mom.isNonnull() && mom->numberOfMothers()>0 && mom->mother(0) && mom->motherRef(0)->fromHardProcessDecayed()) { mom = mom->motherRef(0); }
    const auto lineOfFlight = cand.daughter(0)->vertex() - mom->vertex();
    dl3D = lineOfFlight.r();
    dl2D = lineOfFlight.rho();
    angle3D = angle(lineOfFlight.x(), lineOfFlight.y(), lineOfFlight.z(), cand.px(), cand.py(), cand.pz());
    angle2D =  angle(lineOfFlight.x(), lineOfFlight.y(), 0.0, cand.px(), cand.py(), 0.0);
  }
  info.add("angle3D", angle3D);
  info.add("angle2D", angle2D);
  info.add("decayLength3D", dl3D);
  info.add("decayLength2D", dl2D);

  // initialize daughter and mother information
  info.add("dauIdx", std::vector<UShort_t>());
  info.add("momIdx", std::vector<UShort_t>());

  // push data
  info.pushData(cand);

  // add daughter information
  if (cand.numberOfDaughters()>0)
  {
    for (size_t iDau=0; iDau<cand.numberOfDaughters(); iDau++)
    {
      const auto& dau = findGenDaughter(candR, iDau);
      const auto& dauIdx = fillGenParticleInfo(dau, idx);
      if (dauIdx!=USHRT_MAX) { info.push(idx, "dauIdx", dauIdx); }
    }
  }

  // mother information
  if (cand.numberOfMothers()>0)
  {
    for (size_t iMom=0; iMom<cand.numberOfMothers(); iMom++)
    {
      const auto& mom = findGenMother(candR);
      const auto& momIdx = fillGenParticleInfo(mom, idx);
      if (idx!=USHRT_MAX) { info.push(idx, "momIdx", momIdx); }
    }
  }

  // return index
  return idx;
}


bool
ParticleAnalyzer::addGenParticle(pat::GenericParticle& cand, const reco::GenParticleRefVector& genColl)
{
  if (genColl.empty()) return false;

  // check that it is in input collection if already matched
  if (cand.genParticle()) { return std::find(genColl.begin(), genColl.end(), cand.genParticleRef())!=genColl.end(); }

  // initialize information
  bool isSwap = false;
  reco::GenParticleRef genCand;
  cand.addUserInt("isGenMatched", 1);

  // do gen-reco matching
  double dPtRel = maxGenDeltaPtRel_;
  for (const auto& genPar : genColl)
  {
    // check if compatible
    if (!isCompatible(*genPar, cand)) continue;
    // case: final state particle (direct matching)
    if (cand.status()==1)
    {
      if (!isMatched(genPar->p4(), cand.p4(), maxGenDeltaR_, dPtRel)) continue;
      dPtRel = deltaPt(genPar->pt(), cand.pt());
      genCand = genPar;
    }
    // case: decayed particle (daughter matching)
    else if (cand.status()>1)
    {
      if (genPar->numberOfDaughters()==0) continue;
      // loop over daughters
      genCand = genPar;
      const auto& dauColl = *cand.userData<pat::GenericParticleCollection>("daughters");
      for (const auto& d: dauColl)
      {
        auto& dau = *const_cast<pat::GenericParticle*>(&d);
        reco::GenParticleRefVector genDauColl;
        findGenDaughters(genDauColl, genPar, dau, maxGenIter_);
        if (!addGenParticle(dau, genDauColl)) { genCand = reco::GenParticleRef(); break; }
        if (dau.userInt("isSwap")) { isSwap = true; }
      }
      if (genCand.isNonnull()) break;
    }
  }
  if (genCand.isNull()) return false;

  // add generated particle information
  isSwap = isSwap || (std::abs(genCand->mass()-cand.mass()) > 0.01);
  cand.addGenParticleRef(genCand);
  cand.addUserInt("isSwap", isSwap);
  return true;
}


void
ParticleAnalyzer::findGenDaughters(reco::GenParticleRefVector& genColl, const reco::GenParticleRef& genPar, const pat::GenericParticle& cand, const short& iter)
{
  if (iter>=0 || genPar.isNull()) return;
  for (size_t i=0; i<genPar->numberOfDaughters(); i++)
  {
    const auto& dau = findGenDaughter(genPar, i);
    if (isCompatible(*dau, cand)) { genColl.push_back(dau); }
    findGenDaughters(genColl, dau, cand, iter-1);
  }
}


reco::GenParticleRef
ParticleAnalyzer::findGenDaughter(const reco::GenParticleRef& par, const size_t& idx)
{
  if (par.isNull() || par->numberOfDaughters()<=idx) return reco::GenParticleRef();
  auto dau = par->daughterRef(idx);
  const auto pdgId = dau->pdgId();
  while (dau->numberOfDaughters()>0 && dau->daughterRef(0)->pdgId()==pdgId)
  {
    dau = dau->daughterRef(0);
  }
  return dau;
}


reco::GenParticleRef
ParticleAnalyzer::findGenMother(const reco::GenParticleRef& par)
{
  if (par.isNull() || par->numberOfMothers()==0) return reco::GenParticleRef();
  auto mom = par->motherRef(0);
  const auto pdgId = par->pdgId();
  while (mom->numberOfMothers()>0 && mom->pdgId()==pdgId)
  {
    mom = mom->motherRef(0);
  }
  if (mom->pdgId()==pdgId) return reco::GenParticleRef();
  return mom;
}


void
ParticleAnalyzer::initTree()
{
  if (tree_) return;
  tree_ = fileService_->make<TTree>("ParticleTree","ParticleTree");

  // add event branches
  eventInfo_.initTree(*tree_);

  // add particle branches
  for (const auto& p : std::vector<std::string>({"cand", "trk", "elec", "muon", "tau", "pho", "jet", "pf", "gen", "trig"}))
  {
    if (particleInfo_.find(p)!=particleInfo_.end())
    {
      particleInfo_.at(p).initTree(*tree_, (p+"_"));
    }
  }
}


// ------------ method called once each job just before starting event
//loop  ------------
void
ParticleAnalyzer::beginJob()
{
}


//--------------------------------------------------------------------------------------------------
void 
ParticleAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  bool changed = true;
  EDConsumerBase::Labels triggerResultsLabel;
  EDConsumerBase::labelsForToken(tok_triggerResults_, triggerResultsLabel);
  if (!hltPrescaleProvider_.init(iRun, iSetup, triggerResultsLabel.process, changed)) { throw(cms::Exception("Trigger")<<"HLT provider failed init"); }
  auto& l1tGlobalUtil = *const_cast<l1t::L1TGlobalUtil*>(&hltPrescaleProvider_.l1tGlobalUtil());
  l1tGlobalUtil.setUnprescaledUnmasked(false, true);
  l1tGlobalUtil.retrieveL1Setup(iSetup);
}


// ------------ method called once each job just after ending the event
//loop  ------------
void 
ParticleAnalyzer::endJob()
{
}


//define this as a plug-in
DEFINE_FWK_MODULE(ParticleAnalyzer);

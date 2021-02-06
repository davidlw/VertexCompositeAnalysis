// system include files
#include <string>
#include <vector>
#include <map>
#include <TTree.h>

// user include files
#include "DataFormats/Candidate/interface/Candidate.h"

//
// class decleration for containers
//

typedef std::map<std::string, std::vector<size_t> > TriggerIndexMap;


class Container
{
 public:
  Container() {};
  ~Container() {};

  // getters
  const std::map<std::string, bool     >& boolM()   const { return boolM_;   };
  const std::map<std::string, Char_t   >& charM()   const { return charM_;   };
  const std::map<std::string, short    >& shortM()  const { return shortM_;  };
  const std::map<std::string, int      >& intM()    const { return intM_;    };
  const std::map<std::string, UChar_t  >& ucharM()  const { return ucharM_;  };
  const std::map<std::string, UShort_t >& ushortM() const { return ushortM_; };
  const std::map<std::string, UInt_t   >& uintM()   const { return uintM_;   };
  const std::map<std::string, float    >& floatM()  const { return floatM_;  };
  const std::map<std::string, std::vector<bool>     >& boolVM()   const { return boolVM_;   };
  const std::map<std::string, std::vector<Char_t>   >& charVM()   const { return charVM_;   };
  const std::map<std::string, std::vector<short>    >& shortVM()  const { return shortVM_;  };
  const std::map<std::string, std::vector<int>      >& intVM()    const { return intVM_;    };
  const std::map<std::string, std::vector<UChar_t>  >& ucharVM()  const { return ucharVM_;  };
  const std::map<std::string, std::vector<UShort_t> >& ushortVM() const { return ushortVM_; };
  const std::map<std::string, std::vector<UInt_t>   >& uintVM()   const { return uintVM_;   };
  const std::map<std::string, std::vector<float>    >& floatVM()  const { return floatVM_;  };

  // setters
  template <class T>
  void add(const std::string& n, const T&        v) = delete; // avoid implicit conversion
  void add(const std::string& n, const bool&     v) { boolM_[n]   = v; };
  void add(const std::string& n, const Char_t&   v) { charM_[n]   = v; };
  void add(const std::string& n, const short&    v) { shortM_[n]  = v; };
  void add(const std::string& n, const int&      v) { intM_[n]  = v;   };
  void add(const std::string& n, const UChar_t&  v) { ucharM_[n]  = v; };
  void add(const std::string& n, const UShort_t& v) { ushortM_[n] = v; };
  void add(const std::string& n, const UInt_t&   v) { uintM_[n]   = v; };
  void add(const std::string& n, const float&    v) { floatM_[n]  = v; };
  void add(const std::string& n, const double&   v) { floatM_[n]  = v; };
  void add(const std::string& n, const std::vector<bool>&     v) { boolVM_[n]   = v; };
  void add(const std::string& n, const std::vector<Char_t>&   v) { charVM_[n]   = v; };
  void add(const std::string& n, const std::vector<short>&    v) { shortVM_[n]  = v; };
  void add(const std::string& n, const std::vector<int>&      v) { intVM_[n]  = v;   };
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
  void push(const std::string& n, const double&   v, const bool& c=0) { push(n, float(v), c); }

  void copyData(Container& data, const std::string& n="") const
  {
    for (const auto& d : boolM_   ) { data.add(n+d.first, d.second); }
    for (const auto& d : charM_   ) { data.add(n+d.first, d.second); }
    for (const auto& d : shortM_  ) { data.add(n+d.first, d.second); }
    for (const auto& d : intM_    ) { data.add(n+d.first, d.second); }
    for (const auto& d : ucharM_  ) { data.add(n+d.first, d.second); }
    for (const auto& d : ushortM_ ) { data.add(n+d.first, d.second); }
    for (const auto& d : uintM_   ) { data.add(n+d.first, d.second); }
    for (const auto& d : floatM_  ) { data.add(n+d.first, d.second); }
    for (const auto& d : boolVM_  ) { for (uint i=0; i<d.second.size(); i++) { data.add(n+d.first+Form("%u",i), d.second[i]); } }
    for (const auto& d : charVM_  ) { for (uint i=0; i<d.second.size(); i++) { data.add(n+d.first+Form("%u",i), d.second[i]); } }
    for (const auto& d : shortVM_ ) { for (uint i=0; i<d.second.size(); i++) { data.add(n+d.first+Form("%u",i), d.second[i]); } }
    for (const auto& d : intVM_   ) { for (uint i=0; i<d.second.size(); i++) { data.add(n+d.first+Form("%u",i), d.second[i]); } }
    for (const auto& d : ucharVM_ ) { for (uint i=0; i<d.second.size(); i++) { data.add(n+d.first+Form("%u",i), d.second[i]); } }
    for (const auto& d : ushortVM_) { for (uint i=0; i<d.second.size(); i++) { data.add(n+d.first+Form("%u",i), d.second[i]); } }
    for (const auto& d : uintVM_  ) { for (uint i=0; i<d.second.size(); i++) { data.add(n+d.first+Form("%u",i), d.second[i]); } }
    for (const auto& d : floatVM_ ) { for (uint i=0; i<d.second.size(); i++) { data.add(n+d.first+Form("%u",i), d.second[i]); } }
  };

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

  template <class T>
  void erase(std::map<std::string, T>& c, const std::string& n)
  {
    for (const auto& it : std::map<std::string, T>(c))
    {
      const bool has = it.first.find(n)!=std::string::npos;
      if (has) { c.erase(it.first); }
    }
  };

  void erase(const std::string& n, const std::vector<std::string>& tv)
  {
    for (const auto& t : tv)
    {
      if      (t=="bool"   ) erase(boolM_, n);
      else if (t=="uchar"  ) erase(ucharM_, n);
      else if (t=="char"   ) erase(charM_, n);
      else if (t=="ushort" ) erase(ushortM_, n);
      else if (t=="short"  ) erase(shortM_, n);
      else if (t=="uint"   ) erase(uintM_, n);
      else if (t=="int"    ) erase(intM_, n);
      else if (t=="float"  ) erase(floatM_, n);
      else if (t=="boolV"  ) erase(boolVM_, n);
      else if (t=="ucharV" ) erase(ucharVM_, n);
      else if (t=="charV"  ) erase(charVM_, n);
      else if (t=="ushortV") erase(ushortVM_, n);
      else if (t=="shortV" ) erase(shortVM_, n);
      else if (t=="uintV"  ) erase(uintVM_, n);
      else if (t=="intV"   ) erase(intVM_, n);
      else if (t=="floatV" ) erase(floatVM_, n);
    }
  };

  // tree initializer
  void initTree(TTree& tree, const std::string& n="")
  {
    for (auto& d : boolM_   ) { tree.Branch((n+d.first).c_str(), &d.second, (n+d.first+"/O").c_str()); }
    for (auto& d : ucharM_  ) { tree.Branch((n+d.first).c_str(), &d.second, (n+d.first+"/b").c_str()); }
    for (auto& d : charM_   ) { tree.Branch((n+d.first).c_str(), &d.second, (n+d.first+"/B").c_str()); }
    for (auto& d : ushortM_ ) { tree.Branch((n+d.first).c_str(), &d.second, (n+d.first+"/s").c_str()); }
    for (auto& d : shortM_  ) { tree.Branch((n+d.first).c_str(), &d.second, (n+d.first+"/S").c_str()); }
    for (auto& d : uintM_   ) { tree.Branch((n+d.first).c_str(), &d.second, (n+d.first+"/i").c_str()); }
    for (auto& d : intM_    ) { tree.Branch((n+d.first).c_str(), &d.second, (n+d.first+"/I").c_str()); }
    for (auto& d : floatM_  ) { tree.Branch((n+d.first).c_str(), &d.second, (n+d.first+"/F").c_str()); }
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
  const size_t& size() const { return size_; };

  template <class T>
  bool getIndex(size_t& index, const T& cand) const
  {
    index = size_;
    const auto& parIt = parM_.find(getPar(cand));
    if (parIt==parM_.end()) return false;
    index = parIt->second;
    return true;
  };

  template <class T>
  const T        get(const size_t& i, const std::vector<T>& v, const T& d) const { return (i < size_ ? v[i] : d); };
  template <class T>
  const T        get(const size_t& i, const std::string& n, const T&        d) const = delete; // avoid implicit conversion
  const bool     get(const size_t& i, const std::string& n, const bool&     d) const { return (i < size_ ? boolVM_.at(n)[i]   : d); };
  const int      get(const size_t& i, const std::string& n, const Int_t&    d) const { return (i < size_ ? intVM_.at(n)[i]    : d); };
  const UChar_t  get(const size_t& i, const std::string& n, const UChar_t&  d) const { return (i < size_ ? ucharVM_.at(n)[i]  : d); };
  const UShort_t get(const size_t& i, const std::string& n, const UShort_t& d) const { return (i < size_ ? ushortVM_.at(n)[i] : d); };
  const std::vector<UShort_t> get(const size_t& i, const std::string& n, const std::vector<UShort_t>& d) const { return (i < size_ ? ushortVVM_.at(n)[i] : d); };
  const std::vector<UInt_t  > get(const size_t& i, const std::string& n, const std::vector<UInt_t  >& d) const { return (i < size_ ? uintVVM_.at(n)[i]   : d); };

  // setters
  template <class T>
  void add(const std::string& n, const T& v) { data_.add(n, v); };

  template <class T>
  void add(const size_t& i, const std::string& n, const T&        v) = delete; // avoid implicit conversion
  void add(const size_t& i, const std::string& n, const bool&     v) { if (i < size_) { boolVM_[n][i]   = v; } else { add(n, v); } };
  void add(const size_t& i, const std::string& n, const UShort_t& v) { if (i < size_) { ushortVM_[n][i] = v; } else { add(n, v); } };
  void add(const size_t& i, const std::string& n, const UInt_t&   v) { if (i < size_) { uintVM_[n][i]   = v; } else { add(n, v); } };

  template <class T>
  void push(const std::string& n, const T& v, const bool& c=0) { data_.push(n, v, c); };

  template <class T>
  void pushBack(std::vector<std::vector<T> >& c, const size_t& i, const std::string& n, const T& v) {
    if (i>=c.size()) throw std::logic_error(Form("Invalid index %lu for %s", i, n.c_str()));
    c[i].push_back(v);
  };

  template <class T>
  void pushV(const size_t& i, const std::string& n, const T&        v, const bool& c=0) = delete; // avoid implicit conversion
  void pushV(const size_t& i, const std::string& n, const UChar_t&  v, const bool& c=0) { if(!c || v!=UCHAR_MAX) pushBack(ucharVVM_[n], i, n, v);  };
  void pushV(const size_t& i, const std::string& n, const UShort_t& v, const bool& c=0) { if(!c || v!=USHRT_MAX) pushBack(ushortVVM_[n], i, n, v); };
  void pushV(const size_t& i, const std::string& n, const UInt_t&   v, const bool& c=0) { if(!c || v!=UINT_MAX ) pushBack(uintVVM_[n], i, n, v);   };
  void pushV(const size_t& i, const std::string& n, const float&    v, const bool& c=0) { pushBack(floatVVM_[n], i, n, v); };
  void pushV(const size_t& i, const std::string& n, const double&   v, const bool& c=0) { pushV(i, n, float(v), c); };
  template <class T>
  void push(const size_t& i, const std::string& n, const T& v, const bool& c=0) { (i < size_ ? pushV(i, n, v, c) : push(n, v, c)); }

  template <class T>
  void pushData(const T& cand)
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
    for (const auto& d : data_.floatVM() ) { floatVVM_[d.first].push_back(d.second);  }
    parM_[par] = size_++;
    data_.clear();
  };

  void copyData(Container& data, const size_t& i, const std::string& n="") const
  {
    for (const auto& d : boolVM_   ) { data.add(n+d.first, get(i, d.second, false));        }
    for (const auto& d : charVM_   ) { data.add(n+d.first, get(i, d.second, char(-99)));    }
    for (const auto& d : shortVM_  ) { data.add(n+d.first, get(i, d.second, short(-99)));   }
    for (const auto& d : intVM_    ) { data.add(n+d.first, get(i, d.second, int(-99)));     }
    for (const auto& d : ucharVM_  ) { data.add(n+d.first, get(i, d.second, UChar_t(0)));   }
    for (const auto& d : ushortVM_ ) { data.add(n+d.first, get(i, d.second, UShort_t(0)));  }
    for (const auto& d : uintVM_   ) { data.add(n+d.first, get(i, d.second, UInt_t(0)));    }
    for (const auto& d : floatVM_  ) { data.add(n+d.first, get(i, d.second, float(-99.9))); }
    for (const auto& d : floatVVM_ ) { for (uint j=0; i<size_ && j<d.second[i].size(); j++) { data.add(n+d.first+Form("%u",j), d.second[i][j]); } }
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
    for (auto& d : floatVVM_ ) { d.second.clear(); }
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
    for (auto& d : floatVVM_ ) { tree.Branch((n+d.first).c_str(), &d.second); }
  };

 private:
  typedef std::tuple<float, float, float, float, char, int> Particle;

  Particle getPar(const reco::Candidate& cand) const
  {
    return Particle(cand.pt(), cand.eta(), cand.phi(), cand.mass(), cand.charge(), cand.pdgId());
  };

  Particle getPar(const reco::Track& track) const
  {
    return Particle(track.pt(), track.eta(), track.phi(), 0.0, track.charge(), 0);
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
  std::map<std::string, std::vector<std::vector<float> >    > floatVVM_;
};


class TriggerInfo
{
 public:
  TriggerInfo() : triggerIndex_(-1), filterIndex_(-1), minN_(0), hltPDs_(0),
                  validPrescale_(0), validLumi_(0), hltPrescale_(0), l1Prescale_(0),
                  triggerName_(""), filterName_(""), triggerBit_({}), filterObjects_({}),
                  recordLumi_(0), totalLumi_(0)
  {
  };
  ~TriggerInfo() {};

  // getters
  const int&                  triggerIndex()  const { return triggerIndex_;  }
  const int&                  filterIndex()   const { return filterIndex_;   }
  const int&                  minN()          const { return minN_;          }
  const bool&                 validPrescale() const { return validPrescale_; }
  const UShort_t&             hltPrescale()   const { return hltPrescale_;   }
  const UShort_t&             l1Prescale()    const { return l1Prescale_;    }
  const UChar_t&              hltPDs()        const { return hltPDs_;        }
  const std::string&          triggerName()   const { return triggerName_;   }
  const std::string&          filterName()    const { return filterName_;    }
  const std::array<bool,4>&   triggerBit()    const { return triggerBit_;    }
  const bool                  triggerBit(const size_t& i) const { return (i<4 ? triggerBit_[i] : false); }
  const TriggerIndexMap&      filterObjects() const { return filterObjects_; }
  const bool&                 validLumi()     const { return validLumi_;     }
  const float&                recordLumi()    const { return recordLumi_;    }
  const float&                totalLumi()     const { return totalLumi_;     }

  // setters
  void setInfo(const int& triggerIndex, const int& filterIndex,
               const std::string& triggerName, const std::string& filterName, const int& minN,
               const bool& validPrescale, const UShort_t& hltPrescale, const UShort_t& l1Prescale,
               const UChar_t& hltPDs, const std::array<bool,4>& bit, const TriggerIndexMap& filterObjects)
  {
    triggerIndex_  = triggerIndex;
    filterIndex_   = filterIndex;
    triggerName_   = triggerName;
    filterName_    = filterName;
    minN_          = minN;
    validPrescale_ = validPrescale;
    hltPrescale_   = hltPrescale;
    l1Prescale_    = l1Prescale;
    hltPDs_        = hltPDs;
    triggerBit_    = bit;
    filterObjects_ = filterObjects;
  }

  void setLumiInfo(const double& recordLumi, const double& totalLumi)
  {
    recordLumi_ = recordLumi;
    totalLumi_  = totalLumi;
    validLumi_  = true;
  };

 private:
  int triggerIndex_, filterIndex_, minN_;
  UChar_t hltPDs_;
  bool validPrescale_, validLumi_;
  UShort_t hltPrescale_, l1Prescale_;
  std::string triggerName_, filterName_;
  std::array<bool, 4> triggerBit_;
  TriggerIndexMap filterObjects_;
  float recordLumi_, totalLumi_;
};


class MatchInfo
{
 public:
  MatchInfo() : collection_(""), maxDeltaR_(0), maxDeltaPtRel_(0),
                maxDeltaEta_(0), maxDeltaPhi_(0)
  {
  };
  ~MatchInfo() {};

  // getters
  const std::string& collection()    const { return collection_;    }
  const double&      maxDeltaR()     const { return maxDeltaR_;     }
  const double&      maxDeltaPtRel() const { return maxDeltaPtRel_; }
  const double&      maxDeltaEta()   const { return maxDeltaEta_;   }
  const double&      maxDeltaPhi()   const { return maxDeltaPhi_;   }

  // setters
  void setInfo(const std::string& collection, const double& maxDeltaR, const double& maxDeltaPtRel,
               const double& maxDeltaEta, const double& maxDeltaPhi)
  {
    collection_    = collection;
    maxDeltaR_     = maxDeltaR;
    maxDeltaPtRel_ = maxDeltaPtRel;
    maxDeltaEta_   = maxDeltaEta;
    maxDeltaPhi_   = maxDeltaPhi;
  };

 private:
  std::string collection_;
  double maxDeltaR_, maxDeltaPtRel_, maxDeltaEta_, maxDeltaPhi_;
};


//
//// constants, enums and typedefs
//

typedef std::map<std::string, ParticleContainer> ParticleContainerMap;
typedef std::vector<TriggerInfo> TriggerContainer;
typedef std::map<std::string, MatchInfo> MatchContainer;

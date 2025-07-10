#include "DataFormats/PatCandidates/interface/UserHolder.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include <vector>

namespace {
  pat::UserHolder<
    edm::Ref<std::vector<reco::Track>, reco::Track,
    edm::refhelper::FindUsingAdvance<std::vector<reco::Track>, reco::Track> > > dummy_phi_ref_holder;
}

#include <DataFormats/PatCandidates/interface/UserData.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h>


pat::UserHolder<pat::Electron> dummy1;
pat::UserHolder<pat::Muon> dummy2;
pat::UserHolder<pat::Tau> dummy3;
pat::UserHolder<pat::Jet> dummy4;
pat::UserHolder<pat::Photon> dummy5;
pat::UserHolder<std::vector<math::XYZTLorentzVector> > dummy6;
pat::UserHolder<pat::GenericParticleRefVector> dummy7;
KinematicParametersError dummy8;
pat::UserHolder<KinematicParametersError> dummy9;
pat::UserHolder<reco::TrackRef> dummy10;
pat::UserHolder<reco::Conversion> dummy11;

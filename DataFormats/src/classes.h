#include <DataFormats/PatCandidates/interface/UserData.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParameters.h>
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h>
#include <DataFormats/TrackReco/interface/Track.h>


pat::UserHolder<pat::Electron> dummy1;
pat::UserHolder<pat::Muon> dummy2;
pat::UserHolder<pat::Tau> dummy3;
pat::UserHolder<pat::Jet> dummy4;
pat::UserHolder<std::vector<reco::Candidate::LorentzVector> > dummy5;
pat::UserHolder<pat::GenericParticleRefVector> dummy6;
KinematicParameters dummy7;
pat::UserHolder<KinematicParameters> dummy8;
KinematicParametersError dummy9;
pat::UserHolder<KinematicParametersError> dummy10;
pat::UserHolder<reco::TrackRef> dummy11;

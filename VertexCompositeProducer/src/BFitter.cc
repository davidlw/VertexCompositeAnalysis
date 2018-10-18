// -*- C++ -*-
//
// Package:    VertexCompositeProducer
// Class:      BFitter
// 
/**\class BFitter BFitter.cc VertexCompositeAnalysis/VertexCompositeProducer/src/BFitter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
//
//

#include "VertexCompositeAnalysis/VertexCompositeProducer/interface/BFitter.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include <TMath.h>
#include <TVector3.h>
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"

const float piMassB = 0.13957018;
const float piMassBSquared = piMassB*piMassB;
const float kaonMassB = 0.493677;
const float kaonMassBSquared = kaonMassB*kaonMassB;
const float d0MassB = 1.86484;
const float bMassB = 5.27929;
float piMassB_sigma = 3.5E-7f;
float kaonMassB_sigma = 1.6E-5f;
float d0MassB_sigma = d0MassB*1.e-6;

// Constructor and (empty) destructor
BFitter::BFitter(const edm::ParameterSet& theParameters,  edm::ConsumesCollector && iC)
{
  using std::string;

  token_beamSpot = iC.consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  token_tracks = iC.consumes<reco::TrackCollection>(theParameters.getParameter<edm::InputTag>("trackRecoAlgorithm"));
  token_vertices = iC.consumes<reco::VertexCollection>(theParameters.getParameter<edm::InputTag>("vertexRecoAlgorithm"));
  token_d0s = iC.consumes<reco::VertexCompositeCandidateCollection>(theParameters.getParameter<edm::InputTag>("d0RecoAlgorithm"));
  token_dedx = iC.consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));

  // Second, initialize post-fit cuts
  mPiDCutMin = theParameters.getParameter<double>(string("mPiDCutMin"));
  mPiDCutMax = theParameters.getParameter<double>(string("mPiDCutMax"));
  batTkDCACut = theParameters.getParameter<double>(string("batTkDCACut"));
  batTkChi2Cut = theParameters.getParameter<double>(string("batTkChi2Cut"));
  batTkNhitsCut = theParameters.getParameter<int>(string("batTkNhitsCut"));
  batTkPtCut = theParameters.getParameter<double>(string("batTkPtCut"));
  batTkPtErrCut = theParameters.getParameter<double>(string("batTkPtErrCut"));
  batTkEtaCut = theParameters.getParameter<double>(string("batTkEtaCut"));
  bVtxChi2Cut = theParameters.getParameter<double>(string("bVtxChi2Cut"));
  bRVtxCut = theParameters.getParameter<double>(string("bVtx2DCut"));
  bRVtxSigCut = theParameters.getParameter<double>(string("bVtxSignificance2DCut"));
  bLVtxCut = theParameters.getParameter<double>(string("bVtx3DCut"));
  bLVtxSigCut = theParameters.getParameter<double>(string("bVtxSignificance3DCut"));
  bCollinCut2D = theParameters.getParameter<double>(string("bCollinCut2D"));
  bCollinCut3D = theParameters.getParameter<double>(string("bCollinCut3D"));
  bMassCut = theParameters.getParameter<double>(string("bMassCut"));
  batDauTransImpactSigCut = theParameters.getParameter<double>(string("batDauTransImpactSigCut"));
  batDauLongImpactSigCut = theParameters.getParameter<double>(string("batDauLongImpactSigCut"));
  bVtxChiProbCut = theParameters.getParameter<double>(string("bVtxChiProbCut"));
  bPtCut = theParameters.getParameter<double>(string("bPtCut"));
  bAlphaCut = theParameters.getParameter<double>(string("bAlphaCut"));
  bAlpha2DCut = theParameters.getParameter<double>(string("bAlpha2DCut"));
  isWrongSignB = theParameters.getParameter<bool>(string("isWrongSignB"));
}

BFitter::~BFitter() {
}

// Method containing the algorithm for vertex reconstruction
void BFitter::fitAll(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using std::vector;
  using std::cout;
  using std::endl;
  using namespace reco;
  using namespace edm;
  using namespace std; 

  typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
  typedef ROOT::Math::SVector<double, 3> SVector3;

  // Create std::vectors for Tracks and TrackRefs (required for
  //  passing to the KalmanVertexFitter)
  std::vector<TrackRef> theTrackRefs;
  std::vector<TransientTrack> theTransTracks;

  // Handles for tracks, B-field, and tracker geometry
  Handle<reco::TrackCollection> theTrackHandle;
  Handle<reco::VertexCollection> theVertexHandle;
  Handle<reco::VertexCompositeCandidateCollection> theD0Handle;
  Handle<reco::BeamSpot> theBeamSpotHandle;
  ESHandle<MagneticField> bFieldHandle;
  Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle;

  // Get the tracks, vertices from the event, and get the B-field record
  //  from the EventSetup
  iEvent.getByToken(token_tracks, theTrackHandle); 
  iEvent.getByToken(token_vertices, theVertexHandle);
  iEvent.getByToken(token_d0s, theD0Handle);
  iEvent.getByToken(token_beamSpot, theBeamSpotHandle);  
  iEvent.getByToken(token_dedx, dEdxHandle);

  if( !theTrackHandle->size() ) return;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  magField = bFieldHandle.product();

  // Setup TMVA
//  mvaValValueMap = auto_ptr<edm::ValueMap<float> >(new edm::ValueMap<float>);
//  edm::ValueMap<float>::Filler mvaFiller(*mvaValValueMap);

  bool isVtxPV = 0;
  double xVtx=-99999.0;
  double yVtx=-99999.0;
  double zVtx=-99999.0;
  double xVtxError=-999.0;
  double yVtxError=-999.0;
  double zVtxError=-999.0;
  const reco::VertexCollection vtxCollection = *(theVertexHandle.product());
  reco::VertexCollection::const_iterator vtxPrimary = vtxCollection.begin();
  if(vtxCollection.size()>0 && !vtxPrimary->isFake() && vtxPrimary->tracksSize()>=2)
  {
    isVtxPV = 1;
    xVtx = vtxPrimary->x();
    yVtx = vtxPrimary->y();
    zVtx = vtxPrimary->z();
    xVtxError = vtxPrimary->xError();
    yVtxError = vtxPrimary->yError();
    zVtxError = vtxPrimary->zError();
  }
  else {
    isVtxPV = 0;
    xVtx = theBeamSpotHandle->position().x();
    yVtx = theBeamSpotHandle->position().y();
    zVtx = 0.0;
    xVtxError = theBeamSpotHandle->BeamWidthX();
    yVtxError = theBeamSpotHandle->BeamWidthY();
    zVtxError = 0.0;
  }
  math::XYZPoint bestvtx(xVtx,yVtx,zVtx);

  // Fill vectors of TransientTracks and TrackRefs after applying preselection cuts.
  for(unsigned int indx = 0; indx < theTrackHandle->size(); indx++) {
    TrackRef tmpRef( theTrackHandle, indx );
    bool quality_ok = true;
    if (qualities.size()!=0) {
      quality_ok = false;
      for (unsigned int ndx_ = 0; ndx_ < qualities.size(); ndx_++) {
	if (tmpRef->quality(qualities[ndx_])){
	  quality_ok = true;
	  break;          
	}
      }
    }
    if( !quality_ok ) continue;

    if( tmpRef->normalizedChi2() < batTkChi2Cut &&
        tmpRef->numberOfValidHits() >= batTkNhitsCut &&
        tmpRef->ptError() / tmpRef->pt() < batTkPtErrCut &&
        tmpRef->pt() > batTkPtCut && fabs(tmpRef->eta()) < batTkEtaCut ) {
//      TransientTrack tmpTk( *tmpRef, &(*bFieldHandle), globTkGeomHandle );
      TransientTrack tmpTk( *tmpRef, magField );

      double dzvtx = tmpRef->dz(bestvtx);
      double dxyvtx = tmpRef->dxy(bestvtx);      
      double dzerror = sqrt(tmpRef->dzError()*tmpRef->dzError()+zVtxError*zVtxError);
      double dxyerror = sqrt(tmpRef->d0Error()*tmpRef->d0Error()+xVtxError*yVtxError);

      double dauLongImpactSig = dzvtx/dzerror;
      double dauTransImpactSig = dxyvtx/dxyerror;

      if( fabs(dauTransImpactSig) > batDauTransImpactSigCut && fabs(dauLongImpactSig) > batDauLongImpactSigCut ) {
        theTrackRefs.push_back( tmpRef );
        theTransTracks.push_back( tmpTk );
      }
    }
  }

  const reco::VertexCompositeCandidateCollection theD0s = *(theD0Handle.product());
  for(unsigned it=0; it<theD0s.size(); ++it){

    const reco::VertexCompositeCandidate & theD0 = theD0s[it];

    float massWindow = 0.040;
    if(theD0.mass() > d0MassB + massWindow || theD0.mass() < d0MassB - massWindow) continue;

    vector<RecoChargedCandidate> d0daughters;
    vector<TrackRef> theDaughterTracks;
    
    d0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                            (theD0.daughter(0))) );
    d0daughters.push_back( *(dynamic_cast<const reco::RecoChargedCandidate *>
                            (theD0.daughter(1))) );
   
    for(unsigned int j = 0; j < d0daughters.size(); ++j) {
      theDaughterTracks.push_back(d0daughters[j].track());
    }

    vector<TransientTrack> tracksForKalmanFit;
    for (unsigned int ndx = 0; ndx < theDaughterTracks.size(); ndx++) {
        tracksForKalmanFit.push_back(TransientTrack(theDaughterTracks[ndx], &(*bFieldHandle)));
    }

    vector<float> d0DauMasses;
    vector<float> d0DauMasses_sigma;
    if (theD0.pdgId()>0) {
      d0DauMasses.push_back(piMassB);
      d0DauMasses.push_back(kaonMassB);
      d0DauMasses_sigma.push_back(piMassB_sigma);
      d0DauMasses_sigma.push_back(kaonMassB_sigma);
    } else if (theD0.pdgId()<0) {
      d0DauMasses.push_back(kaonMassB);
      d0DauMasses.push_back(piMassB);
      d0DauMasses_sigma.push_back(kaonMassB_sigma);
      d0DauMasses_sigma.push_back(piMassB_sigma);
    }

    for(unsigned int trdx = 0; trdx < theTrackRefs.size(); trdx++) {

       if ( theTrackRefs[trdx].isNull() ) continue;

       if(isWrongSignB && theTrackRefs[trdx]->charge()*theD0.pdgId()<0) continue;
       if(!isWrongSignB && theTrackRefs[trdx]->charge()*theD0.pdgId()>0) continue;

       bool match = false;

       for(unsigned int j = 0; j < theDaughterTracks.size(); ++j) {
         if (theTrackRefs[trdx]->charge() == theDaughterTracks[j]->charge() &&
             theTrackRefs[trdx]->momentum() == theDaughterTracks[j]->momentum() ) match = true;
         if (match) break;
       }
       if (match) continue; // Track is already used in making the D0

       // pre-selections on B invariant mass and pT to save time
       double ETotPre = sqrt(theTrackRefs[trdx]->momentum().mag2()+piMassBSquared) + theD0.energy(); 
       double PTotPreSq = (theTrackRefs[trdx]->momentum() + theD0.momentum()).mag2(); 
       double totalPt = sqrt((theTrackRefs[trdx]->momentum() + theD0.momentum()).perp2());
       double massPre = sqrt( ETotPre*ETotPre - PTotPreSq);

       if(massPre > mPiDCutMax || massPre < mPiDCutMin) continue;
       if(totalPt < bPtCut ) continue;

       TransientTrack dauPos(theDaughterTracks[0], &(*bFieldHandle) );
       TransientTrack dauNeg(theDaughterTracks[1], &(*bFieldHandle) );
       TransientTrack batDau(theTrackRefs[trdx], &(*bFieldHandle) );

       if (!dauPos.isValid()) continue;
       if (!dauNeg.isValid()) continue;
       if (!batDau.isValid()) continue;

       //Creating a KinematicParticleFactory
       KinematicParticleFactoryFromTransientTrack pFactory;

       float chi = 0.;
       float ndf = 0.;
       vector<RefCountedKinematicParticle> d0Particles;
       d0Particles.push_back(pFactory.particle(dauPos,d0DauMasses[0],chi,ndf,d0DauMasses_sigma[0]));
       d0Particles.push_back(pFactory.particle(dauNeg,d0DauMasses[1],chi,ndf,d0DauMasses_sigma[1]));

       KinematicParticleVertexFitter fitter;
       RefCountedKinematicTree d0VertexFitTree;
       d0VertexFitTree = fitter.fit(d0Particles);
       if (!d0VertexFitTree->isValid()) continue;

       d0VertexFitTree->movePointerToTheTop();

       RefCountedKinematicParticle d0_vFit = d0VertexFitTree->currentParticle();
       RefCountedKinematicVertex d0_vFit_vertex = d0VertexFitTree->currentDecayVertex();

       d0VertexFitTree->movePointerToTheFirstChild();
       RefCountedKinematicParticle d0Pi1 = d0VertexFitTree->currentParticle();
       d0VertexFitTree->movePointerToTheNextChild();
       RefCountedKinematicParticle d0Pi2 = d0VertexFitTree->currentParticle();

       // now apply D0 mass constraint
       KinematicParticleFitter csFitterD0;
       KinematicConstraint * bmeson = new MassKinematicConstraint(d0MassB,d0MassB_sigma);

       d0VertexFitTree->movePointerToTheTop();
       d0VertexFitTree = csFitterD0.fit(bmeson,d0VertexFitTree);
       if (!d0VertexFitTree->isValid()) continue;
       d0VertexFitTree->movePointerToTheTop();
       RefCountedKinematicParticle d0_vFit_withMC = d0VertexFitTree->currentParticle();

       if (!d0_vFit_withMC->currentState().isValid()) continue;

       vector<RefCountedKinematicParticle> bFitParticles;

       bFitParticles.push_back(pFactory.particle(batDau,piMassB,chi,ndf,piMassB_sigma));
       bFitParticles.push_back(d0_vFit_withMC);

       //fit B
       RefCountedKinematicTree bFitTree = fitter.fit(bFitParticles);
       if (!bFitTree->isValid()) continue;

       bFitTree->movePointerToTheTop();
       RefCountedKinematicParticle bCand = bFitTree->currentParticle();
       if (!bCand->currentState().isValid()) continue;

       RefCountedKinematicVertex bDecayVertex = bFitTree->currentDecayVertex();
       if (!bDecayVertex->vertexIsValid()) continue;

       // get children from final B fit
       bFitTree->movePointerToTheFirstChild();
       RefCountedKinematicParticle batPionCand = bFitTree->currentParticle();
       bFitTree->movePointerToTheNextChild();
       RefCountedKinematicParticle d0CandMC = bFitTree->currentParticle();

       if(!batPionCand->currentState().isValid() || !d0CandMC->currentState().isValid()) continue;

       // get batchlor pion and D0 parameters from B fit
       KinematicParameters batPionKP = batPionCand->currentState().kinematicParameters();
       KinematicParameters d0CandKP = d0CandMC->currentState().kinematicParameters();

       GlobalVector bTotalP = GlobalVector (bCand->currentState().globalMomentum().x(),
                                            bCand->currentState().globalMomentum().y(),
                                            bCand->currentState().globalMomentum().z());

       GlobalVector batPionTotalP = GlobalVector(batPionKP.momentum().x(),batPionKP.momentum().y(),batPionKP.momentum().z());
       GlobalVector d0TotalP = GlobalVector(d0CandKP.momentum().x(),d0CandKP.momentum().y(),d0CandKP.momentum().z());

       double batPionTotalE = sqrt( batPionTotalP.mag2() + piMassBSquared );
       double d0TotalE = sqrt( d0TotalP.mag2() + d0MassB*d0MassB );
       double bTotalE = batPionTotalE + d0TotalE;

       const Particle::LorentzVector bP4(bTotalP.x(), bTotalP.y(), bTotalP.z(), bTotalE);

       Particle::Point bVtx((*bDecayVertex).position().x(), (*bDecayVertex).position().y(), (*bDecayVertex).position().z());
       std::vector<double> bVtxEVec;
       bVtxEVec.push_back( d0_vFit_vertex->error().cxx() );
       bVtxEVec.push_back( d0_vFit_vertex->error().cyx() );
       bVtxEVec.push_back( d0_vFit_vertex->error().cyy() );
       bVtxEVec.push_back( d0_vFit_vertex->error().czx() );
       bVtxEVec.push_back( d0_vFit_vertex->error().czy() );
       bVtxEVec.push_back( d0_vFit_vertex->error().czz() );
       SMatrixSym3D bVtxCovMatrix(bVtxEVec.begin(), bVtxEVec.end());
       const Vertex::CovarianceMatrix bVtxCov(bVtxCovMatrix);
       double bVtxChi2(bDecayVertex->chiSquared());
       double bVtxNdof(bDecayVertex->degreesOfFreedom());
       double bNormalizedChi2 = bVtxChi2/bVtxNdof;

       double bRVtxMag = 99999.0;
       double bLVtxMag = 99999.0;
       double bSigmaRvtxMag = 999.0;
       double bSigmaLvtxMag = 999.0;

       GlobalVector bLineOfFlight = GlobalVector (bVtx.x() - xVtx,
                                                    bVtx.y() - yVtx,
                                                    bVtx.z() - zVtx);

       SMatrixSym3D bTotalCov;
       if(isVtxPV) bTotalCov = bVtxCovMatrix + vtxPrimary->covariance();
       else bTotalCov = bVtxCovMatrix + theBeamSpotHandle->rotatedCovariance3D();

       SVector3 distanceVector3D(bLineOfFlight.x(), bLineOfFlight.y(), bLineOfFlight.z());
       SVector3 distanceVector2D(bLineOfFlight.x(), bLineOfFlight.y(), 0.0);

       double bAngle3D = angle(bLineOfFlight.x(), bLineOfFlight.y(), bLineOfFlight.z(),
                         bTotalP.x(), bTotalP.y(), bTotalP.z());
       double bAngle2D = angle(bLineOfFlight.x(), bLineOfFlight.y(), (float)0.0,
                        bTotalP.x(), bTotalP.y(), (float)0.0);
       bLVtxMag = bLineOfFlight.mag();
       bRVtxMag = bLineOfFlight.perp();
       bSigmaLvtxMag = sqrt(ROOT::Math::Similarity(bTotalCov, distanceVector3D)) / bLVtxMag;
       bSigmaRvtxMag = sqrt(ROOT::Math::Similarity(bTotalCov, distanceVector2D)) / bRVtxMag;

       if( bNormalizedChi2 > bVtxChi2Cut ||
           bRVtxMag < bRVtxCut ||
           bRVtxMag / bSigmaRvtxMag < bRVtxSigCut ||
           bLVtxMag < bLVtxCut ||
           bLVtxMag / bSigmaLvtxMag < bLVtxSigCut ||
           cos(bAngle3D) < bCollinCut3D || cos(bAngle2D) < bCollinCut2D || bAngle3D > bAlphaCut || bAngle2D > bAlpha2DCut
       ) continue;

       RecoChargedCandidate PionCand(theTrackRefs[trdx]->charge(), Particle::LorentzVector(batPionTotalP.x(),
                                               batPionTotalP.y(), batPionTotalP.z(),
                                               batPionTotalE), bVtx);
       PionCand.setTrack(theTrackRefs[trdx]);

       AddFourMomenta addp4;

       VertexCompositeCandidate* theB = 0;
       theB = new VertexCompositeCandidate(0, bP4, bVtx, bVtxCov, bVtxChi2, bVtxNdof);
       theB->addDaughter(theD0);
       theB->addDaughter(PionCand);

       if(theD0.pdgId()<0) {theB->setPdgId(521); theB->setCharge(theTrackRefs[trdx]->charge());}
       else if(theD0.pdgId()>0) {theB->setPdgId(-521); theB->setCharge(theTrackRefs[trdx]->charge());}

       addp4.set( *theB );

       if( theB->mass() < bMassB + bMassCut &&
           theB->mass() > bMassB - bMassCut ) theBs.push_back( *theB );
       if(theB) delete theB;
          theB = 0;
    }
  }
}
// Get methods

const reco::VertexCompositeCandidateCollection& BFitter::getB() const {
  return theBs;
}

/*
const std::vector<float>& BFitter::getMVAVals() const {
  return mvaVals_;
}
*/
/*
auto_ptr<edm::ValueMap<float> > BFitter::getMVAMap() const {
  return mvaValValueMap;
}
*/

void BFitter::resetAll() {
    theBs.clear();
//    mvaVals_.clear();
}

#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description='Do a test')
parser.add_argument('-i', dest='inFile', type=str, nargs=1, help='input file name')
parser.add_argument('-o', dest='outDir', type=str, nargs=1, help='output directory')
args = parser.parse_args()

outDir = args.outDir[0]
if '/' != outDir[-1]:
    outDir = outDir+'/'

import ROOT as r
f = r.TFile(args.inFile[0], "READ")
t = f.Get("kshortana/ParticleTree")
tp = t.GetPlayer()

tp.SetScanRedirect(r.kTRUE)
tp.SetScanFileName(outDir+"cand.txt")
t.Scan("cand_charge:cand_mass:cand_pT:cand_eta:cand_phi:cand_charge:cand_decayLength3D:cand_decayLength2D:cand_angle3D:cand_angle2D:cand_decayLengthError3D:cand_decayLengthError2D:cand_pseudoDecayLengthError3D:cand_pseudoDecayLengthError2D:cand_dca:cand_massDau:cand_pTDau:cand_etaDau:cand_phiDau:gen_pdgId", "cand_status == 3")
tp.SetScanRedirect(r.kTRUE)
tp.SetScanFileName(outDir+"dau.txt")
t.Scan("cand_charge:cand_mass:cand_pT:cand_eta:cand_phi", "cand_status == 1")
tp.SetScanRedirect(r.kTRUE)
tp.SetScanFileName(outDir+"trk.txt")
t.Scan("trk_xyDCASignificance:trk_zDCASignificance:trk_nHit:trk_isHP", "")

tp.SetScanRedirect(r.kTRUE)
tp.SetScanFileName(outDir+"trkIdx.txt")
t.Scan("trk_candIdx", "")

#include "jetreader/lib/assert.h"
#include "jetreader/reader/reader.h"
#include "jetreader/lib/parse_csv.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // for getenv, atof, atoi
#include <vector>
#include <assert.h>
#include <cmath>
#include <climits>
#include <sstream>
#include <numeric>
#include <random>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "jetreader/reader/centrality_def.h"
#include "jetreader/reader/centrality.h"
#include "jetreader/reader/bemc_helper.h"
#include "jetreader/reader/track_selector.h"

using namespace std;
using namespace fastjet;
using namespace jetreader;

int main(int argc, const char **argv)
{

    Int_t nEvents = atoi(argv[1]); //! total number of events to run over
    Int_t nFiles = atoi(argv[2]);
    string InputFileList = argv[3];
    TString JobIdName = argv[4];
    TString OutputFileName = "/star/u/youqi/results/";
    OutputFileName.Append(JobIdName);
    OutputFileName.Append(".root");
    
    double trackDCACut = 1.0;       //! primary track DCA cut
    double trackNhitCut = 20;       //! primary track number of hits in tpc
    double trackNhitFracCut = 0.52; //! primary track nhits fraction to remove split tracks

    bool printDebug = false;
    bool removebadruns = true;
    bool removebadtowers = false;

    // set up input and output files
    cout << "Input file: " << InputFileList << endl;
    cout << "Output file: " << OutputFileName << endl;

    jetreader::Reader reader(InputFileList);
    TFile out(OutputFileName, "RECREATE");
    
    // setup events
    reader.eventSelector()->setVzRange(-30, 30);
    reader.eventSelector()->setVrMax(0.5);
    reader.eventSelector()->addTriggerId(480004); // VPDMB-novtx
    reader.eventSelector()->addTriggerId(480007); // VPDMB-5-trgonly
    reader.eventSelector()->addTriggerId(480201); // BHT0*VPDMB-5
    reader.eventSelector()->addTriggerId(480204); // BHT1*BBCMB
    reader.eventSelector()->addTriggerId(480205); // BHT2*BBCMB
    std::string badrunlist = "";
    if (removebadruns)
    {
        badrunlist = "bad.csv";
        reader.eventSelector()->addBadRuns(badrunlist);
    }

    // setup tracks
    reader.trackSelector()->setPtMax(30.0);
    reader.trackSelector()->setPtMin(2.0);
    reader.trackSelector()->rejectEventOnPtFailure(true);
    reader.trackSelector()->setDcaMax(trackDCACut);
    reader.trackSelector()->setNHitsMin(trackNhitCut);
    reader.trackSelector()->setNHitsFracMin(trackNhitFracCut);

    // set up TH1F and TTree for the output file
    TH1::SetDefaultSumw2(true);
    TH2::SetDefaultSumw2(true);

    long event = 0;

    Double_t vz = 0.0, vy = 0.0, vx = 0.0, vzvpd = 0.0, vzdiff = 0.0;
    int refmult = 0, btofmult = 0, ntrack = 0;

    TH1F *h_ntrack = new TH1F("h_ntrack", "number of global tracks", 3000, 0, 3000);
    TH1F *h_refmult = new TH1F("h_refmult", "refmult", 100, 0, 100);
    TH1F *h_btofmult = new TH1F("h_btofmult", "bTOF mult", 1000, 0, 1000);
    TH1F *h_vz = new TH1F("h_vz", "v_{z} [cm]", 400, -200, 200);
    TH1F *h_vy = new TH1F("h_vy", "v_{y} [cm]", 10, -5, 5);
    TH1F *h_vx = new TH1F("h_vx", "v_{x} [cm]", 10, -5, 5);
    TH1F *h_vzvpd = new TH1F("h_vzvpd", "v_{z,VPD} [cm]", 400, -200, 200);
    TH1F *h_vzdiff = new TH1F("h_vzdiff", "v_{z} - v_{z,VPD} [cm]", 50, -25, 25);
    TH2F *h_refmult_btofmult = new TH2F("h_refmult_btofmult", ";refmult;bTOF mult", 6000, 0, 6000, 6000, 0, 6000);

    Double_t pt, eta, phi, nspi, nsk, nsp, dedx, m2;
    Short_t q;

    TTree *t = new TTree("tracks", "tracks");
    t->Branch("pt", &pt, "pt/D");
    t->Branch("eta", &eta, "eta/D");
    t->Branch("phi", &phi, "phi/D");
    t->Branch("nspi", &nspi, "nspi/D");
    t->Branch("nsk", &nsk, "nsk/D");
    t->Branch("nsp", &nsp, "nsp/D");
    t->Branch("dedx", &dedx, "dedx/D");
    t->Branch("m2", &m2, "m2/D");
    t->Branch("q", &q, "q/S");

    out.ls(); // for debugging

    // initialize the reader
    reader.Init();
    cout << "number of events in chain: " << reader.tree()->GetEntries() << endl;
 
    while (reader.next())
    {
        event++;

        h_ntrack->Fill(reader.picoDst()->event()->numberOfGlobalTracks());
        h_refmult->Fill(reader.picoDst()->event()->refMult());
        h_btofmult->Fill(reader.picoDst()->event()->btofTrayMultiplicity());
        h_vz->Fill(reader.picoDst()->event()->primaryVertex().Z());
        h_vy->Fill(reader.picoDst()->event()->primaryVertex().Y());
        h_vx->Fill(reader.picoDst()->event()->primaryVertex().X());
        h_vzvpd->Fill(reader.picoDst()->event()->vzVpd());
        h_vzdiff->Fill(reader.picoDst()->event()->primaryVertex().Z() - reader.picoDst()->event()->vzVpd());
        h_refmult_btofmult->Fill(reader.picoDst()->event()->refMult(), reader.picoDst()->event()->btofTrayMultiplicity());

        // get StPicoTracks
        bool use_primary_tracks_ = true;

        TVector3 vertex = reader.picoDst()->event()->primaryVertex();

        for (int track_id = 0; track_id < reader.picoDst()->numberOfTracks(); ++track_id)
        {
            StPicoTrack *track = (StPicoTrack *)reader.picoDst()->track(track_id);
            jetreader::TrackStatus track_status = reader.trackSelector()->select(track, vertex, use_primary_tracks_);
            if (track_status != jetreader::TrackStatus::acceptTrack)
                continue;
            if (track->pMom().Eta() > 1. || track->pMom().Eta() < -1.)
                continue;
                
            pt = track->pPt();
            eta = track->pMom().Eta();
            phi = track->pMom().Phi();
            nspi = track->nSigmaPion();
            nsk = track->nSigmaKaon();
            nsp = track->nSigmaProton();
            dedx = track->dEdx();
            q = track->charge();
            m2 = -9.;
            Int_t bTofPidTraitsIndex = track->bTofPidTraitsIndex();
            if (bTofPidTraitsIndex >= 0)
            {
                StPicoBTofPidTraits *btofPidTraits = reader.picoDst()->btofPidTraits(bTofPidTraitsIndex);
                Float_t beta = btofPidTraits->btofBeta();
                m2 = track->pMom().Dot(track->pMom()) * (1.0 / (beta * beta) - 1.0);
            }
    
            t->Fill();

        } // track loop
    } // event loop
    out.cd();
    out.Write();
    out.ls(); // for debugging
    out.Close();
    cout << event << " good events out of " << reader.tree()->GetReadEntry() << endl;

    return 0;
}

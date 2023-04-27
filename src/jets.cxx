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
    reader.trackSelector()->rejectEventOnPtFailure(false);
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

    Double_t pt_j = 0.0, eta_j = 0.0, phi_j = 0.0, m_j = 0.0;
    Int_t event_j, mult_j = 0;

    TTree *t = new TTree("jets", "jets");
    t->Branch("pt_j", &pt_j, "pt_j/D");
    t->Branch("eta_j", &eta_j, "eta_j/D");
    t->Branch("phi_j", &phi_j, "phi_j/D");
    t->Branch("m_j", &m_j, "m_j/D");
    t->Branch("mult_j", &mult_j, "mult_j/I");

    Double_t pt[100], eta[100], phi[100], nspi[100], nsk[100], nsp[100], m2[100];
    Short_t q[100];

    t->Branch("pt", pt, "pt[mult_j]/D");
    t->Branch("eta", eta, "eta[mult_j]/D");
    t->Branch("phi", phi, "phi[mult_j]/D");
    t->Branch("nspi", nspi, "nspi[mult_j]/D");
    t->Branch("nsk", nsk, "nsk[mult_j]/D");
    t->Branch("nsp", nsp, "nsp[mult_j]/D");
    t->Branch("m2", m2, "m2[mult_j]/D");
    t->Branch("q", q, "q[mult_j]/S");

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
        vector<fastjet::PseudoJet> eventTracks;

        TVector3 vertex = reader.picoDst()->event()->primaryVertex();

        for (int track_id = 0; track_id < reader.picoDst()->numberOfTracks(); ++track_id)
        {
            StPicoTrack *track = (StPicoTrack *)reader.picoDst()->track(track_id);
            jetreader::TrackStatus track_status = reader.trackSelector()->select(track, vertex, use_primary_tracks_);
            if (track_status != jetreader::TrackStatus::acceptTrack)
                continue;
            fastjet::PseudoJet j;
            double px = track->pMom().Px();
            double py = track->pMom().Py();
            double pz = track->pMom().Pz();
            double E = TMath::Sqrt(track->pMom().Mag2() + 0.13957039 * 0.13957039);
            j.set_user_index(track_id);
            j.reset_momentum(px, py, pz, E);
            eventTracks.push_back(j);
        }

        // run jet clustering
        fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, 0.4);
        fastjet::Selector jetrap = fastjet::SelectorRapRange(-0.6, 0.6);

        fastjet::ClusterSequence cs(eventTracks, jetdef);
        vector<fastjet::PseudoJet> jets = sorted_by_pt(jetrap(cs.inclusive_jets(15.0)));

        for (unsigned int i = 0; i < jets.size(); i++)
        {
            event_j = event;
            pt_j = jets[i].perp();
            eta_j = jets[i].eta();
            phi_j = jets[i].phi();
            m_j = jets[i].m();
            mult_j = jets[i].constituents().size();

            // jet constituent loop
            for (unsigned int i_track = 0; i_track < jets[i].constituents().size(); i_track++)
            {
                pt[i_track] = 0.;
                eta[i_track] = 0.;
                phi[i_track] = 0.;
                nspi[i_track] = 0.;
                nsk[i_track] = 0.;
                nsp[i_track] = 0.;
                m2[i_track] = 0.;
                q[i_track] = 0.;

                int track_id = jets[i].constituents().at(i_track).user_index();
                StPicoTrack *track = (StPicoTrack *)reader.picoDst()->track(track_id);
                pt[i_track] = track->pPt();
                eta[i_track] = track->pMom().Eta();
                phi[i_track] = track->pMom().Phi();
                nspi[i_track] = track->nSigmaPion();
                nsk[i_track] = track->nSigmaKaon();
                nsp[i_track] = track->nSigmaProton();
                q[i_track] = track->charge();
                m2[i_track] = -9.;
                Int_t bTofPidTraitsIndex = track->bTofPidTraitsIndex();
                if (bTofPidTraitsIndex >= 0)
                {
                    StPicoBTofPidTraits *btofPidTraits = reader.picoDst()->btofPidTraits(bTofPidTraitsIndex);
                    Float_t beta = btofPidTraits->btofBeta();
                    m2[i_track] = track->pMom().Dot(track->pMom()) * (1.0 / (beta * beta) - 1.0);
                }
            } // jet constituent loop
            
            t->Fill();
            cout << "filled t" << endl;

        } // jet loop
    } // event loop
    out.cd();
    out.Write();
    out.ls(); // for debugging
    out.Close();
    cout << event << " good events out of " << reader.tree()->GetReadEntry() << endl;

    return 0;
}

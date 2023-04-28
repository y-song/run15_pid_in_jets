void analyze_tree()
{
    // create pT bins
    const int nbin = 16;
    const double pt_min = 2.0;
    const double pt_max = 10.0;
    const double pt_interval = 0.5;
    double pt_bin[nbin + 1];
    for (int ibin = 0; ibin < nbin + 1; ibin++)
    {
        pt_bin[ibin] = pt_min + ibin * pt_interval;
    }

    // set up histograms
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gStyle->SetOptStat(11);

    TObjArray trackm2_tracknspi;
    TObjArray tracknspi;
    TObjArray tracknspi_fits;

    for (int ibin = 0; ibin < nbin; ibin++)
    {
        tracknspi.Add(new TH1F(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), ";track nsigmapi", 50, -5, 5));
        trackm2_tracknspi.Add(new TH2F(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), ";track M^{2} [GeV^{2}];track nsigmapi", 50, 0, 2, 100, -5, 5));
        tracknspi_fits.Add(new TF1(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), "[0]/sqrt(2*TMath::Pi()*[1]^2)*exp((-0.5*x/[1])**2)+[2]/sqrt(2*TMath::Pi()*[3]^2)*exp((-0.5*(x-[3])/[4])**2)", -3, 3));
    }

    TH1F *h_trackm2 = new TH1F("h_trackm2", ";track M^{2} [GeV^{2}]", 50, -0.5, 1.5);
    TH1F *h_tracknspi = new TH1F("h_tracknspi", ";track nsigmapi", 50, -5, 5);
    TH2F *h_pt_nsknspidiff = new TH2F("h_pt_nsknspidiff", ";track p_{T} [GeV];track (nsigmaK - nsigmapi)", nbin, pt_min, pt_max, 60, 0, 3);
    TH2F *h_pt_m2 = new TH2F("h_pt_m2", ";track p_{T} [GeV];track M^{2} [GeV^{2}]", 2*nbin, 0, pt_max/2, 100, -0.5, 1.5);
    TH2F *h_trackm2_tracknspi = new TH2F("h_trackm2_tracknspi", ";track M^{2} [GeV^{2}];track nsigmapi", 75, -1, 2, 100, -5, 5);

    TFile *f = new TFile("result/jets_0428.root");
    TTree *t = (TTree *)f->Get("jets");

    Double_t pt_j = 0.0, eta_j = 0.0, phi_j = 0.0, m_j = 0.0;
    Int_t event_j, mult_j = 0;
    Double_t pt[100], eta[100], phi[100], nspi[100], nsk[100], nsp[100], m2[100];
    Short_t q[100];

    t->SetBranchAddress("pt_j", &pt_j);
    t->SetBranchAddress("eta_j", &eta_j);
    t->SetBranchAddress("phi_j", &phi_j);
    t->SetBranchAddress("m_j", &m_j);
    t->SetBranchAddress("mult_j", &mult_j);
    t->SetBranchAddress("pt", pt);
    t->SetBranchAddress("eta", eta);
    t->SetBranchAddress("phi", phi);
    t->SetBranchAddress("nspi", nspi);
    t->SetBranchAddress("nsk", nsk);
    t->SetBranchAddress("nsp", nsp);
    t->SetBranchAddress("m2", m2);
    t->SetBranchAddress("q", q);

    int njet = t->GetEntries();
    for (int ijet = 0; ijet < njet; ijet++)
    {
        t->GetEntry(ijet);

        for (int itrack = 0; itrack < mult_j; itrack++)
        {
            for (int ibin = 0; ibin < nbin; ibin++)
            {
                if (pt_bin[ibin] < pt[itrack] && pt[itrack] < pt_bin[ibin + 1])
                {
                    ((TH1F *)tracknspi.At(ibin))->Fill(nspi[itrack]);
                    ((TH2F *)trackm2_tracknspi.At(ibin))->Fill(m2[itrack], nspi[itrack]);
                }
            }
            h_trackm2->Fill(m2[itrack]);
            h_pt_nsknspidiff->Fill(pt[itrack], nsk[itrack] - nspi[itrack]);
            h_pt_m2->Fill(pt[itrack], m2[itrack]);
            h_trackm2_tracknspi->Fill(m2[itrack], nspi[itrack]);
        }
    }

    //TCanvas *c = new TCanvas("c", "c");

    //gPad->SetLogy(1);
    //h_trackm2->Draw();
    //h_pt_nsknspidiff->Draw("colz");
    //h_pt_nsknspidiff->ProfileX()->Draw("same");
    h_pt_m2->Draw("colz");
    //h_trackm2_tracknspi->Draw("colz");

    /*vector<Double_t> nsknspidiff;
    for (int ibin = 0; ibin < nbin; ibin++)
    {
        nsknspidiff.push_back( h_pt_nsknspidiff->ProfileX()->GetBinContent(ibin+1) );
    }

    for (int ibin = 0; ibin < nbin; ibin++)
    {
        dist = ((TH1F *)tracknspi.At(ibin));
        fit = ((TF1 *)tracknspi_fits.At(ibin));
        fit->SetParLimits(3, -nsknspidiff.at(ibin)-0.01, -nsknspidiff.at(ibin)+0.01);
        dist->Fit(fit, "rll");
        dist->Draw();
        //trackm2_tracknspi.At(ibin)->Draw("colz");
        c->SaveAs(Form("plots/trackm_tracknspi_%0.1f_pT_%0.1f.pdf", pt_bin[ibin], pt_bin[ibin + 1]), "pdf");
    }*/
}
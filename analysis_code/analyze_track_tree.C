void analyze_track_tree()
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
    TObjArray tracknsk;
    TObjArray tracknsp;
    TObjArray trackm2;
    TObjArray tracknspi_fits;

    for (int ibin = 0; ibin < nbin; ibin++)
    {
        tracknspi.Add(new TH1F(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), ";track nsigmapi", 50, -5, 5));
        tracknsk.Add(new TH1F(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), ";track nsigmak", 50, -5, 5));
        tracknsp.Add(new TH1F(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), ";track nsigmap", 50, -5, 5));
        trackm2.Add(new TH1F(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), ";track M^{2} [GeV^{2}]", 50, -0.5, 1.5));
        trackm2_tracknspi.Add(new TH2F(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), ";track M^{2} [GeV^{2}];track nsigmapi", 50, -0.5, 1.5, 100, -5, 5));
        tracknspi_fits.Add(new TF1(Form("%0.1f<pT<%0.1f GeV", pt_bin[ibin], pt_bin[ibin + 1]), "[0]/sqrt(2*TMath::Pi()*[1]^2)*exp((-0.5*x/[1])**2)+[2]/sqrt(2*TMath::Pi()*[3]^2)*exp((-0.5*(x-[3])/[4])**2)", -3, 3));
    }

    TH1F *h_trackm2 = new TH1F("h_trackm2", ";track M^{2} [GeV^{2}]", 50, -0.5, 1.5);
    TH1F *h_nspi = new TH1F("h_tracknspi", ";track nsigmapi", 50, -5, 5);
    TH1F *h_nsk = new TH1F("h_tracknsk", ";track nsigmak", 50, -5, 5);
    TH1F *h_nsp = new TH1F("h_tracknsp", ";track nsigmap", 50, -5, 5);

    TH2F *h_pt_nsknspidiff = new TH2F("h_pt_nsknspidiff", ";track p_{T} [GeV];track (nsigmaK - nsigmapi)", nbin, pt_min, pt_max, 60, 0, 3);
    TH2F *h_pt_m2 = new TH2F("h_pt_m2", ";track p_{T} [GeV];track M^{2} [GeV^{2}]", 2 * nbin, pt_min, pt_max / 2, 100, -0.5, 1.5);
    TH2F *h_trackm2_tracknspi = new TH2F("h_trackm2_tracknspi", ";track M^{2} [GeV^{2}];track nsigmapi", 75, -1, 2, 100, -5, 5);
    TH2F *h_pt_dedx = new TH2F("h_pt_dedx", ";p_{T} [GeV];dE/dx [keV/cm]", 300, 2, 17, 100, 0, 10);

    // TChain *t = new TChain("tracks");
    // t->Add("/gpfs01/star/pwg/youqi/run15/result/tracks_0428/E4*.root");
    TFile *f = new TFile("result/tracks_0428.root");
    TTree *t = (TTree *)f->Get("tracks");

    Double_t pt, eta, phi, nspi, nsk, nsp, dedx, m2;
    Short_t q;

    t->SetBranchAddress("pt", &pt);
    t->SetBranchAddress("eta", &eta);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("nspi", &nspi);
    t->SetBranchAddress("nsk", &nsk);
    t->SetBranchAddress("nsp", &nsp);
    t->SetBranchAddress("dedx", &dedx);
    t->SetBranchAddress("m2", &m2);
    t->SetBranchAddress("q", &q);

    int ntrack = t->GetEntries();
    for (int itrack = 0; itrack < ntrack; itrack++)
    {
        t->GetEntry(itrack);
        /*for (int ibin = 0; ibin < nbin; ibin++)
        {
            if (pt_bin[ibin] < pt && pt < pt_bin[ibin + 1])
            {
                ((TH1F *)tracknspi.At(ibin))->Fill(nspi);
                ((TH1F *)tracknsk.At(ibin))->Fill(nsk);
                ((TH1F *)tracknsp.At(ibin))->Fill(nsp);
                ((TH1F *)trackm2.At(ibin))->Fill(m2);
                ((TH2F *)trackm2_tracknspi.At(ibin))->Fill(m2, nspi);
            }
        }
        
        h_nspi->Fill(nspi);
        h_nsk->Fill(nsk);
        h_nsp->Fill(nsp);
        h_trackm2->Fill(m2);

        h_pt_nsknspidiff->Fill(pt, nsk - nspi);
        h_trackm2_tracknspi->Fill(m2, nspi);
        h_pt_dedx->Fill(pt, dedx);*/
        h_pt_m2->Fill(pt, m2);
        
    }

    TCanvas *c = new TCanvas("c", "c");
    gPad->SetLogz(1);
    h_pt_m2->Draw("colz");
    c->SaveAs("plots/pt_m2.png", "png");

    /*h_nspi->Draw();
    h_nsk->Draw();
    h_nsp->Draw();
    h_trackm2->Draw();
    h_pt_nsknspidiff->Draw("colz");
    h_pt_nsknspidiff->ProfileX()->Draw("same");
    h_trackm2_tracknspi->Draw("colz");
    h_pt_m2->Draw("colz");

    vector<Double_t> nsknspidiff;
    for (int ibin = 0; ibin < nbin; ibin++)
    {
        nsknspidiff.push_back( h_pt_nsknspidiff->ProfileX()->GetBinContent(ibin+1) );
    }

    for (int ibin = 0; ibin < nbin; ibin++)
    {
        dist = ((TH2F *)trackm2.At(ibin));
        //fit = ((TF1 *)tracknspi_fits.At(ibin));
        //fit->SetParLimits(3, -nsknspidiff.at(ibin)-0.01, -nsknspidiff.at(ibin)+0.01);
        //dist->Fit(fit, "rll");
        dist->Draw();
        trackm2_tracknspi.At(ibin)->Draw("colz");
        c->SaveAs(Form("plots/nspi_m2_%0.1f_pT_%0.1f.png", pt_bin[ibin], pt_bin[ibin + 1]), "png");
    }

    for (int ibin = 0; ibin < nbin; ibin++)
    {
        dist = ((TH1F *)tracknsp.At(ibin));
        dist->Draw();
        c->SaveAs(Form("plots/tracknsp_%0.1f_pT_%0.1f.png", pt_bin[ibin], pt_bin[ibin + 1]), "png");
    }*/
}
#include <fstream>

void photon_luminosity() {
    TCanvas *c1 = new TCanvas("c1", "", 900, 900);
    c1->SetGrid();

    Double_t z[100], lumi[100], lumi_p[100], lumi_m[100];
    std::ifstream file("../data/photon_luminosity.dat");
    for (int i = 0; i < 100; ++i) {
        file >> z[i] >> lumi[i] >> lumi_p[i] >> lumi_m[i];
    }
    file.close();

    TGraph *gr = new TGraph(100, z, lumi);
    gr->SetLineWidth(3);
    gr->SetLineColor(kRed);
    gr->SetLineStyle(10);
    gr->SetTitle("");
    gr->GetXaxis()->SetTitle("z = #sqrt{#tau}");
    gr->GetYaxis()->SetTitle("1/L_{ee} dL_{#gamma #gamma }^{++}/dz");
    gr->GetXaxis()->SetRangeUser(0.0, 1.0);
    gr->GetYaxis()->SetRangeUser(0.0, 10.5);
    gr->Draw("ALP");

    TGraph *gr_p = new TGraph(100, z, lumi_p);
    gr_p->SetLineWidth(3);
    gr_p->Draw("LP same");

    TGraph *gr_m = new TGraph(100, z, lumi_m);
    gr_m->SetLineWidth(3);
    gr_m->SetLineColor(kRed);
    gr_m->SetLineStyle(9);
    gr_m->Draw("LP same");

    c1->Update();
    c1->GetFrame()->SetBorderSize(12);
    c1->Modified();

    c1->SaveAs("photon_luminosity_pp.pdf");
}

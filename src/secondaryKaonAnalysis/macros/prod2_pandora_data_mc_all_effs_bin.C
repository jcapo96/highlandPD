void prod2_pandora_data_mc_all_effs_bin()
{
//=========Macro generated from canvas: bigcan/
//=========  (Thu Feb  3 10:54:12 2022) by ROOT version 6.16/00
   TCanvas *bigcan = new TCanvas("bigcan", "",0,45,1600,1000);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   bigcan->Range(0,0,1,1);
   bigcan->SetFillColor(10);
   bigcan->SetBorderMode(0);
   bigcan->SetBorderSize(0);
   bigcan->SetTicky(1);
   bigcan->SetTopMargin(0.075);
   bigcan->SetBottomMargin(0.125);
   bigcan->SetFrameBorderMode(0);
   bigcan->SetFrameBorderSize(0);
  
// ------------>Primitives in pad: padPion
   TPad *padPion = new TPad("padPion", "",0,0.5,0.5,1);
   padPion->Draw();
   padPion->cd();
   padPion->Range(0,0,1,1);
   padPion->SetFillColor(10);
   padPion->SetBorderMode(0);
   padPion->SetBorderSize(2);
   padPion->SetTicky(1);
   padPion->SetTopMargin(0.075);
   padPion->SetBottomMargin(0.125);
   padPion->SetFrameBorderMode(0);
   padPion->SetFrameBorderSize(0);
   
   TH2D *dummyLo__1 = new TH2D("dummyLo__1","",3,0,3,10,0,1.2);
   dummyLo__1->SetLineWidth(3);
   dummyLo__1->SetMarkerStyle(20);
   dummyLo__1->SetMarkerSize(0.9);
   dummyLo__1->GetXaxis()->SetTitle("Beam Momentum [GeV/c]");
   dummyLo__1->GetXaxis()->SetBinLabel(1,"1");
   dummyLo__1->GetXaxis()->SetBinLabel(2,"2");
   dummyLo__1->GetXaxis()->SetBinLabel(3,"3");
   dummyLo__1->GetXaxis()->CenterTitle(true);
   dummyLo__1->GetXaxis()->SetNdivisions(506);
   dummyLo__1->GetXaxis()->SetLabelFont(42);
   dummyLo__1->GetXaxis()->SetLabelOffset(0.007);
   dummyLo__1->GetXaxis()->SetLabelSize(0.07);
   dummyLo__1->GetXaxis()->SetTitleSize(0.05);
   dummyLo__1->GetXaxis()->SetTitleOffset(1);
   dummyLo__1->GetXaxis()->SetTitleFont(42);
   dummyLo__1->GetYaxis()->SetTitle(" Beam Particle Identification Efficiency");
   dummyLo__1->GetYaxis()->CenterTitle(true);
   dummyLo__1->GetYaxis()->SetNdivisions(506);
   dummyLo__1->GetYaxis()->SetLabelFont(42);
   dummyLo__1->GetYaxis()->SetLabelOffset(0.015);
   dummyLo__1->GetYaxis()->SetLabelSize(0.045);
   dummyLo__1->GetYaxis()->SetTitleSize(0.05);
   dummyLo__1->GetYaxis()->SetTitleOffset(1);
   dummyLo__1->GetYaxis()->SetTitleFont(42);
   dummyLo__1->GetZaxis()->SetLabelFont(42);
   dummyLo__1->GetZaxis()->SetLabelOffset(0.01);
   dummyLo__1->GetZaxis()->SetLabelSize(0.035);
   dummyLo__1->GetZaxis()->SetTitleSize(0.05);
   dummyLo__1->GetZaxis()->SetTitleOffset(1);
   dummyLo__1->GetZaxis()->SetTitleFont(42);
   dummyLo__1->Draw("");
   
   Double_t Graph0_fx3001[3] = {
   0.5,
   1.5,
   2.5};
   Double_t Graph0_fy3001[3] = {
   0.9620662,
   0.9769737,
   0.9857776};
   Double_t Graph0_felx3001[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph0_fely3001[3] = {
   0.01456916,
   0.008971519,
   0.002639887};
   Double_t Graph0_fehx3001[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph0_fehy3001[3] = {
   0.0150287,
   0.01126117,
   0.002937306};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(3,Graph0_fx3001,Graph0_fy3001,Graph0_felx3001,Graph0_fehx3001,Graph0_fely3001,Graph0_fehy3001);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff9999");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",2");
   
   Double_t Graph1_fx3002[3] = {
   0.5,
   1.5,
   2.5};
   Double_t Graph1_fy3002[3] = {
   0.9620662,
   0.9769737,
   0.9857776};
   Double_t Graph1_felx3002[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph1_fely3002[3] = {
   0.005436248,
   0.004949059,
   0.001463713};
   Double_t Graph1_fehx3002[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph1_fehy3002[3] = {
   0.006237267,
   0.006109992,
   0.001620001};
   grae = new TGraphAsymmErrors(3,Graph1_fx3002,Graph1_fy3002,Graph1_felx3002,Graph1_fehx3002,Graph1_fely3002,Graph1_fehy3002);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff0000");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",2");
   
   Double_t Graph2_fx3003[3] = {
   0.5,
   1.5,
   2.5};
   Double_t Graph2_fy3003[3] = {
   0.9324269,
   0.9407643,
   0.9482006};
   Double_t Graph2_felx3003[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph2_fely3003[3] = {
   0.002283642,
   0.001735979,
   0.001291094};
   Double_t Graph2_fehx3003[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph2_fehy3003[3] = {
   0.002356208,
   0.001784253,
   0.001321884};
   grae = new TGraphAsymmErrors(3,Graph2_fx3003,Graph2_fy3003,Graph2_felx3003,Graph2_fehx3003,Graph2_fely3003,Graph2_fehy3003);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",p");
   
   TH2D *dummyLo__2 = new TH2D("dummyLo__2","",3,0,3,10,0,1.2);
   dummyLo__2->SetLineWidth(3);
   dummyLo__2->SetMarkerStyle(20);
   dummyLo__2->SetMarkerSize(0.9);
   dummyLo__2->GetXaxis()->SetTitle("Beam Momentum [GeV/c]");
   dummyLo__2->GetXaxis()->SetBinLabel(1,"1");
   dummyLo__2->GetXaxis()->SetBinLabel(2,"2");
   dummyLo__2->GetXaxis()->SetBinLabel(3,"3");
   dummyLo__2->GetXaxis()->CenterTitle(true);
   dummyLo__2->GetXaxis()->SetNdivisions(506);
   dummyLo__2->GetXaxis()->SetLabelFont(42);
   dummyLo__2->GetXaxis()->SetLabelOffset(0.007);
   dummyLo__2->GetXaxis()->SetLabelSize(0.07);
   dummyLo__2->GetXaxis()->SetTitleSize(0.05);
   dummyLo__2->GetXaxis()->SetTitleOffset(1);
   dummyLo__2->GetXaxis()->SetTitleFont(42);
   dummyLo__2->GetYaxis()->SetTitle(" Beam Particle Identification Efficiency");
   dummyLo__2->GetYaxis()->CenterTitle(true);
   dummyLo__2->GetYaxis()->SetNdivisions(506);
   dummyLo__2->GetYaxis()->SetLabelFont(42);
   dummyLo__2->GetYaxis()->SetLabelOffset(0.015);
   dummyLo__2->GetYaxis()->SetLabelSize(0.045);
   dummyLo__2->GetYaxis()->SetTitleSize(0.05);
   dummyLo__2->GetYaxis()->SetTitleOffset(1);
   dummyLo__2->GetYaxis()->SetTitleFont(42);
   dummyLo__2->GetZaxis()->SetLabelFont(42);
   dummyLo__2->GetZaxis()->SetLabelOffset(0.01);
   dummyLo__2->GetZaxis()->SetLabelSize(0.035);
   dummyLo__2->GetZaxis()->SetTitleSize(0.05);
   dummyLo__2->GetZaxis()->SetTitleOffset(1);
   dummyLo__2->GetZaxis()->SetTitleFont(42);
   dummyLo__2->Draw("sameaxis");
   
   TLegend *leg = new TLegend(0,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph0","Simulation","lf");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph2","Data","lp");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
   TLatex *   tex = new TLatex(0.1,0.94,"#bf{DUNE:ProtoDUNE-SP}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.9,0.94,"#pi^{+}/#mu^{+} Beam");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   padPion->Modified();
   bigcan->cd();
  
// ------------>Primitives in pad: padPosi
   TPad *padPosi = new TPad("padPosi", "",0.5,0.5,1,1);
   padPosi->Draw();
   padPosi->cd();
   padPosi->Range(0,0,1,1);
   padPosi->SetFillColor(10);
   padPosi->SetBorderMode(0);
   padPosi->SetBorderSize(2);
   padPosi->SetTicky(1);
   padPosi->SetTopMargin(0.075);
   padPosi->SetBottomMargin(0.125);
   padPosi->SetFrameBorderMode(0);
   padPosi->SetFrameBorderSize(0);
   
   TH2D *dummyLo__3 = new TH2D("dummyLo__3","",3,0,3,10,0,1.2);
   dummyLo__3->SetLineWidth(3);
   dummyLo__3->SetMarkerStyle(20);
   dummyLo__3->SetMarkerSize(0.9);
   dummyLo__3->GetXaxis()->SetTitle("Beam Momentum [GeV/c]");
   dummyLo__3->GetXaxis()->SetBinLabel(1,"1");
   dummyLo__3->GetXaxis()->SetBinLabel(2,"2");
   dummyLo__3->GetXaxis()->SetBinLabel(3,"3");
   dummyLo__3->GetXaxis()->CenterTitle(true);
   dummyLo__3->GetXaxis()->SetNdivisions(506);
   dummyLo__3->GetXaxis()->SetLabelFont(42);
   dummyLo__3->GetXaxis()->SetLabelOffset(0.007);
   dummyLo__3->GetXaxis()->SetLabelSize(0.07);
   dummyLo__3->GetXaxis()->SetTitleSize(0.05);
   dummyLo__3->GetXaxis()->SetTitleOffset(1);
   dummyLo__3->GetXaxis()->SetTitleFont(42);
   dummyLo__3->GetYaxis()->SetTitle(" Beam Particle Identification Efficiency");
   dummyLo__3->GetYaxis()->CenterTitle(true);
   dummyLo__3->GetYaxis()->SetNdivisions(506);
   dummyLo__3->GetYaxis()->SetLabelFont(42);
   dummyLo__3->GetYaxis()->SetLabelOffset(0.015);
   dummyLo__3->GetYaxis()->SetLabelSize(0.045);
   dummyLo__3->GetYaxis()->SetTitleSize(0.05);
   dummyLo__3->GetYaxis()->SetTitleOffset(1);
   dummyLo__3->GetYaxis()->SetTitleFont(42);
   dummyLo__3->GetZaxis()->SetLabelFont(42);
   dummyLo__3->GetZaxis()->SetLabelOffset(0.01);
   dummyLo__3->GetZaxis()->SetLabelSize(0.035);
   dummyLo__3->GetZaxis()->SetTitleSize(0.05);
   dummyLo__3->GetZaxis()->SetTitleOffset(1);
   dummyLo__3->GetZaxis()->SetTitleFont(42);
   dummyLo__3->Draw("");
   
   Double_t Graph0_fx3004[3] = {
   0.5,
   1.5,
   2.5};
   Double_t Graph0_fy3004[3] = {
   0.9611998,
   0.9776786,
   0.9737073};
   Double_t Graph0_felx3004[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph0_fely3004[3] = {
   0.02218954,
   0.009397987,
   0.008447128};
   Double_t Graph0_fehx3004[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph0_fehy3004[3] = {
   0.02219984,
   0.01088098,
   0.01016761};
   grae = new TGraphAsymmErrors(3,Graph0_fx3004,Graph0_fy3004,Graph0_felx3004,Graph0_fehx3004,Graph0_fely3004,Graph0_fehy3004);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff9999");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(10);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",2");
   
   Double_t Graph1_fx3005[3] = {
   0.5,
   1.5,
   2.5};
   Double_t Graph1_fy3005[3] = {
   0.9611998,
   0.9776786,
   0.9737073};
   Double_t Graph1_felx3005[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph1_fely3005[3] = {
   0.00236324,
   0.004915584,
   0.004731183};
   Double_t Graph1_fehx3005[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph1_fehy3005[3] = {
   0.002505236,
   0.006101834,
   0.005637967};
   grae = new TGraphAsymmErrors(3,Graph1_fx3005,Graph1_fy3005,Graph1_felx3005,Graph1_fehx3005,Graph1_fely3005,Graph1_fehy3005);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff0000");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",2");
   
   Double_t Graph2_fx3006[3] = {
   0.5,
   1.5,
   2.5};
   Double_t Graph2_fy3006[3] = {
   0.9205632,
   0.9615006,
   0.9668118};
   Double_t Graph2_felx3006[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph2_fely3006[3] = {
   0.001802761,
   0.001616257,
   0.002092022};
   Double_t Graph2_fehx3006[3] = {
   0.5,
   0.5,
   0.5};
   Double_t Graph2_fehy3006[3] = {
   0.001840443,
   0.00168268,
   0.002223188};
   grae = new TGraphAsymmErrors(3,Graph2_fx3006,Graph2_fy3006,Graph2_felx3006,Graph2_fehx3006,Graph2_fely3006,Graph2_fehy3006);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",p");
   
   TH2D *dummyLo__4 = new TH2D("dummyLo__4","",3,0,3,10,0,1.2);
   dummyLo__4->SetLineWidth(3);
   dummyLo__4->SetMarkerStyle(20);
   dummyLo__4->SetMarkerSize(0.9);
   dummyLo__4->GetXaxis()->SetTitle("Beam Momentum [GeV/c]");
   dummyLo__4->GetXaxis()->SetBinLabel(1,"1");
   dummyLo__4->GetXaxis()->SetBinLabel(2,"2");
   dummyLo__4->GetXaxis()->SetBinLabel(3,"3");
   dummyLo__4->GetXaxis()->CenterTitle(true);
   dummyLo__4->GetXaxis()->SetNdivisions(506);
   dummyLo__4->GetXaxis()->SetLabelFont(42);
   dummyLo__4->GetXaxis()->SetLabelOffset(0.007);
   dummyLo__4->GetXaxis()->SetLabelSize(0.07);
   dummyLo__4->GetXaxis()->SetTitleSize(0.05);
   dummyLo__4->GetXaxis()->SetTitleOffset(1);
   dummyLo__4->GetXaxis()->SetTitleFont(42);
   dummyLo__4->GetYaxis()->SetTitle(" Beam Particle Identification Efficiency");
   dummyLo__4->GetYaxis()->CenterTitle(true);
   dummyLo__4->GetYaxis()->SetNdivisions(506);
   dummyLo__4->GetYaxis()->SetLabelFont(42);
   dummyLo__4->GetYaxis()->SetLabelOffset(0.015);
   dummyLo__4->GetYaxis()->SetLabelSize(0.045);
   dummyLo__4->GetYaxis()->SetTitleSize(0.05);
   dummyLo__4->GetYaxis()->SetTitleOffset(1);
   dummyLo__4->GetYaxis()->SetTitleFont(42);
   dummyLo__4->GetZaxis()->SetLabelFont(42);
   dummyLo__4->GetZaxis()->SetLabelOffset(0.01);
   dummyLo__4->GetZaxis()->SetLabelSize(0.035);
   dummyLo__4->GetZaxis()->SetTitleSize(0.05);
   dummyLo__4->GetZaxis()->SetTitleOffset(1);
   dummyLo__4->GetZaxis()->SetTitleFont(42);
   dummyLo__4->Draw("sameaxis");
   
   leg = new TLegend(0,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph0","Simulation","lf");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph2","Data","lp");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
      tex = new TLatex(0.1,0.94,"#bf{DUNE:ProtoDUNE-SP}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.9,0.94,"e^{+} Beam");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   padPosi->Modified();
   bigcan->cd();
  
// ------------>Primitives in pad: padProt
   TPad *padProt = new TPad("padProt", "",0,0,0.6,0.5);
   padProt->Draw();
   padProt->cd();
   padProt->Range(0,0,1,1);
   padProt->SetFillColor(10);
   padProt->SetBorderMode(0);
   padProt->SetBorderSize(2);
   padProt->SetTicky(1);
   padProt->SetLeftMargin(0.083);
   padProt->SetTopMargin(0.075);
   padProt->SetBottomMargin(0.125);
   padProt->SetFrameBorderMode(0);
   padProt->SetFrameBorderSize(0);
   
   TH2D *dummy__5 = new TH2D("dummy__5","",7,0,7,10,0,1.2);
   dummy__5->SetLineWidth(3);
   dummy__5->SetMarkerStyle(20);
   dummy__5->SetMarkerSize(0.9);
   dummy__5->GetXaxis()->SetTitle("Beam Momentum [GeV/c]");
   dummy__5->GetXaxis()->SetBinLabel(1,"1");
   dummy__5->GetXaxis()->SetBinLabel(2,"2");
   dummy__5->GetXaxis()->SetBinLabel(3,"3");
   dummy__5->GetXaxis()->SetBinLabel(4,"4");
   dummy__5->GetXaxis()->SetBinLabel(5,"5");
   dummy__5->GetXaxis()->SetBinLabel(6,"6");
   dummy__5->GetXaxis()->SetBinLabel(7,"7");
   dummy__5->GetXaxis()->CenterTitle(true);
   dummy__5->GetXaxis()->SetNdivisions(506);
   dummy__5->GetXaxis()->SetLabelFont(42);
   dummy__5->GetXaxis()->SetLabelOffset(0.007);
   dummy__5->GetXaxis()->SetLabelSize(0.07);
   dummy__5->GetXaxis()->SetTitleSize(0.05);
   dummy__5->GetXaxis()->SetTitleOffset(1);
   dummy__5->GetXaxis()->SetTitleFont(42);
   dummy__5->GetYaxis()->SetTitle(" Beam Particle Identification Efficiency");
   dummy__5->GetYaxis()->CenterTitle(true);
   dummy__5->GetYaxis()->SetNdivisions(506);
   dummy__5->GetYaxis()->SetLabelFont(42);
   dummy__5->GetYaxis()->SetLabelOffset(0.015);
   dummy__5->GetYaxis()->SetLabelSize(0.045);
   dummy__5->GetYaxis()->SetTitleSize(0.05);
   dummy__5->GetYaxis()->SetTitleOffset(0.83);
   dummy__5->GetYaxis()->SetTitleFont(42);
   dummy__5->GetZaxis()->SetLabelFont(42);
   dummy__5->GetZaxis()->SetLabelOffset(0.01);
   dummy__5->GetZaxis()->SetLabelSize(0.035);
   dummy__5->GetZaxis()->SetTitleSize(0.05);
   dummy__5->GetZaxis()->SetTitleOffset(1);
   dummy__5->GetZaxis()->SetTitleFont(42);
   dummy__5->Draw("");
   
   Double_t Graph0_fx3007[5] = {
   0.5,
   1.5,
   2.5,
   5.5,
   6.5};
   Double_t Graph0_fy3007[5] = {
   0.9103079,
   0.9677966,
   0.9790819,
   0.9796334,
   0.9710425};
   Double_t Graph0_felx3007[5] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t Graph0_fely3007[5] = {
   0.04513554,
   0.01615347,
   0.006126828,
   0.01664117,
   0.01567071};
   Double_t Graph0_fehx3007[5] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t Graph0_fehy3007[5] = {
   0.04552726,
   0.01785432,
   0.007147237,
   0.01726852,
   0.01806563};
   grae = new TGraphAsymmErrors(5,Graph0_fx3007,Graph0_fy3007,Graph0_felx3007,Graph0_fehx3007,Graph0_fely3007,Graph0_fehy3007);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff9999");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(10);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",2");
   
   Double_t Graph1_fx3008[5] = {
   0.5,
   1.5,
   2.5,
   5.5,
   6.5};
   Double_t Graph1_fy3008[5] = {
   0.9103079,
   0.9677966,
   0.9790819,
   0.9796334,
   0.9710425};
   Double_t Graph1_felx3008[5] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t Graph1_fely3008[5] = {
   0.007440986,
   0.007251598,
   0.003444698,
   0.004487829,
   0.007333139};
   Double_t Graph1_fehx3008[5] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t Graph1_fehy3008[5] = {
   0.008014892,
   0.009020668,
   0.00404822,
   0.005574515,
   0.009387836};
   grae = new TGraphAsymmErrors(5,Graph1_fx3008,Graph1_fy3008,Graph1_felx3008,Graph1_fehx3008,Graph1_fely3008,Graph1_fehy3008);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff0000");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",2");
   
   Double_t Graph2_fx3009[5] = {
   0.5,
   1.5,
   2.5,
   5.5,
   6.5};
   Double_t Graph2_fy3009[5] = {
   0.8908764,
   0.9430156,
   0.9590263,
   0.9625942,
   0.9639017};
   Double_t Graph2_felx3009[5] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t Graph2_fely3009[5] = {
   0.003218628,
   0.002884074,
   0.002083862,
   0.003118689,
   0.003662067};
   Double_t Graph2_fehx3009[5] = {
   0.5,
   0.5,
   0.5,
   0.5,
   0.5};
   Double_t Graph2_fehy3009[5] = {
   0.003302685,
   0.003024217,
   0.002187677,
   0.003378447,
   0.004037044};
   grae = new TGraphAsymmErrors(5,Graph2_fx3009,Graph2_fy3009,Graph2_felx3009,Graph2_fehx3009,Graph2_fely3009,Graph2_fehy3009);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",p");
   
   TH2D *dummy__6 = new TH2D("dummy__6","",7,0,7,10,0,1.2);
   dummy__6->SetLineWidth(3);
   dummy__6->SetMarkerStyle(20);
   dummy__6->SetMarkerSize(0.9);
   dummy__6->GetXaxis()->SetTitle("Beam Momentum [GeV/c]");
   dummy__6->GetXaxis()->SetBinLabel(1,"1");
   dummy__6->GetXaxis()->SetBinLabel(2,"2");
   dummy__6->GetXaxis()->SetBinLabel(3,"3");
   dummy__6->GetXaxis()->SetBinLabel(4,"4");
   dummy__6->GetXaxis()->SetBinLabel(5,"5");
   dummy__6->GetXaxis()->SetBinLabel(6,"6");
   dummy__6->GetXaxis()->SetBinLabel(7,"7");
   dummy__6->GetXaxis()->CenterTitle(true);
   dummy__6->GetXaxis()->SetNdivisions(506);
   dummy__6->GetXaxis()->SetLabelFont(42);
   dummy__6->GetXaxis()->SetLabelOffset(0.007);
   dummy__6->GetXaxis()->SetLabelSize(0.07);
   dummy__6->GetXaxis()->SetTitleSize(0.05);
   dummy__6->GetXaxis()->SetTitleOffset(1);
   dummy__6->GetXaxis()->SetTitleFont(42);
   dummy__6->GetYaxis()->SetTitle(" Beam Particle Identification Efficiency");
   dummy__6->GetYaxis()->CenterTitle(true);
   dummy__6->GetYaxis()->SetNdivisions(506);
   dummy__6->GetYaxis()->SetLabelFont(42);
   dummy__6->GetYaxis()->SetLabelOffset(0.015);
   dummy__6->GetYaxis()->SetLabelSize(0.045);
   dummy__6->GetYaxis()->SetTitleSize(0.05);
   dummy__6->GetYaxis()->SetTitleOffset(0.83);
   dummy__6->GetYaxis()->SetTitleFont(42);
   dummy__6->GetZaxis()->SetLabelFont(42);
   dummy__6->GetZaxis()->SetLabelOffset(0.01);
   dummy__6->GetZaxis()->SetLabelSize(0.035);
   dummy__6->GetZaxis()->SetTitleSize(0.05);
   dummy__6->GetZaxis()->SetTitleOffset(1);
   dummy__6->GetZaxis()->SetTitleFont(42);
   dummy__6->Draw("sameaxis");
   
   leg = new TLegend(0,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph0","Simulation","lf");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph2","Data","lp");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
      tex = new TLatex(0.083,0.94,"#bf{DUNE:ProtoDUNE-SP}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.9,0.94,"p Beam");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   padProt->Modified();
   bigcan->cd();
  
// ------------>Primitives in pad: padkaon
   TPad *padkaon = new TPad("padkaon", "",0.6,0,1,0.5);
   padkaon->Draw();
   padkaon->cd();
   padkaon->Range(0,0,1,1);
   padkaon->SetFillColor(10);
   padkaon->SetBorderMode(0);
   padkaon->SetBorderSize(2);
   padkaon->SetTicky(1);
   padkaon->SetLeftMargin(0.125);
   padkaon->SetRightMargin(0.125);
   padkaon->SetTopMargin(0.075);
   padkaon->SetBottomMargin(0.125);
   padkaon->SetFrameBorderMode(0);
   padkaon->SetFrameBorderSize(0);
   
   TH2D *dummyHi__7 = new TH2D("dummyHi__7","",2,5,7,10,0,1.2);
   dummyHi__7->SetLineWidth(3);
   dummyHi__7->SetMarkerStyle(20);
   dummyHi__7->SetMarkerSize(0.9);
   dummyHi__7->GetXaxis()->SetTitle("Beam Momentum [GeV/c]");
   dummyHi__7->GetXaxis()->SetBinLabel(1,"6");
   dummyHi__7->GetXaxis()->SetBinLabel(2,"7");
   dummyHi__7->GetXaxis()->CenterTitle(true);
   dummyHi__7->GetXaxis()->SetNdivisions(506);
   dummyHi__7->GetXaxis()->SetLabelFont(42);
   dummyHi__7->GetXaxis()->SetLabelOffset(0.007);
   dummyHi__7->GetXaxis()->SetLabelSize(0.07);
   dummyHi__7->GetXaxis()->SetTitleSize(0.05);
   dummyHi__7->GetXaxis()->SetTitleOffset(1);
   dummyHi__7->GetXaxis()->SetTitleFont(42);
   dummyHi__7->GetYaxis()->SetTitle(" Beam Particle Identification Efficiency");
   dummyHi__7->GetYaxis()->CenterTitle(true);
   dummyHi__7->GetYaxis()->SetNdivisions(506);
   dummyHi__7->GetYaxis()->SetLabelFont(42);
   dummyHi__7->GetYaxis()->SetLabelOffset(0.015);
   dummyHi__7->GetYaxis()->SetLabelSize(0.045);
   dummyHi__7->GetYaxis()->SetTitleSize(0.05);
   dummyHi__7->GetYaxis()->SetTitleOffset(1.2);
   dummyHi__7->GetYaxis()->SetTitleFont(42);
   dummyHi__7->GetZaxis()->SetLabelFont(42);
   dummyHi__7->GetZaxis()->SetLabelOffset(0.01);
   dummyHi__7->GetZaxis()->SetLabelSize(0.035);
   dummyHi__7->GetZaxis()->SetTitleSize(0.05);
   dummyHi__7->GetZaxis()->SetTitleOffset(1);
   dummyHi__7->GetZaxis()->SetTitleFont(42);
   dummyHi__7->Draw("");
   
   Double_t Graph0_fx3010[2] = {
   5.5,
   6.5};
   Double_t Graph0_fy3010[2] = {
   0.984252,
   0.9555556};
   Double_t Graph0_felx3010[2] = {
   0.5,
   0.5};
   Double_t Graph0_fely3010[2] = {
   0.03224714,
   0.0251206};
   Double_t Graph0_fehx3010[2] = {
   0.5,
   0.5};
   Double_t Graph0_fehy3010[2] = {
   0.03379985,
   0.03303224};
   grae = new TGraphAsymmErrors(2,Graph0_fx3010,Graph0_fy3010,Graph0_felx3010,Graph0_fehx3010,Graph0_fely3010,Graph0_fehy3010);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff9999");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(10);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",2");
   
   Double_t Graph1_fx3011[2] = {
   5.5,
   6.5};
   Double_t Graph1_fy3011[2] = {
   0.984252,
   0.9555556};
   Double_t Graph1_felx3011[2] = {
   0.5,
   0.5};
   Double_t Graph1_fely3011[2] = {
   0.006229219,
   0.0136705};
   Double_t Graph1_fehx3011[2] = {
   0.5,
   0.5};
   Double_t Graph1_fehy3011[2] = {
   0.009287584,
   0.01836106};
   grae = new TGraphAsymmErrors(2,Graph1_fx3011,Graph1_fy3011,Graph1_felx3011,Graph1_fehx3011,Graph1_fely3011,Graph1_fehy3011);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ff0000");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",2");
   
   Double_t Graph2_fx3012[2] = {
   5.5,
   6.5};
   Double_t Graph2_fy3012[2] = {
   0.9461988,
   0.9504804};
   Double_t Graph2_felx3012[2] = {
   0.5,
   0.5};
   Double_t Graph2_fely3012[2] = {
   0.005476122,
   0.005917735};
   Double_t Graph2_fehx3012[2] = {
   0.5,
   0.5};
   Double_t Graph2_fehy3012[2] = {
   0.006025787,
   0.00662508};
   grae = new TGraphAsymmErrors(2,Graph2_fx3012,Graph2_fy3012,Graph2_felx3012,Graph2_fehx3012,Graph2_fely3012,Graph2_fehy3012);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillStyle(1000);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0.9);
   grae->Draw(",p");
   
   TH2D *dummyHi__8 = new TH2D("dummyHi__8","",2,5,7,10,0,1.2);
   dummyHi__8->SetLineWidth(3);
   dummyHi__8->SetMarkerStyle(20);
   dummyHi__8->SetMarkerSize(0.9);
   dummyHi__8->GetXaxis()->SetTitle("Beam Momentum [GeV/c]");
   dummyHi__8->GetXaxis()->SetBinLabel(1,"6");
   dummyHi__8->GetXaxis()->SetBinLabel(2,"7");
   dummyHi__8->GetXaxis()->CenterTitle(true);
   dummyHi__8->GetXaxis()->SetNdivisions(506);
   dummyHi__8->GetXaxis()->SetLabelFont(42);
   dummyHi__8->GetXaxis()->SetLabelOffset(0.007);
   dummyHi__8->GetXaxis()->SetLabelSize(0.07);
   dummyHi__8->GetXaxis()->SetTitleSize(0.05);
   dummyHi__8->GetXaxis()->SetTitleOffset(1);
   dummyHi__8->GetXaxis()->SetTitleFont(42);
   dummyHi__8->GetYaxis()->SetTitle(" Beam Particle Identification Efficiency");
   dummyHi__8->GetYaxis()->CenterTitle(true);
   dummyHi__8->GetYaxis()->SetNdivisions(506);
   dummyHi__8->GetYaxis()->SetLabelFont(42);
   dummyHi__8->GetYaxis()->SetLabelOffset(0.015);
   dummyHi__8->GetYaxis()->SetLabelSize(0.045);
   dummyHi__8->GetYaxis()->SetTitleSize(0.05);
   dummyHi__8->GetYaxis()->SetTitleOffset(1.2);
   dummyHi__8->GetYaxis()->SetTitleFont(42);
   dummyHi__8->GetZaxis()->SetLabelFont(42);
   dummyHi__8->GetZaxis()->SetLabelOffset(0.01);
   dummyHi__8->GetZaxis()->SetLabelSize(0.035);
   dummyHi__8->GetZaxis()->SetTitleSize(0.05);
   dummyHi__8->GetZaxis()->SetTitleOffset(1);
   dummyHi__8->GetZaxis()->SetTitleFont(42);
   dummyHi__8->Draw("sameaxis");
   
   leg = new TLegend(0,0,0,0,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   entry=leg->AddEntry("Graph0","Simulation","lf");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph2","Data","lp");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   leg->Draw();
      tex = new TLatex(0.125,0.94,"#bf{DUNE:ProtoDUNE-SP}");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.875,0.94,"K^{+} Beam");
tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   padkaon->Modified();
   bigcan->cd();
   bigcan->Modified();
   bigcan->cd();
   bigcan->SetSelected(bigcan);
}

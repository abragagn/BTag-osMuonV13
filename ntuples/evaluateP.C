#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "Math/MinimizerOptions.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include "TEfficiency.h"

using namespace std;

TString process_;
TString dir_ = "./";
double min_, max_, x1_, x2_;
double mean_, sigma1_, sigma2_;
int nBins_=50;
pair<double, double> CountEventsWithFit(TH1 *hist);

void evaluateP(TString file = "./ntuBsMC2018.root",  TString cutEvt_ = "", TString cut_ = "")
{
    gErrorIgnoreLevel = kWarning;
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 20000 );
    gStyle->SetOptStat(10); //ksiourmen
    gStyle->SetOptFit(0111); //pcev
    TFile *f = new TFile(file);
    if(f->IsZombie()) return false;
    TTree *t = (TTree*)f->Get("PDsecondTree");
    if(file.Contains("Bs")){
        process_ = "BsJPsiPhi";
        min_ = 5.25;
        max_ = 5.60;
    }
    if(file.Contains("Bu")){
        process_ = "BuJPsiK";
        min_ = 5.10;
        max_ = 5.45;
        x2_ = 5.499;
        x1_= 5.18;
    }

    if(file.Contains("MC")) process_ = process_ + "MC";
    if(file.Contains("Data")) process_ = process_ + "Data";

    TH1F *ssB       = new TH1F( "ssB", "ssB", nBins_, min_, max_ );
    TH1F *ssB_NT    = new TH1F( "ssB_NT", "ssB_NT", nBins_, min_, max_ );
    TH1F *ssB_RT    = new TH1F( "ssB_RT", "ssB_RT", nBins_, min_, max_ );
    TH1F *ssB_WT    = new TH1F( "ssB_WT", "ssB_WT", nBins_, min_, max_ );

    TString cut = "1";
    TString cutEvt = "ssbIsTight";

    if(cutEvt_ != "") cutEvt = cutEvt + "&&" + cutEvt_;
    if(cut_ != "")    cut = cut + "&&" + cut_;

    TString base =  "evtWeight*((" + cutEvt;
    TString cutNT = base + ")&&!osMuon)";
    TString cutRT = base + ")&&(" + cut + ")&&osMuon&&osMuonTag==1)";
    TString cutWT = base + ")&&(" + cut + ")&&osMuon&&osMuonTag==0)";

    t->Project("ssB", "ssbMass", base + "))" );
    t->Project("ssB_RT", "ssbMass", cutRT );
    t->Project("ssB_WT", "ssbMass", cutWT );
    t->Project("ssB_NT", "ssbMass", cutNT );

    pair<double, double> nTot;
    pair<double, double> nNT;
    pair<double, double> nRT;
    pair<double, double> nWT;
    pair<double, double> eff;
    pair<double, double> w;
    pair<double, double> power;

    nTot.first = ssB->Integral();
    nRT.first = ssB_RT->Integral();
    nWT.first = ssB_WT->Integral();
    nNT.first = ssB_NT->Integral();

    nTot.second = sqrt(nTot.first);
    nRT.second = sqrt(nRT.first);
    nWT.second = sqrt(nWT.first);
    nNT.second = sqrt(nNT.first);

    if(process_.Contains("Data")){
        nTot = CountEventsWithFit(ssB);
        nRT  = CountEventsWithFit(ssB_RT);
        nWT  = CountEventsWithFit(ssB_WT);        
        nNT  = CountEventsWithFit(ssB_NT);        
    }

    double tagged = nRT.first+nWT.first;

    eff.first = tagged/(tagged+nNT.first);
    w.first = nWT.first/tagged;
    power.first = eff.first*pow((1-2*w.first), 2);

    // ERRORS
    if(process_.Contains("Data")){
        double RT2 = pow(nRT.first,2);
        double WT2 = pow(nWT.first,2);
        double NT2 = pow(nNT.first,2);
        double Tot2 = pow(nTot.first,2);

        double sRT2 = pow(nRT.second,2);
        double sWT2 = pow(nWT.second,2);
        double sNT2 = pow(nNT.second,2);
        double sTot2 = pow(nTot.second,2);

        eff.second = NT2*(sRT2+sWT2) + sNT2*pow(nRT.first+nWT.first,2);
        eff.second /= pow(nRT.first+nWT.first+nNT.first,4);
        eff.second = sqrt(eff.second);

        w.second = sqrt(WT2*sRT2 + RT2*sWT2)/pow(nRT.first+nWT.first,2);

        power.second = sWT2*pow(nRT.first-nWT.first,2)
                        *pow(4*RT2+4*nRT.first*nWT.first+3*nRT.first*nNT.first+nWT.first*nWT.first,2);
        power.second += sRT2*pow(nRT.first-nWT.first,2)
                        *pow(nRT.first*(4*nWT.first+nNT.first)+nWT.first*(4*nWT.first+3*nNT.first),2);
        power.second /= pow(nRT.first+nWT.first,4)*pow(nRT.first+nWT.first+nNT.first,4);
        power.second = sqrt(power.second);
    }else{
        double cl = 0.6827;
        eff.second = (TEfficiency::AgrestiCoull(nTot.first,tagged,cl,1)
                    -TEfficiency::AgrestiCoull(nTot.first,tagged,cl,0))/2;
        w.second = (TEfficiency::AgrestiCoull(tagged,nWT.first,cl,1)
                    -TEfficiency::AgrestiCoull(tagged,nWT.first,cl,0))/2;
        power.second = sqrt(16*pow(eff.first,2)*pow(1-2*w.first,2)*pow(w.second,2) 
                    + pow(1-2*w.first,4)*pow(eff.second,2));
    }

    cout<<endl;

    cout<<"Bs = "<<(int)nTot.first<<" +- "<<(int)nTot.second<<"   --   "<<100*nTot.second/nTot.first<<"  --  "<<ssB->Integral()<<"  --  "<<100*1/sqrt(nTot.first)<<endl;
    cout<<"NT = "<<(int)nNT.first<<" +- "<<(int)nNT.second<<"   --   "<<100*nNT.second/nNT.first<<"  --  "<<ssB_NT->Integral()<<"  --  "<<100*1/sqrt(nNT.first)<<endl;
    cout<<"RT = "<<(int)nRT.first<<" +- "<<(int)nRT.second<<"   --   "<<100*nRT.second/nRT.first<<"  --  "<<ssB_RT->Integral()<<"  --  "<<100*1/sqrt(nRT.first)<<endl;
    cout<<"WT = "<<(int)nWT.first<<" +- "<<(int)nWT.second<<"   --   "<<100*nWT.second/nWT.first<<"  --  "<<ssB_WT->Integral()<<"  --  "<<100*1/sqrt(nWT.first)<<endl;
    cout<<endl;
    cout<<"Eff = "<<100*eff.first<<" +- "<<100*eff.second<<" %"<<endl;
    cout<<"Mistag = "<<100*w.first<<" +- "<<100*w.second<<" %"<<endl;
    cout<<"Power = "<<100*power.first<<" +- "<<100*power.second<<" %"<<endl;

    cout<<endl;

    return;
}

pair<double, double> CountEventsWithFit(TH1 *hist)
{
    TString title = hist->GetTitle();
    cout<<" ---  now fitting "<<title<<endl;
    bool isTot = title == "ssB" ? true : false;

    TRandom3 *r3 = new TRandom3();
    double mean = 5.3663;
    double sigma = 0.015;
    if(process_.Contains("BsJPsiPhi")) mean = 5.3663;
    if(process_.Contains("BuJPsiK"))   mean = 5.2793;
    if(process_.Contains("BdJPsiKx"))  mean = 5.2796;

    TString sgnDef = "[1]*TMath::Gaus(x, [0], [2], true)";
    sgnDef +=       "+[3]*TMath::Gaus(x, [0], [4], true)";

    TString bkgDef = "[5]+[6]*x";
    bkgDef += "+[7]*TMath::Erfc([8]*(x-[9]))";

    TString funcDef = sgnDef + "+" + bkgDef;

    TF1 *func = new TF1("func", funcDef, min_, max_);

    func->SetParName(0, "mean");
    func->SetParName(1, "A1");
    func->SetParName(2, "sigma1");
    func->SetParName(3, "A2");
    func->SetParName(4, "sigma2");

    //SIGNAL
    double limit = hist->GetEntries()*hist->GetBinWidth(1);

    func->SetParameter(0, mean);
    func->SetParameter(1, limit/2*r3->Gaus(1.,0.01));
    func->SetParameter(3, limit/2*r3->Gaus(1.,0.01));
    func->SetParameter(2, sigma*r3->Gaus(1.,0.01));
    func->SetParameter(4, sigma*r3->Gaus(1.,0.01));
    func->SetParLimits(1, 0, 2*limit);
    func->SetParLimits(3, 0, 2*limit);
    func->SetParLimits(2, 0.001, 0.1);
    func->SetParLimits(4, 0.001, 0.1);

    //BKG
    TAxis *xaxis = hist->GetXaxis();
    int binx1 = xaxis->FindBin(x1_);
    int binx2 = xaxis->FindBin(x2_);
    double x2 = hist->GetBinContent(binx2);
    double x1 = hist->GetBinContent(binx1);
    double y2 = hist->GetBinCenter(binx2);
    double y1 = hist->GetBinCenter(binx1);
    double m = (y2-y1)/(x2-x1);
    double q = y1 - m*x1;
    func->SetParameter(5, 10);
    func->SetParameter(6, m);

    func->SetParameter(7, hist->GetBinContent(2)/2);
    func->SetParameter(8, 10);
    func->SetParameter(9, 5);
    func->SetParLimits(7, 0, hist->GetBinContent(2)*1.5);

    if(!isTot){
        func->FixParameter(0, mean_);
        func->FixParameter(2, sigma1_);
        func->FixParameter(4, sigma2_);
    }

    if(process_.Contains("MC")){
        func->FixParameter(5, 0);
        func->FixParameter(6, 0);
        func->FixParameter(7, 0);
        func->FixParameter(8, 0);
        func->FixParameter(9, 0);
    }

    func->SetNpx(2000);

    auto c5 = new TCanvas();
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(.75);
    TFitResultPtr r = hist->Fit("func","MRLS");
    hist->Draw("PE");
    hist->SetMinimum(0);
    TF1 *fit = hist->GetFunction("func");
    if(isTot){
        mean_ = fit->GetParameter(0);
        sigma1_ = fit->GetParameter(2);
        sigma2_ = fit->GetParameter(4);
    }

    TF1 *f1 = new TF1("f1","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f2 = new TF1("f2","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f4 = new TF1("f4","[0]+[1]*x+[2]*TMath::Erfc([3]*(x-[4]))", min_, max_);
    TF1 *f5 = new TF1("f5","[0]+[1]*x", min_, max_);   

    f1->SetParameters(fit->GetParameter(1),fit->GetParameter(0),fit->GetParameter(2));
    f2->SetParameters(fit->GetParameter(3),fit->GetParameter(0),fit->GetParameter(4));
    f4->SetParameters(fit->GetParameter(5),fit->GetParameter(6)
                    ,fit->GetParameter(7),fit->GetParameter(8),fit->GetParameter(9));
    f5->SetParameters(fit->GetParameter(5),fit->GetParameter(6));

    f1->SetLineColor(kBlue);
    f2->SetLineColor(kViolet);
    f4->SetLineColor(kGreen);
    f5->SetLineColor(kOrange);
    f1->SetLineStyle(2);
    f2->SetLineStyle(2);
    f4->SetLineStyle(2);
    f5->SetLineStyle(2);

    f1->Draw("same");
    f2->Draw("same");
    f4->Draw("same");
    f5->Draw("same");

    c5->Print(title + ".png");

    double nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(3);
    nEvt/=hist->GetBinWidth(1);

    TMatrixDSym cov = r->GetCovarianceMatrix();
    double errN = sqrt(cov(1,1)+cov(3,3)+2*cov(1,3))/hist->GetBinWidth(1);

    return make_pair(nEvt, errN);
}

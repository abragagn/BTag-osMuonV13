#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/PyMethodBase.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TMath.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "Math/MinimizerOptions.h"
#include "TFitResult.h"
#include "TPaveStats.h"

using namespace std;

TString process_;
TString dirPath_;
double min_, max_, x1_, x2_;
double mean_, sigma1_, sigma2_;
int nBins_ = 100;

pair<double, double> CountEventsWithFit(TH1 *hist);

int evaluateCalibration(TString file_ = "ntuples/ntuBuData2018.root"
    , int nEvents_ = -1
    , TString method_ = "DNNOsMuonHLTJpsiMu"
    , bool useTightSelection_ = true
    , float muonIDwp_ = 0.21
    )
{

    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(10); //ksiourmen
    gStyle->SetOptFit(1); //pcev
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );

    cout<<"----- Parameters -----"<<endl;
    cout<<"file_ = "<<file_<<endl;
    cout<<"method_ = "<<method_<<endl;
    cout<<"useTightSelection_ = "<<useTightSelection_<<endl;
    cout<<"nEvents_ = "<<nEvents_<<endl;

    cout<<endl<<"----- BEGIN CODE"<<endl;
    TString path = "/lustre/cmswork/abragagn/BPH/BTag/osMuonV13/src/PDAnalysis/Ntu/bin/";
    TString methodpath = "/lustre/cmswork/abragagn/mvaWeights/OsMuonTag/TMVAClassification_";

    auto *f = new TFile(path + file_);
    auto *t = (TTree*)f->Get("PDsecondTree");
    if(f->IsZombie()) return 0;

    bool isData = false;

    if(file_.Contains("Bs")){
        process_ = "BsJPsiPhi";
        min_ = 5.20;
        max_ = 5.65;
        dirPath_ = "./Bs";
        if(file_.Contains("DG0")) process_ += "DG0";
    }
    if(file_.Contains("Bu")){
        process_ = "BuJPsiK";
        min_ = 5.10;
        max_ = 5.65;
        x2_ = 5.64;
        x1_= 5.44;
        dirPath_ = "./Bu";
    }
    if(file_.Contains("MC")){
        process_ = process_ + "MC";
        dirPath_ += "MC";
    }
    if(file_.Contains("Data")){
        process_ = process_ + "Data";
        dirPath_ += "Data";
        isData = true;
    }

    if(file_.Contains("2017")){
        process_ = process_ + "2017";
        dirPath_ += "2017";
    }

    if(file_.Contains("2018")){
        process_ = process_ + "2018";
        dirPath_ += "2018";
    }

    cout<<"process_ = "<<process_<<endl;

    // CALIBRATION
    auto *f2 = new TFile("OSMuonTaggerCalibration" + process_ + ".root");
    f2->cd();
    TF1 *fCal = (TF1*)f2->Get("osMuonCal");
    f2->Close();
    delete f2;
    f->cd();

    cout<<"----- FILE OPEN"<<endl;

    double totalP = 0.;
    double totalPerr = 0;

    //MVA VARIABLES
    float muoPt;
    float muoEta;
    float muoDxy;
    float muoExy;
    float muoDz;
    float muoEz;
    float muoSoftMvaValue;
    float muoDrB;
    float muoPFIso;
    float muoConeCleanPt;
    float muoConeCleanPtRel;
    float muoConeCleanDr;
    float muoConeCleanEnergyRatio;
    float muoConeCleanQ;
    //TAGGING VARIABLES
    int osMuon, osMuonTag, osMuonCharge, ssbLund;
    //EVENT VARIABLES
    float ssbMass;
    int evtWeight, hltJpsiMu, ssbIsTight;

    //BOOKING
    t->SetBranchAddress("muoPt", &muoPt);
    t->SetBranchAddress("muoEta", &muoEta);
    t->SetBranchAddress("muoDxy", &muoDxy);
    t->SetBranchAddress("muoExy", &muoExy);
    t->SetBranchAddress("muoDz", &muoDz);
    t->SetBranchAddress("muoEz", &muoEz);
    t->SetBranchAddress("muoSoftMvaValue", &muoSoftMvaValue);
    t->SetBranchAddress("muoDrB", &muoDrB);
    t->SetBranchAddress("muoPFIso", &muoPFIso);
    t->SetBranchAddress("muoConeCleanPt", &muoConeCleanPt);
    t->SetBranchAddress("muoConeCleanPtRel", &muoConeCleanPtRel);
    t->SetBranchAddress("muoConeCleanDr", &muoConeCleanDr);
    t->SetBranchAddress("muoConeCleanEnergyRatio", &muoConeCleanEnergyRatio);
    t->SetBranchAddress("muoConeCleanQ", &muoConeCleanQ);
    t->SetBranchAddress("osMuon", &osMuon);
    t->SetBranchAddress("osMuonTag", &osMuonTag);
    t->SetBranchAddress("evtWeight", &evtWeight);
    t->SetBranchAddress("muoCharge", &osMuonCharge);
    t->SetBranchAddress("ssbLund", &ssbLund);
    t->SetBranchAddress("ssbMass", &ssbMass);
    t->SetBranchAddress("hltJpsiMu", &hltJpsiMu);
    t->SetBranchAddress("ssbIsTight", &ssbIsTight);

    //COMPUTE MVA
    TMVA::PyMethodBase::PyInitialize();
    TMVA::Reader reader("!Color:Silent");
    double mvaValue = -1.;

    reader.AddVariable("muoPt", &muoPt);
    reader.AddVariable("muoEta", &muoEta);
    reader.AddVariable("muoDxy", &muoDxy);
    reader.AddVariable("muoExy", &muoExy);
    reader.AddVariable("muoDz", &muoDz);
    reader.AddVariable("muoEz", &muoEz);
    reader.AddVariable("muoSoftMvaValue", &muoSoftMvaValue);
    reader.AddVariable("muoDrB", &muoDrB);
    reader.AddVariable("muoPFIso", &muoPFIso);
    reader.AddVariable("muoConeCleanPt", &muoConeCleanPt);
    reader.AddVariable("muoConeCleanPtRel", &muoConeCleanPtRel);
    reader.AddVariable("muoConeCleanDr", &muoConeCleanDr);
    reader.AddVariable("muoConeCleanEnergyRatio", &muoConeCleanEnergyRatio);
    reader.AddVariable("muoConeCleanQ", &muoConeCleanQ);
    reader.BookMVA( method_, methodpath + method_ + ".weights.xml" );

    //HISTOGRAMS BOOKING
    int nBinsMva = 200;
    auto *mistag    = new TH1F( "mistag",    "mistag",    nBinsMva, 0.0, 1.0 );
    auto *mistag_RT = new TH1F( "mistag_RT", "mistag_RT", nBinsMva, 0.0, 1.0 );
    auto *mistag_WT = new TH1F( "mistag_WT", "mistag_WT", nBinsMva, 0.0, 1.0 );
    auto *hMassTot    = new TH1F( "hMassTot","hMassTot", nBins_, min_, max_ );
    auto *hMassRT  = new TH1F( "hMassRT","hMassRT", nBins_, min_, max_ );
    auto *hMassWT  = new TH1F( "hMassWT","hMassWT", nBins_, min_, max_ );
    auto *hMassNT  = new TH1F( "hMassNT","hMassNT", nBins_, min_, max_ );

    cout<<"----- BOOKING COMPLETED"<<endl;

    //EVENT LOOP
    cout<<"----- BEGIN LOOP"<<endl;

    if(nEvents_ == -1) nEvents_ = t->GetEntries();

    for(int i=0; i<nEvents_; ++i){
        if(i%100000==0) cout<<"----- at event "<<i<<endl;

        t->GetEntry(i);

        //EVENT SELECTION
        if(!hltJpsiMu) continue;
        if(useTightSelection_ && !ssbIsTight) continue;
        hMassTot->Fill(ssbMass, evtWeight);

        //MUON SELECTION
        if(!osMuon){
            hMassNT->Fill(ssbMass, evtWeight);
            continue;
        }
        if(muoSoftMvaValue <= muonIDwp_) continue;

        //TAGGING
        mvaValue = reader.EvaluateMVA(method_);
        double evtWpred = 1. - mvaValue;
        double evtW = fCal->Eval(evtWpred); // calibration
        double evtWerr = sqrt(pow(fCal->GetParError(0),2)+pow((fCal->GetParError(1))*(evtWpred),2));
        totalP += pow(1.-2.*evtW, 2)*evtWeight;
        double errp = sqrt( pow(8*evtW-4,2)*pow(evtWerr,2) );
        totalPerr += pow(evtWeight*errp,2);

        //cout<<evtWpred<<"   "<<evtW<<" +- "<<evtWerr<<endl;

        mistag->Fill(evtW, evtWeight);

        int evtTag = -1*osMuonCharge;
        bool isTagRight = TMath::Sign(1, ssbLund) == evtTag;
        if(isTagRight!=osMuonTag) cout<<"!!?? isTagRight != osMuonTag"<<endl;

        if(isTagRight){
            mistag_RT->Fill(evtW, evtWeight);
            hMassRT->Fill(ssbMass, evtWeight);
        }else{
            mistag_WT->Fill(evtW, evtWeight);
            hMassWT->Fill(ssbMass, evtWeight);
        }
    }

    cout<<endl<<"----- EVENTS LOOP ENDED"<<endl;

    // PERFORMANCE OUTPUT
    double nRT = hMassRT->Integral(); //Integral() takes in consideration event weights
    double nWT = hMassWT->Integral();
    double nNT = hMassNT->Integral();
    double nTot = hMassTot->Integral();

    if(isData){ //for data fit mass
        nTot = CountEventsWithFit(hMassTot).first; //Fit of the total histogram need to be called first
        nRT = CountEventsWithFit(hMassRT).first;
        nWT = CountEventsWithFit(hMassWT).first;
        nNT = CountEventsWithFit(hMassNT).first;
    }

    double effBase = (double)(nRT+nWT)/(nRT+nWT+nNT);
    double wBase = (double)nWT/(nRT+nWT);
    double pBase = effBase*pow(1.-2.*wBase,2);

    cout<<"nTot "<<nTot<<endl;
    cout<<"nRT "<<nRT<<endl;
    cout<<"nWT "<<nWT<<endl;
    cout<<"nNT "<<nNT<<endl;

    cout<<endl;
    cout<<"Base efficiency = "<<100*effBase<<"%"<<endl;
    cout<<"Base mistag = "<<100*wBase<<"%"<<endl;
    cout<<"Base power = "<<100*pBase<<"%"<<endl;

    double park = totalP;
    totalP /= nTot;
    totalPerr /= pow((double)nTot,2);
    totalPerr += nTot*pow(park,2)/pow((double)nTot,4);
    totalPerr = sqrt(totalPerr);

    cout<<endl;
    cout<<"Per-event-mistag power = "<<100.*totalP<<" +- "<<100.*totalPerr<<" % (+"<<100*(totalP - pBase)/pBase<<"%)"<<endl;
    cout<<endl;


    auto *c2 = new TCanvas();
    gPad->SetGrid();
    // mistag->SetMarkerStyle(20);
    // mistag->SetMarkerSize(.75);
    mistag->SetLineWidth(2);
    mistag->SetTitle("mistagDistribution " + process_);
    mistag->GetXaxis()->SetTitle("mistag probability");
    mistag->GetXaxis()->SetNdivisions(10+100*(int)(nBins_/10), kFALSE);
    gStyle->SetOptStat(0);
    mistag->Draw("HIST L");
    c2->Print("mistagDistribution " + process_ + ".pdf");

    f->Close();
    delete f;
    return 1;

}

// ----------------------------------------------
pair<double, double> CountEventsWithFit(TH1 *hist)
{
    TString title = hist->GetTitle();
    // cout<<" ---  now fitting "<<title<<endl;

    bool isTot = title == "hMassTot" ? true : false;
    bool lowStat = hist->GetEntries()<=250 ? true : false;
    bool highStat = hist->GetEntries()>1500 ? true : false;

    TRandom3 *r3 = new TRandom3();
    double mean = 5.3663;
    double sigma = 0.015;
    if(process_.Contains("BsJPsiPhi")) mean = 5.3663;
    if(process_.Contains("BuJPsiK"))   mean = 5.2793;
    if(process_.Contains("BdJPsiKx"))  mean = 5.2796;

    TString sgnDef = "[1]*TMath::Gaus(x, [0], [2], true)";
    sgnDef +=       "+[3]*TMath::Gaus(x, [0], [4], true)";

    TString bkgDef = "[5]";
    if(!lowStat)  bkgDef += "+[6]*x";
    if(highStat) bkgDef += "+[7]*TMath::Erfc([8]*(x-[9]))";

    bkgDef = "(" + bkgDef + ">= 0 ? " + bkgDef + " : 0 )";

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
    if(m>0) m = -1;

    func->SetParameter(5, 10);
    func->SetParameter(6, m);

    if(lowStat){
        func->SetParameter(5, 1);
        func->SetParLimits(5, 0, 1e3);
    }

    if(highStat){
        func->SetParameter(7, hist->GetBinContent(2)/2);
        func->SetParameter(8, 10);
        func->SetParameter(9, 5);
        func->SetParLimits(7, 0, hist->GetBinContent(2)*1.5);
    }

    //FIXING PARAMETERS
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
    hist->SetMarkerSize(.5);
    TFitResultPtr r = hist->Fit("func","RLSQ");
    int fitstatus = r;
    int covstatus = r->CovMatrixStatus();
    if(fitstatus != 0) cout<<"STATUS of "<<title<<" --> "<<fitstatus<<endl;
    if(covstatus != 3) cout<<"COV STATUS of "<<title<<" --> "<<covstatus<<endl;
    TF1 *fit = hist->GetFunction("func");
    if(isTot){
        mean_ = fit->GetParameter("mean");
        sigma1_ = fit->GetParameter("sigma1");
        sigma2_ = fit->GetParameter("sigma2");
    }

    double nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(3);
    nEvt/=hist->GetBinWidth(1);

    TMatrixDSym cov = r->GetCovarianceMatrix();
    double errN = sqrt(cov(1,1)+cov(3,3)+2*cov(1,3))/hist->GetBinWidth(1);

    return make_pair(nEvt, errN);

}
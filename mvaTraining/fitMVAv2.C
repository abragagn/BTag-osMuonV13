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
int nBins_ = 50;
bool output_ = true;
double systematics_ = 0.;

pair<double, double> CountEventsWithFit(TH1 *hist);

int fitMVAv2(TString file_ = "ntuples/ntuBsDG0MC2018.root"
    , int nEvents_ = -1
    , int nBinCal_ = 100 // number of bin for calibration
    , bool output_par = true
    , bool useSyst_ = false
    , TString method_ = "DNNOsMuonHLTJpsiMu"
    , bool useTightSelection_ = true
    , float muonIDwp_ = 0.21
    )
{

    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(10); //ksiourmen
    gStyle->SetOptFit(1); //pcev
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 15000 );
    output_ = output_par;

    cout<<"----- Parameters -----"<<endl;
    cout<<"file_ = "<<file_<<endl;
    cout<<"method_ = "<<method_<<endl;
    cout<<"useTightSelection_ = "<<useTightSelection_<<endl;
    cout<<"nEvents_ = "<<nEvents_<<endl;
    cout<<"nBinCal_ = "<<nBinCal_<<endl;

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
        if(useSyst_) systematics_ = 0.2;
    }

    if(file_.Contains("2018")){
        process_ = process_ + "2018";
        dirPath_ += "2018";
        if(useSyst_) systematics_ = 0.35;
    }

    cout<<"----- FILE OPEN"<<endl;

    // PER-EVENT VARIABLES
    double evtW[2] = {-1., -1.}; // per-event mistag rate
    double totalP = 0.; // total tagging power
    double totalPbinned = 0.;
    double totalPbinned_err = 0.;

    double pass = (1.-0.)/nBinCal_; // (1-0) -> mva score range

    TH1F *hMassCalRT[nBinCal_];
    TH1F *hMassCalWT[nBinCal_];
    double *wCalc = new double[nBinCal_]; // measured mistag rate
    double *wCalcEdgeL = new double[nBinCal_];
    double *wCalcEdgeH = new double[nBinCal_];

    for(int i=0;i<nBinCal_;++i){
        hMassCalRT[i] = new TH1F(TString::Format("hMassCalRT%i", i),TString::Format("hMassCalRT%i", i),nBins_, min_, max_);
        hMassCalWT[i] = new TH1F(TString::Format("hMassCalWT%i", i),TString::Format("hMassCalWT%i", i),nBins_, min_, max_);
        wCalc[i] = 0.;
        wCalcEdgeL[i] = 0.;
        wCalcEdgeH[i] = 0;
    }

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
    int nBinsMva = 100;
    auto *mva    = new TH1F( "mva",    "mva",    nBinsMva, 0.0, 1.0 );
    auto *mva_RT = new TH1F( "mva_RT", "mva_RT", nBinsMva, 0.0, 1.0 );
    auto *mva_WT = new TH1F( "mva_WT", "mva_WT", nBinsMva, 0.0, 1.0 );
    auto *hMassTot    = new TH1F( "hMassTot","hMassTot", nBins_, min_, max_ );
    auto *hMassRT  = new TH1F( "hMassRT","hMassRT", nBins_, min_, max_ );
    auto *hMassWT  = new TH1F( "hMassWT","hMassWT", nBins_, min_, max_ );
    auto *hMassNT  = new TH1F( "hMassNT","hMassNT", nBins_, min_, max_ );

    cout<<"----- BOOKING COMPLETED"<<endl;

    //EVENT LOOP
    cout<<"----- BEGIN LOOP"<<endl;

    if(nEvents_ == -1) nEvents_ = t->GetEntries();

    for(int i=0; i<nEvents_; ++i){
        if(i%1000000==0) cout<<"----- at event "<<i+1<<"/"<<nEvents_<<endl;

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
        mva->Fill(mvaValue, evtWeight);
        evtW[0] = 1-mvaValue;
        totalP += pow(1.-2.*evtW[0], 2)*evtWeight;

        int evtTag = -1*osMuonCharge;
        bool isTagRight = TMath::Sign(1, ssbLund) == evtTag;
        if(isTagRight!=osMuonTag) cout<<"!!?? isTagRight != osMuonTag"<<endl;

        for(int j=0;j<nBinCal_;++j){
            if( (evtW[0]>=(double)j*pass) && (evtW[0]<((double)j*pass+pass)) ){
                if(isTagRight) hMassCalRT[j]->Fill(ssbMass, evtWeight);
                else           hMassCalWT[j]->Fill(ssbMass, evtWeight);
                wCalc[j] += evtW[0]*evtWeight;
                break;
            }
        }

        if(isTagRight){
            mva_RT->Fill(mvaValue, evtWeight);
            hMassRT->Fill(ssbMass, evtWeight);
        }else{
            mva_WT->Fill(mvaValue, evtWeight);
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

    totalP /= (double)(nRT+nWT+nNT);
    cout<<endl;
    cout<<"Per-event-mistag power (not calibrated) = "<<100.*totalP<<"% (+"<<100*(totalP - pBase)/pBase<<"%)"<<endl;
    cout<<endl;

    // CALIBRATION
    for(int j=0;j<nBinCal_;++j){
        wCalc[j] /= (hMassCalRT[j]->Integral() + hMassCalWT[j]->Integral());
    }

    vector<double> vX;
    vector<double> vY;
    vector<double> vEXL;
    vector<double> vEXH;
    vector<double> vEYL;
    vector<double> vEYH;

    int minEntries = 0;
    if(isData) minEntries = 20;
    int rebinThr = 1500;

    for(int j=0;j<nBinCal_;++j){
        pair<double, double> calRT; // .first = nEvt; .second = sigma(nEvt)
        pair<double, double> calWT;
        double wMeas;
        double wMeasErr;
        double wMeasErrL;
        double wMeasErrH;

        calRT.first = hMassCalRT[j]->Integral();
        calWT.first = hMassCalWT[j]->Integral();
        calRT.second = sqrt(calRT.first);
        calWT.second = sqrt(calWT.first);

        if( calRT.first<=minEntries && calWT.first<=minEntries ) continue;

        if(isData){
            if(calRT.first<rebinThr) hMassCalRT[j]->Rebin();
            if(calWT.first<rebinThr) hMassCalWT[j]->Rebin();
            if(calRT.first >= minEntries)
                calRT = CountEventsWithFit(hMassCalRT[j]);
            if(calWT.first >= minEntries)
                calWT = CountEventsWithFit(hMassCalWT[j]);
            if(calRT.second < 1 ) calRT.second = 2*sqrt(calRT.first);
            if(calWT.second < 1 ) calWT.second = 2*sqrt(calWT.first);
            wMeas = calWT.first/(calWT.first + calRT.first);
            wMeasErr = sqrt(pow(calWT.first,2)*pow(calRT.second,2) 
                          + pow(calRT.first,2)*pow(calWT.second,2))/pow(calRT.first+calWT.first,2);

            wMeasErr += wMeasErr*systematics_;

            wMeasErrH = wMeasErr;
            wMeasErrL = wMeasErr;
            if(wMeas + wMeasErrH > 1.) wMeasErrH = 1. - wMeas;
            if(wMeas - wMeasErrL < 0.) wMeasErrL = wMeas - 0.;
        }else{
            wMeas = calWT.first/(calWT.first + calRT.first);
            wMeasErrH = TEfficiency::AgrestiCoull(calWT.first+calRT.first, calWT.first,0.6827,true) - wMeas;
            wMeasErrL = wMeas - TEfficiency::AgrestiCoull(calWT.first+calRT.first, calWT.first,0.6827,false);
            wMeasErr = max(wMeasErrH, wMeasErrL);
        }

        vX.push_back( wCalc[j] );
        vEXL.push_back( wCalc[j] - ((double)j*pass) );
        vEXH.push_back( ((double)j*pass+pass) - wCalc[j] );
        vY.push_back( wMeas ); // measured mistag
        vEYL.push_back(wMeasErrL);
        vEYH.push_back(wMeasErrH);

        totalPbinned += (calRT.first+calWT.first)*pow(1-2*wMeas,2);
        // totalPbinned_err //TODO
        cout<<"BIN "<<j<<" -- wCalc "<<wCalc[j];
        cout<<" -- nRT "<<(int)calRT.first<<" +- "<<(int)calRT.second<<" -- nWT "<<(int)calWT.first<<" +- "<<(int)calWT.second;
        cout<<" -- wMeas "<<wMeas <<" +- "<<wMeasErr<<endl<<endl;
    }

    totalPbinned /= (double)nTot;

    cout<<endl;
    cout<<"Per-event-mistag power (calibrated bins) = "<<100.*totalPbinned<<"% +- "<<100.*totalPbinned_err<<" (+"<<100*(totalPbinned - pBase)/pBase<<"%)"<<endl;
    cout<<endl;

    cout<<endl;

    auto *c1 = new TCanvas("c1","c1",1000,1600);
    TPad *pad1 = new TPad("pad1", "",0.0,0.3,1.0,1.0);
    TPad *pad2 = new TPad("pad2", "",0.0,0.0,1.0,0.3);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    gPad->SetGrid();

    auto *gCal = new TGraphAsymmErrors(vX.size(),&vX[0],&vY[0],0,0,&vEYL[0],&vEYH[0]);
    auto *gCalErr = new TGraphAsymmErrors(vX.size(),&vX[0],&vY[0],&vEXL[0],&vEXH[0],&vEYL[0],&vEYH[0]);
    auto *fCal = new TF1("osMuonCal","[0]+[1]*x",0.,1.);
    gCal->SetName("osMuonCalGraph");
    gCalErr->SetName("osMuonCalGraphBins");

    TFitResultPtr fitresultCal = gCal->Fit("osMuonCal","S");
    fCal = gCal->GetFunction("osMuonCal");

    double q = fCal->GetParameter(0);
    double m = fCal->GetParameter(1);

    cout<<endl;
    cout<<"q = "<<q<<" +- "<<fCal->GetParError(0)<<" ["<<abs(q)/fCal->GetParError(0)<<" s.d.]"<<endl;
    cout<<"m = "<<m<<" +- "<<fCal->GetParError(1)<<" ["<<abs(m-1)/fCal->GetParError(1)<<" s.d.]"<<endl;

    vector<double> wResY;
    vector<double> wResEY;
    vector<double> wResEYH;
    vector<double> wResEYL;
    vector<double> wRatioY;
    vector<double> wRatioEY;
    vector<double> wRatioEYH;
    vector<double> wRatioEYL;

    for (unsigned int j=0;j<vX.size();++j){
        double dev = vY[j] - fCal->Eval(vX[j]);
        double sigma;
        if(dev>=0) sigma = vEYL[j];
        else       sigma = vEYH[j];

        wResEYH.push_back(vEYH[j]/sigma);
        wResEYL.push_back(vEYL[j]/sigma);
        wResY.push_back(dev/sigma);
        wResEY.push_back(1.);
    }

    auto *gCalRes = new TGraphAsymmErrors(vX.size(),&vX[0],&wResY[0],&vEXL[0],&vEXH[0],&wResEYL[0],&wResEYH[0]);

    gCal->SetMarkerStyle(20);
    gCal->SetMarkerSize(1);
    gCal->SetMaximum(1.);
    gCal->SetMinimum(0.);
    gCal->GetXaxis()->SetLimits(0.,1.);
    gCal->SetTitle("calibration " + process_);
    gCal->GetXaxis()->SetTitle("mistag calc.");
    gCal->GetYaxis()->SetTitle("mistag meas.");
    gCal->Draw("AP");
    gCalErr->Draw("EZ same");
    fCal->Draw("same");
    gPad->Modified(); gPad->Update();
    TPaveStats *st = (TPaveStats*)gCal->FindObject("stats");
    st->SetX1NDC(0.65);
    st->SetX2NDC(0.95);
    st->SetY1NDC(0.15);
    st->SetY2NDC(0.30);
    gPad->Modified(); gPad->Update();

    pad2->cd();
    gPad->SetGrid();
    gCalRes->SetMarkerStyle(20);
    gCalRes->SetMarkerSize(1);
    gCalRes->GetXaxis()->SetLimits(0.0,1.02);
    gCalRes->Draw("APZ");
    gCalRes->SetTitle("");
    gCalRes->GetYaxis()->SetTitle("# s.d.");
    auto *y0_ = new TF1("","0.",0.,1.02);
    y0_->SetLineColor(kBlack);
    y0_->Draw("SAME");

    if(output_) c1->Print("calibration" + process_ + ".pdf");

    auto *c2 = new TCanvas();
    gPad->SetGrid();
    mva->SetMarkerStyle(20);
    mva->SetMarkerSize(.75);
    mva->SetTitle("dnnDistribution " + process_);
    mva->GetXaxis()->SetTitle("dnn score (right tag probability)");
    mva->GetXaxis()->SetNdivisions(10+100*(int)(nBins_/10), kFALSE);
    gStyle->SetOptStat(0);
    mva->Draw("HIST PL");
    if(output_) c2->Print("dnnDistribution_" + process_ + ".pdf");
    gStyle->SetOptStat(10);

    //FUNCTIONS
    if(output_){
        auto *fo = new TFile("OSMuonTaggerCalibration" + process_ + ".root", "RECREATE");
        fo->cd();
        fCal->Write();
        gCal->Write();
        // gCalErr->Write();
        fitresultCal->Write();
        fo->Close();
        delete fo;
    }
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
    double y2 = hist->GetBinContent(binx2);
    double y1 = hist->GetBinContent(binx1);
    double x2 = hist->GetBinCenter(binx2);
    double x1 = hist->GetBinCenter(binx1);
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
    hist->SetMarkerSize(.75);
    TFitResultPtr r = hist->Fit("func","LRSQ");
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

    hist->Draw("PE");
    hist->SetMinimum(0);

    // PLOTTING
    TF1 *f1 = new TF1("f1","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f2 = new TF1("f2","[0]*TMath::Gaus(x, [1], [2], true)", min_, max_);
    TF1 *f4;
    TF1 *f5 = new TF1("f5","[0]+[1]*x", min_, max_);

    f1->SetParameters(fit->GetParameter("A1"),fit->GetParameter("mean"),fit->GetParameter("sigma1"));
    f2->SetParameters(fit->GetParameter("A2"),fit->GetParameter("mean"),fit->GetParameter("sigma2"));

    if(lowStat){
        f4 = new TF1("f4","[0]", min_, max_);
        f4->SetParameter(0, fit->GetParameter(5));
    }

    if(highStat){
        f4 = new TF1("f4","[0]+[1]*x+[2]*TMath::Erfc([3]*(x-[4]))", min_, max_);
        f4->SetParameters(fit->GetParameter(5),fit->GetParameter(6)
                        ,fit->GetParameter(7),fit->GetParameter(8),fit->GetParameter(9));
        f5->SetParameters(fit->GetParameter(5),fit->GetParameter(6));
    }

    if(!lowStat && !highStat){
        f4 = new TF1("f4","[0]+[1]*x", min_, max_);
        f4->SetParameters(fit->GetParameter(5),fit->GetParameter(6));
    }

    f1->SetLineColor(kBlue);
    f2->SetLineColor(kViolet);
    f4->SetLineColor(kGreen);
    if(highStat) f5->SetLineColor(kOrange);
    f1->SetLineStyle(2);
    f2->SetLineStyle(2);
    f4->SetLineStyle(2);
    if(highStat) f5->SetLineStyle(2);

    f1->Draw("same");
    f2->Draw("same");
    f4->Draw("same");
    if(highStat) f5->Draw("same");

    if(output_) c5->Print(dirPath_ + "/" + title + ".pdf");

    double nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(3);
    nEvt/=hist->GetBinWidth(1);

    TMatrixDSym cov = r->GetCovarianceMatrix();
    double errN = sqrt(cov(1,1)+cov(3,3)+2*cov(1,3))/hist->GetBinWidth(1);

    return make_pair(nEvt, errN);

}
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
#include "TCut.h"
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

void customFOM( TString filename = "./2017/ntuBuMC2017.root"
    , TString var = "muoPt"
    , TString cutDir = "down"
    , TString cutEvt_ = ""
    , float min = -1
    , float max = -1
    , bool verbose = false )
{
    if((cutDir != "down") && (cutDir != "up")) return;

    gStyle->SetOptTitle(1); gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);   gStyle->SetStatBorderSize(0);
    gStyle->SetStatX(.49);  gStyle->SetStatY(.89);

    TFile *f = new TFile(filename);
    TTree *t = (TTree*)f ->Get("PDsecondTree");

    int nBins=20;
    if(min == -1) min=t->GetMinimum(var);
    if(max == -1) max=t->GetMaximum(var);

    TH1F *Tot = new TH1F( "Tot", "Tot", 100, 5.2, 5.65 );
    TH1F *Sgn = new TH1F( "Sgn", "Sgn", nBins, min, max );
    TH1F *Bkg = new TH1F( "Bkg", "Bkg", nBins, min, max );

    TString cutEvt = "ssbIsTight";

    if(cutEvt_ != "") cutEvt = cutEvt + "&&" + cutEvt_;

    TString base =  "evtWeight*((" + cutEvt;
    TString cutTT = base + ")&&osMuon)";
    TString sgnDef = base + ")&&osMuon&&osMuonTag==1)";
    TString bkgDef = base + ")&&osMuon&&osMuonTag==0)";


    t->Project("Tot", "ssbMass", base + "))");
    t->Project("Sgn", var, sgnDef);
    t->Project("Bkg", var, bkgDef);

    int nTot = Tot->Integral();
    float nSgn=0, nBkg=0;

    float *x = new float[nBins];
    float *yE = new float[nBins];
    float *yW = new float[nBins];
    float *yP = new float[nBins];

    for(int i=0; i<nBins; ++i){
        x[i]=0;
        yE[i]=0;
        yW[i]=0;
        yP[i]=0;
    }

    float bestCut = 0.;
    float bestP = 0., bestE = 0., bestW = 0.5;

    if(verbose) cout<<endl<<"nTot = "<<nTot<<endl;
    if(verbose) cout<<"x[bin]  nSgn  nBkg  yE[bin]  yW[bin]  yP[bin]"<<endl<<endl;

    if(cutDir == "down"){
        for(int bin=nBins; bin>0; --bin){

            float varCut = Sgn->GetBinLowEdge(bin);
            nSgn += Sgn->GetBinContent(bin);
            nBkg += Bkg->GetBinContent(bin);
            if(nSgn + nBkg == 0) continue;
            float eff = (nSgn + nBkg)/ nTot;
            float w = nBkg / (nSgn + nBkg);
            float p = eff*pow(1-2*w,2);
            x[bin-1] = varCut;
            yE[bin-1] = 100*eff ;
            yW[bin-1] = 100*w;
            yP[bin-1] = 100*p;

            if(verbose) cout<<x[bin-1]<<"    "<<nSgn<<"    "<<nBkg<<"    ";
            if(verbose) cout<<yE[bin-1]<<"    "<<yW[bin-1]<<"    "<<yP[bin-1]<<endl;

            if(p > bestP){
                bestCut = varCut;
                bestP = p;
                bestE = eff;
                bestW = w;
            }
        }       
    }else{
        for(int bin=0; bin<=nBins; ++bin){

            float varCut = Sgn->GetBinLowEdge(bin+1);
            nSgn += Sgn->GetBinContent(bin);
            nBkg += Bkg->GetBinContent(bin);
            if(nSgn + nBkg == 0) continue;
            float eff = (nSgn + nBkg)/ nTot;
            float w = nBkg / (nSgn + nBkg);
            float p = eff*pow(1-2*w,2);
            x[bin] = varCut;
            yE[bin] = 100*eff ;
            yW[bin] = 100*w;
            yP[bin] = 100*p;

            if(verbose) cout<<x[bin]<<"    "<<nSgn<<"    "<<nBkg<<"    ";
            if(verbose) cout<<yE[bin]<<"    "<<yW[bin]<<"    "<<yP[bin]<<endl;

            if(p > bestP){
                bestCut = varCut;
                bestP = p;
                bestE = eff;
                bestW = w;
            }
        }           
    }

    TGraph* grP = new TGraph(nBins, x, yP);
    TGraph* grE = new TGraph(nBins, x, yE);
    TGraph* grW = new TGraph(nBins, x, yW);

    TCanvas *c1 = new TCanvas("c1", "c1", 600, 1000);
    c1->Divide(1,3);
    TString opt = "AL";
    c1->cd(1); grP->SetTitle("Power"); grP->Draw(opt); 
    c1->cd(2); grE->SetTitle("Efficiency"); grE->Draw(opt);
    c1->cd(3); grW->SetTitle("Mistag"); grW->Draw(opt);

    c1->Print("FOM_" + var + ".png");
    c1->DrawClone();

    TCanvas *c2 = new TCanvas();
    grW->Draw(opt);
    c1->DrawClone();

    cout<<"best cut = "<<bestCut<<endl;
    cout<<"(P="<<100*bestP<<"%, eff="<<100*bestE<<"%, w="<<100*bestW<<"%)"<<endl;

}

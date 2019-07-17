#ifndef PDSecondNtupleData_h
#define PDSecondNtupleData_h
#include <vector>
#include "NtuTool/Common/interface/TreeWrapper.h"
using namespace std;

class PDSecondNtupleData: public virtual TreeWrapper {

public:

void Reset() { autoReset(); }

PDSecondNtupleData() {


}
virtual ~PDSecondNtupleData() {
}

void initTree() {
    treeName = "PDsecondTree";

    setBranch( "ssbPt", &ssbPt, "ssbPt/F", &b_ssbPt );
    setBranch( "ssbEta", &ssbEta, "ssbEta/F", &b_ssbEta );
    setBranch( "ssbPhi", &ssbPhi, "ssbPhi/F", &b_ssbPhi );
//    setBranch( "ssbDxy", &ssbDxy, "ssbDxy/F", &b_ssbDxy );
//    setBranch( "ssbExy", &ssbExy, "ssbExy/F", &b_ssbExy );
//    setBranch( "ssbDz", &ssbDz, "ssbDz/F", &b_ssbDz );
//    setBranch( "ssbEz", &ssbEz, "ssbEz/F", &b_ssbEz );

    setBranch( "ssbMass", &ssbMass, "ssbMass/F", &b_ssbMass );
//    setBranch( "jpsiMass", &jpsiMass, "jpsiMass/F", &b_jpsiMass );
    setBranch( "ssbSVT", &ssbSVT, "ssbSVT/I", &b_ssbSVT );
    setBranch( "ssbPVT", &ssbPVT, "ssbPVT/I", &b_ssbPVT );
    setBranch( "ssbLund", &ssbLund, "ssbLund/I", &b_ssbLund );
    setBranch( "ssbIsTight", &ssbIsTight, "ssbIsTight/I", &b_ssbIsTight );

    setBranch( "ssbLxy", &ssbLxy , "ssbLxy/F" , &b_ssbLxy );
    setBranch( "ssbCt2D", &ssbCt2D , "ssbCt2D/F" , &b_ssbCt2D );
    setBranch( "ssbCt2DErr", &ssbCt2DErr , "ssbCt2DErr/F" , &b_ssbCt2DErr );
    setBranch( "ssbCt2DSigmaUnit", &ssbCt2DSigmaUnit , "ssbCt2DSigmaUnit/F" , &b_ssbCt2DSigmaUnit );
    setBranch( "ssbCt3D", &ssbCt3D , "ssbCt3D/F" , &b_ssbCt3D );
    setBranch( "ssbCt3DErr", &ssbCt3DErr , "ssbCt3DErr/F" , &b_ssbCt3DErr );
    setBranch( "ssbCt3DSigmaUnit", &ssbCt3DSigmaUnit , "ssbCt3DSigmaUnit/F" , &b_ssbCt3DSigmaUnit );


    setBranch( "osMuon", &osMuon, "osMuon/I", &b_osMuon );
    setBranch( "osMuonJet", &osMuonJet, "osMuonJet/I", &b_osMuonJet );
    setBranch( "osMuonTag", &osMuonTag, "osMuonTag/I", &b_osMuonTag );
    setBranch( "osMuonChargeInfo", &osMuonChargeInfo, "osMuonChargeInfo/I", &b_osMuonChargeInfo );
    setBranch( "osMuonTagMvaValue", &osMuonTagMvaValue, "osMuonTagMvaValue/F", &b_osMuonTagMvaValue );
    setBranch( "osMuonTagMistag", &osMuonTagMistag, "osMuonTagMistag/F", &b_osMuonTagMistag );

    setBranch( "muoPt", &muoPt, "muoPt/F", &b_muoPt );
    setBranch( "muoEta", &muoEta, "muoEta/F", &b_muoEta );
    setBranch( "muoPhi", &muoPhi, "muoPhi/F", &b_muoPhi );
    setBranch( "muoCharge", &muoCharge, "muoCharge/I", &b_muoCharge );
    setBranch( "muoDxy", &muoDxy, "muoDxy/F", &b_muoDxy );
    setBranch( "muoExy", &muoExy, "muoExy/F", &b_muoExy );
    setBranch( "muoDz", &muoDz, "muoDz/F", &b_muoDz );
    setBranch( "muoEz", &muoEz, "muoEz/F", &b_muoEz );
    setBranch( "muoLund", &muoLund, "muoLund/I", &b_muoLund );
    setBranch( "muoAncestor", &muoAncestor, "muoAncestor/I", &b_muoAncestor );
    setBranch( "muoSoftMvaValue", &muoSoftMvaValue, "muoSoftMvaValue/F", &b_muoSoftMvaValue );

    setBranch( "muoDrB", &muoDrB, "muoDrB/F", &b_muoDrB );
    setBranch( "muoPFIso", &muoPFIso, "muoPFIso/F", &b_muoPFIso );

    setBranch( "muoJetPt", &muoJetPt, "muoJetPt/F", &b_muoJetPt );
    setBranch( "muoJetPtRel", &muoJetPtRel, "muoJetPtRel/F", &b_muoJetPtRel );
    setBranch( "muoJetDr", &muoJetDr, "muoJetDr/F", &b_muoJetDr );
    setBranch( "muoJetEnergyRatio", &muoJetEnergyRatio, "muoJetEnergyRatio/F", &b_muoJetEnergyRatio );
    setBranch( "muoJetQ", &muoJetQ, "muoJetQ/F", &b_muoJetQ );
    setBranch( "muoJetCSV", &muoJetCSV, "muoJetCSV/F", &b_muoJetCSV );
    setBranch( "muoJetDFprob", &muoJetDFprob, "muoJetDFprob/F", &b_muoJetDFprob );
    setBranch( "muoJetSize", &muoJetSize, "muoJetSize/I", &b_muoJetSize );
    setBranch( "muoJetNDau", &muoJetNDau, "muoJetNDau/I", &b_muoJetNDau );
    setBranch( "muoJetNF", &muoJetNF, "muoJetNF/F", &b_muoJetNF );
    setBranch( "muoJetCF", &muoJetCF, "muoJetCF/F", &b_muoJetCF );
    setBranch( "muoJetNCH", &muoJetNCH, "muoJetNCH/I", &b_muoJetNCH );

    setBranch( "muoConePt", &muoConePt, "muoConePt/F", &b_muoConePt );
    setBranch( "muoConePtRel", &muoConePtRel, "muoConePtRel/F", &b_muoConePtRel );
    setBranch( "muoConeDr", &muoConeDr, "muoConeDr/F", &b_muoConeDr );
    setBranch( "muoConeEnergyRatio", &muoConeEnergyRatio, "muoConeEnergyRatio/F", &b_muoConeEnergyRatio );
    setBranch( "muoConeQ", &muoConeQ, "muoConeQ/F", &b_muoConeQ );
    setBranch( "muoConeSize", &muoConeSize, "muoConeSize/I", &b_muoConeSize );
    setBranch( "muoConeNF", &muoConeNF, "muoConeNF/F", &b_muoConeNF );
    setBranch( "muoConeCF", &muoConeCF, "muoConeCF/F", &b_muoConeCF );
    setBranch( "muoConeNCH", &muoConeNCH, "muoConeNCH/I", &b_muoConeNCH );

    setBranch( "muoConeCleanPt", &muoConeCleanPt, "muoConeCleanPt/F", &b_muoConeCleanPt );
    setBranch( "muoConeCleanPtRel", &muoConeCleanPtRel, "muoConeCleanPtRel/F", &b_muoConeCleanPtRel );
    setBranch( "muoConeCleanDr", &muoConeCleanDr, "muoConeCleanDr/F", &b_muoConeCleanDr );
    setBranch( "muoConeCleanEnergyRatio", &muoConeCleanEnergyRatio, "muoConeCleanEnergyRatio/F", &b_muoConeCleanEnergyRatio );
    setBranch( "muoConeCleanQ", &muoConeCleanQ, "muoConeCleanQ/F", &b_muoConeCleanQ );
    setBranch( "muoConeCleanSize", &muoConeCleanSize, "muoConeCleanSize/I", &b_muoConeCleanSize );
    setBranch( "muoConeCleanNF", &muoConeCleanNF, "muoConeCleanNF/F", &b_muoConeCleanNF );
    setBranch( "muoConeCleanCF", &muoConeCleanCF, "muoConeCleanCF/F", &b_muoConeCleanCF );
    setBranch( "muoConeCleanNCH", &muoConeCleanNCH, "muoConeCleanNCH/I", &b_muoConeCleanNCH );

    setBranch( "muoHowMany", &muoHowMany, "muoHowMany/I", &b_muoHowMany );

    setBranch( "evtNumber", &evtNumber, "evtNumber/I", &b_evtNumber );
    setBranch( "evtWeight", &evtWeight, "evtWeight/I", &b_evtWeight );
    setBranch( "evtNb", &evtNb, "evtNb/I", &b_evtNb );

    setBranch( "hltJpsiMu", &hltJpsiMu , "hltJpsiMu/I" , &b_hltJpsiMu );
    setBranch( "hltJpsiTrkTrk", &hltJpsiTrkTrk , "hltJpsiTrkTrk/I" , &b_hltJpsiTrkTrk );
    setBranch( "hltJpsiTrk", &hltJpsiTrk , "hltJpsiTrk/I" , &b_hltJpsiTrk );
}

float ssbPt, ssbEta, ssbPhi, ssbMass, jpsiMass, ssbDxy, ssbExy, ssbDz, ssbEz;
float ssbLxy, ssbCt2D, ssbCt2DErr, ssbCt2DSigmaUnit, ssbCt3D, ssbCt3DErr, ssbCt3DSigmaUnit;
int ssbSVT, ssbPVT, ssbLund, evtNumber, hltJpsiMu, hltJpsiTrkTrk, hltJpsiTrk, ssbIsTight, evtNb, evtWeight;

TBranch *b_ssbPt, *b_ssbEta, *b_ssbPhi, *b_ssbMass, *b_jpsiMass, *b_ssbDxy, *b_ssbExy, *b_ssbDz, *b_ssbEz;
TBranch *b_ssbLxy, *b_ssbCt2D, *b_ssbCt2DErr, *b_ssbCt2DSigmaUnit, *b_ssbCt3D, *b_ssbCt3DErr, *b_ssbCt3DSigmaUnit;
TBranch *b_ssbSVT, *b_ssbPVT, *b_ssbLund, *b_evtNumber, *b_hltJpsiMu, *b_hltJpsiTrkTrk, *b_hltJpsiTrk, *b_ssbIsTight, *b_evtWeight, *b_evtNb;

int muoLund, muoAncestor, muoCharge;
float muoSoftMvaValue;
float muoPt, muoEta, muoPhi, muoDxy, muoExy, muoDz, muoEz;

TBranch *b_muoPt, *b_muoEta, *b_muoPhi, *b_muoDxy, *b_muoExy, *b_muoDz, *b_muoEz, *b_muoCharge;
TBranch *b_muoLund, *b_muoAncestor, *b_muoSoftMvaValue;

int osMuon, osMuonJet, osMuonTag, osMuonChargeInfo;
float osMuonTagMvaValue, osMuonTagMistag;
float muoDrB, muoDzPV, muoPFIso, muoQCone;
int muoHowMany;

TBranch *b_osMuon, *b_osMuonJet, *b_osMuonTag, *b_osMuonChargeInfo, *b_muoDrB, *b_muoPFIso, *b_muoQCone, *b_muoDzPV, *b_muoHowMany;
TBranch *b_osMuonTagMvaValue, *b_osMuonTagMistag;

int muoJetSize, muoJetNDau, muoJetNCH;
float muoJetPtRel, muoJetDr, muoJetEnergyRatio, muoJetQ, muoJetQtight, muoJetCSV, muoJetDFprob, muoJetPt, muoJetNF, muoJetCF;

TBranch *b_muoJetSize, *b_muoJetPtRel, *b_muoJetDr, *b_muoJetEnergyRatio, *b_muoJetQ, *b_muoJetQtight, *b_muoJetCSV, *b_muoJetDFprob, *b_muoJetPt;
TBranch *b_muoJetNDau, *b_muoJetNF, *b_muoJetCF, *b_muoJetNCH;

int muoConeSize, muoConeNCH;
float muoConePtRel, muoConeDr, muoConeEnergyRatio, muoConeQ, muoConeCSV, muoConeDFprob, muoConePt, muoConeNF, muoConeCF;

TBranch *b_muoConeNF, *b_muoConeCF, *b_muoConeNCH;
TBranch *b_muoConeSize, *b_muoConePtRel, *b_muoConeDr, *b_muoConeEnergyRatio, *b_muoConeQ, *b_muoConeCSV, *b_muoConeDFprob, *b_muoConePt;

int muoConeCleanSize, muoConeCleanNCH;
float muoConeCleanPtRel, muoConeCleanDr, muoConeCleanEnergyRatio, muoConeCleanQ, muoConeCleanCSV, muoConeCleanDFprob, muoConeCleanPt, muoConeCleanNF, muoConeCleanCF;

TBranch *b_muoConeCleanNF, *b_muoConeCleanCF, *b_muoConeCleanNCH;
TBranch *b_muoConeCleanSize, *b_muoConeCleanPtRel, *b_muoConeCleanDr, *b_muoConeCleanEnergyRatio, *b_muoConeCleanQ, *b_muoConeCleanCSV, *b_muoConeCleanDFprob, *b_muoConeCleanPt;

private:

PDSecondNtupleData ( const PDSecondNtupleData& a );
PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif


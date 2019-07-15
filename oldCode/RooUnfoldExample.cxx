//=====================================================================-*-C++-*-
//
//
//
//==============================================================================

//TODO SCALE hJtUnfBg correctly before subtracting from data
//TODO SCALE All histograms before unfolding

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TFile.h"
#include <TStopwatch.h>

#include "AliJHistManager.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
#include "RooUnfold/src/RooUnfoldBayes.h"
#include "RooUnfold/src/RooUnfoldSvd.h"
//#include "RooUnfold/src/RooUnfoldTUnfold.h"
#endif

TString dirnameData = "AliJJetJtTask_kEMCEJE/AliJJetJtHistManager";
//TString dirnameData = "AliJJetJtTask/AliJJetJtHistManager";
TString dirname = "AliJJetJtTask/AliJJetJtHistManager";

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

//AliJHistManager *hmg;
const int Njet = 9;
const int rebin = 1;
const int doWeight = 1;
const int nR = 1;
const int Nsets = 2*nR;
const bool doSignal = kFALSE;
//TString option = "width";
TString option = "";
bool useRndmBg = false;
bool doData = true;
bool do2D = false;
bool scaleResponse = false; //Don't use without option "width"
double n_jet_meas[Nsets][Njet];
double n_jet_true[Nsets][Njet];
double n_jet_pythia[Njet];
int  n_jet_data[Nsets][Njet];
double n_bg_meas[Nsets][Njet];
double n_bg_true[Nsets][Njet];
int n_bg_data[Nsets][Njet];

TH2D *responseMatrix[Nsets][Njet];
TH2D *responseMatrixConst[Nsets][Njet];
TH2D *responseMatrixPythia[Njet];
TH2D *responseMatrixJetPt[Nsets];
TH2D *responseMatrixJetPtCoarse[Nsets];
TH2D *responseCoarse;
TH1D *hJtMeas[Nsets][Njet];
TH1D *hJtMeas_unscaled[Nsets][Njet];
TH1D *hJtMeasBg[Nsets][Njet];
TH1D *hJtMeasBgRndm[Nsets][Njet];
TH1D *hJtPythia[Njet];
TH1D *hJtMeasBgConst[Nsets][Njet];
TH1D *hJtMeasConst[Nsets][Njet];
TH1D *hJtMeasConst_unscaled[Nsets][Njet];
TH1D *hJtTrue[Nsets][Njet];
TH1D *hJtTrueBg[Nsets][Njet];
TH1D *hJtTrueBgRndm[Nsets][Njet];
TH1D *hJtTrueConst[Nsets][Njet];
TH1D *hJtTrueBgConst[Nsets][Njet];
TH1D *hJtUnfBg[Nsets][Njet];
TH1D *hJtUnfBg_unscaled[Nsets][Njet];
TH1D *hJtUnfBgData[Nsets][Njet];

TH1D *hJtData[Nsets][Njet];
TH1D *hJtDataSub[Nsets][Njet];
TH1D *hJtDataConst[Nsets][Njet];
TH1D *hJtDataConstSub[Nsets][Njet];
TH1D *hJetPtData[Nsets];
TH1D *hJetPtBinData[Nsets][Njet];
TH1D* hRecoBayesData[Nsets][Njet];
TH1D* hRecoBayesDataSub[Nsets][Njet];
TH1D *hRecoSVDData[Nsets][Njet];
TH1D *hRecoSVDDataSub[Nsets][Njet];
TH1D *hRecoBayesDataConst[Nsets][Njet];
TH1D *hRecoSVDDataConst[Nsets][Njet];
TH1D *hRecoBayesDataConstSub[Nsets][Njet];
TH1D *hRecoSVDDataConstSub[Nsets][Njet];
TH1D *hJetPtReco[Nsets];
TH1D *hJtBgData[Nsets][Njet];
TH1D *hJtBgRndmData[Nsets][Njet];

TH1D *hJetPtMeas[Nsets];
TH1D *hJetPtBinMeas[Nsets][Njet];
TH1D *hJetPtBinPythia[Njet];
TH1D *hJetPtTrue[Nsets];
TH1D *hJetPtBinTrue[Nsets][Njet];
TH1D *hRecoBayes[Nsets][Njet];
TH1D *hRecoBayes2D[Nsets][Njet];
TH1D *hRecoBayes2DSub[Nsets][Njet];
TH1D *hRecoBayesConst2D[Nsets][Njet];
TH1D *hRecoBayesConst2DSub[Nsets][Njet];
TH1D *hRecoBayesConst[Nsets][Njet];
TH1D *hRecoBayes1D[Nsets][Njet];
TH1D *hRecoBayesConst1D[Nsets][Njet];
TH1D *hRecoBayes1DSub[Nsets][Njet];
TH1D *hRecoBayesConst1DSub[Nsets][Njet];
TH1D *hRecoBayesConst1D[Nsets][Njet];
TH1D *hRecoSVD[Nsets][Njet];
TH1D *hRecoSVD2D[Nsets][Njet];
TH1D *hRecoSVDConst2D[Nsets][Njet];
TH1D *hRecoSVDSub[Nsets][Njet];
TH1D *hRecoSVD2DSub[Nsets][Njet];
TH1D *hRecoSVDConst2DSub[Nsets][Njet];
TH1D *hRecoSVDConst[Nsets][Njet];
TH1D *hRecoSVDConstSub[Nsets][Njet];
TH1D *hJtUnscaled[Nsets][Njet];
TH2D *hJtResponse[Nsets][Njet];
TH2D *hJetPtResponse[Nsets];
TH1D *hJetPtReco[Nsets];
TH1D *hJetPtRecoSVD[Nsets];
TH1D *hJetPtReco2[Nsets];

TH1D *hRecoBayes_unscaled[Nsets][Njet];
TH1D *hRecoBayesConst_unscaled[Nsets][Njet];
TH1D *hRecoSVD_unscaled[Nsets][Njet];
TH1D *hRecoSVDConst_unscaled[Nsets][Njet];
TH1D *hRecoBayesSub_unscaled[Nsets][Njet];
TH1D *hRecoBayesConstSub_unscaled[Nsets][Njet];
TH1D *hRecoSVDSub_unscaled[Nsets][Njet];
TH1D *hRecoSVDConstSub_unscaled[Nsets][Njet];

TH1D *hRecoBayesSub2D[Nsets][Njet];
TH1D *hRecoSVDSub2D[Nsets][Njet];

TH1D *hJtMeasSub[Nsets][Njet];
TH1D *hJtMeasConstSub[Nsets][Njet];
TH1D *hJtMeasSub_unscaled[Nsets][Njet];
TH1D *hJtMeasConstSub_unscaled[Nsets][Njet];
TH1D *hRecoBayesSub[Nsets][Njet];
TH1D *hRecoBayesConstSub[Nsets][Njet];

TH1D *hBgTrkNumberBinMeas[Nsets][Njet];
TH1D *hBgTrkNumberBinTrue[Nsets][Njet];
TH1D *hBgTrkNumberBinData[Nsets][Njet];


RooUnfoldResponse *jTresponse[Nsets][Njet];
RooUnfoldResponse *jTresponseConst[Nsets][Njet];

TFile *inFile, *outFile,*inFile2,*toyFile,*dataFile;


//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= gRandom->Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= gRandom->Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}

void loadBg(){
  for(int is = 0; is < Nsets; is++){
    for(int ij = 0; ij < Njet; ij++){
      if(doWeight){
        name = Form("%s/BgJtWeightBin/BgJtWeightBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      }else{
        name = Form("%s/BgJtBin/BgJtBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      }
      cout << name << endl;
      hJtMeasBg[is][ij] = (TH1D*)inFile->Get(name);
      if(!hJtMeasBg[is][ij]) continue;
      hJtMeasBg[is][ij]->Rebin(rebin);
      hJtMeasBg[is][ij]->Print();

      name = Form("%s/BgTrkNumberBin/BgTrkNumberBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      hBgTrkNumberBinMeas[is][ij] = (TH1D*)inFile->Get(name);
      name = Form("%s/BgTrkNumberBin/BgTrkNumberBinNFin%02dJetPt%02d",dirname.Data(),is+nR*2,ij);
      hBgTrkNumberBinTrue[is][ij] = (TH1D*)inFile->Get(name);
      writeToFile("Pythia/Measured","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinMeas[is][ij]);
      writeToFile("Pythia/True","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinTrue[is][ij]);
      writeToFile("Pythia/BayesUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinMeas[is][ij]);
      writeToFile("Pythia/BayesSubUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinMeas[is][ij]);
      writeToFile("Pythia/SVDUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinMeas[is][ij]);
      writeToFile("Pythia/SVDSubUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinMeas[is][ij]);
      writeToFile("PythiaUnscaled/BayesSubUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinMeas[is][ij]);
      writeToFile("PythiaUnscaled/SVDSubUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinMeas[is][ij]);
      writeToFile("PythiaUnscaled/BayesUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinMeas[is][ij]);
      writeToFile("PythiaUnscaled/SVDUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinMeas[is][ij]);


      if(doWeight){
        name = Form("%s/BgJtWeightBin/BgJtWeightBinNFin%02dJetPt%02d",dirname.Data(),is+nR*2,ij);
      }else{
        name = Form("%s/BgJtBin/BgJtBinNFin%02dJetPt%02d",dirname.Data(),is+nR*2,ij);
      }
      hJtTrueBg[is][ij] = (TH1D*)inFile->Get(name);
      if(!hJtTrueBg[is][ij]) continue;
      hJtTrueBg[is][ij]->Rebin(rebin);


      /*if(doWeight){
        name = Form("%s/BgRndmJtBin/BgRndmJtBinNFin%02dJetPt%02d",dirname.Data(),is,ij); //FIXME
      }else{
        name = Form("%s/BgRndmJtBin/BgRndmJtBinNFin%02dJetPt%02d",dirname.Data(),is,ij); //FIXME
      }
      cout << name << endl;
      hJtMeasBgRndm[is][ij] = (TH1D*)inFile->Get(name);
      hJtMeasBgRndm[is][ij]->Print();
      hJtMeasBgRndm[is][ij]->Rebin(rebin);


      name = Form("%s/hNumber/hNumberNFin%02d",dirname.Data(),is);
      cout << name << endl;
      TH1D *hNumberMeas = (TH1D*)inFile->Get(name);
      hNumberMeas->Print();

      if(doWeight){
        name = Form("%s/BgRndmJtWeightBin/BgRndmJtWeightBinNFin%02dJetPt%02d",dirname.Data(),is+nR*2,ij);
      }else{
        name = Form("%s/BgRndmJtBin/BgRndmJtBinNFin%02dJetPt%02d",dirname.Data(),is+nR*2,ij);
      }
      hJtTrueBgRndm[is][ij] = (TH1D*)inFile->Get(name);
      hJtTrueBgRndm[is][ij]->Print();
      hJtTrueBgRndm[is][ij]->Rebin(rebin);
      name = Form("%s/hNumber/hNumberNFin%02d",dirname.Data(),is+nR*2);
      cout << name << endl;
      TH1D *hNumberTrue = (TH1D*)inFile->Get(name);
      hNumberTrue->Print();

      int N_rnmdBg_meas = hNumberMeas->GetBinContent(7+ij);
      int N_rnmdBg_true = hNumberTrue->GetBinContent(7+ij);

      hJtMeasBgRndm[is][ij]->Scale(1.0/N_rnmdBg_meas,"width");
      hJtMeasBgTrue[is][ij]->Scale(1.0/N_rnmdBg_true,"width");

      if(doWeight){
        writeToFile("Pythia/Measured","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/True","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtTrueBgRndm[is][ij]);
        writeToFile("Pythia/BayesUnfolding","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/SVDUnfolding","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding/WeightCorrected","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding/WeightCorrected","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
      }else{
        writeToFile("Pythia/Measured","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/True","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtTrueBgRndm[is][ij]);
        writeToFile("Pythia/BayesUnfolding","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/SVDUnfolding","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding/WeightCorrected","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding/WeightCorrected","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBgRndm[is][ij]);
      }*/
      

      n_jet_meas[is][ij] = hJetPtBinMeas[is][ij]->Integral();
      n_jet_true[is][ij] = hJetPtBinTrue[is][ij]->Integral();
      n_bg_meas[is][ij] = hBgTrkNumberBinMeas[is][ij]->Integral();
      n_bg_true[is][ij] = hBgTrkNumberBinTrue[is][ij]->Integral();

      if(is == 0){
        n_jet_pythia[ij] = hJetPtBinPythia[ij]->Integral();
      }

      hJtMeasBg[is][ij]->Scale(1.0/n_bg_meas[is][ij],option);
      hJtTrueBg[is][ij]->Scale(1.0/n_bg_true[is][ij],option);
      /*if(doWeight){
        writeToFile("Pythia/Measured","BgJtWeightBin",hJtMeasBg[is][ij]);
        writeToFile("Pythia/True","BgJtWeightBin",hJtTrueBg[is][ij]);
      }else{
        writeToFile("Pythia/Measured","BgJtBin",hJtMeasBg[is][ij]);
        writeToFile("Pythia/True","BgJtBin",hJtTrueBg[is][ij]);
      }*/

      cout << "is: " << is << " ij: " << ij << endl;
      cout << "N_jets, true: " << n_jet_true[is][ij] << " meas: " << n_jet_meas[is][ij] << endl;
      cout << "N_bg, true: " << n_bg_true[is][ij] << " meas: " << n_bg_meas[is][ij] << endl;

      if(doWeight){
        writeToFile("Pythia/Measured","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/True","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtTrueBg[is][ij]);
        writeToFile("Pythia/BayesUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/2DUnfoldingBayes","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/2DUnfoldingBayesSub","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVD","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVDSub","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/SVDUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("PythiaUnscaled/BayesSubUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("PythiaUnscaled/SVDSubUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("PythiaUnscaled/BayesUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("PythiaUnscaled/SVDUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
      }else{
        writeToFile("Pythia/Measured","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/True","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtTrueBg[is][ij]);
        writeToFile("Pythia/BayesUnfolding","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/2DUnfoldingBayes","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/2DUnfoldingBayesSub","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVD","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVDSub","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/SVDUnfolding","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding/WeightCorrected","BgJtBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding/WeightCorrected","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtMeasBg[is][ij]);
      }
    }
  }
}

void loadData(){
  for(int is = 0; is < Nsets; is++){
    name = Form("%s/JetPt/JetPtNFin%02d",dirnameData.Data(),is);
    hJetPtData[is]= (TH1D*)dataFile->Get(name);
    for(int ij = 0; ij < Njet; ij++){
      cout << "is: " << is << " ij: " << ij << endl;
      name = Form("%s/JetPtBin/JetPtBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij);
      hJetPtBinData[is][ij]= (TH1D*)dataFile->Get(name);
      hJetPtBinData[is][ij]->Print();
      if(!hJetPtBinData[is][ij]) continue;
      double N_jetData = hJetPtBinData[is][ij]->GetEntries();
      cout << "Data jets: " << N_jetData << endl;
      writeToFile("Data/Measured","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinData[is][ij]);
      writeToFile("Data/BayesUnfolding","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinData[is][ij]);
      writeToFile("Data/BayesSubUnfolding","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinData[is][ij]);
      writeToFile("Data/SVDSubUnfolding","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinData[is][ij]);
      writeToFile("Data/SVDUnfolding","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinData[is][ij]);
      writeToFile("Data/BayesSubUnfolding/WeightCorrected","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinData[is][ij]);
      writeToFile("Data/SVDSubUnfolding/WeightCorrected","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinData[is][ij]);

      if(doWeight){
        name = Form("%s/BgJtWeightBin/BgJtWeightBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij);
      }else{
        name = Form("%s/BgJtBin/BgJtBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij);
      }
      cout << name << endl;
      hJtBgData[is][ij] = (TH1D*)dataFile->Get(name);
      if(!hJtBgData[is][ij]) continue;
      hJtBgData[is][ij]->Print();
      hJtBgData[is][ij]->Rebin(rebin);

      if(doWeight){
        //name = Form("%s/BgRndmJtWeightBin/BgRndmJtWeightBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij); //FIXME
        name = Form("%s/BgRndmJtBin/BgRndmJtBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij); //FIXME
      }else{
        name = Form("%s/BgRndmJtBin/BgRndmJtBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij); //FIXME
      }
      cout << name << endl;
      hJtBgRndmData[is][ij] = (TH1D*)dataFile->Get(name);
      hJtBgRndmData[is][ij]->Print();
      hJtBgRndmData[is][ij]->Rebin(rebin);
      
      name = Form("%s/hNumber/hNumberNFin%02d",dirnameData.Data(),is);
      cout << name << endl;
      TH1D *hNumber = (TH1D*)dataFile->Get(name);
      hNumber->Print();
      int N_rnmdBg = hNumber->GetBinContent(7+ij);
      hJtBgRndmData[is][ij]->Scale(1.0/N_rnmdBg);

      if(doWeight){
        writeToFile("Data/Measured","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/BayesUnfolding","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/BayesSubUnfolding","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/SVDUnfolding","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/SVDSubUnfolding","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/BayesSubUnfolding/WeightCorrected","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/SVDSubUnfolding/WeightCorrected","BgRndmJtWeightBin",Form("BgRndmJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
      }else{
        writeToFile("Data/Measured","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/BayesUnfolding","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/BayesSubUnfolding","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/SVDUnfolding","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/SVDSubUnfolding","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/BayesSubUnfolding/WeightCorrected","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
        writeToFile("Data/SVDSubUnfolding/WeightCorrected","BgRndmJtBin",Form("BgRndmJtBinNFin%02dJetPt%02d",is,ij),hJtBgRndmData[is][ij]);
      }

      name = Form("%s/BgTrkNumberBin/BgTrkNumberBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij);
      hBgTrkNumberBinData[is][ij] = (TH1D*)dataFile->Get(name);
      writeToFile("Data/Measured","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinData[is][ij]);
      writeToFile("Data/BayesUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinData[is][ij]);
      writeToFile("Data/BayesSubUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinData[is][ij]);
      writeToFile("Data/SVDUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinData[is][ij]);
      writeToFile("Data/SVDSubUnfolding","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinData[is][ij]);
      writeToFile("Data/BayesSubUnfolding/WeightCorrected","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinData[is][ij]);
      writeToFile("Data/SVDSubUnfolding/WeightCorrected","BgTrkNumberBin",Form("BgTrkNumberBinNFin%02dJetPt%02d",is,ij),hBgTrkNumberBinData[is][ij]);

      n_jet_data[is][ij] = hJetPtBinData[is][ij]->GetEntries();
      n_bg_data[is][ij] = hBgTrkNumberBinData[is][ij]->GetEntries();

      hJtBgData[is][ij]->Scale(1.0/n_bg_data[is][ij],option);
      if(doWeight){
        writeToFile("Data/Measured","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/BayesUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/BayesSubUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/SVDUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/SVDSubUnfolding","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/BayesSubUnfolding/WeightCorrected","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/SVDSubUnfolding/WeightCorrected","BgJtWeightBin",Form("BgJtWeightBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
      }else{
        writeToFile("Data/Measured","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/BayesUnfolding","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/BayesSubUnfolding","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/SVDUnfolding","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/SVDSubUnfolding","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/BayesSubUnfolding/WeightCorrected","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
        writeToFile("Data/SVDSubUnfolding/WeightCorrected","BgJtBin",Form("BgJtBinNFin%02dJetPt%02d",is,ij),hJtBgData[is][ij]);
      }

      if(doWeight){
        name = Form("%s/JtWeightBin/JtWeightBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij);
      }else{
        name = Form("%s/JtBin/JtBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij);
      }
      hJtDataConst[is][ij] = (TH1D*)dataFile->Get(name);
      if(!hJtDataConst[is][ij]) continue;
      hJtDataConst[is][ij]->Rebin(rebin);

      if(doWeight){
        name = Form("%s/JetConeJtWeightBin/JetConeJtWeightBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij);
      }else{
        name = Form("%s/JetConeJtBin/JetConeJtBinNFin%02dJetPt%02d",dirnameData.Data(),is,ij);
      }
      hJtData[is][ij] = (TH1D*)dataFile->Get(name);
      if(!hJtData[is][ij]) continue;
      hJtData[is][ij]->Rebin(rebin);

      hJtData[is][ij]->Scale(1.0/n_jet_data[is][ij],option);
      hJtDataConst[is][ij]->Scale(1.0/n_jet_data[is][ij],option);

      hJtDataSub[is][ij] = (TH1D*)hJtData[is][ij]->Clone();
      hJtDataConstSub[is][ij] = (TH1D*)hJtDataConst[is][ij]->Clone();

      hJtUnfBgData[is][ij] =  (TH1D*)hJtUnfBg[is][ij]->Clone();
      //hJtUnfBgData[is][ij]->Scale(hJtData[is][ij]->GetEntries()/hJtMeas[is][ij]->GetEntries());
      //hJtUnfBgData[is][ij]->Scale(n_jet_data[is][ij]/n_jet_meas[is][ij]);

      hJtDataSub[is][ij]->Add(hJtUnfBgData[is][ij],-1);
      hJtDataConstSub[is][ij]->Add(hJtUnfBgData[is][ij],-1);

      if(doSignal){ //FIXME if signal used
        TH1D* hJtDataIncl = (TH1D*)hJtData[is][ij]->Clone();
        hJtDataIncl->SetName(Form("DataInclNFin%02dJetPt%02d",is,ij));
        hJtDataIncl->Write();

        TH1D* hJtDataConstIncl = (TH1D*)hJtDataConst[is][ij]->Clone();
        hJtDataConstIncl->SetName(Form("DataConstInclNFin%02dJetPt%02d",is,ij));
        hJtDataConstIncl->Write();
        cout << "is: " << is << " ij: " << ij << " hJtData integral before substraction: " << hJtData[is][ij]->Integral() << endl;
        //hJtData[is][ij]->Add(hJtDataBg[is][ij],-1); //TODO 
        cout << "is: " << is << " ij: " << ij << " hJtData integral after substraction: " << hJtData[is][ij]->Integral() << endl;
        //hJtDataConst[is][ij]->Add(hJtDataBg[is][ij],-1); //TODO
      }

      outFile->cd();
      if(doWeight){
        writeToFile("Data/Measured","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),hJtDataConst[is][ij]);
        writeToFile("Data/Measured","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hJtData[is][ij]);
        writeToFile("Data/Sub","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),hJtDataConstSub[is][ij]);
        writeToFile("Data/Sub","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hJtDataSub[is][ij]);
      }else{
        writeToFile("Data/Measured","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hJtDataConst[is][ij]);
        writeToFile("Data/Measured","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hJtData[is][ij]);
        writeToFile("Data/Sub","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hJtDataConstSub[is][ij]);
        writeToFile("Data/Sub","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hJtDataSub[is][ij]);
      }
    }
  }
}

void UnfoldToy(){
  TH1D *hJtMeasToy;
  TH1D *hJtTrueToy;
  TH1D *hJtMeasBinToy[Njet];
  TH1D *hJtTrueBinToy[Njet];
  TH2D *responseMatrixBinToy[Njet];
  TH2D *responseMatrixToy;
  TH1D *hJetPtBinMeasToy[Njet];
  TH1D *hJetPtBinTrueToy[Njet];
  TH1D *hRecoBayes1DBinToy[Njet];
  TH1D *hRecoBayes2DBinToy[Njet];
  TH1D *hRecoSVD2DBinToy[Njet];


  TString name = "events";
  TH1D *hEvents= (TH1D*)toyFile->Get(name);
  Int_t nEvents = hEvents->GetBinContent(1);
  cout << "Number of events: " << nEvents << endl;
  Int_t number_jets = hEvents->GetBinContent(2);
  cout << "Number of jets: " << hEvents->GetBinContent(2) << endl;
  cout << "Number of jets: " << number_jets << endl;

  name = Form("MultiTrueNFin00");
  TH1D* hMultiTrue = (TH1D*)toyFile->Get(name);
  cout << "Number of true jets: " << hMultiTrue->GetEntries() << endl;
  int is = 0;
  name = Form("MeasuredConstNFin00");
  hJtMeasToy = (TH1D*)toyFile->Get(name);

  name = Form("TrueConstNFin00");
  hJtTrueToy = (TH1D*)toyFile->Get(name);

  name = Form("MeasuredCorrConstNFin00");
  hJtMeasCorrToy = (TH1D*)toyFile->Get(name);

  name = Form("TrueCorrConstNFin00");
  hJtTrueCorrToy = (TH1D*)toyFile->Get(name);

  name = Form("ResponseMatrixNFin00");
  responseMatrixToy = (TH2D*)toyFile->Get(name);
  name = Form("ResponseMatrixCorrNFin00");
  TH2D *responseMatrixCorrToy = (TH2D*)toyFile->Get(name);

  hJtMeasToy->Scale(1.0/number_jets);
  hJtTrueToy->Scale(1.0/number_jets);
  hJtMeasCorrToy->Scale(1.0/number_jets);
  hJtTrueCorrToy->Scale(1.0/number_jets);

  RooUnfoldResponse *response  = CreateResponse(hJtMeasToy,responseMatrixToy);
  RooUnfoldBayes *unfoldBayesToy = new RooUnfoldBayes(response, hJtMeasToy, 4,false); 
  RooUnfoldSvd *unfoldSVDToy = new RooUnfoldSvd(response, hJtMeasToy, 20);  

  RooUnfoldResponse *responseCorr  = CreateResponse(hJtMeasCorrToy,responseMatrixCorrToy);
  RooUnfoldBayes *unfoldBayesCorrToy = new RooUnfoldBayes(responseCorr, hJtMeasCorrToy, 4,false); 
  RooUnfoldSvd *unfoldSVDCorrToy = new RooUnfoldSvd(responseCorr, hJtMeasCorrToy, 20);  

  TH1D *hRecoBayesToy = (TH1D*) unfoldBayesToy->Hreco();
  TH1D *hRecoSVDToy = (TH1D*) unfoldSVDToy->Hreco();

  TH1D *hRecoBayesCorrToy = (TH1D*) unfoldBayesCorrToy->Hreco();
  TH1D *hRecoSVDCorrToy = (TH1D*) unfoldSVDCorrToy->Hreco();

  hRecoBayesToy->SetName(Form("UnfoldedBayesToyNFin%02d",is));
  hRecoBayesToy->Print();
  hRecoSVDToy->SetName(Form("UnfoldedSVDToyNFin%02d",is));
  hRecoSVDToy->Print();

  hRecoBayesCorrToy->SetName(Form("UnfoldedBayesCorrToyNFin%02d",is));
  hRecoBayesCorrToy->Print();
  hRecoSVDCorrToy->SetName(Form("UnfoldedSVDCorrToyNFin%02d",is));
  hRecoSVDCorrToy->Print();

  responseMatrixToy->SetName(Form("responseMatrixToyNFin%02d",is));
  responseMatrixToy->Print();
  TH1D *hJtFoldedToy = foldDist(hJtTrueToy,responseMatrixToy,1);
  hJtFoldedToy->SetName(Form("FoldedToyNFin%02d",is));
  TH1D *hJtFoldedCorrToy = foldDist(hJtTrueCorrToy,responseMatrixCorrToy,1);
  hJtFoldedCorrToy->SetName(Form("FoldedCorrToyNFin%02d",is));
  outFile->cd();
  hJtFoldedToy->Write();
  hJtFoldedCorrToy->Write();
  hJtMeasToy->Write();
  hJtTrueToy->Write();
  hJtMeasCorrToy->Write();
  hJtTrueCorrToy->Write();
  hRecoBayesToy->Write();
  hRecoSVDToy->Write();
  hRecoBayesCorrToy->Write();
  hRecoSVDCorrToy->Write();
  responseMatrixToy->Write();

  for(int ij = 0; ij < Njet; ij++){
    name = Form("MeasuredConstBinNFin00JetPt%02d",ij);
    hJtMeasBinToy[ij] = (TH1D*)toyFile->Get(name);
    hJtMeasBinToy[ij]->Print();
    name = Form("TrueConstBinNFin00JetPt%02d",ij);
    hJtTrueBinToy[ij] = (TH1D*)toyFile->Get(name);
    name = Form("ResponseMatrixBinNFin00JetPt%02d",ij);
    responseMatrixBinToy[ij] = (TH2D*)toyFile->Get(name);

    name = Form("JetPtMeasBinNFin%02dJetPt%02d",is,ij);
    hJetPtBinMeasToy[ij]= (TH1D*)toyFile->Get(name);

    name = Form("JetPtTrueBinNFin%02dJetPt%02d",is,ij);
    hJetPtBinTrueToy[ij]= (TH1D*)toyFile->Get(name);
    hJetPtBinTrueToy[ij]->Print();

    double N_jetTrue = hJetPtBinTrueToy[ij]->Integral();
    double N_jetMeas = hJetPtBinMeasToy[ij]->Integral();
    hJtMeasBinToy[ij]->Scale(1.0/N_jetTrue);
    hJtTrueBinToy[ij]->Scale(1.0/N_jetMeas);
    RooUnfoldResponse *response  = CreateResponse(hJtMeasBinToy[ij],responseMatrixBinToy[ij]);
    RooUnfoldBayes *unfoldBayesToy = new RooUnfoldBayes(response, hJtMeasBinToy[ij], 4,false); 
    RooUnfoldSvd *unfoldSVDToy = new RooUnfoldSvd(response, hJtMeasBinToy[ij], 20);  

    TH1D *hRecoBayesToy = (TH1D*) unfoldBayesToy->Hreco();
    TH1D *hRecoSVDToy = (TH1D*) unfoldSVDToy->Hreco();
    cout << "Measured integral: " << hJtMeasBinToy[ij]->Integral() << " True integral: " << hJtTrueBinToy[ij]->Integral() << " Unfolded integral: " << hRecoBayesToy->Integral() << endl;


    hRecoBayesToy->SetName(Form("UnfoldedBayesToyNFin%02dJetPt%02d",is,ij));
    hRecoBayesToy->Print();
    hRecoSVDToy->SetName(Form("UnfoldedSVDToyNFin%02dJetPt%02d",is,ij));
    hRecoSVDToy->Print();
    responseMatrixBinToy[ij]->SetName(Form("responseMatrixBinToyNFin%02dJetPt%02d",is,ij));
    responseMatrixBinToy[ij]->Print();
    TH1D *hJtFoldedToy = foldDist(hJtTrueBinToy[ij],responseMatrixBinToy[ij],1);
    hJtFoldedToy->SetName(Form("FoldedBinToyNFin%02dJetPt%02d",is,ij));
    outFile->cd();
    hJtFoldedToy->Write();
    hJtMeasBinToy[ij]->Write();
    hJtTrueBinToy[ij]->Write();
    hRecoBayesToy->Write();
    hRecoSVDToy->Write();
    responseMatrixBinToy[ij]->Write();
  }

  name = Form("JetPtMeasNFin%02d",is);
  TH1D *hJetPtMeasToy= (TH1D*)toyFile->Get(name);
  name = Form("JetPtTrueNFin%02d",is);
  TH1D *hJetPtTrueToy= (TH1D*)toyFile->Get(name);

  outFile->cd();
  hJetPtMeasToy->Write();
  hJetPtTrueToy->Write();

  name = Form("JetPtCorrCoarseNFin%02d",is);
  TH2D *responseMatrixJetPtCoarseToy = (TH2D*)toyFile->Get(name);

  cout << "==================================== RESPONSE ================================" << endl;
  TH2D *responseCoarseToy = (TH2D*)responseMatrixJetPtCoarseToy->Clone();
  outFile->cd();
  double ptTrue, ptObs;
  double NTrue, NObs, Ntot;
  double binW;
  int ib,ib2;
  TH1D *projX = (TH1D*)responseCoarseToy->ProjectionX();
  double eff[10];
  double corr[10];
  TH1D *hJetPtUnfoldedToy = (TH1D*)projX->Clone();
  TH1D *hJetPtMeasCoarseToy = (TH1D*)projX->Clone();
  TH1D *hJetPtTrueCoarseToy = (TH1D*)projX->Clone();
  hJetPtMeasCoarseToy->Reset();
  hJetPtTrueCoarseToy->Reset();
  hJetPtUnfoldedToy->Reset();
  for(int ij = 0; ij < Njet; ij++){
    hJetPtMeasCoarseToy->SetBinContent(ij+1,hJetPtBinMeasToy[ij]->Integral());
    hJetPtTrueCoarseToy->SetBinContent(ij+1,hJetPtBinTrueToy[ij]->Integral());
  }
  outFile->cd();
  hJetPtMeasCoarseToy->SetName(Form("hJetPtMeasCoarseToy%02d",is));
  hJetPtTrueCoarseToy->SetName(Form("hJetPtTrueCoarseToy%02d",is));
  hJetPtMeasCoarseToy->Write();
  hJetPtTrueCoarseToy->Write();
  double factor;
  double sumP;
  //responseCoarseToy->Write();

  //First create P(E|C) i.e. the actual normalized response matrix from correlation matrix, also the efficiency, i.e. N(True jets in bin)/N(Observed jets in any bin corresponding to true jets)
  //corr = 1/efficiency
  for(int ibx = 1 ; ibx <= responseCoarseToy->GetNbinsX() ; ibx++){
    hJetPtUnfoldedToy->SetBinContent(ibx,1); //First create uniform matrix for initial guess
    ptTrue = responseCoarseToy->GetXaxis()->GetBinCenter(ibx); //True pT is on X axis
    binW = responseCoarseToy->GetXaxis()->GetBinWidth(ibx);
    cout << "True pT: " << ptTrue - binW/2 <<"-" << ptTrue + binW/2;
    ib = responseCoarseToy->GetBin(ibx,0);
    Ntrue  = projX->GetBinContent(ibx);
    cout << " True N: " << Ntrue << endl;
    cout << " Missed N: " << responseCoarseToy->GetBinContent(ib) << endl;
    Ntot = 0;
    if(Ntrue > 0){
      cout << " fraction missed: " << responseCoarseToy->GetBinContent(ib)/Ntrue << endl;
      for(int iby = 1 ; iby <= responseCoarseToy->GetNbinsY(); iby++){
        ptObs = responseCoarseToy->GetYaxis()->GetBinCenter(iby);
        binW = responseCoarseToy->GetYaxis()->GetBinWidth(iby);
        cout << "Observed pT: " << ptObs - binW/2 <<"-" << ptObs + binW/2;
        ib = responseCoarseToy->GetBin(ibx,iby);
        NObs = responseCoarseToy->GetBinContent(ib);
        Ntot += NObs;
        cout << ", N: " << NObs <<  ", P(E|C): " << 1.0*NObs/Ntrue <<endl;
        //if(N > 0) cout << "Fill bin (" << jtobs << "," << jt << ") with " << N << endl;
        //response->Fill(jtobs,jt,N);
        responseCoarseToy->SetBinContent(ib,1.0*NObs/Ntrue);
      }
    }else{
      break;
    }
    cout << "Total N: " << Ntot;
    cout << ", Efficiency: " << Ntot/Ntrue;
    eff[ibx-1] = Ntot/Ntrue;
    corr[ibx-1] = Ntrue/Ntot;
    cout << ", Correction: " << corr[ibx-1] << endl;
  }
  responseCoarseToy->SetTitle(Form("Probability of effect Y given Cause X is=%02d",is));
  responseCoarseToy->SetName("responseCoarseToy");
  responseCoarseToy->Write();
  outFile->cd();
  TH2D *inverseP = (TH2D*)responseCoarseToy->Clone();
  inverseP->Reset();
  //inverseP->SetName(Form("inverseP%02d",is));
  inverseP->SetTitle("Probability of cause X given observed Y");


  //Start the iteration loop
  int ib3;
  double sum = 0;
  for(int it = 0; it < 10 ;it++){
    hJetPtUnfoldedToy->Scale(1.0/hJetPtUnfoldedToy->Integral()); //Scale the unfolded histogram so that it gives the probability
    //First create P(C|E), labeled inverseP
    for(int ij = 0; ij < Njet; ij++){
      for(int ij2 = 0; ij2 < Njet; ij2++){
        ib = inverseP->GetBin(ij+1,ij2+1);  //Get bin x,y
        ib2 = responseCoarseToy->GetBin(ij+1,ij2+1); //Get bin x,y
        if(ib != ib2) cout << "ib != ib2" << endl << endl;
        sum = 0; //Sum over z P(E_y | C_z)*P( C_z )
        for(int ij3 = 0; ij3 < Njet; ij3++){
          ib3 = responseCoarseToy->GetBin(ij3+1,ij2+1); //Bin z,y
          sum += responseCoarseToy->GetBinContent(ib3)*hJetPtUnfoldedToy->GetBinContent(ij3+1); // P(E_y | C_z) * P( C_z)
        }
        cout << "Sum is: " << sum << endl;
        double PEy_Cx = responseCoarseToy->GetBinContent(ib2); //Get P(E_y | C_x)
        double PCx = hJetPtUnfoldedToy->GetBinContent(ij+1); //Get P( C_x)
        cout << "responseCoarseToy->GetBinContent(ib2) gives " << PEy_Cx << endl; 
        cout << "hJetPtUnfoldedToy->GetBinContent(ij+1) gives " << PCx << endl;
        factor = PEy_Cx * PCx / sum; //Factor is P(E_y | C_x) * P( C_x) / sum
        cout << "Set P(C" << ij << "| E" << ij2 << ") to " << factor << endl;
        inverseP->SetBinContent(ib,factor); // Set bin content
      }
    } 
    outFile->cd();
    inverseP->SetName(Form("inversePToy%02d",it));
    inverseP->Write();
    double sumP = 0;
    for(int ij2 = 0; ij2 < Njet; ij2++){
      sumP = 0;
      for(int ij = 0; ij < Njet; ij++){
        ib = inverseP->GetBin(ij+1,ij2+1);
        sumP +=   inverseP->GetBinContent(ib);
      }
      cout << "ij2: " << ij2 << " sum of P is: " << sumP << endl;
    }

    //Then calculate new distribution from P(C | E)
    hJetPtUnfoldedToy->Reset();
    for(int ij = 0; ij < Njet ; ij++){
      cout << "Unfolding ij: " << ij << endl;
      hRecoBayes1DBinToy[ij] = (TH1D*)hJtMeasBinToy[ij]->Clone(); //jT unfolding
      hRecoBayes1DBinToy[ij]->Reset(); //jT unfolding
      hRecoBayes1DBinToy[ij]->SetName(Form("hRecoBayes1DBinToyNFin%02dJetPt%02dIt%02d",is,ij,it)); //jT unfolding
      hRecoBayes1DBinToy[ij]->Print(); //jT unfolding
      for(int ij2 = 0 ; ij2 < Njet; ij2++){
        ib =  inverseP->GetBin(1+ij,1+ij2);  //Get bin x,y
        hRecoBayes1DBinToy[ij]->Add(hJtMeasBinToy[ij2],corr[ij]*inverseP->GetBinContent(ib)); 
        //cout << "Add content to hJetPtUnfolded: " << corr[ij] << " * " << inverseP->GetBinContent(ib) << " * " << hJetPtBinMeas[is][ij2]->GetEntries() << " = " << corr[ij]*inverseP->GetBinContent(ib)*hJetPtBinMeas[is][ij2]->GetEntries() << endl;
        hJetPtUnfoldedToy->AddBinContent(ij+1,corr[ij]*inverseP->GetBinContent(ib)*hJetPtMeasCoarseToy->GetBinContent(ij2+1)); //Bin of new histogram is sum over old bins j 1/eff * n(E_j) * P(C_x | E_j)

        //cout << "Histogram ij2: " << ij2 << " was scaled by " << corr[ij] << " * " << inverseP->GetBinContent(ib) << " = " << corr[ij]*inverseP->GetBinContent(ib) << endl;
      }
      cout << "Total content added: " << hJetPtUnfoldedToy->GetBinContent(ij+1) << endl;
      hRecoBayes1DBinToy[ij]->Write();
    }
    hJetPtUnfoldedToy->SetName(Form("hJetPtUnfolded2DToyIt%02d",it));
    outFile->cd();
    hJetPtUnfoldedToy->Write();
  }


  for(int ij = 0 ; ij < Njet; ij++){
    RooUnfoldResponse *response  = CreateResponse(hRecoBayes1DBinToy[ij],responseMatrixBinToy[ij]);
    RooUnfoldBayes *unfoldBayesJetCone = new RooUnfoldBayes(response, hRecoBayes1DBinToy[ij], 4,false); 
    RooUnfoldSvd *unfoldSVDJetCone = new RooUnfoldSvd(response, hRecoBayes1DBinToy[ij], 20);  

    cout << "is: " << is  << " ij: " << ij << endl;
    hRecoBayes2DBinToy[ij] = (TH1D*) unfoldBayesJetCone->Hreco();
    hRecoSVD2DBinToy[ij] = (TH1D*) unfoldSVDJetCone->Hreco();

    hRecoBayes2DBinToy[ij]->SetName(Form("UnfoldedBayes2DToyNFin%02dJetPt%02d",is,ij));
    hRecoBayes2DBinToy[ij]->Print();
    hRecoSVD2DBinToy[ij]->SetName(Form("UnfoldedSVD2ToyDNFin%02dJetPt%02d",is,ij));
    hRecoSVD2DBinToy[ij]->Print();
    outFile->cd();
    hRecoBayes2DBinToy[ij]->Write();
    hRecoSVD2DBinToy[ij]->Write();
  }
}

void loadJt(){
  TString name;
  for(int is = 0; is < Nsets; is++){
    for(int ij = 0; ij < Njet; ij++){
      if(is == 0){
        if(doWeight){
          name = Form("AliJJetJtTask/AliJJetJtMCHistManager/JtWeightBinPythia/JtWeightBinPythiaJetPt%02d",ij);
        }else{
          name = Form("AliJJetJtTask/AliJJetJtMCHistManager/JtBinPythia/JtBinPythiaJetPt%02d",ij);
        }
        hJtPythia[ij]= (TH1D*)inFile2->Get(name);
        hJtPythia[ij]->Rebin(rebin);
        outFile->cd();
        hJtPythia[ij]->SetName(Form("TruePythiaJetPt%02d",ij));
        hJtPythia[ij]->Print();
        hJtPythia[ij]->Scale(1.0/n_jet_pythia[ij],option);
        if(doWeight){
          writeToFile("Pythia/Pythia","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",0,ij),hJtPythia[ij]);
        }else{
          writeToFile("Pythia/Pythia","JtWeightBin",Form("JtBinNFin%02dJetPt%02d",0,ij),hJtPythia[ij]);
        }
      }
      cout << "is: " << is << " ij: " << ij << endl;
      if(doWeight){
        name = Form("%s/JtWeightBin/JtWeightBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      }else{
        name = Form("%s/JtBin/JtBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      }
      hJtMeasConst[is][ij] = (TH1D*)inFile->Get(name);
      if(!hJtMeasConst[is][ij]) continue;
      hJtMeasConst[is][ij]->Rebin(rebin);
      name = Form("%s/JtBin/JtBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      cout << name << endl;
      hJtMeasConst_unscaled[is][ij] = (TH1D*)inFile->Get(name);
      hJtMeasConst_unscaled[is][ij]->Rebin(rebin);
      hJtMeasConst_unscaled[is][ij]->Print();

      if(doWeight){
        name = Form("%s/JetConeJtWeightBin/JetConeJtWeightBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      }else{
        name = Form("%s/JetConeJtBin/JetConeJtBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      }
      cout << name << endl;
      hJtMeas[is][ij] = (TH1D*)inFile->Get(name);
      if(!hJtMeas[is][ij]) continue;
      hJtMeas[is][ij]->Rebin(rebin);
      name = Form("%s/JetConeJtBin/JetConeJtBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      hJtMeas_unscaled[is][ij] = (TH1D*)inFile->Get(name);
      hJtMeas_unscaled[is][ij]->Rebin(rebin);

      if(doWeight){
        name = Form("%s/JtWeightBin/JtWeightBinNFin%02dJetPt%02d",dirname.Data(),is+nR*2,ij);
      }else{
        name = Form("%s/JtBin/JtBinNFin%02dJetPt%02d",dirname.Data(),is+nR*2,ij);
      }
      hJtTrueConst[is][ij]= (TH1D*)inFile->Get(name);
      //hJtTrueConst[is][ij]->Print();
      hJtTrueConst[is][ij]->Rebin(rebin);

      if(doWeight){
        name = Form("%s/JetConeJtWeightBin/JetConeJtWeightBinNFin%02dJetPt%02d",dirname.Data(),is+2*nR,ij);
      }else{
        name = Form("%s/JetConeJtBin/JetConeJtBinNFin%02dJetPt%02d",dirname.Data(),is+2*nR,ij);
      }
      hJtTrue[is][ij]= (TH1D*)inFile->Get(name);
      hJtTrue[is][ij]->Rebin(rebin);
      if(!hJtTrue[is][ij]) continue;
      //hJtTrue[is][ij]->Print();

      cout << "hJtMeas Integral: " << hJtMeas[is][ij]->Integral() << endl;
      cout << "n_jet_meas: " << n_jet_meas[is][ij] << endl;
      hJtMeas[is][ij]->Scale(1.0/n_jet_meas[is][ij],option);
      hJtMeasConst[is][ij]->Scale(1.0/n_jet_meas[is][ij],option);
      hJtTrue[is][ij]->Scale(1.0/n_jet_true[is][ij],option);
      hJtTrueConst[is][ij]->Scale(1.0/n_jet_true[is][ij],option);
      hJtMeas_unscaled[is][ij]->Scale(1.0/n_jet_meas[is][ij],option);
      hJtMeasConst_unscaled[is][ij]->Scale(1.0/n_jet_meas[is][ij],option);

      hJtMeasSub[is][ij] = (TH1D*)hJtMeas[is][ij]->Clone();
      hJtMeasConstSub[is][ij] = (TH1D*)hJtMeasConst[is][ij]->Clone();
      hJtMeasSub_unscaled[is][ij] = (TH1D*)hJtMeas_unscaled[is][ij]->Clone();
      hJtMeasConstSub_unscaled[is][ij] = (TH1D*)hJtMeasConst_unscaled[is][ij]->Clone();

      if(doWeight){
        name = Form("AliJJetJtTask/AliJJetJtMCHistManager/JetConeJtWeightBinUnfBg/JetConeJtWeightBinUnfBgNFin%02dJetPt%02d",is,ij);
      }else{
        name = Form("AliJJetJtTask/AliJJetJtMCHistManager/JetConeJtBinUnfBg/JetConeJtBinUnfBgNFin%02dJetPt%02d",is,ij);
      }
      hJtUnfBg[is][ij]= (TH1D*)inFile->Get(name);
      name = Form("AliJJetJtTask/AliJJetJtMCHistManager/JetConeJtBinUnfBg/JetConeJtBinUnfBgNFin%02dJetPt%02d",is,ij);
      hJtUnfBg_unscaled[is][ij]= (TH1D*)inFile->Get(name);
      if(hJtUnfBg[is][ij]){
        hJtUnfBg[is][ij]->Rebin(rebin);
        hJtUnfBg[is][ij]->Scale(1.0/n_jet_meas[is][ij],option);
        hJtUnfBg_unscaled[is][ij]->Rebin(rebin);
        hJtUnfBg_unscaled[is][ij]->Scale(1.0/n_jet_meas[is][ij],option);
        hJtMeasConstSub[is][ij]->Add(hJtUnfBg[is][ij],-1);
        hJtMeasSub[is][ij]->Add(hJtUnfBg[is][ij],-1);
        hJtMeasSub_unscaled[is][ij]->Add(hJtUnfBg_unscaled[is][ij],-1);
        hJtMeasConstSub_unscaled[is][ij]->Add(hJtUnfBg_unscaled[is][ij],-1);
      }else{
        cout << name << " Not found " << endl;
      }

      cout << "ij: " << ij << " hJtMeas Entries: " << hJtMeas[is][ij]->GetEntries() << " hJtUnfBg Entries: " << hJtUnfBg[is][ij]->GetEntries() << endl;
      cout << "ij: " << ij << " hJtMeas Integral: " << hJtMeas[is][ij]->Integral();
      cout << " hJtUnfBg Integral: " << hJtUnfBg[is][ij]->Integral() << endl;
      if(doSignal){ //FIXME IF signal used
        TH1D* hJtMeasIncl = (TH1D*)hJtMeas[is][ij]->Clone();
        hJtMeasIncl->SetName(Form("MeasuredInclNFin%02dJetPt%02d",is,ij));
        hJtMeasIncl->Write();
        TH1D* hJtTrueIncl = (TH1D*)hJtTrue[is][ij]->Clone();
        hJtTrueIncl->SetName(Form("TrueInclNFin%02dJetPt%02d",is,ij));
        hJtTrueIncl->Write();
        TH1D* hJtMeasConstIncl = (TH1D*)hJtMeasConst[is][ij]->Clone();
        hJtMeasConstIncl->SetName(Form("MeasuredConstInclNFin%02dJetPt%02d",is,ij));
        hJtMeasConstIncl->Write();
        TH1D* hJtTrueConstIncl = (TH1D*)hJtTrueConst[is][ij]->Clone();
        hJtTrueConstIncl->SetName(Form("TrueConstInclNFin%02dJetPt%02d",is,ij));
        hJtTrueConstIncl->Write();
        cout << "is: " << is << " ij: " << ij << " hJtMeas integral before substraction: " << hJtMeas[is][ij]->Integral() << endl;
        hJtMeas[is][ij]->Add(hJtMeasBg[is][ij],-1);
        cout << "is: " << is << " ij: " << ij << " hJtMeas integral after substraction: " << hJtMeas[is][ij]->Integral() << endl;
        hJtTrue[is][ij]->Add(hJtTrueBg[is][ij],-1);
        hJtMeasConst[is][ij]->Add(hJtMeasBg[is][ij],-1);
        hJtTrueConst[is][ij]->Add(hJtTrueBg[is][ij],-1);
      }

      outFile->cd();
      if(doWeight){
        writeToFile("Pythia/Measured","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasConst[is][ij]);
        writeToFile("Pythia/Measured","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeas[is][ij]);
        writeToFile("PythiaUnscaled/Measured","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hJtMeasConst_unscaled[is][ij])); 
        writeToFile("PythiaUnscaled/Measured","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hJtMeas_unscaled[is][ij])); 
        writeToFile("Pythia/True","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),hJtTrueConst[is][ij]);
        writeToFile("Pythia/True","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hJtTrue[is][ij]);
        writeToFile("Pythia/Sub","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasConstSub[is][ij]);
        writeToFile("Pythia/Sub","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hJtMeasSub[is][ij]);
      }else{
        writeToFile("Pythia/Measured","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hJtMeasConst[is][ij]);
        writeToFile("Pythia/Measured","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hJtMeas[is][ij]);
        writeToFile("PythiaUnscaled/Measured","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hJtMeasConst[is][ij]); 
        writeToFile("PythiaUnscaled/Measured","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hJtMeas[is][ij]); 
        writeToFile("Pythia/True","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hJtTrueConst[is][ij]);
        writeToFile("Pythia/True","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hJtTrue[is][ij]);
        writeToFile("Pythia/Sub","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hJtMeasConstSub[is][ij]);
        writeToFile("Pythia/Sub","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hJtMeasSub[is][ij]);
      }
      if(hJtUnfBg[is][ij]){
        writeToFile("Pythia/Measured","jTUnfBg",Form("jTUnfBgNFin%02d,JetPt%02d",is,ij),hJtUnfBg[is][ij]);
        writeToFile("PythiaUnscaled/Measured","jTUnfBg",Form("jTUnfBgNFin%02d,JetPt%02d",is,ij),addWeight(hJtUnfBg_unscaled[is][ij])); //FIXME
      }
    }
  }
}

void UnfoldJetConst(){
  cout << "==================================== JET CONST ================================" << endl;
  for(int is = 0; is < Nsets; is++){
    for(int ij = 0; ij < Njet; ij++){
      cout << "is: " << is << " ij: " << ij << endl;
      name = Form("AliJJetJtTask/AliJJetJtMCHistManager/ConstMatchSuccess/ConstMatchSuccessNFin%02dJetPt%02d",is,ij);
      TH1D *ConstMatch = (TH1D*)inFile->Get(name);
      name = Form("AliJJetJtTask/AliJJetJtMCHistManager/ConstJtCorrBin/ConstJtCorrBinNFin%02dJetPt%02d",is,ij);
      responseMatrixConst[is][ij] = (TH2D*)inFile->Get(name);
      responseMatrixConst[is][ij]->Print();
      responseMatrixConst[is][ij]->Rebin(rebin);
      if(scaleResponse) responseMatrixConst[is][ij]->Scale(1.0,"width");

      cout << "==================================== RESPONSE ================================" << endl;

      RooUnfoldResponse *response  = CreateResponse(hJtMeasConst[is][ij],responseMatrixConst[is][ij]);
      RooUnfoldBayes *unfoldBayesJetConst = new RooUnfoldBayes(response, hJtMeasConst[is][ij], 4,false); 
      RooUnfoldBayes *unfoldBayesJetConstSub = new RooUnfoldBayes(response, hJtMeasConstSub[is][ij], 4,false); 
      RooUnfoldSvd *unfoldSVDJetConst = new RooUnfoldSvd(response, hJtMeasConst[is][ij], 20);  
      RooUnfoldSvd *unfoldSVDJetConstSub = new RooUnfoldSvd(response, hJtMeasConstSub[is][ij], 20);  

      RooUnfoldBayes *unfoldBayesJetConst_unscaled = new RooUnfoldBayes(response, hJtMeasConst_unscaled[is][ij], 4,false); 
      RooUnfoldBayes *unfoldBayesJetConstSub_unscaled = new RooUnfoldBayes(response, hJtMeasConstSub_unscaled[is][ij], 4,false); 
      RooUnfoldSvd *unfoldSVDJetConst_unscaled = new RooUnfoldSvd(response, hJtMeasConst_unscaled[is][ij], 20);  
      RooUnfoldSvd *unfoldSVDJetConstSub_unscaled = new RooUnfoldSvd(response, hJtMeasConstSub_unscaled[is][ij], 20);  

      cout << "is: " << is  << " ij: " << ij << endl;
      hRecoBayesConst[is][ij] = (TH1D*) unfoldBayesJetConst->Hreco();
      hRecoBayesConstSub[is][ij] = (TH1D*) unfoldBayesJetConstSub->Hreco();
      hRecoSVDConst[is][ij] = (TH1D*) unfoldSVDJetConst->Hreco();
      hRecoSVDConstSub[is][ij] = (TH1D*) unfoldSVDJetConstSub->Hreco();

      hRecoBayesConst_unscaled[is][ij] = (TH1D*) unfoldBayesJetConst_unscaled->Hreco();
      hRecoBayesConstSub_unscaled[is][ij] = (TH1D*) unfoldBayesJetConstSub_unscaled->Hreco();
      hRecoSVDConst_unscaled[is][ij] = (TH1D*) unfoldSVDJetConst_unscaled->Hreco();
      hRecoSVDConstSub_unscaled[is][ij] = (TH1D*) unfoldSVDJetConstSub_unscaled->Hreco();

      cout << "Measured integral: " << hJtMeasConst[is][ij]->Integral() << " True integral: " << hJtTrueConst[is][ij]->Integral() << " Unfolded integral: " << hRecoBayes[is][ij]->Integral() << endl;


      responseMatrixConst[is][ij]->SetName(Form("responseMatrixConstNFin%02dJetPt%02d",is,ij));
      responseMatrixConst[is][ij]->Print();
      outFile->cd();
      writeToFile("Pythia","ConstMatch",ConstMatch);
      writeToFile("Pythia","ResponseMatrixConst",responseMatrixConst[is][ij]);

      RooUnfoldBayes *unfoldBayesData = new RooUnfoldBayes(response,hJtDataConst[is][ij],4,false);
      RooUnfoldSvd *unfoldSVDData = new RooUnfoldSvd(response, hJtDataConst[is][ij],20);
      hRecoBayesDataConst[is][ij] = (TH1D*)unfoldBayesData->Hreco();
      hRecoSVDDataConst[is][ij] = (TH1D*)unfoldSVDData->Hreco();

      RooUnfoldBayes *unfoldBayesDataSub = new RooUnfoldBayes(response,hJtDataConstSub[is][ij],4,false);
      RooUnfoldSvd *unfoldSVDDataSub = new RooUnfoldSvd(response,hJtDataConstSub[is][ij],20);
      hRecoBayesDataConstSub[is][ij] = (TH1D*)unfoldBayesDataSub->Hreco();
      hRecoSVDDataConstSub[is][ij] = (TH1D*)unfoldSVDDataSub->Hreco();


      TH1D *hFolded = foldDist(hJtTrueConst[is][ij],responseMatrixConst[is][ij],1);
      outFile->cd();
      if(doWeight){
        TString sname = Form("JtWeightBinNFin%02dJetPt%02d",is,ij);
        hRecoBayesData[is][ij]->Print();
        writeToFile("Data/BayesUnfolding","JtWeightBin" ,sname, hRecoBayesDataConst[is][ij]);
        writeToFile("Data/SVDUnfolding","JtWeightBin" ,sname, hRecoSVDDataConst[is][ij]);
        writeToFile("Data/BayesSubUnfolding","JtWeightBin" ,sname, hRecoBayesDataConstSub[is][ij]);
        writeToFile("Data/SVDSubUnfolding","JtWeightBin" ,sname, hRecoSVDDataConstSub[is][ij]);
        writeToFile("Pythia/BayesUnfolding","JtWeightBin",sname,hRecoBayesConst[is][ij]);
        writeToFile("Pythia/SVDUnfolding","JtWeightBin",sname,hRecoSVDConst[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding","JtWeightBin",sname,hRecoBayesConstSub[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding","JtWeightBin",sname,hRecoSVDConstSub[is][ij]);
        writeToFile("Pythia/Folded","JtWeightBin",sname,hFolded);
        /*writeToFile("Data/BayesSubUnfolding/WeightCorrected","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),removeWeight(hRecoBayesDataConstSub[is][ij]));
        writeToFile("Data/SVDSubUnfolding/WeightCorrected","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),removeWeight(hRecoSVDDataConstSub[is][ij]));
        writeToFile("Pythia/BayesSubUnfolding/WeightCorrected","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),removeWeight(hRecoBayesConstSub[is][ij]));
        writeToFile("Pythia/SVDSubUnfolding/WeightCorrected","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),removeWeight(hRecoSVDConstSub[is][ij]));*/
      }else{
        TString sname = Form("JtBinNFin%02dJetPt%02d",is,ij);
        hRecoBayesData[is][ij]->Print();
        writeToFile("Data/BayesUnfolding","JtBin" ,sname, hRecoBayesDataConst[is][ij]);
        writeToFile("Data/SVDUnfolding","JtBin" ,sname, hRecoSVDDataConst[is][ij]);
        writeToFile("Data/BayesSubUnfolding","JtBin" ,sname, hRecoBayesDataConstSub[is][ij]);
        writeToFile("Data/SVDSubUnfolding","JtBin" ,sname, hRecoSVDDataConstSub[is][ij]);
        writeToFile("Pythia/BayesUnfolding","JtBin",sname,hRecoBayesConst[is][ij]);
        writeToFile("Pythia/SVDUnfolding","JtBin",sname,hRecoSVDConst[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding","JtBin",sname,hRecoBayesConstSub[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding","JtBin",sname,hRecoSVDConstSub[is][ij]);
        writeToFile("Pythia/Folded","JtBin",sname,hFolded);
        /*writeToFile("Data/BayesSubUnfolding/WeightCorrected","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoBayesDataConstSub[is][ij]));
        writeToFile("Data/SVDSubUnfolding/WeightCorrected","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoSVDDataConstSub[is][ij]));
        writeToFile("Pythia/BayesSubUnfolding/WeightCorrected","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoBayesConstSub[is][ij]));
        writeToFile("Pythia/SVDSubUnfolding/WeightCorrected","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoSVDConstSub[is][ij]));*/
      }
      writeToFile("PythiaUnscaled/BayesSubUnfolding","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoBayesConstSub_unscaled[is][ij]));
      writeToFile("PythiaUnscaled/SVDSubUnfolding","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoSVDConstSub_unscaled[is][ij]));
      writeToFile("PythiaUnscaled/BayesSubUnfolding","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoBayesConst_unscaled[is][ij]));
      writeToFile("PythiaUnscaled/SVDSubUnfolding","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoSVDConst_unscaled[is][ij]));
    }
  }
}


void UnfoldPythia(){
  cout << "==================================== PYTHIA ================================" << endl;
  for(int ij = 0; ij < Njet; ij++){
    cout << "ij: " << ij << endl;

    name = Form("AliJJetJtTask/AliJJetJtMCHistManager/TrackJtCorrBinPythia/TrackJtCorrBinPythiaJetPt%02d",ij);
    responseMatrixPythia[ij] = (TH2D*)inFile2->Get(name);
    responseMatrixPythia[ij]->Print();
    responseMatrixPythia[ij]->Rebin(rebin);
    if(scaleResponse) responseMatrixPythia[ij]->Scale(1.0,"width");
    responseMatrixPythia[ij]->SetName(Form("responseMatrixPythiaJetPt%02d",ij));
    writeToFile("Pythia/Pythia","response",responseMatrixPythia[ij]);


    for(int is = 0; is < Nsets; is++){

      /*name = Form("AliJJetJtTask/AliJJetJtMCHistManager/JetConeJtWeightBinUnfBg/JetConeJtWeightBinUnfBgNFin%02dJetPt%02d",is,ij);
        hJtUnfBg[is][ij]= (TH1D*)inFile->Get(name);
        if(hJtUnfBg[is][ij]){
        hJtUnfBg[is][ij]->Print();
        }else{
        cout << name << " Not found " << endl;
        }

        hJtMeasSub[is][ij] = (TH1D*)hJtMeas[is][ij]->Clone();

        if(hJtUnfBg[is][ij]){
        hJtMeasSub[is][ij]->Add(hJtUnfBg[is][ij],-1);
        }*/ //TODO Uncomment when UnfBg done

      cout << "==================================== RESPONSE ================================" << endl;

      RooUnfoldResponse *response  = CreateResponse(hJtMeas[is][ij],responseMatrixPythia[ij]);
      RooUnfoldBayes *unfoldBayesPythia = new RooUnfoldBayes(response, hJtMeas[is][ij], 4,false); 
      RooUnfoldBayes *unfoldBayesPythiaSub = new RooUnfoldBayes(response, hJtMeasSub[is][ij], 4,false); 
      RooUnfoldSvd *unfoldSVDPythia = new RooUnfoldSvd(response, hJtMeas[is][ij], 20);  
      //RooUnfoldSvd *unfoldSVDJetConstSub = new RooUnfoldSvd(response, hJtMeasSub[is][ij], 20);  

      cout << "is: " << is  << " ij: " << ij << endl;
      hRecoBayes[is][ij] = (TH1D*) unfoldBayesPythia->Hreco();
      hRecoBayesSub[is][ij] = (TH1D*) unfoldBayesPythiaSub->Hreco();
      hRecoSVD[is][ij] = (TH1D*) unfoldSVDPythia->Hreco();
      //hRecoSVDSub[is][ij] = (TH1D*) unfoldSVDJetConeSub->Hreco();
      cout << "Measured integral: " << hJtMeas[is][ij]->Integral() << " True integral: " << hJtPythia[ij]->Integral() << " Unfolded integral: " << hRecoBayes[is][ij]->Integral() << endl;


      hRecoBayes[is][ij]->SetName(Form("UnfoldedBayesPythiaNFin%02dJetPt%02d",is,ij));
      hRecoBayes[is][ij]->Print();
      hRecoBayesSub[is][ij]->SetName(Form("UnfoldedBayesPythiaSubNFin%02dJetPt%02d",is,ij));
      hRecoBayesSub[is][ij]->Print();
      hRecoSVD[is][ij]->SetName(Form("UnfoldedSVDPythiaNFin%02dJetPt%02d",is,ij));
      hRecoSVD[is][ij]->Print();
      //hRecoBayesSub[is][ij]->SetName(Form("UnfoldedBayesSubNFin%02dJetPt%02d",is,ij));
      //hRecoBayesSub[is][ij]->Print();
      //hRecoSVDSub[is][ij]->SetName(Form("UnfoldedSVDSubNFin%02dJetPt%02d",is,ij));
      //hRecoSVDSub[is][ij]->Print();
      //hRecoTUnfold[is]->SetName(Form("UnfoldedTUnfoldNFin%02dJetPt%02d",is,ij));
      //hJtMeas[is][ij]->SetName(Form("MeasuredConstNFin%02dJetPt%02d",is,ij));
      //hJtMeas[is][ij]->Print();
      //hJtMeasSub[is][ij]->SetName(Form("MeasuredConstSubNFin%02dJetPt%02d",is,ij));
      //hJtMeasSub[is][ij]->Print();
      outFile->cd();
      if(doWeight){
        TSTring sname = Form("JtWeightBinNFin%02dJetPt%02d",is,ij);
        writeToFile("Pythia/Pythia","JtWeightBin",sname,hRecoBayes[is][ij]);
        writeToFile("Pythia/Pythia","JtWeightBin",sname,hRecoBayesSub[is][ij]);
        writeToFile("Pythia/Pythia","JtWeightBin",sname,hRecoSVD[is][ij]);
      }else{
        TSTring sname = Form("JtBinNFin%02dJetPt%02d",is,ij);
        writeToFile("Pythia/Pythia","JtBin",sname,hRecoBayes[is][ij]);
        writeToFile("Pythia/Pythia","JtBin",sname,hRecoBayesSub[is][ij]);
        writeToFile("Pythia/Pythia","JtBin",sname,hRecoSVD[is][ij]);
      }
    }
  }
}

void UnfoldJetPt(){
  for(int is = 0; is < Nsets; is++){
    for(int ij = 0; ij < Njet; ij++){
      cout << "is: " << is << " ij: " << ij << endl;
      name = Form("%s/JetPtBin/JetPtBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      cout << name << endl;
      hJetPtBinMeas[is][ij]= (TH1D*)inFile->Get(name);
      hJetPtBinMeas[is][ij]->Print();
      if(!hJetPtBinMeas[is][ij]) continue;

      name = Form("%s/JetPtBin/JetPtBinNFin%02dJetPt%02d",dirname.Data(),is+2*nR,ij);
      hJetPtBinTrue[is][ij]= (TH1D*)inFile->Get(name);
      hJetPtBinTrue[is][ij]->Print();

      double N_jetTrue = hJetPtBinTrue[is][ij]->GetEntries();
      double N_jetMeas = hJetPtBinMeas[is][ij]->GetEntries();
      cout << "True jets: " << N_jetTrue << endl;
      cout << "Measured jets: " << N_jetMeas << endl;
      outFile->cd();
      writeToFile("Pythia/True","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinTrue[is][ij]);
      writeToFile("Pythia/Measured","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinMeas[is][ij]);
      writeToFile("Pythia/BayesUnfolding","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinMeas[is][ij]);
      writeToFile("Pythia/BayesSubUnfolding","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinMeas[is][ij]);
      writeToFile("Pythia/SVDUnfolding","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinMeas[is][ij]);
      writeToFile("Pythia/SVDSubUnfolding","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinMeas[is][ij]);
      writeToFile("Pythia/BayesSubUnfolding/WeightCorrected","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinMeas[is][ij]);
      writeToFile("Pythia/SVDSubUnfolding/WeightCorrected","JetPtBin",Form("JetPtBinNFin%02dJetPt%02d",is,ij),hJetPtBinMeas[is][ij]);
      if(is == 0){
        name = Form("AliJJetJtTask/AliJJetJtMCHistManager/JetPtBinPythia/JetPtBinPythiaJetPt%02d",ij);
        hJetPtBinPythia[ij] = (TH1D*)inFile->Get(name);
      }
    }
    name = Form("%s/JetPt/JetPtNFin%02d",dirname.Data(),is);
    hJetPtMeas[is]= (TH1D*)inFile2->Get(name);
    if(!hJetPtMeas[is]) continue;

    name = Form("%s/JetPt/JetPtNFin%02d",dirname.Data(),is+2*nR);
    hJetPtTrue[is]= (TH1D*)inFile2->Get(name);

    name = Form("AliJJetJtTask/AliJJetJtMCHistManager/JetPtCorr/JetPtCorrNFin%02d",is);
    cout << "is: " << is << endl;
    cout << name << endl;
    responseMatrixJetPt[is] = (TH2D*)inFile->Get(name);
    if(scaleResponse) responseMatrixJetPt[is]->Scale(1.0,"width");
    writeToFile("Pythia","JetPtCorr",Form("JetPtCorrNFin%02d",is),responseMatrixJetPt[is]);


    //RooUnfoldResponse *response  = CreateResponse(hJetPtMeas[is],responseMatrixJetPt[is]);
    RooUnfoldResponse *response  = CreateResponseInverse(hJetPtMeas[is],responseMatrixJetPt[is]);

    RooUnfoldBayes *unfoldBayesJetPt = new RooUnfoldBayes(response, hJetPtMeas[is], 4,false); 
    RooUnfoldSvd   *unfoldSVDJetPt = new RooUnfoldSvd(response, hJetPtMeas[is],20);

    outFile->cd();
    hJetPtReco[is] = (TH1D*) unfoldBayesJetPt.Hreco();
    writeToFile("Pythia/BayesUnfolding","JetPt",Form("JetPtNFin%02d",is),hJetPtReco[is]);
    writeToFile("Pythia/BayesSubUnfolding","JetPt",Form("JetPtNFin%02d",is),hJetPtReco[is]);
    writeToFile("Pythia/BayesSubUnfolding/WeightCorrected","JetPt",Form("JetPtNFin%02d",is),hJetPtReco[is]);

    hJetPtRecoSVD[is] = (TH1D*) unfoldSVDJetPt.Hreco();
    writeToFile("Pythia/SVDUnfolding","JetPt",Form("JetPtNFin%02d",is),hJetPtRecoSVD[is]);
    writeToFile("Pythia/SVDSubUnfolding","JetPt",Form("JetPtNFin%02d",is),hJetPtRecoSVD[is]);
    writeToFile("Pythia/SVDSubUnfolding/WeightCorrected","JetPt",Form("JetPtNFin%02d",is),hJetPtRecoSVD[is]);

    writeToFile("Pythia/Measured","JetPt",Form("JetPtNFin%02d",is),hJetPtMeas[is]);
    writeToFile("Pythia/True","JetPt",Form("JetPtNFin%02d",is),hJetPtTrue[is]);
    TH1D *hJetPtFolded = foldDist(hJetPtTrue[is],responseMatrixJetPt[is],0);
    writeToFile("Pythia/Folded","JetPt",Form("JetPtNFin%02d",is),hJetPtFolded);
  }
}

void UnfoldJetCone(){
  cout << "==================================== JET CONE ================================" << endl;
  for(int is = 0; is < Nsets; is++){
    for(int ij = 0; ij < Njet; ij++){


      name = Form("AliJJetJtTask/AliJJetJtMCHistManager/TrackMatchSuccess/TrackMatchSuccessNFin%02dJetPt%02d",is,ij);
      TH1D *TrackMatch = (TH1D*)inFile->Get(name);

      name = Form("AliJJetJtTask/AliJJetJtMCHistManager/TrackJtCorrBin/TrackJtCorrBinNFin%02dJetPt%02d",is,ij);
      responseMatrix[is][ij] = (TH2D*)inFile->Get(name);
      if(scaleResponse) responseMatrix[is][ij]->Scale(1.0,"width");
      responseMatrix[is][ij]->Rebin(rebin);

      cout << "==================================== RESPONSE ================================" << endl;

      RooUnfoldResponse *response  = CreateResponse(hJtMeas[is][ij],responseMatrix[is][ij]);
      RooUnfoldBayes *unfoldBayesJetCone = new RooUnfoldBayes(response, hJtMeas[is][ij], 4,false); 
      RooUnfoldBayes *unfoldBayesJetConeSub = new RooUnfoldBayes(response, hJtMeasSub[is][ij], 4,false); 
      RooUnfoldBayes *unfoldBayesJetConeSub2D = new RooUnfoldBayes(response, hRecoBayes1D[is][ij], 4,false); 
      RooUnfoldSvd *unfoldSVDJetCone = new RooUnfoldSvd(response, hJtMeas[is][ij], 20);  
      RooUnfoldSvd *unfoldSVDJetConeSub = new RooUnfoldSvd(response, hJtMeasSub[is][ij], 20);  
      RooUnfoldSvd *unfoldSVDJetConeSub2D = new RooUnfoldSvd(response, hRecoBayes1D[is][ij], 20);  
      cout << "==================================== RESPONSE CREATED ========================" << endl;

      RooUnfoldBayes *unfoldBayesJetCone_unscaled = new RooUnfoldBayes(response, hJtMeas_unscaled[is][ij], 4,false); 
      RooUnfoldBayes *unfoldBayesJetConeSub_unscaled = new RooUnfoldBayes(response, hJtMeasSub_unscaled[is][ij], 4,false); 
      RooUnfoldSvd *unfoldSVDJetCone_unscaled = new RooUnfoldSvd(response, hJtMeas_unscaled[is][ij], 20);  
      RooUnfoldSvd *unfoldSVDJetConeSub_unscaled = new RooUnfoldSvd(response, hJtMeasSub_unscaled[is][ij], 20);  

      outFile->cd();
      hRecoBayes[is][ij] = (TH1D*) unfoldBayesJetCone->Hreco();
      hRecoBayesSub[is][ij] = (TH1D*) unfoldBayesJetConeSub->Hreco();
      hRecoSVD[is][ij] = (TH1D*) unfoldSVDJetCone->Hreco();
      hRecoSVDSub[is][ij] = (TH1D*) unfoldSVDJetConeSub->Hreco();
      if(do2D){
        hRecoBayesSub2D[is][ij] = (TH1D*)unfoldBayesJetConeSub2D->Hreco();
        hRecoSVDSub2D[is][ij] = (TH1D*)unfoldSVDJetConeSub2D->Hreco();
      }


      hRecoBayes_unscaled[is][ij] = (TH1D*) unfoldBayesJetCone_unscaled->Hreco();
      hRecoBayesSub_unscaled[is][ij] = (TH1D*) unfoldBayesJetConeSub_unscaled->Hreco();
      hRecoSVD_unscaled[is][ij] = (TH1D*) unfoldSVDJetCone_unscaled->Hreco();
      hRecoSVDSub_unscaled[is][ij] = (TH1D*) unfoldSVDJetConeSub_unscaled->Hreco();

      cout << "Measured integral: " << hJtMeas[is][ij]->Integral() << " True integral: " << hJtTrue[is][ij]->Integral() << " Unfolded integral: " << hRecoBayes[is][ij]->Integral() << endl;

      responseMatrix[is][ij]->SetName(Form("responseMatrixNFin%02dJetPt%02d",is,ij));
      writeToFile("Pythia","TrackMatch",TrackMatch);

      RooUnfoldBayes *unfoldBayesData = new RooUnfoldBayes(response,hJtData[is][ij],4,false);
      RooUnfoldSvd *unfoldSVDData = new RooUnfoldSvd(response, hJtData[is][ij],20);
      hRecoBayesData[is][ij] = (TH1D*)unfoldBayesData->Hreco();
      hRecoSVDData[is][ij] = (TH1D*)unfoldSVDData->Hreco();

      RooUnfoldBayes *unfoldBayesDataSub = new RooUnfoldBayes(response,hJtDataSub[is][ij],4,false);
      RooUnfoldSvd *unfoldSVDDataSub = new RooUnfoldSvd(response,hJtDataSub[is][ij],20);
      hRecoBayesDataSub[is][ij] = (TH1D*)unfoldBayesDataSub->Hreco();
      hRecoSVDDataSub[is][ij] = (TH1D*)unfoldSVDDataSub->Hreco();

      writeToFile("Pythia","responseMatrix",responseMatrix[is][ij]);
      if(hJtUnfBg[is][ij]){
        writeToFile("Pythia/Measured","jTUnfBgNFin",Form("jTUnfBgNFin%02dJetPt%02d",is,ij),hJtUnfBg[is][ij]);
      }

      TH1D *hFolded = foldDist(hJtTrue[is][ij],responseMatrix[is][ij],1);
      outFile->cd();
      if(doWeight){
        TString sname = Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij);
        hRecoBayesData[is][ij]->Print();
        writeToFile("Data/BayesUnfolding","JetConeJtWeightBin" ,sname, hRecoBayesData[is][ij]);
        writeToFile("Data/SVDUnfolding","JetConeJtWeightBin" ,sname, hRecoSVDData[is][ij]);
        writeToFile("Data/BayesSubUnfolding","JetConeJtWeightBin" ,sname, hRecoBayesDataSub[is][ij]);
        writeToFile("Data/SVDSubUnfolding","JetConeJtWeightBin" ,sname, hRecoSVDDataSub[is][ij]);
        writeToFile("Pythia/BayesUnfolding","JetConeJtWeightBin",sname,hRecoBayes[is][ij]);
        writeToFile("Pythia/SVDUnfolding","JetConeJtWeightBin",sname,hRecoSVD[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding","JetConeJtWeightBin",sname,hRecoBayesSub[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding","JetConeJtWeightBin",sname,hRecoSVDSub[is][ij]);
        writeToFile("Pythia/Folded","JetConeJtWeightBin",sname,hFolded);
        writeToFile("Data/BayesSubUnfolding/WeightCorrected","JetConeJtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),removeWeight(hRecoBayesDataSub[is][ij]));
        writeToFile("Data/SVDSubUnfolding/WeightCorrected","JetConeJtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),removeWeight(hRecoSVDDataSub[is][ij]));
        writeToFile("Pythia/BayesSubUnfolding/WeightCorrected","JetConeJtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),removeWeight(hRecoBayesSub[is][ij]));
        writeToFile("Pythia/SVDSubUnfolding/WeightCorrected","JetConeJtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),removeWeight(hRecoSVDSub[is][ij]));
      }else{
        TString sname = Form("JetConeJtBinNFin%02dJetPt%02d",is,ij);
        hRecoBayesData[is][ij]->Print();
        writeToFile("Data/BayesUnfolding","JetConeJtBin" ,sname, hRecoBayesData[is][ij]);
        writeToFile("Data/SVDUnfolding","JetConeJtBin" ,sname, hRecoSVDData[is][ij]);
        writeToFile("Data/BayesSubUnfolding","JetConeJtBin" ,sname, hRecoBayesDataSub[is][ij]);
        writeToFile("Data/SVDSubUnfolding","JetConeJtBin" ,sname, hRecoSVDDataSub[is][ij]);
        writeToFile("Pythia/BayesUnfolding","JetConeJtBin",sname,hRecoBayes[is][ij]);
        writeToFile("Pythia/SVDUnfolding","JetConeJtBin",sname,hRecoSVD[is][ij]);
        writeToFile("Pythia/BayesSubUnfolding","JetConeJtBin",sname,hRecoBayesSub[is][ij]);
        writeToFile("Pythia/SVDSubUnfolding","JetConeJtBin",sname,hRecoSVDSub[is][ij]);
        writeToFile("Pythia/Folded","JetConeJtBin",sname,hFolded);
        writeToFile("Data/BayesSubUnfolding/WeightCorrected","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoBayesDataSub[is][ij]));
        writeToFile("Data/SVDSubUnfolding/WeightCorrected","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoSVDDataSub[is][ij]));
        writeToFile("Pythia/BayesSubUnfolding/WeightCorrected","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoBayesSub[is][ij]));
        writeToFile("Pythia/SVDSubUnfolding/WeightCorrected","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoSVDSub[is][ij]));
      }

      writeToFile("PythiaUnscaled/BayesSubUnfolding","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoBayesSub_unscaled[is][ij]));
      writeToFile("PythiaUnscaled/SVDSubUnfolding","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoSVDSub_unscaled[is][ij]));
      writeToFile("PythiaUnscaled/BayesSubUnfolding","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoBayes_unscaled[is][ij]));
      writeToFile("PythiaUnscaled/SVDSubUnfolding","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),addWeight(hRecoSVD_unscaled[is][ij]));
    }
  }
}

void Unfold2D(){
  cout << "==================================== JET CONE 2D ================================" << endl;
  for(int is = 0; is < Nsets; is++){
    for(int ij = 0; ij < Njet; ij++){
      cout << "is: " << is << " ij: " << ij << endl;
      name = Form("%s/JetPtBin/JetPtBinNFin%02dJetPt%02d",dirname.Data(),is,ij);
      cout << name << endl;
      hJetPtBinMeas[is][ij]= (TH1D*)inFile->Get(name);
      cout << "hJetPtBinMeas: ";
      hJetPtBinMeas[is][ij]->Print();
      if(!hJetPtBinMeas[is][ij]) continue;

      name = Form("%s/JetPtBin/JetPtBinNFin%02dJetPt%02d",dirname.Data(),is+2*nR,ij);
      cout << name << endl;
      hJetPtBinTrue[is][ij]= (TH1D*)inFile->Get(name);
      cout << "hJetPtBinTrue: ";
      hJetPtBinTrue[is][ij]->Print();


      //double N_jetTrue = hJetPtBinTrue[is][ij]->GetEntries();
      //double N_jetMeas = hJetPtBinMeas[is][ij]->GetEntries();
      double N_jetTrue = hJetPtBinTrue[is][ij]->Integral();
      double N_jetMeas = hJetPtBinMeas[is][ij]->Integral();
      cout << "True jets: " << N_jetTrue << endl;
      cout << "Measured jets: " << N_jetMeas << endl;
      outFile->cd();

    }

    name = Form("%s/JetPt/JetPtNFin%02d",dirname.Data(),is);
    hJetPtMeas[is]= (TH1D*)inFile2->Get(name);
    if(!hJetPtMeas[is]) continue;

    name = Form("%s/JetPt/JetPtNFin%02d",dirname.Data(),is+2*nR);
    hJetPtTrue[is]= (TH1D*)inFile2->Get(name);
    outFile->cd();
    hJetPtMeas[is]->Write();
    hJetPtTrue[is]->Write();

    name = Form("AliJJetJtTask/AliJJetJtMCHistManager/JetPtCorrCoarse/JetPtCorrCoarseNFin%02d",is);
    responseMatrixJetPtCoarse[is] = (TH2D*)inFile->Get(name);

    cout << "==================================== RESPONSE ================================" << endl;
    responseCoarse = (TH2D*)responseMatrixJetPtCoarse[is]->Clone();
    outFile->cd();
    //responseCoarse->Draw();
    double ptTrue, ptObs;
    double NTrue, NObs, Ntot;
    double binW;
    int ib,ib2;
    TH1D *projX = (TH1D*)responseCoarse->ProjectionX();
    double eff[10];
    double corr[10];
    TH1D *hJetPtUnfolded = (TH1D*)projX->Clone();
    TH1D *hJetPtMeasCoarse = (TH1D*)projX->Clone();
    TH1D *hJetPtTrueCoarse = (TH1D*)projX->Clone();
    hJetPtMeasCoarse->Reset();
    hJetPtTrueCoarse->Reset();
    hJetPtUnfolded->Reset();
    for(int ij = 0; ij < Njet; ij++){
      //hJetPtMeasCoarse->SetBinContent(ij+1,hJetPtBinMeas[is][ij]->GetEntries());
      //hJetPtTrueCoarse->SetBinContent(ij+1,hJetPtBinTrue[is][ij]->GetEntries());
      hJetPtMeasCoarse->SetBinContent(ij+1,hJetPtBinMeas[is][ij]->Integral());
      hJetPtTrueCoarse->SetBinContent(ij+1,hJetPtBinTrue[is][ij]->Integral());
    }
    outFile->cd();
    hJetPtMeasCoarse->SetName(Form("hJetPtMeasCoarse%02d",is));
    hJetPtTrueCoarse->SetName(Form("hJetPtTrueCoarse%02d",is));
    writeToFile("Pythia", "hJetPtMeasCoarse", hJetPtMeasCoarse);
    writeToFile("Pythia", "hJetPtTrueCoarse", hJetPtTrueCoarse);
    double factor;
    double sumP;
    writeToFile("Pythia","responseCoarse",responseCoarse);

    //First create P(E|C) i.e. the actual normalized response matrix from correlation matrix, also the efficiency, i.e. N(True jets in bin)/N(Observed jets in any bin corresponding to true jets)
    //corr = 1/efficiency
    for(int ibx = 1 ; ibx <= responseCoarse->GetNbinsX() ; ibx++){
      hJetPtUnfolded->SetBinContent(ibx,1); //First create uniform matrix for initial guess
      ptTrue = responseCoarse->GetXaxis()->GetBinCenter(ibx); //True pT is on X axis
      binW = responseCoarse->GetXaxis()->GetBinWidth(ibx);
      cout << "True pT: " << ptTrue - binW/2 <<"-" << ptTrue + binW/2;
      ib = responseCoarse->GetBin(ibx,0);
      Ntrue  = projX->GetBinContent(ibx);
      cout << " True N: " << Ntrue << endl;
      cout << " Missed N: " << responseCoarse->GetBinContent(ib) << endl;
      Ntot = 0;
      if(Ntrue > 0){
        cout << " fraction missed: " << responseCoarse->GetBinContent(ib)/Ntrue << endl;
        for(int iby = 1 ; iby <= responseCoarse->GetNbinsY(); iby++){
          ptObs = responseCoarse->GetYaxis()->GetBinCenter(iby);
          binW = responseCoarse->GetYaxis()->GetBinWidth(iby);
          cout << "Observed pT: " << ptObs - binW/2 <<"-" << ptObs + binW/2;
          ib = responseCoarse->GetBin(ibx,iby);
          NObs = responseCoarse->GetBinContent(ib);
          Ntot += NObs;
          cout << ", N: " << NObs <<  ", P(E|C): " << 1.0*NObs/Ntrue <<endl;
          //if(N > 0) cout << "Fill bin (" << jtobs << "," << jt << ") with " << N << endl;
          //response->Fill(jtobs,jt,N);
          responseCoarse->SetBinContent(ib,1.0*NObs/Ntrue);
        }
      }else{
        break;
      }
      cout << "Total N: " << Ntot;
      cout << ", Efficiency: " << Ntot/Ntrue;
      eff[ibx-1] = Ntot/Ntrue;
      corr[ibx-1] = Ntrue/Ntot;
      cout << ", Correction: " << corr[ibx-1] << endl;
    }
    responseCoarse->SetTitle(Form("Probability of effect Y given Cause X is=%02d",is));
    responseCoarse->Write();
    writeToFile("Pythia","responseCoarse",responseCoarse);
    outFile->cd();
    TH2D *inverseP = (TH2D*)responseCoarse->Clone();
    inverseP->Reset();
    //inverseP->SetName(Form("inverseP%02d",is));
    inverseP->SetTitle("Probability of cause X given observed Y");


    //Start the iteration loop
    int ib3;
    double sum = 0;
    for(int it = 0; it < 10 ;it++){
      hJetPtUnfolded->Scale(1.0/hJetPtUnfolded->Integral()); //Scale the unfolded histogram so that it gives the probability
      //First create P(C|E), labeled inverseP
      for(int ij = 0; ij < Njet; ij++){
        for(int ij2 = 0; ij2 < Njet; ij2++){
          ib = inverseP->GetBin(ij+1,ij2+1);  //Get bin x,y
          ib2 = responseCoarse->GetBin(ij+1,ij2+1); //Get bin x,y
          if(ib != ib2) cout << "ib != ib2" << endl << endl;
          sum = 0; //Sum over z P(E_y | C_z)*P( C_z )
          for(int ij3 = 0; ij3 < Njet; ij3++){
            ib3 = responseCoarse->GetBin(ij3+1,ij2+1); //Bin z,y
            sum += responseCoarse->GetBinContent(ib3)*hJetPtUnfolded->GetBinContent(ij3+1); // P(E_y | C_z) * P( C_z)
          }
          cout << "Sum is: " << sum << endl;
          double PEy_Cx = responseCoarse->GetBinContent(ib2); //Get P(E_y | C_x)
          double PCx = hJetPtUnfolded->GetBinContent(ij+1); //Get P( C_x)
          cout << "responseCoarse->GetBinContent(ib2) gives " << PEy_Cx << endl; 
          cout << "hJetPtUnfolded->GetBinContent(ij+1) gives " << PCx << endl;
          factor = PEy_Cx * PCx / sum; //Factor is P(E_y | C_x) * P( C_x) / sum
          cout << "Set P(C" << ij << "| E" << ij2 << ") to " << factor << endl;
          inverseP->SetBinContent(ib,factor); // Set bin content
        }
      } 
      outFile->cd();
      inverseP->SetName(Form("inverseP%02d",it));
      inverseP->Write();
      double sumP = 0;
      for(int ij2 = 0; ij2 < Njet; ij2++){
        sumP = 0;
        for(int ij = 0; ij < Njet; ij++){
          ib = inverseP->GetBin(ij+1,ij2+1);
          sumP +=   inverseP->GetBinContent(ib);
        }
        cout << "ij2: " << ij2 << " sum of P is: " << sumP << endl;
      }

      //Then calculate new distribution from P(C | E)
      hJetPtUnfolded->Reset();
      for(int ij = 0; ij < Njet ; ij++){
        cout << "Unfolding ij: " << ij << endl;
        //hRecoBayesSub2D[is][ij] = (TH1D*)hRecoBayesSub[is][ij]->Clone();
        //hRecoBayesSub2D[is][ij]->Reset();
        //hRecoBayesSub2D[is][ij]->SetName(Form("hRecoBayesSub2DNFin%02dJetPt%02dIt%02d",is,ij,it));

        hRecoBayes1D[is][ij] = (TH1D*)hJtMeas[is][ij]->Clone(); //jT unfolding
        hRecoBayes1D[is][ij]->Reset(); //jT unfolding
        hRecoBayes1D[is][ij]->SetName(Form("hRecoBayes1DNFin%02dJetPt%02dIt%02d",is,ij,it)); //jT unfolding
        hRecoBayes1D[is][ij]->Print(); //jT unfolding
        hRecoBayesConst1D[is][ij] = (TH1D*)hJtMeasConst[is][ij]->Clone();
        hRecoBayesConst1D[is][ij]->Reset();
        hRecoBayesConst1D[is][ij]->SetName(Form("hRecoBayesConst1DNFin%02dJetPt%02dIt%02d",is,ij,it));
        hRecoBayesConst1D[is][ij]->Print();
        for(int ij2 = 0 ; ij2 < Njet; ij2++){
          ib =  inverseP->GetBin(1+ij,1+ij2);  //Get bin x,y
          hRecoBayes1D[is][ij]->Add(hJtMeasSub[is][ij2],corr[ij]*inverseP->GetBinContent(ib)); 
          hRecoBayesConst1D[is][ij]->Add(hJtMeasConstSub[is][ij2],corr[ij]*inverseP->GetBinContent(ib)); 
          //hRecoBayesSub2D[is][ij]->Add(hRecoBayesSub[is][ij2],corr[ij]*inverseP->GetBinContent(ib));
          cout << "Add content to hJetPtUnfolded: " << corr[ij] << " * " << inverseP->GetBinContent(ib) << " * " << hJetPtBinMeas[is][ij2]->GetEntries() << " = " << corr[ij]*inverseP->GetBinContent(ib)*hJetPtBinMeas[is][ij2]->GetEntries() << endl;
          hJetPtUnfolded->AddBinContent(ij+1,corr[ij]*inverseP->GetBinContent(ib)*hJetPtMeasCoarse->GetBinContent(ij2+1)); //Bin of new histogram is sum over old bins j 1/eff * n(E_j) * P(C_x | E_j)

          //cout << "Histogram ij2: " << ij2 << " was scaled by " << corr[ij] << " * " << inverseP->GetBinContent(ib) << " = " << corr[ij]*inverseP->GetBinContent(ib) << endl;
        }
        cout << "Total content added: " << hJetPtUnfolded->GetBinContent(ij+1) << endl;
        writeToFile("Pythia/JetPtUnfolding","hJetConeJtWeightBin",Form("hJetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hRecoBayes1D[is][ij]);
        writeToFile("Pythia/JetPtUnfolding","hJtWeightBin",Form("hJtWeightBinNFin%02dJetPt%02d",is,ij),hRecoBayesConst1D[is][ij]);
      }
      hJetPtUnfolded->SetName(Form("hJetPtUnfolded2DIt%02d",it));
      writeToFile("Pythia/2DUnfolding","hJetPt",Form("hJetPtIt%02d",it),hJetPtUnfolded);
    }
    cout << "Debug 0" << endl;


    for(int ij = 0 ; ij < Njet; ij++){
      cout << "Debug 0.1" << endl;
      hRecoBayes1DSub[is][ij] = (TH1D*)hRecoBayes1D[is][ij]->Clone();
      hRecoBayes1DSub[is][ij]->Add(hJtUnfBg[is][ij],-1);
      hRecoBayesConst1DSub[is][ij] = (TH1D*)hRecoBayesConst1D[is][ij]->Clone();
      hRecoBayesConst1DSub[is][ij]->Add(hJtUnfBg[is][ij],-1);

      name = Form("AliJJetJtTask/AliJJetJtMCHistManager/TrackJtCorrBin/TrackJtCorrBinNFin%02dJetPt%02d",is,ij);
      responseMatrix[is][ij] = (TH2D*)inFile->Get(name);
      responseMatrix[is][ij]->Rebin(rebin);
      cout << "Debug 0.2" << endl;

      RooUnfoldResponse *response  = CreateResponse(hRecoBayes1D[is][ij],responseMatrix[is][ij]);
      cout << "Debug  0.25" << endl;
      RooUnfoldBayes *unfoldBayesJetCone = new RooUnfoldBayes(response, hRecoBayes1D[is][ij], 4,false); 
      RooUnfoldBayes *unfoldBayesJetConeSub = new RooUnfoldBayes(response, hRecoBayes1DSub[is][ij], 4,false); 
      cout << "Debug  0.3" << endl;
      RooUnfoldSvd *unfoldSVDJetCone = new RooUnfoldSvd(response, hRecoBayes1D[is][ij], 20);  
      RooUnfoldSvd *unfoldSVDJetConeSub = new RooUnfoldSvd(response, hRecoBayes1DSub[is][ij], 20);  

      cout << "Debug " << endl;

      RooUnfoldBayes *unfoldBayes= new RooUnfoldBayes(response, hRecoBayesConst1D[is][ij], 4,false); 
      RooUnfoldBayes *unfoldBayesSub = new RooUnfoldBayes(response, hRecoBayesConst1DSub[is][ij], 4,false); 
      cout << "Debug 1" << endl;
      RooUnfoldSvd *unfoldSVD = new RooUnfoldSvd(response, hRecoBayesConst1D[is][ij], 20);  
      RooUnfoldSvd *unfoldSVDSub = new RooUnfoldSvd(response, hRecoBayesConst1DSub[is][ij], 20);  

      cout << "Debug 2" << endl;


      cout << "is: " << is  << " ij: " << ij << endl;
      hRecoBayes2D[is][ij] = (TH1D*) unfoldBayesJetCone->Hreco();
      hRecoBayes2DSub[is][ij] = (TH1D*) unfoldBayesJetConeSub->Hreco();
      hRecoSVD2D[is][ij] = (TH1D*) unfoldSVDJetCone->Hreco();
      hRecoSVD2DSub[is][ij] = (TH1D*) unfoldSVDJetConeSub->Hreco();

      hRecoBayesConst2D[is][ij] = (TH1D*) unfoldBayes->Hreco();
      hRecoBayesConst2DSub[is][ij] = (TH1D*) unfoldBayesSub->Hreco();
      hRecoSVDConst2D[is][ij] = (TH1D*) unfoldSVD->Hreco();
      hRecoSVDConst2DSub[is][ij] = (TH1D*) unfoldSVDSub->Hreco();

      hRecoBayes2D[is][ij]->SetName(Form("UnfoldedBayes2DNFin%02dJetPt%02d",is,ij));
      if(doWeight){
        writeToFile("Pythia/2DUnfoldingBayes","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hRecoBayes2D[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVD","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hRecoBayes2DSub[is][ij]);
        writeToFile("Pythia/2DUnfoldingBayesSub","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hRecoSVD2D[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVDSub","JetConeJtWeightBin",Form("JetConeJtWeightBinNFin%02dJetPt%02d",is,ij),hRecoSVD2DSub[is][ij]);

        writeToFile("Pythia/2DUnfoldingBayes","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),hRecoBayesConst2D[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVD","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),hRecoBayesConst2DSub[is][ij]);
        writeToFile("Pythia/2DUnfoldingBayesSub","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),hRecoSVDConst2D[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVDSub","JtWeightBin",Form("JtWeightBinNFin%02dJetPt%02d",is,ij),hRecoSVDConst2DSub[is][ij]);
      }else{
        writeToFile("Pythia/2DUnfoldingBayes","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hRecoBayes2D[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVD","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hRecoBayes2DSub[is][ij]);
        writeToFile("Pythia/2DUnfoldingBayesSub","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hRecoSVD2D[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVDSub","JetConeJtBin",Form("JetConeJtBinNFin%02dJetPt%02d",is,ij),hRecoSVD2DSub[is][ij]);

        writeToFile("Pythia/2DUnfoldingBayes","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hRecoBayesConst2D[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVD","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hRecoBayesConst2DSub[is][ij]);
        writeToFile("Pythia/2DUnfoldingBayesSub","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hRecoSVDConst2D[is][ij]);
        writeToFile("Pythia/2DUnfoldingSVDSub","JtBin",Form("JtBinNFin%02dJetPt%02d",is,ij),hRecoSVDConst2DSub[is][ij]);
      }
      /*hRecoBayes2D[is][ij]->Print();
        hRecoSVD2D[is][ij]->SetName(Form("UnfoldedSVD2DNFin%02dJetPt%02d",is,ij));
        hRecoSVD2D[is][ij]->Print();
        hRecoBayes2DSub[is][ij]->SetName(Form("UnfoldedBayes2DSubNFin%02dJetPt%02d",is,ij));
        hRecoBayes2DSub[is][ij]->Print();
        hRecoSVD2DSub[is][ij]->SetName(Form("UnfoldedSVD2DSubNFin%02dJetPt%02d",is,ij));
        hRecoSVD2DSub[is][ij]->Print();
        hRecoBayes2D[is][ij]->Write();
        hRecoBayes2DSub[is][ij]->Write();
        hRecoSVD2D[is][ij]->Write();
        hRecoSVD2DSub[is][ij]->Write();*/
    }
  }
}

RooUnfoldResponse *CreateResponse(TH1D *hMeas, TH2D *responseMatrix){
  cout << "==================================== CREATE RESPONSE===================================" << endl;
  double jtobs = 1;
  double jt= 1;
  double N = 0;
  int ib = 0;
  RooUnfoldResponse *response = new RooUnfoldResponse(hMeas,hMeas);
  for(int iby = 1 ; iby <= responseMatrix->GetNbinsY(); iby++){
    jt = responseMatrix->GetYaxis()->GetBinCenter(iby);
    ib = responseMatrix->GetBin(0,iby);
    N  = responseMatrix->GetBinContent(ib);
    response->Miss(jt,N); //FIXME uncomment
    for(int ibx = 1 ; ibx <= responseMatrix->GetNbinsX() ; ibx++){
      jtobs = responseMatrix->GetXaxis()->GetBinCenter(ibx);
      ib = responseMatrix->GetBin(ibx,iby);
      N = responseMatrix->GetBinContent(ib);
      //if(N > 0) cout << "Fill bin (" << jtobs << "," << jt << ") with " << N << endl;
      response->Fill(jtobs,jt,N);
    }
  }
  return response;
}

void writeToFile(TString folder, TString subFolder, TString name, TH1 *histo){
  //cout << "Write to " << folder.Data() << "/" << subFolder.Data() << " with name " << name << endl;
  histo->SetName(name);
  writeToFile(folder,subFolder,histo);
}
void writeToFile(TString folder, TString subFolder, TH1 *histo){
  //histo->Print();
  outFile->cd();
  TDirectory *newdir2 = outFile->GetDirectory(Form("%s/%s",folder.Data(),subFolder.Data()));
  if(!newdir2){
    TDirectory *newdir = outFile->mkdir(Form("%s/%s",folder.Data(),subFolder.Data()));
    TDirectory *newdir2 = outFile->GetDirectory(Form("%s/%s",folder.Data(),subFolder.Data()));
  }
  newdir2->cd();
  histo->Write();
}

RooUnfoldResponse *CreateResponseInverse(TH1D *hMeas, TH2D *responseMatrix){
  cout << "==================================== CREATE RESPONSE===================================" << endl;
  double jtobs = 1;
  double jt= 1;
  double N = 0;
  int ib = 0;
  RooUnfoldResponse *response = new RooUnfoldResponse(hMeas,hMeas);
  for(int ibx = 1 ; ibx <= responseMatrix->GetNbinsX(); ibx++){
    jt = responseMatrix->GetXaxis()->GetBinCenter(ibx);
    ib = responseMatrix->GetBin(ibx,0);
    N  = responseMatrix->GetBinContent(ib);
    response->Miss(jt,N); //FIXME uncomment
    for(int iby = 1 ; iby <= responseMatrix->GetNbinsY() ; iby++){
      jtobs = responseMatrix->GetYaxis()->GetBinCenter(iby);
      ib = responseMatrix->GetBin(ibx,iby);
      N = responseMatrix->GetBinContent(ib);
      //if(N > 0) cout << "Fill bin (" << jtobs << "," << jt << ") with " << N << endl;
      response->Fill(jtobs,jt,N);
    }
  }
  return response;
}

TH1D *foldDist(TH1D *hTrue, TH2D *responseMatrix, int inverse){
  //cout << "hTrue Integral: " << hTrue->Integral() << endl;
  TH1D *hFolded = (TH1D*)hTrue->Clone();
  hFolded->Reset();
  int ib = 0;
  double jt,jtobs,N;
  double Ntrue;
  double effx[150];
  double corrx[150];
  double NTrue, NObs, Ntot;
  double binW;
  if(inverse){
    TH1D *projX = (TH1D*)responseMatrix->ProjectionY();
  }else{
    TH1D *projX = (TH1D*)responseMatrix->ProjectionX();
  }
  for(int ibx = 1 ; ibx <= responseMatrix->GetNbinsX() ; ibx++){
    if(inverse) {
      ib = responseMatrix->GetBin(0,ibx);
    }
    else{
      ib = responseMatrix->GetBin(ibx,0);
    }
    Ntrue  = projX->GetBinContent(ibx);
    jt = projX->GetBinCenter(ibx);
    //cout << "ibx: " << ibx <<  " bin center: " << jt << " Ntrue: " << Ntrue << endl;
    Ntot = 0;
    if(Ntrue > 0){
      //cout << " fraction missed: " << responseMatrix->GetBinContent(ib)/Ntrue << endl;
      for(int iby = 1 ; iby <= responseMatrix->GetNbinsY(); iby++){
        if(inverse){
          jtobs = responseMatrix->GetXaxis()->GetBinCenter(iby);
          binW = responseMatrix->GetXaxis()->GetBinWidth(iby);
        }else{
          jtobs = responseMatrix->GetYaxis()->GetBinCenter(iby);
          binW = responseMatrix->GetYaxis()->GetBinWidth(iby);
        }
        //cout << "Observed pT: " << jtobs - binW/2 <<"-" << jtobs + binW/2;
        if(inverse){
          ib = responseMatrix->GetBin(iby,ibx);
        }else{
          ib = responseMatrix->GetBin(ibx,iby);
        }
        NObs = responseMatrix->GetBinContent(ib);
        Ntot += NObs;
        //cout << ", N: " << NObs <<  ", P(E|C): " << 1.0*NObs/Ntrue <<endl;
        responseMatrix->SetBinContent(ib,1.0*NObs/Ntrue);
        responseMatrix->SetBinContent(ib,1.0*NObs/Ntrue);
      }
      effx[ibx-1] = Ntot/Ntrue;
      corrx[ibx-1] = Ntrue/Ntot;
    }else{
      effx[ibx-1] = 1;
      corrx[ibx-1] = 1;

    }
    //cout << "Total N: " << Ntot;
    //cout << ", Efficiency: " << Ntot/Ntrue;
    //cout << ", Correction: " << corrx[ibx-1] << endl;
  }
  /*projX = (TH1D*)responseMatrix->ProjectionX();
    for(int ibx = 1; ibx < projX->GetNbinsX(); ibx++){
    cout << "ibx: " << ibx << " projection content: " << projX->GetBinContent(ibx) << endl; 
    }*/
  double P = 0;
  for(int ibx = 1; ibx <= responseMatrix->GetNbinsX(); ibx++){
    if(inverse){
      jt = hFolded->GetYaxis()->GetBinCenter(ibx);
    }else{
      jt = hFolded->GetXaxis()->GetBinCenter(ibx);
    }
    for(int iby = 1; iby <= responseMatrix->GetNbinsY(); iby++){
      if(inverse){
        jtobs = responseMatrix->GetXaxis()->GetBinCenter(iby);
        ib = responseMatrix->GetBin(iby,ibx);
      }else{
        jtobs = responseMatrix->GetYaxis()->GetBinCenter(iby);
        ib = responseMatrix->GetBin(ibx,iby);
      }
      P = responseMatrix->GetBinContent(ib);
      //cout << "hTrue->GetBinContent: " << hTrue->GetBinContent(ibx) << " P: " << P << ", hFolded filled with " << (P*hTrue->GetBinContent(ibx) )  << endl;
      hFolded->Fill(jtobs,hTrue->GetBinContent(ibx)*P);
    }
  }
  //cout << "hFolded Integral: " << hFolded->Integral() << endl;
  return hFolded;
}

TH1D* removeWeight(TH1D* histo){
  TH1D *histo2 = (TH1D*)histo->Clone();
  int N = histo2->GetNbinsX();
  double content = 0;
  double binc = 0;
  for(int i = 1; i <= N ; i++){
    content = histo2->GetBinContent(i); 
    binc = histo2->GetBinCenter(i); //j_T value
    histo2->SetBinContent(i, content*binc);
  }
  return histo2;
}

TH1D* addWeight(TH1D* histo){
  TH1D *histo2 = (TH1D*)histo->Clone();
  int N = histo2->GetNbinsX();
  double content = 0;
  double binc = 0;
  for(int i = 1; i <= N ; i++){
    content = histo2->GetBinContent(i); 
    binc = histo2->GetBinCenter(i); //j_T value
    histo2->SetBinContent(i, content/binc);
  }
  return histo2;
}

void RooUnfoldExample()
{
  //gROOT->LoadMacro("RooUnfold/src/RooUnfoldResponse.cxx");
  //gROOT->LoadMacro("RooUnfold/src/RooUnfoldBayes.cxx");
  //RooUnfoldResponse *jTresponse, jTresponse2;
  //gROOT->LoadMacro("AliJHistManager.cxx+");
  //hmg = new AliJHistManager();
  //TString inputname = "AnalysisResults1409.root";
  //TString inputname2 = "legotrain_318_20160526-1921_LHCb4_fix_CF_pPb_MC_ptHardMerged.root";
  //TString inputname = "legotrain_345_20160922-1851_LHCb4_fix_CF_pPb_MC_ptHardMerged.root";
  //TString inputname = "legotrain_349_20161031-1658_LHCb4_fix_CF_pPb_MC_ptHardMerged.root";
  TString inputname = "CF_pPb_MC_legotrain/legotrain_397_20170804-0029_LHCb4_fix_CF_pPb_MC_ptHardMerged.root";
  //TString dataName = "legotrain_CF_pPb-1053_20170223-2002_LHC13bc.root"; //Minimum Bias
  //TString dataName = "legotrain_CF_pPb-1053_20170223-2002_LHC13bcde.root"; // Combined
  TString dataName= "CF_pPb_legotrain/legotrain_CF_pPb_CF_pPb-1209_20170807_LHC13bcde.root";
  //TString dataName = "legotrain_CF_pPb-1053_20170223-2002_LHC13e_pass2.root"; //Triggered
  //TString dataName = "legotrain_CF_pPb-1053_20170223-2002_LHC13e_pass2.root";
  //TString inputname = "LHC13b4_fix_pthard5.root";
  //TString inputname = "legotrain_349_20161031-1658_LHCb4_fix_CF_pPb_MC_ptHard5.root";
  //TString inputname = "legotrain_318_20160526-1921_LHCb4_fix_CF_pPb_MC_ptHardMerged3.root";
  //TString inputname = "MC_data/legotrain_318_20160526-1921_LHCb4_fix_CF_pPb_MC/legotrain_318_20160526-1921_LHCb4_fix_CF_pPb_MC_ptHard5.root";
  //TString inputname = "legotrain_318_20160526-1921_LHCb4_fix_CF_pPb_MC_ptHard5.root";
  //TString inputname = "legotrain_336_20160624-1842_LHCb4_fix_CF_pPb_MC_ptHard5.root";
  inFile = new TFile(inputname, "update" );
  inFile2 = new TFile(inputname, "READ" );
  dataFile = new TFile(dataName, "READ" );
  inFile->Print();
  TString name;
  if(doWeight){
    outFile = new TFile("output_weight.root","recreate");
  }else{
    outFile = new TFile("output_noweight.root","recreate");
  }
  toyFile = new TFile("toyMC.root","read");

  //UnfoldToy();
  UnfoldJetPt();
  loadBg();
  cout << "Load jT:" << endl;
  loadJt();
  cout << "Load Data:" << endl;
  loadData();
  //UnfoldPythia();
  if(do2D){
    Unfold2D();
  }
  cout << "UnfoldJetCone" << endl;
  UnfoldJetCone();
  cout << "UnfoldJetConst" << endl;
  UnfoldJetConst();

  outFile->Close();
  inFile->Close();
  inFile2->Close();
  toyFile->Close();
  dataFile->Close();

}


#ifndef __CINT__
int main () { RooUnfoldExample(); return 0; }  // Main program when run stand-alone
#endif

//#include "include/common.h"
#include "include/Filipad2.h"
//#include "drawMacros.h"

TF1 *fitJtHisto(TH1D *histo,TString method, double cut, int ij, int iF);
TGraphErrors *makeErrorGraph(TH1D *histo, TH1D *errors);
TGraphErrors *makeErrorGraph(TGraphErrors *histo, TGraphErrors *errors);
TH1D *makeSystError(TH1D *histo1, TH1D *histo2);
TGraphErrors *makeSystError(TGraphErrors *histo1, TGraphErrors *histo2);
TGraphErrors *makeSystError(TGraphErrors *histo1, TGraphErrors *histo2,bool abs);
TGraphErrors *combineErrors(TGraphErrors *histo1, TGraphErrors *histo2);
TGraphErrors *combineErrors(TGraphErrors *histo1, TGraphErrors *histo2,TGraphErrors *histo3);
TH1D *combineErrors(TH1D *errors1, TH1D *errors2);
TH1D *combineErrors(TH1D *errors1, TH1D *errors2,TH1D *errors3);
int findStart(TH1D *histo, double start);
double getChi2(TH1D *histo, TF1 *fit,double lowX, double highX, int debug);
TGraphErrors *makeSignalGraph(TH1D *histo);
void drawFit(TH1D *histo, TF1 *fit,TString method,double xlow, double xhigh, int xlog, double ylow, double yhigh, int ylog, TString title,TString title2,TString comment,TString file, int is, int ij, int iF);
void drawFitScan(TH2D *histo,TString title,TString comment, TString file, int is, int ij,TString proj);
TH2D *scanFit(TH1D *histo, double cut, int ij, int iF,int is,TString proj,bool draw);
void scanFit2(TH1D *histo, double cut, int ij, int iF);
void drawComparison(TH1D **histo, double xlow, double xhigh,int xlog, double ylow, double yhigh, int ylog,TString title, TString file, int is, int ij);
void drawComparisonChi2(TGraph **histo, double xlow, double xhigh,int xlog, double ylow, double yhigh, int ylog,TString title, TString file);
void drawErrors(TGraphErrors *histo, TGraphErrors *errBg, TGraphErrors *errUnf, TGraphErrors *errFit, TGraphErrors *errTot,double xlow, double xhigh, int xlog, double ylow, double yhigh, int ylog, TString title,TString title2,TString comment,TString file, int is, int ij);
void drawErrors2(TGraphErrors *histo, TGraphErrors *errBg,double xlow, double xhigh, int xlog, double ylow, double yhigh, int ylog, TString title,TString title2,TString comment,TString file, int is, int ij,double limit);
void drawErrors(TGraphErrors *histo, TGraphErrors *errBg,double xlow, double xhigh, int xlog, double ylow, double yhigh, int ylog, TString title,TString title2,TString comment,TString file, int is, int ij);
TGraphErrors* smoothSystError(TGraphErrors *histo, bool draw, TString title);
void drawRMSHisto(TH1D *histo,TString title,TString file,int is,int ij, TString comment);

double JetPtBins2[] = {5,10,20,30,40,60,80,100,150,500};
double JetPtCenter[] = {7.5,15,25,35,50,70,90,125,325};
double JetPtError[] = {2.5,5,5,5,10,10,10,25,175};
const int NsetsTot = 5;
const int Nsets = 5;
const int Njets = 8;
int qColor[] = { kBlack, kRed, kBlue, kGreen+2, kOrange-3, kMagenta+1, kViolet, kPink, kGray, 40, 48, 9, kGreen+2, kOrange-3, kMagenta+1 };
int qMarker[] = { 20, 24, 21, 23, 27, 29,30,20, 24, 21, 25, 33, 27, 34, 28, 20, 21, 22, 23, 34, 28, 20, 21, 22 };

double gausRMS[Nsets][Njets];
double gausYield[Nsets][Njets];
double invgRMS[Nsets][Njets];
double invgYield[Nsets][Njets];
double gausRMSe[Nsets][Njets];
double gausYielde[Nsets][Njets];
double invgRMSe[Nsets][Njets];
double invgYielde[Nsets][Njets];

double fitB1[Nsets][Njets];
double fitB2[Nsets][Njets];
double fitB3[Nsets][Njets];
double fitB4[Nsets][Njets];
double fitB5[Nsets][Njets];
double fitPeak[Nsets][Njets];
double fitB1e[Nsets][Njets];
double fitB2e[Nsets][Njets];
double fitB3e[Nsets][Njets];
double fitB4e[Nsets][Njets];
double fitB5e[Nsets][Njets];
double fitPeake[Nsets][Njets];
double chi2[Nsets][Njets];
double chi2n[Nsets][Njets];
double chi2dof[Nsets][Njets];
double bins[Njets][Njets];

const double B1low = 0.1;
const double B1start[Njets] = {0.15,0.15,0.17,0.18,0.15,0.15,0.21,0.15};
const double B1high = 0.5;
const double B2low = 40;
const double B2high = 150;
const double B2start[Njets] = {50,50,66.20,67.90,50,50,75.89,77.62};
const int B3cut = 6;
const double B3low = 0;
const double B3high1 = 10;
const double B3high2 = 25;
const double B3start[Njets] = {7,7,5.06,4.90,10,10,7,12.93};
const double B4low = 2.0;
const double B4high = 15;
const double B4start[Njets] = {3,3,9.91,8.80,3,3,5.62,4.18};
const double B5low = 0.90;
const double B5high = 4.5;
const double B5start[Njets] = {1.5,1.5,4.5,4.5,1.5,1.5,2.88,1.62 };
const double peakstart = 0.5;
const double peaklow = 0.3;
const double peakhigh = 0.8;

//const TString setTitle[NsetsTot] = {"pPb LHC13b","pp LHC12h"};
//const TString taskName[Nsets] = {"AliJJetJtTask/AliJJetJtHistManager","AliJJetJtTask/AliJJetJtHistManager"};
//const TString setTitle[Nsets] = {"Pythia True","Pythia Measured","Pythia Bayes","Pythia SVD"};
//const TString taskName[Nsets] = {"Pythia/True","Pythia/Measured","Pythia/BayesSubUnfolding","Pythia/SVDSubUnfolding"};
//const TString topcomment = "_systematics_Triggered";
const TString topcomment = "_systematics";
const TString datasetTrig = "Triggered";
const TString datasetMB = "Minimum Bias";
const TString setTitle[NsetsTot] = {"Perp. Cone Bg + Bayes","Random Bg + Bayes","Perp. Cone  + SVD","Pythia True","Perp. Cone Bg + Bayes + fit"};
const TString taskName[NsetsTot] = {"Data/BayesSubUnfolding","Data/BayesSubUnfolding","Data/SVDSubUnfolding","Pythia/True","Data/BayesSubUnfolding"};
//const TString fileNames[Nsets] = {"unfoldedMB.root","unfoldedMB.root","unfoldedMB.root"};
const TString fileNames[NsetsTot] = {"unfoldedTriggered.root","unfoldedTriggered.root","unfoldedTriggered.root","unfoldedTriggered.root","unfoldedTriggered.root"};
//const TString fileNamesMB[Nsets] = {"unfoldedTriggered.root","unfoldedTriggered.root","unfoldedTriggered.root","unfoldedTriggered.root"};
const TString fileNamesMB[NsetsTot] = {"unfoldedMB.root","unfoldedMB.root","unfoldedMB.root","unfoldedMB.root","unfoldedMB.root"};
const TString comment[NsetsTot] = {"perconeBgBayes","randomBgBayes","perpconeBgSVD","PythiaTrue","ModFit"};
const int bgMethod[NsetsTot] = {0,1,0,0,0};
const int ij_Start = 4;
const int ij_Results = 4;
const int ij_MB = 4;
const double fitCut[NsetsTot] = {1,1,1,1,0.99};

//const TString setTitle[Nsets] = {"Measured","True","Bayes Sub","Bayes"};
//const TString taskName[Nsets] = {"Pythia/Measured","Pythia/True","Pythia/BayesSubUnfolding","Pythia/BayesUnfolding"};
const TString finderName[2] = {"Full jets, Anti-k_{T} R = 0.4","Charged Jets, Anti-k_{T} R = 0.4"};

int iC = 0;
int logx = 1;
const bool doWeight = 1;
const bool kScale = 0;
const double mSize = 0.3;
const bool saveFigs = true;
bool doAll = false;
const int nRebin = 2;
TFile *fin[Nsets];
TFile *finMB[Nsets];
const TString finderType[2] = {"Full","Charged"};
const int finderR[2] = {4,4};

void makeSystErrors(){
  //const int Nsets = 4;
  for(int iF = 0 ; iF < Nsets; iF++){
    cout << "Open " << setTitle[iF] << " file: " << fileNames[iF] << endl;
    fin[iF] = TFile::Open(fileNames[iF]);
    finMB[iF] = TFile::Open(fileNamesMB[iF]);
  }
  TString name;

  TH1D *hJetMultiplicity[Njets][Nsets];
  TH1D *hJtBg[Njets][Nsets];
  TH1D *BgTrkNumberBin[Njets][Nsets];
  TH1D *hJt[Njets][Nsets];
  TH1D *hJtSignal[Njets][Nsets];
  TH1D *hJtConst[Njets][Nsets];
  TH1D *hJtSignalConst[Njets][Nsets];
  TH1D *hJetPt[Nsets];
  TH1D *hJetPtBin[Nsets][Njets];
  TH1D *errorsBg[Njets];
  TH1D *errorsUnf[Njets];
  TH1D *errorsTot[Njets];
  TH1D *errorsFit[Njets];
  TH2D *hChi2Scans[Njets][Nsets];
  TGraphErrors *errGraph[Njets];
  TGraphErrors *errGraphBg[Njets];
  TGraphErrors *errGraphFit[Njets];
  TGraphErrors *hJtSignalGraph[Njets];
  TGraphErrors *hJtPythia[Njets];

  TGraphErrors *errorsBgGausRMS;
  TGraphErrors *errorsFitGausRMS;
  TGraphErrors *errorsUnfGausRMS;
  TGraphErrors *errorsTotGausRMS;
  TGraphErrors *errorsBgGausYield;
  TGraphErrors *errorsFitGausYield;
  TGraphErrors *errorsUnfGausYield;
  TGraphErrors *errorsTotGausYield;
  TGraphErrors *GausRMSerrGraph;
  TGraphErrors *GausRMSerrGraphBg;
  TGraphErrors *GausYielderrGraph;
  TGraphErrors *GausRMSerrGraphFitted;
  TGraphErrors *GausYielderrGraphFitted;

  TGraphErrors *errorsBgGausRMSFitted;
  TGraphErrors *errorsFitGausRMSFitted;
  TGraphErrors *errorsUnfGausRMSFitted;
  TGraphErrors *errorsTotGausRMSFitted;
  TGraphErrors *errorsBgGausYieldFitted;
  TGraphErrors *errorsFitGausYieldFitted;
  TGraphErrors *errorsUnfGausYieldFitted;
  TGraphErrors *errorsTotGausYieldFitted;
  //TGraphErrors *GausRMSerrGraphBg;

  TGraphErrors *errorsBginvgRMS;
  TGraphErrors *errorsFitinvgRMS;
  TGraphErrors *errorsUnfinvgRMS;
  TGraphErrors *errorsTotinvgRMS;
  TGraphErrors *errorsBginvgYield;
  TGraphErrors *errorsFitinvgYield;
  TGraphErrors *errorsUnfinvgYield;
  TGraphErrors *errorsTotinvgYield;
  TGraphErrors *invgRMSerrGraph;
  TGraphErrors *invgRMSerrGraphFitted;
  TGraphErrors *invgRMSerrGraphBg;
  TGraphErrors *invgYielderrGraph;
  TGraphErrors *invgYielderrGraphFitted;
  TGraphErrors *errorsBginvgRMSFitted;
  TGraphErrors *errorsFitinvgRMSFitted;
  TGraphErrors *errorsUnfinvgRMSFitted;
  TGraphErrors *errorsTotinvgRMSFitted;
  TGraphErrors *errorsBginvgYieldFitted;
  TGraphErrors *errorsFitinvgYieldFitted;
  TGraphErrors *errorsUnfinvgYieldFitted;
  TGraphErrors *errorsTotinvgYieldFitted;

  //TGraphErrors *invgRMSerrGraphBg;

  TF1 *gausfits[Njets][Nsets];
  TGraphErrors *gausRMSg[Nsets];
  TGraphErrors *gausYieldg[Nsets];
  TGraphErrors *invgRMSg[Nsets];
  TGraphErrors *invgYieldg[Nsets];
  TGraph* chi2g[Nsets];

  TGraphErrors *B1g[Nsets];
  TGraphErrors *B2g[Nsets];
  TGraphErrors *B3g[Nsets];
  TGraphErrors *B4g[Nsets];
  TGraphErrors *B5g[Nsets];

  int n_jet[Njets][Nsets];
  int n_bg[Njets][Nsets];

  int is = 0;
  TFile *errorOutput = TFile::Open("errors_test.root","recreate");


  for(int ij = ij_Start; ij < Njets; ij++){
    cout << "ij: " << ij << Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[ij], JetPtBins2[ij+1]) << endl;
    for(int iF = 0; iF < Nsets;iF++){
      cout << "\t" << fileNames[iF] << endl;
      cout << "\t" << taskName[iF] << endl;

      if(doAll){
        name = Form("%s/JetMultiplicityBin/JetMultiplicityBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
        if(ij > ij_MB) {
          hJetMultiplicity[ij][iF]  = (TH1D*)fin[iF]->Get(name);
        }else{
          hJetMultiplicity[ij][iF] = (TH1D*)finMB[iF]->Get(name);
        }
      }

      if(bgMethod[iF] == 0){ 
        if(doWeight){
          name = Form("%s/BgJtWeightBin/BgJtWeightBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
        }else{
          name = Form("%s/BgJtBin/BgJtBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
        }
      }else{
        if(doWeight){
          name = Form("%s/BgRndmJtWeightBin/BgRndmJtWeightBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
        }else{
          name = Form("%s/BgRndmJtBin/BgRndmJtBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
        }
      }
      cout << "\t\t" << name << endl;
      if(ij > ij_MB) {
        hJtBg[ij][iF] = (TH1D*)fin[iF]->Get(name);
      }else{
        hJtBg[ij][iF] = (TH1D*)finMB[iF]->Get(name);
      }
      cout << "\t\t";
      //hJtBg[ij][iF]->Print();


      if(kScale){
        name = Form("%s/BgTrkNumberBin/BgTrkNumberBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
        if(ij > ij_MB){
          BgTrkNumberBin[ij][iF] = (TH1D*)fin[iF]->Get(name);
        }else{
          BgTrkNumberBin[ij][iF] = (TH1D*)finMB[iF]->Get(name);
        }
        n_bg[ij][iF] = BgTrkNumberBin[ij][iF]->GetEntries();
        name = Form("%s/JetPtBin/JetPtBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
        if(ij > ij_MB){
          hJetPtBin[ij][iF] = (TH1D*)fin[iF]->Get(name);
        }else{
          hJetPtBin[ij][iF] = (TH1D*)finMB[iF]->Get(name);
        }
        n_jet[ij][iF] = hJetPtBin[ij][iF]->GetEntries();
      }else{
        n_bg[ij][iF] = 1;
        n_jet[ij][iF] = 1;
      }
      hJtBg[ij][iF]->Scale(1.0/n_bg[ij][iF], "width");
      hJtBg[ij][iF]->Rebin(nRebin);
      hJtBg[ij][iF]->Scale(1.0/nRebin);
    }

    for(int iF = 0; iF < Nsets; iF++){
      cout << "\tiF: " << iF << endl;
      cout << "\t" << fileNames[iF] << endl;
      if(doWeight){
        name = Form("%s/JetConeJtWeightBin/JetConeJtWeightBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
      }else{
        name = Form("%s/JetConeJtBin/JetConeJtBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
      }
      if(ij > ij_MB){
        hJt[ij][iF] = (TH1D*)fin[iF]->Get(name);
      }else{
        hJt[ij][iF] = (TH1D*)finMB[iF]->Get(name);
      }
      cout << "\t\t";
      //hJt[ij][iF]->Print();
      hJt[ij][iF]->Scale(1.0/n_jet[ij][iF],"width");
      hJt[ij][iF]->Rebin(nRebin);
      hJt[ij][iF]->Scale(1.0/nRebin);
      hJtSignal[ij][iF] = (TH1D*)hJt[ij][iF]->Clone();
      hJtSignal[ij][iF]->Add(hJtBg[ij][iF],-1);
      cout << "ij: " << ij << " iF: " << iF << "hJtSignal Integral: " << hJtSignal[ij][iF]->Integral() << endl;
    }
    for(int iF = 0; iF < Nsets; iF++){
      cout << "\tiF: " << iF << "/" << Nsets << endl;
      if(doWeight){
        name = Form("%s/JtWeightBin/JtWeightBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
      }else{
        name = Form("%s/JtBin/JtBinNFin%02dJetPt%02d",taskName[iF].Data(),is,ij);
      }
      cout << "\t\t" << name << endl;
      if(ij > ij_MB){
        hJtConst[ij][iF] = (TH1D*)fin[iF]->Get(name);
      }else{
        hJtConst[ij][iF] = (TH1D*)finMB[iF]->Get(name);
      }
      cout << "\t\t";
      //hJtConst[ij][iF]->Print();
      hJtConst[ij][iF]->Scale(1.0/n_jet[ij][iF],"width");
      hJtConst[ij][iF]->Rebin(nRebin);
      hJtConst[ij][iF]->Scale(1.0/nRebin);
      hJtSignalConst[ij][iF] = (TH1D*)hJtConst[ij][iF]->Clone();
      hJtSignalConst[ij][iF]->Add(hJtBg[ij][iF],-1);
    }
    cout << "\t\tDone" << endl;
    double ylow = 5e-5;
    double yhigh = 1e3;
    double xlow,xhigh;
    if(logx){
      xlow = 0.1;
      xhigh = 10;
    }else{
      xlow = 0;
      xhigh = 5;
    }
    //drawComparison(hJt[ij],xlow,xhigh,1,ylow,yhigh,1,"Jet cone","JetConejT/JetConejT",is,ij);

    //drawComparison(hJtSignal[ij],xlow,xhigh,logx,ylow,yhigh,1,"Jet Cone Subtracted", "JetConejTSignal/JetConejTSignal",is,ij);

    /*drawComparison(hJtSignalConst[ij],xlow,xhigh,logx,ylow,yhigh,1,"Subtracted jet constituent","ConstjTSignal/ConstjTSignal",is,ij);
      drawComparison(hJtConst[ij],xlow,xhigh,logx,ylow,yhigh,1,"Jet Constituents","ConstjT/ConstjT",is,ij);
      drawComparisonBg(hJt[ij],hJtBg[ij],xlow,xhigh,logx,ylow,yhigh,1,"Jet Cone","JetConejTIncl/JetConejTIncl",is,ij);*/
    TString proj = "xz";
    for(int iF = 0; iF < Nsets; iF++){
      cout << "\t\tiF: " << iF << "/" << Nsets << endl;
      gausfits[ij][iF] =  fitJtHisto(hJtSignal[ij][iF],"norm",fitCut[iF],ij,iF);
      //if(iF == 0 && ij == 5) scanFit2(hJtSignal[ij][iF],fitCut[iF],ij,iF);
      //drawFit(hJtSignal[ij][iF],gausfits[ij][iF],"norm",xlow,xhigh,logx,ylow,yhigh,1,"Subtracted jet cone",setTitle[iF],comment[iF],"JetConejTSignalFit/JetConejTSignalFit",is,ij, iF); //TODO Uncomment
      if(iF == 3 ){
        //drawFit(hJtSignal[ij][iF],gausfits[ij][iF],"norm",xlow,xhigh,logx,ylow,yhigh,1,"Subtracted jet cone","Pythia #sqrt{s_{NN}} = 5.02 TeV",comment[iF],"JetConejTSignalFit/JetConejTSignalFit",is,ij, iF); //TODO Uncomment
      }else{
        //drawFit(hJtSignal[ij][iF],gausfits[ij][iF],"norm",xlow,xhigh,logx,ylow,yhigh,1,"Subtracted jet cone","",comment[iF],"JetConejTSignalFit/JetConejTSignalFit",is,ij, iF); //TODO Uncomment

      }
      if(ij > 9){
        hChi2Scans[ij][iF] = scanFit(hJtSignal[ij][iF],fitCut[iF],ij,iF,is,proj,true);
        drawFit(hJtSignal[ij][iF],gausfits[ij][iF],"norm",xlow,xhigh,logx,ylow,yhigh,1,"Subtracted jet cone",setTitle[iF],comment[iF],"JetConejTSignalFit/JetConejTSignalFit",is,ij, iF);
        //drawFit(hJtSignal[ij][iF],gausfits[ij][iF],xlow,xhigh,logx,ylow,yhigh,1,"Subtracted Jet Constituents", setTitle[iF],comment[iF], "JetConstjTSignalFit/JetConstjTSignalFit",is,ij);
        //drawFitScan(hChi2Scans[ij][iF],"Title",comment[iF],"hChi2Scans",is,ij,proj);
      }
    }
    cout << "\t\tiF Loop finished" << endl;
    //drawComparison(hJtBg[ij],xlow,xhigh,logx,ylow,yhigh,1,"Background","SystematicErrors/Background/BackgroundComparison",is,ij);

    //drawComparison(gausfits[ij],xlow,xhigh,logx,ylow,yhigh,1,"Fits","JetConejTSignalFitComparison/JetConejTSignalFitComparison",is,ij);

    //drawErrorComparison(hJtSignal[ij][0], hJtSignal[ij][1], xlow, xhigh,logx, ylow, yhigh, 1,"Background", "Perpendicular Cone", "Random", "SystematicErrors/Background/Background", is, ij);

    for(int iF = 0; iF < Nsets; iF++){
      cout << "iF: " << iF << "/" << Nsets << endl;
      hJtSignal[ij][iF]->Rebin(2);
      hJtSignal[ij][iF]->Scale(1.0/2);
      cout << hJtSignal[ij][iF]->GetNbinsX() << endl;
    }
    errorsBg[ij] = makeSystError(hJtSignal[ij][0],hJtSignal[ij][1]);
    errorsBg[ij]->Print();
    errorsUnf[ij] = makeSystError(hJtSignal[ij][0],hJtSignal[ij][2]);
    errorsBg[ij]->SetName(Form("%s_Bg",hJtSignal[ij][0]->GetName()));
    errorsUnf[ij]->SetName(Form("%s_Unf",hJtSignal[ij][1]->GetName()));
    errorsUnf[ij]->Print();
    errorsFit[ij] = makeSystError(hJtSignal[ij][0],hJtSignal[ij][4]);
    errorsTot[ij] = combineErrors(errorsBg[ij], errorsUnf[ij],errorsFit[ij]);
    errorsTot[ij]->SetName(Form("%s_Tot",hJtSignal[ij][1]->GetName()));
    errorsBg[ij]->Print();
    //errGraphBg[ij] = makeErrorGraph(hJtSignal[ij][1], errorsBg[ij]);
    errGraphBg[ij] = makeErrorGraph(hJtSignal[ij][0],errorsBg[ij]);
    errGraphBg[ij]->SetName(Form("%s_BgSystematics",hJtSignal[ij][1]->GetName()));
    errGraphFit[ij] = makeErrorGraph(hJtSignal[ij][0],errorsFit[ij]);
    errGraphFit[ij]->SetName(Form("%s_FitSystematics",hJtSignal[ij][0]->GetName()));

    errGraph[ij] = makeErrorGraph(hJtSignal[ij][0],errorsTot[ij]);
    errGraph[ij]->SetName(Form("%s_Systematics",hJtSignal[ij][0]->GetName()));
    hJtSignalGraph[ij] = makeSignalGraph(hJtSignal[ij][0]);
    hJtSignalGraph[ij]->SetName(Form("%s_Statistics",hJtSignal[ij][0]->GetName()));
    hJtPythia[ij] = makeSignalGraph(hJtSignal[ij][3]);
    hJtPythia[ij]->SetName(Form("%s_Pythia",hJtSignal[ij][3]->GetName()));


    //errGraph[ij] = makeErrorGraph(hJtSignal[ij][1],errorsBg[ij]);
    //errGraph[ij] = makeErrorGraph(hJtSignal[ij][1],errorsUnf[ij]);

    //drawWithErrors(hJtSignal[ij][1], errGraph[ij],xlow, xhigh, logx, ylow, yhigh, 1, "Jet Cone Subtracted",dataset,"_data","JetConejtWithSystematics/JetConejtWithSystematics", is, ij);
    //drawErrors(hJtSignal[ij][1],errorsBg[ij],errorsUnf[ij],errorsTot[ij], xlow,xhigh,logx,ylow,yhigh,1,"Jet Cone Subtracted",dataset,"_data","SystematicErrors/SystematicErrors",is,ij);
    //drawErrors(hJtSignal[ij][1],errorsBg[ij],errorsUnf[ij],errorsTot[ij], ,double xlow, double xhigh, int xlog, double ylow, double yhigh, int ylog, TString title,TString title2,TString comment,TString file, int is, int ij){

    //cout << "hJtSignal[ij][0] name: " << hJtSignal[ij][0]->GetName() << endl;
    //TH1D *errors = makeSystError(hJtSignal[ij][0],hJtSignal[ij][1]);
    //TH1D *errors = makeSystErrorRelative(hJtSignal[ij][0],hJtSignal[ij][1]);
    //errors->SetName(hJtSignal[ij][0]->GetName());
    errorOutput->cd();
    cout << Form("JetConeJtWeightBinNFin%02dJetPt%02d_FitFunction",0,ij) << endl;
    gausfits[ij][0]->SetName(Form("JetConeJtWeightBinNFin%02dJetPt%02d_FitFunction",0,ij));
    gausfits[ij][0]->Print();
    gausfits[ij][0]->Write();
    errorsBg[ij]->Write();
    errorsUnf[ij]->Write();
    errorsFit[ij]->Write();
    errGraph[ij]->Write();
    errGraphBg[ij]->Write();
    errGraphFit[ij]->Write();
    hJtSignalGraph[ij]->Write();
    hJtPythia[ij]->Write();
    for(int iF = 0; iF < Nsets; iF++){
      hJtSignal[ij][iF]->Delete();
    }
    hJtPythia[ij]->Delete();
    hJtSignalGraph[ij]->Delete();
    errorsBg[ij]->Delete();
    errorsUnf[ij]->Delete();
    errorsFit[ij]->Delete();
    errGraph[ij]->Delete();

  }
  for(int iF = 0; iF < Nsets; iF++){
    //cout << "iF " << iF << endl;
    B1g[iF] = new TGraphErrors(Njets,JetPtCenter,fitB1[iF],JetPtError,fitB1e[iF]);
    B2g[iF] = new TGraphErrors(Njets,JetPtCenter,fitB2[iF],JetPtError,fitB2e[iF]);
    B3g[iF] = new TGraphErrors(Njets,JetPtCenter,fitB3[iF],JetPtError,fitB3e[iF]);
    B4g[iF] = new TGraphErrors(Njets,JetPtCenter,fitB4[iF],JetPtError,fitB4e[iF]);
    B5g[iF] = new TGraphErrors(Njets,JetPtCenter,fitB5[iF],JetPtError,fitB5e[iF]);

    gausRMSg[iF] = new TGraphErrors(Njets,JetPtCenter,gausRMS[iF],JetPtError, gausRMSe[iF]);
    gausYieldg[iF] = new TGraphErrors(Njets,JetPtCenter,gausYield[iF],JetPtError, gausYielde[iF]);
    invgRMSg[iF] = new TGraphErrors(Njets,JetPtCenter,invgRMS[iF],JetPtError, invgRMSe[iF]);
    invgYieldg[iF] = new TGraphErrors(Njets,JetPtCenter,invgYield[iF],JetPtError, invgYielde[iF]);
    chi2g[iF] = new TGraph(Njets,JetPtCenter,chi2dof[iF]);
    gausRMSg[iF]->SetName(Form("gGausRMS%02d",iF));
    gausYieldg[iF]->SetName(Form("gGausYield%02d",iF));
    invgRMSg[iF]->SetName(Form("gGammaRMS%02d",iF));
    invgYieldg[iF]->SetName(Form("gGammaYield%02d",iF));
    chi2g[iF]->SetName(Form("gChi2%02d",iF));
    gausRMSg[iF]->Write();
    gausYieldg[iF]->Write();
    invgRMSg[iF]->Write();
    invgYieldg[iF]->Write();
  }
  TGraphErrors *errorsBgGausRMS2 = makeSystError(gausRMSg[0],gausRMSg[1],false);
  TGraphErrors *errorsUnfGausRMS2 = makeSystError(gausRMSg[0],gausRMSg[2],false);
  TGraphErrors *errorsFitGausRMS2 = makeSystError(gausRMSg[0],gausRMSg[4],false);
  TGraphErrors *errorsBgGausYield2 = makeSystError(gausYieldg[0],gausYieldg[1],false);
  TGraphErrors *errorsUnfGausYield2 = makeSystError(gausYieldg[0],gausYieldg[2],false);
  TGraphErrors *errorsFitGausYield2 = makeSystError(gausYieldg[0],gausYieldg[4],false);
  TGraphErrors *errorsBginvgRMS2 = makeSystError(invgRMSg[0],invgRMSg[1],false);
  TGraphErrors *errorsUnfinvgRMS2 = makeSystError(invgRMSg[0],invgRMSg[2],false);
  TGraphErrors *errorsFitinvgRMS2 = makeSystError(invgRMSg[0],invgRMSg[4],false);
  TGraphErrors *errorsBginvgYield2 = makeSystError(invgYieldg[0],invgYieldg[1],false);
  TGraphErrors *errorsUnfinvgYield2 = makeSystError(invgYieldg[0],invgYieldg[2],false);
  TGraphErrors *errorsFitinvgYiel2 = makeSystError(invgYieldg[0],invgYieldg[4],false);

  TString dataset = datasetTrig;
  drawErrors2(gausRMSg[0],errorsBgGausRMS2,40,150,0,-0.2,0.2,0,"Narrow RMS","Background","_data","SystematicErrors/SystematicErrorsGausRMS_Bg",is,Njets,0.09);
  drawErrors2(gausRMSg[0],errorsUnfGausRMS2,40,150,0,-0.2,0.2,0,"Narrow RMS","Unfolding","_data","SystematicErrors/SystematicErrorsGausRMS_Unf",is,Njets,0.08);
  drawErrors2(invgRMSg[0],errorsBginvgRMS2,40,150,0,-0.2,0.2,0,"Wide RMS","Background","_data","SystematicErrors/SystematicErrorsGammaRMS_Bg",is,Njets,0.05);
  drawErrors2(invgRMSg[0],errorsUnfinvgRMS2,40,150,0,-0.2,0.2,0,"Wide RMS","Unfolding","_data","SystematicErrors/SystematicErrorsGammaRMS_Unf",is,Njets,0.08);
  drawErrors2(gausRMSg[0],errorsBgGausYield2,40,150,0,-0.2,0.2,0,"Narrow RMS","Background","_data","SystematicErrors/SystematicErrorsGausYield_Bg",is,Njets,0.09);
  drawErrors2(gausRMSg[0],errorsUnfGausYield2,40,150,0,-0.2,0.2,0,"Narrow RMS","Unfolding","_data","SystematicErrors/SystematicErrorsGausYield_Unf",is,Njets,0.08);
  drawErrors2(invgRMSg[0],errorsBginvgYield2,40,150,0,-0.2,0.2,0,"Wide RMS","Background","_data","SystematicErrors/SystematicErrorsGammaYield_Bg",is,Njets,0.05);
  drawErrors2(invgRMSg[0],errorsUnfinvgYield2,40,150,0,-0.2,0.2,0,"Wide RMS","Unfolding","_data","SystematicErrors/SystematicErrorsGammaYield_Unf",is,Njets,0.08);
  drawComparisonChi2(chi2g, 0, 200,0, 0, 5, 0,"Chi2/dof", "chi2");
  //return;
  //GausRMS errors:
  bool kDraw = true;
  errorsBgGausRMS = makeSystError(gausRMSg[0],gausRMSg[1]);
  errorsBgGausRMS->SetName(Form("%s_Bg", gausRMSg[0]->GetName()));
  errorsBgGausRMSFitted = smoothSystError(errorsBgGausRMS,kDraw,"Gaus RMS Bg");
  errorsUnfGausRMS = makeSystError(gausRMSg[0],gausRMSg[2]);
  errorsUnfGausRMS->SetName(Form("%s_Unf",gausRMSg[0]->GetName()));
  errorsUnfGausRMSFitted = smoothSystError(errorsUnfGausRMS,kDraw,"Gaus RMS Unfolding");
  errorsFitGausRMS = makeSystError(gausRMSg[0],gausRMSg[4]);
  errorsFitGausRMS->SetName(Form("%s_Fit",gausRMSg[0]->GetName()));
  errorsFitGausRMSFitted = smoothSystError(errorsFitGausRMS, true, "Gaus RMS Fit");
  //errorsTotGausRMS = combineErrors(errorsBgGausRMS,errorsUnfGausRMS,errorsFitGausRMS);
  errorsTotGausRMS = combineErrors(errorsBgGausRMS,errorsUnfGausRMS);
  errorsTotGausRMSFitted = combineErrors(errorsBgGausRMSFitted,errorsUnfGausRMSFitted,errorsFitGausRMSFitted);
  errorsTotGausRMS->SetName(Form("%s_Tot",gausRMSg[0]->GetName()));
  GausRMSerrGraph = makeErrorGraph(gausRMSg[0],errorsTotGausRMS);
  GausRMSerrGraph->SetName(Form("%s_Systematics",gausRMSg[0]->GetName()));
  GausRMSerrGraph->Write();
  GausRMSerrGraphFitted = makeErrorGraph(gausRMSg[0],errorsTotGausRMSFitted);
  GausRMSerrGraphFitted->SetName(Form("%s_FittedSystematics",gausRMSg[0]->GetName()));
  GausRMSerrGraphFitted->Write();

  errorsBgGausRMS->Write();
  errorsUnfGausRMS->Write();
  errorsFitGausRMS->Write();
  errorsTotGausRMS->Write();



  /*GausRMSerrGraphBg = makeErrorGraph(gausRMSg[1],errorsBgGausRMS);
  GausRMSerrGraphBg->SetName(Form("%s_BgSystematics",gausRMSg[1]->GetName()));
  GausRMSerrGraphBg->Write();
  cout << GausRMSerrGraphBg->GetName() << " Written " << endl;*/

  //GausYield errors:
  errorsBgGausYield = makeSystError(gausYieldg[0],gausYieldg[1]);
  errorsBgGausYield->SetName(Form("%s_Bg", gausYieldg[0]->GetName()));
  errorsBgGausYieldFitted = smoothSystError(errorsBgGausYield,kDraw,"Gaus Yield Bg");
  errorsFitGausYield = makeSystError(gausYieldg[0],gausYieldg[4]);
  errorsFitGausYield->SetName(Form("%s_Fit", gausYieldg[0]->GetName()));
  errorsFitGausYieldFitted = smoothSystError(errorsFitGausYield,kDraw,"Gaus Yield Fit");
  errorsUnfGausYield = makeSystError(gausYieldg[0],gausYieldg[2]);
  errorsUnfGausYield->SetName(Form("%s_Unf",gausYieldg[0]->GetName()));
  errorsUnfGausYieldFitted = smoothSystError(errorsUnfGausYield,kDraw,"Gaus Yield Unfolding");
  errorsTotGausYield = combineErrors(errorsBgGausYield,errorsUnfGausYield);
  errorsTotGausYieldFitted = combineErrors(errorsBgGausYieldFitted,errorsUnfGausYieldFitted,errorsFitGausYieldFitted);
  GausYielderrGraph = makeErrorGraph(gausYieldg[0],errorsTotGausYield);
  GausYielderrGraph->SetName(Form("%s_Systematics",gausYieldg[0]->GetName()));
  GausYielderrGraph->Write();
  cout << GausYielderrGraph->GetName() << " Written " << endl;
  GausYielderrGraphFitted = makeErrorGraph(gausYieldg[0],errorsTotGausYieldFitted);
  GausYielderrGraphFitted->SetName(Form("%s_FittedSystematics",gausYieldg[0]->GetName()));
  GausYielderrGraphFitted->Write();

  errorsBgGausYield->Write();
  errorsUnfGausYield->Write();
  errorsFitGausYield->Write();
  errorsTotGausYield->Write();

  //invgRMS errors:
  errorsBginvgRMS = makeSystError(invgRMSg[0],invgRMSg[1]);
  errorsBginvgRMS->SetName(Form("%s_Bg", invgRMSg[0]->GetName()));
  errorsBginvgRMSFitted = smoothSystError(errorsBginvgRMS,kDraw,"Gamma RMS Bg");
  errorsUnfinvgRMS = makeSystError(invgRMSg[0],invgRMSg[2]);
  errorsUnfinvgRMS->SetName(Form("%s_Unf",invgRMSg[0]->GetName()));
  errorsUnfinvgRMSFitted = smoothSystError(errorsUnfinvgRMS,kDraw,"Gamma RMS Unfolding");
  errorsFitinvgRMS = makeSystError(invgRMSg[0],invgRMSg[4]);
  errorsFitinvgRMS->SetName(Form("%s_Fit",invgRMSg[0]->GetName()));
  errorsFitinvgRMSFitted = smoothSystError(errorsFitinvgRMS,kDraw,"Gamma RMS Fit");
  //errorsTotinvgRMS = combineErrors(errorsBginvgRMS,errorsUnfinvgRMS,errorsFitinvgRMS);
  errorsTotinvgRMS = combineErrors(errorsBginvgRMS,errorsUnfinvgRMS);
  errorsTotinvgRMSFitted = combineErrors(errorsBginvgRMSFitted,errorsUnfinvgRMSFitted,errorsFitinvgRMSFitted);
  invgRMSerrGraph = makeErrorGraph(invgRMSg[0],errorsTotinvgRMS);
  invgRMSerrGraph->SetName(Form("%s_Systematics",invgRMSg[0]->GetName()));
  invgRMSerrGraph->Write();
  invgRMSerrGraphFitted = makeErrorGraph(invgRMSg[0],errorsTotinvgRMSFitted);
  invgRMSerrGraphFitted->SetName(Form("%s_FittedSystematics",invgRMSg[0]->GetName()));
  invgRMSerrGraphFitted->Write();

  /*invgRMSerrGraphBg = makeErrorGraph(invgRMSg[0],errorsBginvgRMS);
  invgRMSerrGraphBg->SetName(Form("%s_BgSystematics",invgRMSg[1]->GetName()));
  invgRMSerrGraphBg->Write();
  cout << invgRMSerrGraphBg->GetName() << " Written " << endl;*/

  errorsBginvgRMS->Write();
  errorsUnfinvgRMS->Write();
  errorsFitinvgRMS->Write();
  errorsTotinvgRMS->Write();

  //invgYield errors:
  errorsBginvgYield = makeSystError(invgYieldg[0],invgYieldg[1]);
  errorsBginvgYieldFitted = smoothSystError(errorsBginvgYield,kDraw,"Gamma Yield Bg");
  errorsBginvgYield->SetName(Form("%s_Bg", invgYieldg[0]->GetName()));
  errorsUnfinvgYield = makeSystError(invgYieldg[0],invgYieldg[2]);
  errorsUnfinvgYieldFitted = smoothSystError(errorsUnfinvgYield,kDraw,"Gamma Yield Unfolding");
  errorsUnfinvgYield->SetName(Form("%s_Unf",invgYieldg[0]->GetName()));
  errorsFitinvgYield = makeSystError(invgYieldg[0],invgYieldg[4]);
  errorsFitinvgYield->SetName(Form("%s_Fit",invgYieldg[0]->GetName()));
  errorsFitinvgYieldFitted = smoothSystError(errorsFitinvgYield,kDraw,"Gamma Yield Fit");
  errorsTotinvgYield = combineErrors(errorsBginvgYield,errorsUnfinvgYield);
  errorsTotinvgYieldFitted = combineErrors(errorsBginvgYieldFitted,errorsUnfinvgYieldFitted);
  invgYielderrGraph = makeErrorGraph(invgYieldg[0],errorsTotinvgYield);
  invgYielderrGraph->SetName(Form("%s_Systematics",invgYieldg[0]->GetName()));
  invgYielderrGraph->Write();
  invgYielderrGraphFitted = makeErrorGraph(invgYieldg[0],errorsTotinvgYieldFitted);
  invgYielderrGraphFitted->SetName(Form("%s_FittedSystematics",invgYieldg[0]->GetName()));
  invgYielderrGraphFitted->Write();

  errorsBginvgYield->Write();
  errorsUnfinvgYield->Write();
  errorsFitinvgYield->Write();
  errorsTotinvgYield->Write();

  /*if(ij > ij_MB){
    dataset = datasetTrig;
    }else{
    dataset = datasetMB;
    }*/
  //drawErrors(gausRMSg[0],errorsBgGausRMS,0,150,0,0,0.5,0,"Jet Cone",dataset,"_data","SystematicErrors/SystematicErrorsGausRMSBg",is,ij);
  //drawErrors(gausRMSg[0],errorsBginvgRMS,0,150,0,0,0.5,0,"Jet Cone",dataset,"_data","SystematicErrors/SystematicErrorsinvgRMSBg",is,ij);
  //

  drawErrors(gausRMSg[0],errorsBgGausRMS, errorsUnfGausRMS,errorsFitGausRMS, errorsTotGausRMS,0,150,0,0,0.5,0,"Gaus RMS",dataset,"_data","SystematicErrors/SystematicErrorsGausRMS",is,Njets);
  drawErrors(gausYieldg[0],errorsBgGausYield, errorsUnfGausYield, errorsFitGausYield,errorsTotGausYield,0,150,0,0,5,0,"Gaus Yield",dataset,"_data","SystematicErrors/SystematicErrorsGausYield",is,Njets);
  drawErrors(invgRMSg[0],errorsBginvgRMS, errorsUnfinvgRMS, errorsFitinvgRMS,errorsTotinvgRMS,0,150,0,0,1,0,"Gamma RMS",dataset,"_data","SystematicErrors/SystematicErrorsinvgRMS",is,Njets);
  drawErrors(invgYieldg[0],errorsBginvgYield, errorsUnfinvgYield, errorsFitinvgYield,errorsTotinvgYield,0,150,0,0,5,0,"Gamma Yield",dataset,"_data","SystematicErrors/SystematicErrorsinvgYield",is,Njets);

  /*drawErrors(gausRMSg[0],errorsBgGausRMSFitted, errorsUnfGausRMSFitted,errorsFitGausRMSFitted, errorsTotGausRMSFitted,0,150,0,0,1.5,0,"Gaus RMS Fitted",dataset,"_data","SystematicErrors/SystematicErrorsGausRMSFitted",is,Njets);
    drawErrors(gausYieldg[0],errorsBgGausYieldFitted, errorsUnfGausYieldFitted, errorsFitGausYieldFitted,errorsTotGausYieldFitted,0,150,0,0,7,0,"Gaus Yield Fitted",dataset,"_data","SystematicErrors/SystematicErrorsGausYieldFitted",is,Njets);
    drawErrors(invgRMSg[0],errorsBginvgRMSFitted, errorsUnfinvgRMSFitted, errorsFitinvgRMSFitted,errorsTotinvgRMSFitted,0,150,0,0,1.5,0,"Gamma RMS Fitted",dataset,"_data","SystematicErrors/SystematicErrorsinvgRMSFitted",is,Njets);
    drawErrors(invgYieldg[0],errorsBginvgYieldFitted, errorsUnfinvgYieldFitted, errorsFitinvgYieldFitted,errorsTotinvgYieldFitted,0,150,0,0,7,0,"Gamma Yield Fitted",dataset,"_data","SystematicErrors/SystematicErrorsinvgYieldFitted",is,Njets);*/
  /*
     drawComparison2(gausRMSg,0,150,0,0,1.5,0,"Gaus RMS" ,"SystematicErrors/GausRMS/GausRMS"  ,"RMS",0,0);
     drawComparison2(invgRMSg,0,150,0,0,3,0,"Gamma RMS","SystematicErrors/GammaRMS/GammaRMS","RMS",0,0);
     drawComparison2(gausYieldg,0,150,0,0,20,0,"Gaus Yield" ,"SystematicErrors/GausYield/GausYield"  ,"Yield",0,0);
     drawComparison2(invgYieldg,0,150,0,0,20,0,"Gamma Yield","SystematicErrors/GammaYield/GammaYield","Yield",0,0);

     drawComparison2(B1g,0,150,0,0,0.3,0,"B1","SystematicErrors/B1/B1","B1",0,0);
     drawComparison2(B2g,0,150,0,0,200,0,"B2","SystematicErrors/B2/B2","B2",0,0);
     drawComparison2(B3g,0,150,0,0,30,0,"B3","SystematicErrors/B3/B3","B3",0,0);
     drawComparison2(B4g,0,150,0,0,15,0,"B4","SystematicErrors/B4/B4","B4",0,0);
     drawComparison2(B5g,0,150,0,0,7,0,"B5","SystematicErrors/B5/B5","B5",0,0);
     */
  TGraphErrors *gausRMSg2[1];
  TGraphErrors *invgRMSg2[1];
  gausRMSg2[0] = gausRMSg[0];
  invgRMSg2[0] = invgRMSg2[0];
  //drawComparison2(gausRMSg2,0,150,0,0,1.5,0,"Gaus RMS", "SystematicErrors/GausRMS/GausRMS0", "RMS",0,0);
  //drawComparison2(invgRMSg2,0,150,0,0,1.5,0,"Gaus RMS", "SystematicErrors/GausRMS/GausRMS0", "RMS",0,0);

  for(int iF = 0; iF < Nsets; iF++){
    fin[iF]->Close();
    finMB[iF]->Close();
  }
  errorOutput->Close();

  cout << "Gaus RMS \\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(gausRMS[iF][ij] > 0){
        cout << setprecision(3) <<  "&" << gausRMS[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << gausRMS[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }

  cout << "Gaus Yield\\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(gausYield[iF][ij] > 0){
        cout << setprecision(3) <<  "&" << gausYield[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << gausYield[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }

  cout << "Gamma RMS \\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(invgRMS[iF][ij] > 0){
        cout << setprecision(3) <<  "&" << invgRMS[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << invgRMS[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }

  cout << "Gamma Yield\\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(invgYield[iF][ij] > 0){
        cout << setprecision(3) <<  "&" << invgYield[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << invgYield[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }

  cout << "Chi squared \\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      cout << setprecision(2) <<  "&" << chi2[iF][ij];
    }
    cout << "\\\\" << endl;
  }
  cout << "Nonzero Bins\\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      cout << setprecision(2) <<  "&" << bins[iF][ij];
    }
    cout << "\\\\" << endl;
  }
  cout << "Chi squared/n \\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      cout << setprecision(2) <<  "&" << chi2n[iF][ij];
    }
    cout << "\\\\" << endl;
  }

  cout << "Prob \\\\" << endl;
  for(int iF = 0; iF < Nsets; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      cout << setprecision(5) << "&" << TMath::Prob(chi2[iF][ij],(bins[iF][ij]-5));
    }
    cout << "\\\\" << endl;
  }

  cout << "Chi squared/dof \\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(chi2dof[iF][ij] < 1.0){
        cout << setprecision(2) <<  "&" << chi2dof[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << chi2dof[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }

  cout << "B1\\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(fitB1[iF][ij] > 0 && fitB1[iF][ij] < B1high*0.99 && fitB1[iF][ij] > B1low * 1.01 && chi2dof[iF][ij] < 1.0){
        cout << setprecision(3) <<  "&" << fitB1[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << fitB1[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }


  cout << "B2\\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(fitB2[iF][ij] > 0 && fitB2[iF][ij] < B2high*0.99 && fitB2[iF][ij] > B2low * 1.01 && chi2dof[iF][ij] < 1.0){
        cout << setprecision(3) <<  "&" << fitB2[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << fitB2[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }
  cout << "B3\\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(fitB3[iF][ij] > 0 && fitB3[iF][ij] < B3high2*0.99 && fitB3[iF][ij] > B3low * 1.01 && chi2dof[iF][ij] < 1.0){
        cout << setprecision(3) <<  "&" << fitB3[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << fitB3[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }
  cout << "B4\\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(fitB4[iF][ij] > 0 && fitB4[iF][ij] < B4high*0.99 && fitB4[iF][ij] > B4low * 1.01 && chi2dof[iF][ij] < 1.0){
        cout << setprecision(3) <<  "&" << fitB4[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << fitB4[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }
  cout << "B5\\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(fitB5[iF][ij] > 0 && fitB5[iF][ij] < B5high*0.99 && fitB5[iF][ij] > B5low * 1.01 && chi2dof[iF][ij] < 1.0){
        cout << setprecision(3) <<  "&" << fitB5[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << fitB5[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }

  cout << "Gamma peak\\\\" << endl;
  for(int iF = 0; iF < Nsets ; iF++){
    cout << comment[iF].Data();
    for(int ij = ij_Start; ij < Njets; ij++){
      if(fitPeak[iF][ij] > 0 && fitPeak[iF][ij] < peakhigh*0.99 && fitPeak[iF][ij] > peaklow * 1.01 && chi2dof[iF][ij] < 1.0){
        cout << setprecision(3) <<  "&" << fitPeak[iF][ij];
      }else{
        cout << setprecision(3) <<  "& {\\color{red}" << fitPeak[iF][ij] << "}";
      }
    }
    cout << "\\\\" << endl;
  }


  }

  //MAIN END

  TF1 *fitJtHisto(TH1D *histo,TString method, double cut, int ij, int iF){
    TF1 *gaussfit = new TF1("gaussfit","gausn",0,10);
    gaussfit->FixParameter(1,0);
    gaussfit->SetParLimits(2,B1low,B1high);
    gaussfit->SetParLimits(0,B2low,B2high);
    histo->Fit("gaussfit","QN");
    TF1 *invG;
    if(method == "alt"){
      invG = new TF1("invG","[0]* (pow([1]*([2]+1),[2])/TMath::Gamma([2]))*exp(-[1]*([2]+1)/x) * pow(x, -[2]-1)",0,10);
    }else{
      invG = new TF1("invG","[0]* (pow([1],[2])/TMath::Gamma([2]))*exp(-[1]/x) * pow(x, -[2]-1)",0,10);
    }
    invG->SetParameter(0,B3start[ij]);
    if(ij < B3cut){
      invG->SetParLimits(0,B3low,B3high1);
    }else{
      invG->SetParLimits(0,B3low,B3high2);
    }
    if(method == "alt"){
      invG->SetParameter(1,peakstart);
      invG->SetParLimits(1,peaklow,peakhigh);
    }else{
      invG->SetParameter(1,B5start[ij]);
      invG->SetParLimits(1,B5low,B5high);
    }
    invG->SetParameter(2,B4start[ij]);
    invG->SetParLimits(2,B4low,B4high);
    histo->Fit("invG","QN","",cut,3);
    //invG->Print();
    TF1 *gaussfit3;
    if(method == "alt"){
      gaussfit3= new TF1("gaussfit3","gausn(0) + [3]* (pow([4]*([5]+1),[5])/TMath::Gamma([5]))*exp(-[4]*([5]+1)/x) * pow(x, -[5]-1)",0,10);
    }else{
      gaussfit3= new TF1("gaussfit3","gausn(0) + [3]* (pow([4],[5])/TMath::Gamma([5]))*exp(-[4]/x) * pow(x, -[5]-1)",0,10);
    }
    gaussfit3->SetParameter(0, gaussfit->GetParameter(0));
    gaussfit3->SetParameter(1, gaussfit->GetParameter(1));
    gaussfit3->SetParameter(2, gaussfit->GetParameter(2));
    gaussfit3->SetParameter(3, invG->GetParameter(0));
    gaussfit3->SetParameter(4, invG->GetParameter(1));
    gaussfit3->SetParameter(5, invG->GetParameter(2));
    gaussfit3->FixParameter(1,0);

    gaussfit3->SetParLimits(0,B2low,B2high); //B2
    gaussfit3->SetParLimits(2,B1low,B1high); //B1
    if(ij < B3cut){
      gaussfit3->SetParLimits(3,B3low,B3high1); //B3
    }else{
      gaussfit3->SetParLimits(3,B3low,B3high2); //B3
    }
    if(method == "alt"){
      gaussfit3->SetParLimits(4,peaklow,peakhigh); //Peak
    }else{
      gaussfit3->SetParLimits(4,B5low,B5high); //B5
    }
    gaussfit3->SetParLimits(5,B4low,B4high); //B4

    bins[iF][ij] = histo->FindLastBinAbove(1e-5);
    double end = histo->GetBinCenter(bins[iF][ij]);
    bins[iF][ij] = bins[iF][ij] - findStart(histo,0.1);
    histo->Fit("gaussfit3","QN","",0.1,end);
    //chi2[iF][ij] = gaussfit3->GetChisquare();
    chi2[iF][ij] = getChi2(histo,gaussfit3,0.1,end,0);
    chi2dof[iF][ij] = chi2[iF][ij]/(bins[iF][ij]-5);
    chi2n[iF][ij] = chi2[iF][ij]/bins[iF][ij];
    double B2 = gaussfit3->GetParameter(0);
    double B1 = gaussfit3->GetParameter(2);
    double B2e = gaussfit3->GetParError(0);
    double B1e = gaussfit3->GetParError(2);
    double gaussigma = TMath::Sqrt(2)*B1;
    double gausyield = B2*B1/TMath::Sqrt(2 * TMath::Pi());
    double gausyielde = TMath::Sqrt((B2*B2*B1e*B1e + B1*B1*B2e*B2e)/(2*TMath::Pi()));
    double gaussigmae = TMath::Sqrt(2)*B1e;

    double constant = gaussfit3->GetParameter(3); //B3
    double constante = gaussfit3->GetParError(3);
    double alpha = gaussfit3->GetParameter(5); // B4
    double alphae = gaussfit3->GetParError(5); 
    double peak,peake,beta,betae;
    if(method == "alt"){
      peak = gaussfit3->GetParameter(4); 
      peake = gaussfit3->GetParError(4);
      beta = peak*(alpha+1);//B5
      betae = TMath::Sqrt(pow( (alpha+1) * peake,2) + pow(peak*alphae,2));
    }else{
      beta = gaussfit3->GetParameter(4); // B5
      betae = gaussfit3->GetParError(4);
      peak = beta/(alpha+1);
      peake = TMath::Sqrt(pow(betae/(alpha+1),2)+pow(alphae*beta/pow(alpha+1,2),2));
    }
    double gammaYield = constant * beta /(alpha - 1);
    double gammaRMS = beta /TMath::Sqrt((alpha-2)*(alpha-3));
    double gammaYielde =TMath::Sqrt(pow(beta *  constante/(alpha-1),2) + pow(constant*beta*alphae/pow(alpha-1,2),2) + pow(constant * betae/(alpha-1),2));
    double gammaRMSe = TMath::Sqrt(pow((5-2*alpha)*beta*alphae/pow(2*((alpha-2)*(alpha-3)),1.5),2) + pow(betae/TMath::Sqrt((alpha-2)*(alpha-3)),2));

    /*cout << "Gauss Constant: " << gaussfit3->GetParameter(0) << ", center: " << gaussfit3->GetParameter(1) << ", RMS: " << gaussfit3->GetParameter(2) << endl;
      cout << "RMS: " << gaussigma << endl;
      cout << "Yield :" << gausyield << endl;
      cout << "Inverse Gamma Constant: " < <gaussfit3->GetParameter(3) << ", scale(beta): " << beta << ", Shape(alpha): " << alpha << endl;
      cout << "RMS: " << TMath::Sqrt(beta*beta/((alpha-1)*(alpha-1)*(alpha-2))) << endl;
      cout << "Yield: " << gammaYield << endl;*/
    fitB2[iF][ij] = B2;
    fitB1[iF][ij] = B1;
    fitB3[iF][ij] = constant;
    fitB4[iF][ij] = alpha;
    fitB5[iF][ij] = beta;
    fitPeak[iF][ij] = peak;
    fitB2e[iF][ij] = B2e;
    fitB1e[iF][ij] = B1e;
    fitB3e[iF][ij] = constante;
    fitB4e[iF][ij] = alphae;
    fitB5e[iF][ij] = betae;
    fitPeake[iF][ij] = peake;

    gausRMS[iF][ij] = gaussigma;
    gausRMSe[iF][ij] = gaussigmae;
    gausYield[iF][ij] = gausyield;
    gausYielde[iF][ij] = gausyielde;
    invgRMS[iF][ij] = gammaRMS;
    invgRMSe[iF][ij] = gammaRMSe;
    invgYield[iF][ij] = gammaYield;
    invgYielde[iF][ij] = gammaYielde;
    return gaussfit3;
  }

  void scanFit2(TH1D *histo, double cut, int ij, int iF){
    TRandom *r3 = new TRandom();
    double Chi2 = 10;
    int counter = 0;
    cout << counter << endl;
    counter++;
    double b1low = r3->Uniform(B1low*0.5,B1low*2);
    double b2low = r3->Uniform(B2low*0.5,B2low*2);
    double b3low = r3->Uniform(B3low*0.5,B3low*2);
    double b4low = r3->Uniform(B4low*0.5,B4low*2);
    double b5low = r3->Uniform(B5low*0.5,B5low*2);
    double b1high = r3->Uniform(B1high*0.5,B1high*2);
    double b2high = r3->Uniform(B2high*0.5,B2high*2);
    double b3high;
    if(ij < B3cut){
      b3high = r3->Uniform(B3high1*0.5,B3high1*2);
    }else{
      b3high = r3->Uniform(B3high1*0.5,B3high2*2);
    }
    double b4high = r3->Uniform(B4high*0.5,B4high*2);
    double b5high = r3->Uniform(B5high*0.5,B5high*2);
    double cut_ = r3->Uniform(cut*0.5,cut*2);

    TF1 *gaussfit = new TF1("gaussfit","gausn",0,10);
    gaussfit->FixParameter(1,0);
    gaussfit->SetParLimits(0,b2low,b2high);
    gaussfit->SetParameter(0,B2start[ij]);
    gaussfit->SetParLimits(2,b1low,b1high);
    gaussfit->SetParameter(2,B1start[ij]);
    histo->Fit("gaussfit","QN");
    TF1 *invG = new TF1("invG","[0]* (pow([1],[2])/TMath::Gamma([2]))*exp(-[1]/x) * pow(x, -[2]-1)",0,10);
    invG->SetParameter(0,B3start[ij]);
    if(ij < B3cut){
      invG->SetParLimits(0,b3low,b3high);
    }else{
      invG->SetParLimits(0,b3low,b3high);
    }
    invG->SetParameter(1,B5start[ij]);
    invG->SetParLimits(1,b5low,b5high);
    invG->SetParameter(2,B4start[ij]);
    invG->SetParLimits(2,b4low,b4high);
    histo->Fit("invG","QN","",cut,3);
    TF1 *gaussfit3= new TF1("gaussfit3","gausn(0) + [3]* (pow([4],[5])/TMath::Gamma([5]))*exp(-[4]/x) * pow(x, -[5]-1)",0,10);

    gaussfit3->SetParameter(0, gaussfit->GetParameter(0));
    gaussfit3->SetParameter(2, gaussfit->GetParameter(2));
    gaussfit3->SetParameter(3, invG->GetParameter(0));
    gaussfit3->SetParameter(4, invG->GetParameter(1));
    gaussfit3->SetParameter(5, invG->GetParameter(2));
    gaussfit3->FixParameter(1,0);

    gaussfit3->SetParLimits(0,b2low,b2high); //B2
    gaussfit3->SetParLimits(2,b1low,b1high); //B1
    gaussfit3->SetParLimits(3,b3low,b3high); //B3
    gaussfit3->SetParLimits(4,b5low,b5high); //B5
    gaussfit3->SetParLimits(5,b4low,b4high); //B4

    int lastBin = histo->FindLastBinAbove(1e-5);
    double end = histo->GetBinCenter(lastBin);
    int nBins  = lastBin - findStart(histo,0.1);
    histo->Fit("gaussfit3","QN","",0.1,end);
    double chi2_ = getChi2(histo,gaussfit3,0.1,end,0);
    double chi2dof_ = chi2_/(nBins-5);
    double chi2n_ = chi2_/nBins;
    Chi2 = chi2dof_;
    cout << Form("%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f",b1low,b1high,b2low,b2high,b3low,b3high,b4low,b4high,b5low,b5high,cut,Chi2) << endl;
  }

  TH2D *scanFit(TH1D *histo, double cut, int ij, int iF,int is,TString proj,bool draw){
    int Ndiv = 10;
    TH3D *hChi2;
    TH1D *hnarrowRMS = new TH1D("hnarrowRMS","Narrow component RMS distribution",100,0,2);
    TH1D *hwideRMS= new TH1D("hwideRMS","Wide component RMS distribution",100,0,2);
    TH1D *histChi2 = new TH1D("histChi2","Chi2 distribution",100,0,10);
    if(ij < B3cut){
      hChi2 = new TH3D("hChi2","Chi2 from parameters",Ndiv,B5low*0.99,B5high*1.01,Ndiv,B4low*0.99,B4high*1.01,Ndiv,B3low*0.99,B3high1*1.01);
    }else{
      hChi2 = new TH3D("hChi2","Chi2 from parameters",Ndiv,B5low*0.99,B5high*1.01,Ndiv,B4low*0.99,B4high*1.01,Ndiv,B3low*0.99,B3high2*1.01);
    }
    hChi2->GetXaxis()->SetTitle("B5");
    hChi2->GetYaxis()->SetTitle("B4");
    hChi2->GetZaxis()->SetTitle("B3");
    double par1, par2,par3,par4,par5;

    for(int i1 = 0; i1 < Ndiv; i1++){
      par5 = B5low + i1*(B5high-B5low)/(Ndiv-1);
      for(int i2 = 0; i2 < Ndiv; i2++){
        par4 = B4low + i2*(B4high-B4low)/(Ndiv-1);
        for(int i3 = 0; i3 < Ndiv; i3++){
          if(ij < B3cut){
            par3 = B3low + i3*(B3high1-B3low)/(Ndiv-1);
          }else{
            par3 = B3low + i3*(B3high2-B3low)/(Ndiv-1);
          }
          for(int i4 = 0; i4 < Ndiv; i4++){
            par1 = B1low + i4*(B1high-B1low);
            for(int i5 = 0; i5 < Ndiv; i5++){
              par2 = B2low + i5*(B2high-B2low);
              TF1 *gaussfit = new TF1("gaussfit","gausn",0,10);
              gaussfit->FixParameter(1,0);
              gaussfit->SetParameter(0,par2);
              gaussfit->SetParLimits(0,B2low,B2high);
              gaussfit->SetParameter(2,par1);
              gaussfit->SetParLimits(2,B1low,B1high);
              histo->Fit("gaussfit","QN");
              TF1 *invG = new TF1("invG","[0]* (pow([1],[2])/TMath::Gamma([2]))*exp(-[1]/x) * pow(x, -[2]-1)",0,10);
              invG->SetParameter(0,par3);
              if(ij < B3cut){
                //cout << "Setting parameter B3 to " << B3start1 << " Limits: " << B3low << " to " << B3high1 << endl;;
                invG->SetParLimits(0,B3low,B3high1);
              }else{
                //cout << "Setting parameter B3 to " << B3start2 << " Limits: " << B3low << " to " << B3high2 << endl;;
                invG->SetParLimits(0,B3low,B3high2);
              }
              //cout << "Setting parameter B5 to " << par1 << " Limits: " << B5low << " to " << B5high << endl;;
              invG->SetParameter(1,par5);
              invG->SetParLimits(1,B5low,B5high);
              //cout << "Setting parameter B4 to " << par2 << " Limits: " << B4low << " to " << B4high << endl;;
              invG->SetParameter(2,par4);
              invG->SetParLimits(2,B4low,B4high);
              histo->Fit("invG","QN","",cut,3);
              //invG->Print();
              TF1 *gaussfit3= new TF1("gaussfit3","gausn(0) + [3]* (pow([4],[5])/TMath::Gamma([5]))*exp(-[4]/x) * pow(x, -[5]-1)",0,10);
              if(gaussfit->GetParameter(2) < B1low || gaussfit->GetParameter(2) > B1high) cout << "B1 set to " << gaussfit->GetParameter(2) << " Limits: " << B1low << " to " << B1high << endl; 
              if(gaussfit->GetParameter(0) < B2low || gaussfit->GetParameter(0) > B2high) cout << "B2 set to " << gaussfit->GetParameter(0) << " Limits: " << B2low << " to " << B2high << endl; 
              //if(gaussfit->GetParameter(2) < B1low) cout << "Lower" << endl;
              //if(gaussfit->GetParameter(2) > B1high) cout << "Higher" << endl;
              if(invG->GetParameter(0) < B3low || invG->GetParameter(0) > B3high2) cout << "B3 ERROR" << endl; 
              if(invG->GetParameter(1) < B5low || invG->GetParameter(1) > B5high) cout << "B5 ERROR" << endl; 
              if(invG->GetParameter(2) < B4low || invG->GetParameter(2) > B4high) cout << "B4 ERROR" << endl; 

              gaussfit3->SetParameter(0, gaussfit->GetParameter(0));
              //gaussfit3->SetParameter(1, gaussfit->GetParameter(1));
              gaussfit3->SetParameter(2, gaussfit->GetParameter(2));
              gaussfit3->SetParameter(3, invG->GetParameter(0));
              gaussfit3->SetParameter(4, invG->GetParameter(1));
              gaussfit3->SetParameter(5, invG->GetParameter(2));
              gaussfit3->FixParameter(1,0);

              gaussfit3->SetParLimits(0,B2low,B2high); //B2
              gaussfit3->SetParLimits(2,B1low,B1high); //B1
              if(ij < B3cut){
                gaussfit3->SetParLimits(3,B3low,B3high1); //B3
              }else{
                gaussfit3->SetParLimits(3,B3low,B3high2); //B3
              }
              gaussfit3->SetParLimits(4,B5low,B5high); //B5
              gaussfit3->SetParLimits(5,B4low,B4high); //B4

              int lastBin = histo->FindLastBinAbove(1e-5);
              double end = histo->GetBinCenter(lastBin);
              int nBins  = lastBin - findStart(histo,0.1);
              //cout << "Debug 0";
              histo->Fit("gaussfit3","QN","",0.1,end);
              //cout << "Debug 1";
              //chi2[iF][ij] = gaussfit3->GetChisquare();
              double chi2_ = getChi2(histo,gaussfit3,0.1,end,0);
              double chi2dof_ = chi2_/(nBins-5);
              double chi2n_ = chi2_/nBins;
              histChi2->Fill(chi2dof_);
              hChi2->Fill(par5,par4,par3,chi2dof_);
              hnarrowRMS->Fill(TMath::Sqrt(2)*gaussfit3->GetParameter(2));
              double alpha = gaussfit3->GetParameter(5); // B4
              double beta = gaussfit3->GetParameter(4); // B5
              hwideRMS->Fill(beta / TMath::Sqrt((alpha-2)*(alpha-3)));
              gaussfit3->Delete();
              invG->Delete();
              gaussfit->Delete();
            }//i5
          }//i4
        } //i3
      } //i2
    } //i1
    hChi2->Scale(1.0/10);
    if(draw){
      //drawFitScan((TH2D*)hChi2->Project3D(proj),"Title",comment[iF],"hChi2Scans",iF,ij,proj);
      //drawRMSHisto(hnarrowRMS,"Wide component RMS","NarrowRMSHist",is,ij,comment[iF]);
      //drawRMSHisto(hwideRMS,"Narrow component RMS","WideRMSHist",is,ij,comment[iF]);
      //drawRMSHisto(histChi2,"Chi2","Chi2Hist",is,ij,comment[iF]);
    }
    hnarrowRMS->Delete();
    hwideRMS->Delete();
    return (TH2D*)hChi2->Project3D(proj);
  }

  void drawRMSHisto(TH1D *histo, TString title, TString file,int is,int ij,TString comment){
    TCanvas *can = new TCanvas(Form("c%02d",iC++), Form("c%02d",iC++), 150, 150);
    histo->SetMarkerStyle(qMarker[0]);
    histo->SetMarkerColor(qColor[0]);
    histo->SetMarkerSize(mSize);
    TLegend *leg = new TLegend(0.45,0.45,0.92,0.80, "","brNDC");
    leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.027);
    leg->AddEntry((TObject*)0, title,"");
    histo->Draw();
    leg->Draw();
    gPad->Modified();
    gPad->GetCanvas()->Update();
    TString figname = Form("figs/%sNFin%02dJetPt%02d%s.pdf",file.Data(),is,ij,comment.Data());
    cout << "Save " << figname << endl;
    if(saveFigs) gPad->GetCanvas()->SaveAs(figname);
  }

  TGraphErrors *makeErrorGraph(TH1D *histo, TH1D *errors){
    int n = histo->GetNbinsX();
    double center[150];
    double width[150];
    double y[150];
    double err[150];
    for(int i = 0; i < n; i++){
      center[i] = histo->GetBinCenter(i+1);
      width[i] = histo->GetBinWidth(i+1) / 2;
      y[i] = histo->GetBinContent(i+1);
      err[i] = errors->GetBinContent(i+1);
    }
    TGraphErrors *errGraph = new TGraphErrors(n, center, y, width, err);
    return errGraph;
  }


  TGraphErrors *makeErrorGraph(TGraphErrors *histo, TGraphErrors *errors){
    double xin[2], yin[2], exin[2], eyin[2];
    double x[100], y[100], ex[100]={0}, ey[100];
    int N = histo->GetN();
    for(int i = 0 ; i < N; i++){
      histo->GetPoint(i, xin[0],yin[0]);exin[0] = histo->GetErrorX(i); eyin[0] = histo->GetErrorY(i);
      errors->GetPoint(i, xin[1],yin[1]);exin[1] = errors->GetErrorX(i); eyin[1] = errors->GetErrorY(i);
      x[i] = xin[0];
      y[i] = yin[0];
      ex[i] = exin[0];
      ey[i] = yin[1];
    }
    TGraphErrors *gr = new TGraphErrors(N, x, y, ex, ey);
    return gr;
  }

  TH1D *makeSystError(TH1D *histo1, TH1D *histo2){
    int N = histo1->GetNbinsX();
    TH1D *errors = (TH1D*)histo1->Clone();
    errors->Clear();
    if(histo2->GetNbinsX() != N) return NULL;
    double x1,x2,x3;
    for(int i = 1 ; i <= N; i++){
      x1 = TMath::Max(0.0,histo1->GetBinContent(i));
      x2 = TMath::Max(0.0,histo2->GetBinContent(i));
      x3 = TMath::Abs(x1-x2);
      errors->SetBinContent(i,x3);
    }
    return errors;
  }

  TGraphErrors* smoothSystError(TGraphErrors *histo, bool draw, TString title){
    double xin,yin,exin,eyin;
    double x[100],y[100],ex[100]={0},ey[100];
    double maxY = 0;
    double maxX = 0;
    //TF1 *fit = new TF1("fit","[0]+[1]*x+[2]*x*x+[3]*x*x*x",5,500);
    //TF1 *fit = new TF1("fit","[0]+[1]*x+[2]*x*x",20,500);
    TF1 *fit = new TF1("fit","[0]+[1]*x",20,500);
    histo->Fit(fit, "QN","",20,200);
    int N = histo->GetN();
    for(int i = 0; i < N; i++){
      histo->GetPoint(i,xin,yin);exin = histo->GetErrorX(i); eyin = histo->GetErrorY(i);
      x[i] = xin;
      if(xin > maxX) maxX = xin;
      if(yin > maxY) maxY = yin;
      y[i] = fit->Eval(xin); 
      ey[i] = eyin;
      ex[i] = exin;
      if(y[i] > maxY) maxY = y[i];
    }
    TGraphErrors *gr = new TGraphErrors(N,x,y,ex,ey);
    if(draw){
      Filipad2 *fpad = new Filipad2(iC++, 1.1, 0.4, 150, 50, 0.7, 6,1);
      fpad->Draw();
      TPad *mpad = fpad->GetPad(1); //upper pad

      gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
      gStyle->SetMarkerSize(1.6);
      mpad->SetLogx(0);mpad->SetGridx(0);
      mpad->SetLogy(0);mpad->SetGridy(0);
      mpad->cd();
      gPad->SetLeftMargin(0.17);
      TLegend *leg;
      leg = new TLegend(0.45,0.45,0.92,0.80, "","brNDC");
      leg->AddEntry((TObject*)0, title,"");
      leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.027);

      TH2F *hfrSignal = new TH2F("hfr3"," ", 10, 0, maxX*1.5,10, 0, maxY*1.5);
      hset( *hfrSignal, "p_{T} [GeV]","Error");
      hfrSignal->Draw();

      histo->SetMarkerColor(qColor[0]);
      histo->SetLineColor(qColor[0]);
      histo->SetMarkerStyle(qMarker[0]);
      histo->SetMarkerSize(mSize);
      gr->SetMarkerColor(qColor[1]);
      gr->SetMarkerStyle(qMarker[1]);
      gr->SetLineColor(qColor[1]);
      gr->SetMarkerSize(mSize);

      fit->SetLineColor(qColor[2]);
      gr->Draw("PZ same");
      fit->Draw("l same");
      histo->Draw("PZ same");
      leg->AddEntry(histo,"input","p");
      leg->AddEntry(gr,"output","p");
      leg->AddEntry(fit,"Fit a + bx","l");
      leg->AddEntry((TObject*)0, Form("a: %.4f",fit->GetParameter(0)),"");
      leg->AddEntry((TObject*)0, Form("b: %.4f",fit->GetParameter(1)),"");
      leg->Draw();
      mpad = fpad->GetPad(2); //lower pad
      //gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
      gStyle->SetMarkerSize(1.6);
      mpad->SetLogx(0);mpad->SetGridx(0);
      mpad->SetLogy(0);mpad->SetGridy(0);
      mpad->cd();
      gPad->SetLeftMargin(0.17);

      TH2F *hfrSignal2 = new TH2F("hfr4", "", 10, 0,maxX*1.5,10,0.1,2);
      hset( *hfrSignal2, "p_{T}",Form("Ratio Error/Fit"));
      hfrSignal2->Draw();
      //tempRatio->Divide(fit);
      //tempRatio->Draw("PZ same");
      //gPad->Modified();
      gPad->Modified();
      gPad->GetCanvas()->Update();
    }
    return gr;
  }


  TGraphErrors *makeSystError(TGraphErrors *histo1, TGraphErrors *histo2){
    return makeSystError(histo1, histo2,true);
  }

  TGraphErrors *makeSystError(TGraphErrors *histo1, TGraphErrors *histo2,bool abs){
    double xin[2], yin[2], exin[2], eyin[2];
    double x[100], y[100], ex[100]={0}, ey[100];
    int N = histo1->GetN();
    for(int i = 0 ; i < N; i++){
      histo1->GetPoint(i, xin[0],yin[0]);exin[0] = histo1->GetErrorX(i); eyin[0] = histo1->GetErrorY(i);
      histo2->GetPoint(i, xin[1],yin[1]);exin[1] = histo2->GetErrorX(i); eyin[1] = histo2->GetErrorY(i);
      x[i] = xin[0];
      if(abs){
        y[i] = TMath::Abs(yin[1]-yin[0]);
      }else{
        y[i] = yin[1]-yin[0];
      }
      ex[i] = exin[0];
      ey[i] = TMath::Sqrt(pow(eyin[0],2)+pow(eyin[1],2));
    }
    TGraphErrors *gr = new TGraphErrors(N, x, y, ex, ey);
    return gr;
  }

  TGraphErrors *combineErrors(TGraphErrors *histo1, TGraphErrors *histo2){
    double xin[2], yin[2], exin[2], eyin[2];
    double x[100], y[100], ex[100]={0}, ey[100];
    int N = histo1->GetN();
    for(int i = 0 ; i < N; i++){
      histo1->GetPoint(i, xin[0],yin[0]);exin[0] = histo1->GetErrorX(i); eyin[0] = histo1->GetErrorY(i);
      histo2->GetPoint(i, xin[1],yin[1]);exin[1] = histo2->GetErrorX(i); eyin[1] = histo2->GetErrorY(i);
      x[i] = xin[0];
      y[i] = TMath::Sqrt(pow(yin[1],2)+pow(yin[0],2));
      ex[i] = exin[0];
      ey[i] = TMath::Sqrt(pow(eyin[0],2)+pow(eyin[1],2));
    }
    TGraphErrors *gr = new TGraphErrors(N, x, y, ex, ey);
    return gr;
  }

  TGraphErrors *combineErrors(TGraphErrors *histo1, TGraphErrors *histo2,TGraphErrors *histo3){
    double xin[3], yin[3], exin[3], eyin[3];
    double x[100], y[100], ex[100]={0}, ey[100];
    int N = histo1->GetN();
    for(int i = 0 ; i < N; i++){
      histo1->GetPoint(i, xin[0],yin[0]);exin[0] = histo1->GetErrorX(i); eyin[0] = histo1->GetErrorY(i);
      histo2->GetPoint(i, xin[1],yin[1]);exin[1] = histo2->GetErrorX(i); eyin[1] = histo2->GetErrorY(i);
      histo3->GetPoint(i, xin[2],yin[2]);exin[2] = histo3->GetErrorX(i); eyin[2] = histo3->GetErrorY(i);
      x[i] = xin[0];
      y[i] = TMath::Sqrt(pow(yin[1],2)+pow(yin[0],2)+pow(yin[2],2));
      ex[i] = exin[0];
      ey[i] = TMath::Sqrt(pow(eyin[0],2)+pow(eyin[1],2)+pow(eyin[2],2));
    }
    TGraphErrors *gr = new TGraphErrors(N, x, y, ex, ey);
    return gr;
  }

  TH1D *combineErrors(TH1D *errors1, TH1D *errors2){
    int n = errors1->GetNbinsX();
    TH1D *errors = (TH1D*)errors1->Clone();
    errors->Reset();
    double y1,y2,y;
    for(int i = 1 ; i < n+1; i++){
      y1 = errors1->GetBinContent(i);
      y2 = errors2->GetBinContent(i);
      y = TMath::Sqrt(y1*y1+y2*y2);
      errors->SetBinContent(i,y);
    }
    return errors;
  }

  TH1D *combineErrors(TH1D *errors1, TH1D *errors2,TH1D *errors3){
    int n = errors1->GetNbinsX();
    TH1D *errors = (TH1D*)errors1->Clone();
    errors->Reset();
    double y1,y2,y3,y;
    for(int i = 1 ; i < n+1; i++){
      y1 = errors1->GetBinContent(i);
      y2 = errors2->GetBinContent(i);
      y3 = errors3->GetBinContent(i);
      y = TMath::Sqrt(y1*y1+y2*y2+y3*y3);
      errors->SetBinContent(i,y);
    }
    return errors;
  }

  int findStart(TH1D *histo, double start){
    int N = histo->GetNbinsX();
    double center;
    for(int i = 1; i < N; i++){
      center = histo->GetBinCenter(i);
      if(center > start){
        return i;
      }
    }
    return -5;
  }
  double getChi2(TH1D *histo, TF1 *fit,double lowX, double highX, int debug){
    if(highX < lowX) cout << "WARNING " << endl << endl << endl << endl;
    int N = histo->GetNbinsX();
    double center;
    double content;
    //double error;
    double sum = 0;
    double expected;
    for(int i = 1; i < N+1; i++){
      center = histo->GetBinCenter(i);
      if(debug > 0){
        cout << "i: " << i << "/" << (N+1) << endl;
        cout << "center: " << center << ", lowX: " << lowX << ", highX: " << highX << endl;
      }
      if(center < lowX || center > highX) continue;
      content = histo->GetBinContent(i);
      expected = fit->Eval(center);
      sum += pow(content-expected,2)/expected;
      if(debug > 0){
        cout << "content: " << content << ", Expected: " << expected << ", Sum: " << sum << endl;
      }
      //error = histo->GetBinError(i);
    }
    return sum;
  }
  TGraphErrors *makeSignalGraph(TH1D *histo){
    int n = histo->GetNbinsX();
    double center[150];
    double width[150];
    double y[150];
    double err[150];
    for(int i = 0; i< n; i++){
      center[i] = histo->GetBinCenter(i+1);
      width[i] = histo->GetBinWidth(i+1) / 2;
      y[i] = histo->GetBinContent(i+1);
      err[i] = histo->GetBinError(i+1);
    }
    TGraphErrors *errGraph = new TGraphErrors(n,center,y,width,err);
    return errGraph;
  }

  void drawFitScan(TH2D *histo,TString title, TString comment, TString file, int is, int ij,TString proj){
    TCanvas *can = new TCanvas(Form("c%02d",iC++), Form("c%02d",iC++), 150, 150);

    //gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    //gStyle->SetMarkerSize(1.6);
    //mpad->SetLogx(0);mpad->SetGridx(0);
    //mpad->SetLogy(0);mpad->SetGridy(0);
    //gPad->SetLeftMargin(0.17);

    //leg = new TLegend(0.19,0.15,0.72,0.45, "","brNDC");
    //leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.045);

    //histo->GetXaxis()->SetTitle("B5");
    //histo->GetYaxis()->SetTitle("B4");
    histo->Draw("colz");
    TString figname;
    figname = Form("figs/FitScan%sNFin%02dJetPt%02d%s.pdf",proj.Data(),is,ij,comment.Data());
    if(saveFigs) gPad->GetCanvas()->SaveAs(figname);
  }
  void drawFit(TH1D *histo, TF1 *fit,TString method, double xlow, double xhigh, int xlog, double ylow, double yhigh, int ylog, TString title,TString title2,TString comment,TString file, int is, int ij, int iF){
    Filipad2 *fpadSignal = new Filipad2(iC++, 1.1, 0.4, 150, 50, 0.7, 6,1);
    fpadSignal->Draw();
    TPad *mpad = fpadSignal->GetPad(1); //upper pad

    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(ylog);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.22);
    TLegend *leg;
    if(xlog){
      //TLegend *leg(0.15,0.10,0.62,0.55, "","brNDC");
      leg = new TLegend(0.57,0.55,0.94,0.83, "","brNDC");
    }else{
      leg = new TLegend(0.45,0.45,0.92,0.80, "","brNDC");
    }
    leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.027);

    TH2F *hfrSignal = new TH2F("hfr3"," ", 10, xlow , xhigh,10, ylow, yhigh);
    if(doWeight){
      hset( *hfrSignal, "j_{T} [GeV]","#frac{1}{N_{jets}} #frac{dN}{j_{T}dj_{T}}");
    }else{
      hset( *hfrSignal, "j_{T} [GeV]","#frac{dN}{dj_{T}}");
    }
    hfrSignal->GetYaxis()->SetTitleOffset(1.6);
    //    0.9, 0.9, 0.1,0.09, 0.01,0.001, 0.05,0.05, 505, 505);
    //hset2( iCanvas, 5, 500, 0.5, 1.2 , "p_{T}",
    //    "Ratio to detector level", 10, 00);
    hfrSignal->Draw();

    TF1 *fitGaus = new TF1("gausfit","gausn",0,10);
    TF1 *fitInvG;
    if(method == "alt"){
      fitInvG = new TF1("invGfit","[0]* (pow([1]*([2]+1),[2])/TMath::Gamma([2]))*exp(-[1]*([2]+1)/x) * pow(x, -[2]-1)",0,10);
    }else{
      fitInvG = new TF1("invGfit","[0]* (pow([1],[2])/TMath::Gamma([2]))*exp(-[1]/x) * pow(x, -[2]-1)",0,10);
    }
    fitGaus->SetParameter(0, fit->GetParameter(0));
    fitGaus->SetParameter(1, fit->GetParameter(1));
    fitGaus->SetParameter(2, fit->GetParameter(2));
    fitInvG->SetParameter(0, fit->GetParameter(3));
    fitInvG->SetParameter(1, fit->GetParameter(4));
    fitInvG->SetParameter(2, fit->GetParameter(5));
    fitGaus->SetLineColor(3);
    fitGaus->SetLineStyle(2);
    fitInvG->SetLineColor(1);
    fitInvG->SetLineStyle(2);
    if(title2.Length() > 1){
      leg->AddEntry((TObject*)0, title2, "");
    }else{
      leg->AddEntry((TObject*)0, "p-Pb #sqrt{s_{NN}} = 5.02 TeV","");
    }
    leg->AddEntry((TObject*)0, finderName[is],"");
    if(ij == Njets){
      leg->AddEntry((TObject*)0, Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[0], JetPtBins2[Njets]),"");
    }else{
      leg->AddEntry((TObject*)0, Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[ij], JetPtBins2[ij+1]),"");
    }
    //leg->AddEntry((TObject*)0, title, "");
    histo->SetMarkerColor(1);
    histo->SetLineColor(1);
    histo->SetMarkerStyle(25);
    histo->SetMarkerSize(mSize);
    for(int ib = 1; ib < histo->GetNbinsX() ; ib++){
      if(histo->GetBinCenter(ib) < 0.1){
        histo->SetBinContent(ib,0);
      }
    }
    TH1D *tempRatio = (TH1D*)histo->Clone();
    histo->Draw("PZ same");
    leg->AddEntry(histo,"jT signal","p");
    fit->SetLineColor(2);
    fitGaus->Draw("l same");
    fitInvG->Draw("l same");
    fit->Draw("l same");
    leg->AddEntry(fit,"Fit","l");
    leg->AddEntry(fitGaus,"Gaussian component","l");
    leg->AddEntry(fitInvG,"Inverse Gamma component","l");
    /*leg->AddEntry((TObject*)0, Form("Chi^{2}: %.2f", chi2dof[iF][ij]),"");
    leg->AddEntry((TObject*)0, Form("B_{1}: %.2f", fitB1[iF][ij]),"");
    leg->AddEntry((TObject*)0, Form("B_{2}: %.2f", fitB2[iF][ij]),"");
    leg->AddEntry((TObject*)0, Form("B_{3}: %.2f", fitB3[iF][ij]),"");
    leg->AddEntry((TObject*)0, Form("B_{4}: %.2f", fitB4[iF][ij]),"");
    leg->AddEntry((TObject*)0, Form("B_{5}: %.2f", fitB5[iF][ij]),"");
    leg->AddEntry((TObject*)0, Form("Peak: %.2f", fitPeak[iF][ij]),"");
    */


    leg->Draw();
    mpad = fpadSignal->GetPad(2); //lower pad
    //gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(0);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.22);

    TH2F *hfrSignal2 = new TH2F("hfr4", "", 10, xlow,xhigh,10,0.1,2);
    hset( *hfrSignal2, "j_{T}",Form("#frac{Signal}{Fit}"));
    //hset2(*hfrSignal2, "pT","Ratio");
    //hset2( iCanvas, 5, 500, 0.1, 20.5 , "p_{T}", 
    //    "Ratio to detector level", 10, 00);
    //hset( *hfrSignal2, "p_{T} [GeV/c] ", "Ratio to particle level",
    //    0.9, 0.9, 0.1,0.09, 0.01,0.001, 0.05,0.05, 505, 505);
    hfrSignal2->Draw();
    tempRatio->Divide(fit);
    tempRatio->Draw("PZ same");
    //gPad->Modified();
    gPad->Modified();
    gPad->GetCanvas()->Update();
    //cout << "ppdf: " << file.Data() << endl;
    //ppdf(Form("figs/%sNFin%02dJetPt%02d%s",file.Data(),is,ij,topcomment.Data()),gPad->GetCanvas());
    TString figname;
    if(xlog){
      figname = Form("figs/%sNFin%02dJetPt%02d%s.pdf",file.Data(),is,ij,comment.Data());
    }else{
      figname = Form("figs/%sNFin%02dJetPt%02d_linx%s.pdf",file.Data(),is,ij,comment.Data());
    }
    if(saveFigs) gPad->GetCanvas()->SaveAs(figname);
  }
  void drawComparison(TH1D **histo, double xlow, double xhigh,int xlog, double ylow, double yhigh, int ylog,TString title, TString file, int is, int ij){
    //cout << "DrawComparison" << endl;
    Filipad2 *fpadSignal = new Filipad2(iC++, 1.1, 0.4, 150, 50, 0.7, 6,1);
    fpadSignal->Draw();
    TPad *mpad = fpadSignal->GetPad(1); //upper pad

    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(ylog);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.17);

    TLegend *leg = new TLegend(0.15,0.12,0.62,0.35, "","brNDC");
    leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.045);

    TH2F *hfrSignal = new TH2F("hfr3"," ", 10, xlow , xhigh,10, ylow, yhigh);
    if(doWeight){
      hset( *hfrSignal, "j_{T} [GeV]","#frac{dN}{j_{T}dj_{T}}");
      //hset( *hfrSignal, "p_{T} [GeV/c] ", "#frac{d#sigma}{dp_{T}d#eta}",
    }else{
      hset( *hfrSignal, "j_{T} [GeV]","#frac{dN}{dj_{T}}");
    }
    //    0.9, 0.9, 0.1,0.09, 0.01,0.001, 0.05,0.05, 505, 505);
    //hset2( iCanvas, 5, 500, 0.5, 1.2 , "p_{T}", 
    //    "Ratio to detector level", 10, 00);
    hfrSignal->Draw();

    //leg->AddEntry((TObject*)0, finderName[is],"");
    leg->AddEntry((TObject*)0, Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[ij], JetPtBins2[ij+1]),"");
    leg->AddEntry((TObject*)0, title, "");
    for(int iF = 0; iF < Nsets; iF++){
      //cout << "iF: " << iF << "ij: " << ij << endl;
      //cout << "\t";
      //histo[iF]->Print();
      histo[iF]->Draw("PZ same");
      histo[iF]->SetMarkerColor(iF+1);
      histo[iF]->SetLineColor(iF+1);
      histo[iF]->SetMarkerStyle(25);
      histo[iF]->SetMarkerSize(mSize);
      leg->AddEntry(histo[iF],setTitle[iF].Data(),"p");
    }

    leg->Draw();
    mpad = fpadSignal->GetPad(2); //lower pad
    //gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(0);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.17);

    TH2F *hfrSignal2 = new TH2F("hfr4", "", 10, xlow,xhigh,10,0.1,2);
    //hset( *hfrSignal2, "j_{T}",Form("Ratio to %s",setTitle[iF].Data()));
    hset( *hfrSignal2, "j_{T}","Ratio");
    //hset2(*hfrSignal2, "pT","Ratio");
    //hset2( iCanvas, 5, 500, 0.1, 20.5 , "p_{T}", 
    //    "Ratio to detector level", 10, 00);
    //hset( *hfrSignal2, "p_{T} [GeV/c] ", "Ratio to particle level",
    //    0.9, 0.9, 0.1,0.09, 0.01,0.001, 0.05,0.05, 505, 505);
    hfrSignal2->Draw();
    TH1D *tempRatio[Nsets];
    for(int iF = 0; iF < Nsets; iF++){
      //cout << "iF: " << endl;
      tempRatio[iF] = (TH1D*)histo[iF]->Clone();
      tempRatio[iF]->Divide(histo[0]);
      tempRatio[iF]->Draw("PZ same");
    }
    gPad->Modified();
    gPad->GetCanvas()->Update();
    TString figname;
    if(xlog){
      figname = Form("figs/%sNFin%02dJetPt%02d%s.pdf",file.Data(),is,ij,topcomment.Data());
    }else{
      figname = Form("figs/%sNFin%02dJetPt%02d_linx%s.pdf",file.Data(),is,ij,topcomment.Data());
    }
    if(saveFigs) gPad->GetCanvas()->SaveAs(figname);
  }

  void drawComparisonChi2(TGraph **histo, double xlow, double xhigh,int xlog, double ylow, double yhigh, int ylog,TString title, TString file){
    TCanvas *c = mc(iC++,1.1);
    //cout << "DrawComparison" << endl;

    //Filipad2 *fpadSignal = new Filipad2(iC++, 1.1, 0.4, 150, 50, 0.7, 6,1);
    //fpadSignal->Draw();
    //TPad *mpad = fpadSignal->GetPad(1); //upper pad
    //TPad *mpad = (TPad*)c->GetPadSave();
    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    gPad->SetLogx(xlog);gPad->SetGridx(0);
    gPad->SetLogy(ylog);gPad->SetGridy(0);
    gPad->cd();
    gPad->SetLeftMargin(0.17);

    TLegend *leg = new TLegend(0.15,0.65,0.62,0.85, "","brNDC");
    leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.045);

    TH2F *hfr = new TH2F("hfr3"," ", 10, xlow , xhigh,10, ylow, yhigh);
    hset( *hfr, "p_{T} [GeV]","Chi^{2}");

    hfr->Draw();
    leg->AddEntry((TObject*)0, title, "");
    for(int iF = 0; iF < Nsets; iF++){
      histo[iF]->Draw("PZ same");
      histo[iF]->SetMarkerColor(qColor[iF]);
      histo[iF]->SetLineColor(qColor[iF]);
      histo[iF]->SetMarkerStyle(qMarker[iF]);
      histo[iF]->SetMarkerSize(1);
      leg->AddEntry(histo[iF],setTitle[iF].Data(),"p");
    }

    leg->Draw();

    gPad->Modified();
    gPad->GetCanvas()->Update();
    TString figname;
    figname = Form("figs/%s%s.pdf",file.Data(),topcomment.Data());
    if(saveFigs) gPad->GetCanvas()->SaveAs(figname);
  }

  void drawErrors(TGraphErrors *histo, TGraphErrors *errBg, TGraphErrors *errUnf, TGraphErrors *errFit, TGraphErrors *errTot,double xlow, double xhigh, int xlog, double ylow, double yhigh, int ylog, TString title,TString title2,TString comment,TString file, int is, int ij){
    Filipad2 *fpadSignal = new Filipad2(iC++, 1.1, 0.5, 150, 50, 0.7, 6,1);
    fpadSignal->Draw();
    TPad *mpad = fpadSignal->GetPad(1); //upper pad

    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(ylog);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.17);
    TLegend *leg;
    if(xlog){
      leg = new TLegend(0.15,0.10,0.62,0.55, "","brNDC");
    }else{
      leg = new TLegend(0.45,0.45,0.92,0.80, "","brNDC");
    }
    leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.045);

    TH2F* hfrSignal = new TH2F("hfr3"," ", 10, xlow , xhigh,10, ylow, yhigh);
    hset( *hfrSignal, "j_{T} [GeV]","Absolute Error");
    //    0.9, 0.9, 0.1,0.09, 0.01,0.001, 0.05,0.05, 505, 505);
    //hset2( iCanvas, 5, 500, 0.5, 1.2 , "p_{T}", 
    //    "Ratio to detector level", 10, 00);
    hfrSignal->Draw();

    if(ij == Njets){
      leg->AddEntry((TObject*)0, Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[0], JetPtBins2[Njets]),"");
    }else{
      leg->AddEntry((TObject*)0, Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[ij], JetPtBins2[ij+1]),"");
    }
    leg->AddEntry((TObject*)0, title, "");
    leg->AddEntry((TObject*)0, title2, "");
    errBg->SetLineColor(qColor[1]);
    errUnf->SetLineColor(qColor[2]);
    errFit->SetLineColor(qColor[3]);
    errTot->SetLineColor(qColor[4]);
    errBg->SetMarkerColor(qColor[1]);
    errUnf->SetMarkerColor(qColor[2]);
    errFit->SetMarkerColor(qColor[3]);
    errTot->SetMarkerColor(qColor[4]);
    errBg->SetMarkerSize(mSize);
    errUnf->SetMarkerSize(mSize);
    errFit->SetMarkerSize(mSize);
    errTot->SetMarkerSize(mSize);
    errBg->SetMarkerStyle(qMarker[1]);
    errUnf->SetMarkerStyle(qMarker[2]);
    errFit->SetMarkerStyle(qMarker[3]);
    errTot->SetMarkerStyle(qMarker[3]);

    /*TF1 * bgFit = new TF1("bgFit","pol2",0,5);
      TF1 * unfFit = new TF1("unfFit","pol2",0,5);
      TF1 * totFit = new TF1("totFit","pol2",0,5);*/
    /*TF1 * bgFit = new TF1("bgFit","pol2",0.1,1);
      TF1 * unfFit = new TF1("unfFit","pol2",0.1,1);
      TF1 * FitFit = new TF1("FitFit","pol2",0,1,1);
      TF1 * totFit = new TF1("totFit","pol2",0.1,1);*/
    errBg->Draw("PZ same");
    /*bgFit->SetLineColor(qColor[1]);
      unfFit->SetLineColor(qColor[2]);
      totFit->SetLineColor(qColor[3]);
      bgFit->SetLineStyle(2);
      unfFit->SetLineStyle(2);
      totFit->SetLineStyle(2);*/

    errUnf->Draw("PZ same");
    errFit->Draw("PZ same");
    errTot->Draw("PZ same");

    /*errUnf->Fit("unfFit","QN","",0.1,5);
      errBg->Fit("bgFit","QN","",0.1,5);
      errTot->Fit("totFit","QN","",0.1,5);
      totFit->Draw("l same");


      TH1D *hint = new TH1D(Form("hint%d",1),
      "Fitted gaussian with .95 conf.band", 20, 0.1,5);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
      hint->SetStats(kFALSE);
      hint->SetFillColor(qColor[3]);
      hint->SetFillColorAlpha(qColor[3], 0.35);
      hint->Draw("e3 same");*/


    leg->AddEntry(errBg,"Background method","p");
    leg->AddEntry(errUnf,"Unfolding method","p");
    leg->AddEntry(errFit,"Fitting Cut", "p");
    leg->AddEntry(errTot,"Quadratic Sum","p");
    //leg->AddEntry(totFit,"Fit","l");
    leg->Draw();
    mpad = fpadSignal->GetPad(2); //lower pad
    //gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(0);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.17);

    TH2F *hfrSignal2 = new TH2F("hfr4", "", 10, xlow,xhigh,10,0,1.05);
    hset( *hfrSignal2, "j_{T}",Form("Relative error"));
    //hset2(*hfrSignal2, "pT","Ratio");
    //hset2( iCanvas, 5, 500, 0.1, 20.5 , "p_{T}", 
    //    "Ratio to detector level", 10, 00);
    //hset( *hfrSignal2, "p_{T} [GeV/c] ", "Ratio to particle level",
    //    0.9, 0.9, 0.1,0.09, 0.01,0.001, 0.05,0.05, 505, 505);
    hfrSignal2->Draw();
    TGraphErrors *tempRatioBg = grrDivide(errBg,histo);
    TGraphErrors *tempRatioUnf = grrDivide(errUnf,histo);
    TGraphErrors *tempRatioFit = grrDivide(errFit,histo);
    TGraphErrors *tempRatioTot = grrDivide(errTot,histo);
    tempRatioBg->SetMarkerStyle(qMarker[1]);
    tempRatioUnf->SetMarkerStyle(qMarker[2]);
    tempRatioFit->SetMarkerStyle(qMarker[3]);
    tempRatioTot->SetMarkerStyle(qMarker[4]);
    tempRatioBg->SetLineColor(qColor[1]);
    tempRatioUnf->SetLineColor(qColor[2]);
    tempRatioFit->SetLineColor(qColor[3]);
    tempRatioTot->SetLineColor(qColor[4]);
    tempRatioBg->SetMarkerColor(qColor[1]);
    tempRatioUnf->SetMarkerColor(qColor[2]);
    tempRatioFit->SetMarkerColor(qColor[3]);
    tempRatioTot->SetMarkerColor(qColor[4]);
    tempRatioBg->SetMarkerSize(mSize);
    tempRatioUnf->SetMarkerSize(mSize);
    tempRatioFit->SetMarkerSize(mSize);
    tempRatioTot->SetMarkerSize(mSize);
    /*tempRatioUnf->Fit("unfFit","QN","",0.1,1);
      tempRatioBg->Fit("bgFit","QN","",0.1,1);
      tempRatioTot->Fit("totFit","QN","",0.1,1);*/
    //bgFit->Draw("l same");
    //unfFit->Draw("l same");
    //totFit->Draw("l same");

    /*TH1D *hint = new TH1D(Form("hint%d",1),
      "Fitted gaussian with .95 conf.band", 20, 0.1,1);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
      hint->SetStats(kFALSE);
      hint->SetFillColor(qColor[3]);
      hint->SetFillColorAlpha(qColor[3], 0.35);
      hint->Draw("e3 same");*/

    tempRatioBg->Draw("PZ same");
    tempRatioUnf->Draw("PZ same");
    tempRatioFit->Draw("PZ same");
    tempRatioTot->Draw("PZ same");

    //tempRatio->Divide(fit);
    //tempRatio->Draw("PZ same");
    //gPad->Modified();
    gPad->Modified();
    gPad->GetCanvas()->Update();
    //cout << "ppdf: " << file.Data() << endl;
    //ppdf(Form("figs/%sNFin%02dJetPt%02d%s",file.Data(),is,ij,topcomment.Data()),gPad->GetCanvas());
    TString figname;
    if(xlog){
      figname = Form("figs/%sNFin%02dJetPt%02d%s.pdf",file.Data(),is,ij,comment.Data());
    }else{
      figname = Form("figs/%sNFin%02dJetPt%02d_linx%s.pdf",file.Data(),is,ij,comment.Data());
    }
    if(saveFigs) gPad->GetCanvas()->SaveAs(figname);
  }

  void drawErrors(TGraphErrors *histo, TGraphErrors *errBg,double xlow, double xhigh, int xlog, double ylow, double yhigh, int ylog, TString title,TString title2,TString comment,TString file, int is, int ij){
    Filipad2 *fpadSignal = new Filipad2(iC++, 1.1, 0.5, 150, 50, 0.7, 6,1);
    fpadSignal->Draw();
    TPad *mpad = fpadSignal->GetPad(1); //upper pad

    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(ylog);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.17);

    TLegend *leg;
    if(xlog){
      leg = new TLegend(0.15,0.10,0.62,0.55, "","brNDC");
    }else{
      leg = new TLegend(0.45,0.45,0.92,0.80, "","brNDC");
    }
    leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.045);

    TH2F *hfrSignal = new TH2F("hfr3"," ", 10, xlow , xhigh,10, ylow, yhigh);
    hset( *hfrSignal, "j_{T} [GeV]","Absolute Error");
    //    0.9, 0.9, 0.1,0.09, 0.01,0.001, 0.05,0.05, 505, 505);
    //hset2( iCanvas, 5, 500, 0.5, 1.2 , "p_{T}", 
    //    "Ratio to detector level", 10, 00);
    hfrSignal->Draw();

    if(ij == Njets){
      leg->AddEntry((TObject*)0, Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[0], JetPtBins2[Njets]),"");
    }else{
      leg->AddEntry((TObject*)0, Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[ij], JetPtBins2[ij+1]),"");
    }
    leg->AddEntry((TObject*)0, title, "");
    leg->AddEntry((TObject*)0, title2, "");
    errBg->SetLineColor(qColor[1]);
    errBg->SetMarkerColor(qColor[1]);
    errBg->SetMarkerSize(mSize);
    errBg->SetMarkerStyle(qMarker[1]);

    /*TF1 * bgFit = new TF1("bgFit","pol2",0,5);
      TF1 * unfFit = new TF1("unfFit","pol2",0,5);
      TF1 * totFit = new TF1("totFit","pol2",0,5);*/
    TF1 * bgFit = new TF1("bgFit","pol2",0.1,1);
    errBg->Draw("PZ same");
    bgFit->SetLineColor(qColor[1]);
    bgFit->SetLineStyle(2);


    /*errUnf->Fit("unfFit","QN","",0.1,5);
      errBg->Fit("bgFit","QN","",0.1,5);
      errTot->Fit("totFit","QN","",0.1,5);
      totFit->Draw("l same");


      TH1D *hint = new TH1D(Form("hint%d",1),
      "Fitted gaussian with .95 conf.band", 20, 0.1,5);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
      hint->SetStats(kFALSE);
      hint->SetFillColor(qColor[3]);
      hint->SetFillColorAlpha(qColor[3], 0.35);
      hint->Draw("e3 same");*/


    leg->AddEntry(errBg,"Background method","p");
    //leg->AddEntry(errUnf,"Unfolding method","p");
    //leg->AddEntry(errTot,"Quadratic Sum","p");
    //leg->AddEntry(totFit,"Fit","l");
    leg->Draw();
    mpad = fpadSignal->GetPad(2); //lower pad
    //gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(0);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.17);

    TH2F *hfrSignal2 = new TH2F("hfr4", "", 10, xlow,xhigh,10,0,1.05);
    hset( *hfrSignal2, "p_{T,jet}",Form("Relative error"));
    //hset2(*hfrSignal2, "pT","Ratio");
    //hset2( iCanvas, 5, 500, 0.1, 20.5 , "p_{T}", 
    //    "Ratio to detector level", 10, 00);
    //hset( *hfrSignal2, "p_{T} [GeV/c] ", "Ratio to particle level",
    //    0.9, 0.9, 0.1,0.09, 0.01,0.001, 0.05,0.05, 505, 505);
    hfrSignal2->Draw();
    TGraphErrors *tempRatioBg = grrDivide(errBg,histo);
    //TGraphErrors *tempRatioUnf = grrDivide(errUnf,histo);
    //TGraphErrors *tempRatioTot = grrDivide(errTot,histo);
    tempRatioBg->SetMarkerStyle(qMarker[1]);
    tempRatioBg->SetLineColor(qColor[1]);
    tempRatioBg->SetMarkerColor(qColor[1]);
    tempRatioBg->SetMarkerSize(mSize);
    tempRatioBg->Fit("bgFit","QN","",0.1,1);
    //bgFit->Draw("l same");
    //unfFit->Draw("l same");
    //totFit->Draw("l same");

    /*TH1D *hint = new TH1D(Form("hint%d",1),
      "Fitted gaussian with .95 conf.band", 20, 0.1,1);
      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hint);
      hint->SetStats(kFALSE);
      hint->SetFillColor(qColor[3]);
      hint->SetFillColorAlpha(qColor[3], 0.35);
      hint->Draw("e3 same");*/

    tempRatioBg->Draw("PZ same");

    //tempRatio->Divide(fit);
    //tempRatio->Draw("PZ same");
    //gPad->Modified();
    gPad->Modified();
    gPad->GetCanvas()->Update();
    //cout << "ppdf: " << file.Data() << endl;
    //ppdf(Form("figs/%sNFin%02dJetPt%02d%s",file.Data(),is,ij,topcomment.Data()),gPad->GetCanvas());
    TString figname;
    if(xlog){
      figname = Form("figs/%sNFin%02dJetPt%02d%s.pdf",file.Data(),is,ij,comment.Data());
    }else{
      figname = Form("figs/%sNFin%02dJetPt%02d_linx%s.pdf",file.Data(),is,ij,comment.Data());
    }
    if(saveFigs) gPad->GetCanvas()->SaveAs(figname);
  }

  void drawErrors2(TGraphErrors *histo, TGraphErrors *error,double xlow, double xhigh, int xlog, double ylow, double yhigh, int ylog, TString title,TString title2,TString comment,TString file, int is, int ij,double limit){
    Filipad2 *fpadSignal = new Filipad2(iC++, 1.1, 0.5, 150, 50, 0.7, 6,1);
    fpadSignal->Draw();
    TPad *mpad = fpadSignal->GetPad(1); //upper pad

    gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(ylog);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.17);

    TLegend *leg;
    if(xlog){
      leg = new TLegend(0.15,0.10,0.62,0.55, "","brNDC");
    }else{
      leg = new TLegend(0.45,0.55,0.92,0.80, "","brNDC");
    }
    leg->SetFillStyle(0);leg->SetBorderSize(0);leg->SetTextSize(0.045);

    TH2F *hfrSignal = new TH2F("hfr3"," ", 10, xlow , xhigh,10, ylow, yhigh);
    hset( *hfrSignal, "j_{T} [GeV]","Absolute Error");
    hfrSignal->Draw();

    /*if(ij == Njets){
      leg->AddEntry((TObject*)0, Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[0], JetPtBins2[Njets]),"");
    }else{
      leg->AddEntry((TObject*)0, Form("p_{T,jet}: %3.0f-%3.0f GeV/c",JetPtBins2[ij], JetPtBins2[ij+1]),"");
    }*/
    leg->AddEntry((TObject*)0, title, "");
    //leg->AddEntry((TObject*)0, title2, "");
    leg->AddEntry(error,title2,"p");
    error->SetLineColor(qColor[1]);
    error->SetMarkerColor(qColor[1]);
    error->SetMarkerSize(mSize);
    error->SetMarkerStyle(qMarker[1]);
    TF1 * bgFit = new TF1("bgFit","pol2",0.1,1);
    error->Draw("PZ same");
    bgFit->SetLineColor(qColor[1]);
    bgFit->SetLineStyle(2);

    leg->Draw();
    mpad = fpadSignal->GetPad(2); //lower pad
    //gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(1.6);
    mpad->SetLogx(xlog);mpad->SetGridx(0);
    mpad->SetLogy(0);mpad->SetGridy(0);
    mpad->cd();
    gPad->SetLeftMargin(0.17);

    TH2F *hfrSignal2 = new TH2F("hfr4", "", 10, xlow,xhigh,10,-0.3,0.3);
    hset( *hfrSignal2, "p_{T,jet}",Form("Relative error"));
    hfrSignal2->Draw();
    TGraphErrors *tempRatioBg = grrDivide(error,histo);
    tempRatioBg->SetMarkerStyle(qMarker[1]);
    tempRatioBg->SetLineColor(qColor[1]);
    tempRatioBg->SetMarkerColor(qColor[1]);
    tempRatioBg->SetMarkerSize(mSize);

    tempRatioBg->Draw("PZ same");
    TLine *line = new TLine(40,limit,150,limit);
    TLine *line2 = new TLine(40,-1*limit,150,-1*limit);
    line->SetLineColor(1);
    line->SetLineStyle(2);
    line->Draw();
    line2->SetLineColor(1);
    line2->SetLineStyle(2);
    line2->Draw();

    gPad->Modified();
    gPad->GetCanvas()->Update();
    TString figname;
    if(xlog){
      figname = Form("figs/%sNFin%02dJetPt%02d%s.pdf",file.Data(),is,ij,comment.Data());
    }else{
      figname = Form("figs/%sNFin%02dJetPt%02d_linx%s.pdf",file.Data(),is,ij,comment.Data());
    }
    if(saveFigs) gPad->GetCanvas()->SaveAs(figname);
  }

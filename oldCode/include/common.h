int qColor[] = { kBlack, kRed, kBlue, kGreen+2, kOrange-3, kMagenta+1, kViolet, kPink, kGray, 40, 48, 9, kGreen+2, kOrange-3, kMagenta+1 };
int qMarker[] = { 20, 24, 21, 23, 27, 29,30,20, 24, 21, 25, 33, 27, 34, 28, 20, 21, 22, 23, 34, 28, 20, 21, 22 };
TH2F * hfr = NULL;
TLegend * leg = NULL;
TH1D * h = NULL;
TH2D * h2 = NULL;

enum { kJNRef=2, kJNJetPtBin=10 };
enum { kJIncommingJet, kJOutgoingJet, kJFullJet, kJFullJet08, kJChargedJet, kJChargedJet08, kJNJetType };
enum { kJEtaAll, kJEtaAlice, kJNJetSelection};
enum { kAll, kOpp, kJNDiJetSelection };

TString StrJetType[] = { "Incomming", "Outgoing", "FullJet R=0.4", "FullJet R=0.8", "ChargedJet R=0.4", "ChargedJet R=0.8" } ;
TString StrJetSel[] = { "-5<#eta<5", "-0.4<#eta<0.4" };
TString DiJetSel[] = { "Leading-sub", "Leading-sub,phi>2/#pi", "Marta" };

TString StrJetType2[] = { "Incomming", "Jet_{ideal}", "Jet_{Full,R=0.4}", "Jet_{Full,R=0.8}", "Jet_{Charged,R=0.4}", "Jet_{ChargedJet,R=0.8}" } ;
TString StrJetType3[] = { "Incomming", "Outgoing", "Full,R=0.4", "Full,R=0.8", "Charged,R=0.4", "Charged,R=0.8", "MPI" } ;
TString StrJetSel2[] = { "|#eta|<5", "|#eta|<0.4" };
TString jetBinStr[] ={"M_{jj,ch}","M_{jj,outgoing}", "p_{T,ch}","p_{T,outgoing}", "M_{jj,dijet}","p_{T,each dijet}" };
int SelJetType[]={ 0, 1, 1, 1, 1, 1, 1 };

TH2F * hset2(int ic, double x0, double x1, double y0, double y1, TString xtitle, TString ytitle, int logxy, int gridxy){
  gPad->SetLeftMargin(0.2);
  hfr= new TH2F(Form("hfr%d",ic),"", 1, x0, x1, 1, y0, y1);
  hset( *hfr, xtitle, ytitle );
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  hfr->Draw();
  gPad->SetLogx( (logxy/10)%10 );
  gPad->SetLogy( logxy%10 );
  gPad->SetGridx( gridxy/10%10 );
  gPad->SetGridy( gridxy%10 );

  return hfr;
}

TLegend * setleg( double x0, double x1, double y0, double y1, const char * title=NULL ){
  leg = new TLegend(x0,x1,y0,y1,title,"brNDC");
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.04);
  leg->Draw();
  return leg;
}

TH1D * GetNormalizedTH1D( TDirectory * dir, TString name, double nom = -1 ){
  GetNormalizedTH1D2( dir, name, 0, nom );
}

TH1D * GetNormalizedTH1D2( TDirectory * dir, TString name, double rebin = 0, double nom = -1 ){
  h = (TH1D*) dir->Get(name); 
  if( !h ) {
    cout<< name <<endl;
    gSystem->Exit(-1);
  }
  //h->Sumw2();
  if( rebin >1 ) h->Rebin( 2 );
  if( nom < 0 ) nom = h->GetEntries();
  h->Scale(1./nom, "width");
  cout<<Form("%-30s N=%-20d Mean=%-20.4f RMS=%-20.4f %s", h->GetName(), int(h->GetEntries()),  h->GetMean(), h->GetRMS(), h->GetTitle() )<<endl;
  return h;
}

TH2D * GetNormalizedTH2D( TDirectory * dir, TString name, double nom = -1 ){
  h2 = (TH2D*) dir->Get(name); 
  if( !h2 ) {
    cout<< name <<endl;
    gSystem->Exit(-1);
  }
  //h->Sumw2();
  if( nom < 0 ) nom = h2->GetEntries();
  //h2->Scale(1./nom, "width");
  cout<<Form("%-30s N=%-20d Mean=%-20.4f RMS=%-20.4f %s", h2->GetName(), int(h2->GetEntries()),  h2->GetMean(), h2->GetRMS(), h2->GetTitle() )<<endl;
  return h2;
}

TVector * JetPtBins;
void LoadVector(){
  JetPtBins = (TVector*) gDirectory->Get("JCard/Jet:PtBins");
}

const char * StrJetPtBin( int ipt, TString ptstr="p_{T}" ){
  return Form("%.0f<%s<%0.f", (*JetPtBins)[ipt+1], ptstr.Data(), (*JetPtBins)[ipt+2]);
}


void SaveAs(TString name){
  gPad->GetCanvas()->SaveAs(name+".pdf");
}

AliJHistManager * LoadHistManager(TString name="hmg"){
    gROOT->LoadMacro("AliJHistManager.cxx+");
    AliJHistManager * hmg = new AliJHistManager(name);
    hmg->LoadConfig();
    return hmg;
}

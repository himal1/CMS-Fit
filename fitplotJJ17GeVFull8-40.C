{
//setting the enviroment
  gROOT->SetStyle("Plain");

  gStyle->SetOptStat("kTRUE");
  gStyle->SetOptTitle(kFALSE);

  //gROOT->LoadMacro("tdrstyle.C");
  //setTDRStyle();

  
  // https://root.cern.ch/root/html/tutorials/roofit/rf804_mcstudy_constr.C.html
  // we need to add the RooMC study to this macro for full validation of the fits
  //https://root.cern.ch/root/html/tutorials/roofit/rf801_mcstudy.C.html
  //  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  using namespace RooFit;

  Int_t nbins=64;//128 previously
  Float_t low=8.,high=40.;
  RooRealVar mass("mass","Four Muon Mass (GeV)",low,high);
  //

  //RooRealVar run("run","run",0,9000000000);
  //RooRealVar event("event","event",0,9000000000);
  
  RooRealVar VtxProb("VtxProb","VtxProb",0.05001,1.);
  RooRealVar Psi1To2Significance("Psi1To2Significance","Psi1To2Significance",0.,9.);
  RooRealVar Psi1_Mass("Psi1_Mass","Psi1_Mass",2.8,3.40);
  RooRealVar Psi2_Mass("Psi2_Mass","Psi2_Mass",2.8,3.40);
  RooRealVar Psi1_pT("Psi1_pT","Psi1_pT",5.500001.,100);//applying Jpsi pt Cut
  RooRealVar Psi2_pT("Psi2_pT","Psi2_pT",5.500001.,100);//applying JpsiPt cut
  RooRealVar Psi1_eta("Psi1_eta","Psi1_eta",-4.4,4.4);
  RooRealVar Psi2_eta("Psi2_eta","Psi2_eta",-4.4,4.4);
  RooRealVar Psi1_CTxy("Psi1_CTxy","Psi1_CTxy",-1.,1.);
  RooRealVar Psi2_CTxy("Psi2_CTxy","Psi2_CTxy",-1.,1.);
  RooRealVar Psi1To2_dY("Psi1To2_dY","Psi1To2_dY",0.,1.499999999999999);//applying rapidityCut
  RooRealVar FourMu_Rapidity("FourMu_Rapidity","FourMu_Rapidity",-5,5);
  RooRealVar FourMu_RapidityAbs("FourMu_RapidityAbs","FourMu_RapidityAbs",-5,5);
  RooRealVar FourMu_pT("pt","pt",3.00000000001,100);//applying 4Mu pt cut3.0
  
  RooArgSet args(Psi1_Mass, Psi2_Mass, Psi1_CTxy, Psi2_CTxy, Psi1To2Significance); 
  args.add(RooArgSet(Psi1To2_dY, mass, FourMu_pT, FourMu_Rapidity,FourMu_RapidityAbs));

  RooDataSet *data =  new RooDataSet("data","data",args,0);
  data = RooDataSet::read("ALlJPsiPtGreater5_5.dat",args);

  // global variables
  Float_t before = 0;
  Float_t after = 0;
  Float_t resolution = 0.143;

  
  // --- global fit 8 - 40 GeV
  //   1  a0          -1.36753e+00   1.80275e-02   2.04273e-04  -2.37879e+00
  //   2  a1           6.88030e-01   2.47770e-02   2.32385e-04  -9.66876e-01
  //   3  a2          -4.35819e-01   2.56387e-02   1.85958e-04   1.12705e+00
  //   4  a3           3.22583e-01   2.10228e-02   1.88903e-04   3.20015e+00
  //   5  a4          -1.30328e-01   1.70113e-02   1.89297e-04   3.00583e+00
  // 6  mean         1.77864e+01   1.12818e-01   8.12834e-03  -2.02796e-02
  // 7  nbkg         3.14703e+03   5.93403e+01   2.70403e-03   9.87522e-02
  // 8  nsig         6.21707e+01   2.09396e+01   6.78832e-03   3.11422e-02
  //  9  width        2.80571e-01   1.53744e-01   4.01544e-02  -6.12993e-03
  // ----
  
  // background

  // float----------------------------
  
  RooRealVar a0("a0","a0",0.,-5.,5.);
  RooRealVar a1("a1","a1",0.,-5.,5.);
  RooRealVar a2("a2","a2",0.,-5.,5.);
  RooRealVar a3("a3","a3",0.,-5.,5.);
  RooRealVar a4("a4","a4",0.,-5.,5.);

  

  // fix----------------------

  //RooRealVar a0("a0","a0",0.);
  //RooRealVar a1("a1","a1",0.);
  //RooRealVar a2("a2","a2",0.);
  //RooRealVar a3("a3","a3",0.);
  //RooRealVar a4("a4","a4",0.);  


  //---

  RooChebychev bkg("bkg","Background",mass,RooArgSet(a0,a1,a2,a3,a4));


  before = after = 0;

  // signal

  Float_t i = 177.8;
  resolution = 80.31 + (2.01 * (i/10.)) + (0.1045*(i/10.)*(i/10.)); 
  
  //fixed variable
  RooRealVar mean("mean","{4#mu}_{mass}",i/10.,16,19);
  RooRealVar width("width","Width",0.280571,0.1,0.5);
  RooRealVar sigma("sigma","sigma 4mu mass",resolution/1000.);  
  RooVoigtian Voig("voig","voig PDF",mass,mean,sigma,width);

  //------------------------------------------------------------
  // combined hypothesis

  RooRealVar nsig("nsig","nsignal",100,0,500);
  RooRealVar nbkg("nbkg","nbackground",2500,0,4000);
  RooAddPdf sum("sum","S+B",RooArgList(Voig,bkg),RooArgList(nsig,nbkg));
  sum.fitTo(*data,Extended(kTRUE),Range(low,high));

//-------------------------------

 //setting the enviroment
  gROOT->SetStyle("Plain");

  gStyle->SetOptStat("kTRUE");
  gStyle->SetOptTitle(kFALSE);

  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  //-------------------------------
 
  //--------------------------------------------------
  //--- plotting

  int W = 800;
  int H = 600;

  int H_ref = 600;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;
  TString canvName = "FigExample_";

  canvName += W;
  canvName += "-";
  canvName += H;
  canvName += "_Stefan";  

  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);

  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx();
  canv->SetTicky();
  //canv->SetLogy();

  //-------

  RooPlot* Massframe = mass.frame(nbins);    //initially nbins

  Massframe->SetAxisRange(low,high,"X");
  

  //  Massframe->GetXaxis()->SetNdivisions(0);

  Massframe->GetXaxis()->SetTitle("4 Muon Mass  (GeV/c^{2})");

  Massframe->GetYaxis()->SetTitle("Events / 2 GeV");

  Massframe->GetYaxis()->SetTitleOffset(1.0);
  Massframe->GetYaxis()->SetTitleSize(0.055);
  Massframe->SetMaximum(150);
  //Massframe->SetMinimum(-0.0000001)
  //Massframe->GetYaxis()->SetRangeUser(1,105.);
  
  //Massframe->SetMarkerStyle(21);
  //gr->SetMarkerSize(1.3);


  gStyle->SetEndErrorSize(0);
  Massframe->SetMarkerSize(1.0);
  data->plotOn(Massframe,DataError(RooAbsData::SumW2)) ;
  sum.plotOn(Massframe,Name("Total"));
  //sum.plotOn(Massframe,Name("Background"),Components(backge),LineStyle(kDashed));



  //nsigz=2.;
  sum.plotOn(Massframe,Name("V"),Components(Voig),LineStyle(kDashed),LineColor(kRed));
  
  sum.plotOn(Massframe,Name("Background2"),Components(bkg),LineStyle(kDashed));

  //nsig=2.;
  //sum.plotOn(Massframe,Name("Signal"),Components(gaussH),LineStyle(kDashed),LineColor(kRed)); 
  
  Massframe->SetMarkerStyle(21);
  Massframe->Draw();
  RooPlot* frame1 = mass.frame(nbins);
  data->plotOn(frame1,DataError(RooAbsData::SumW2)) ;
  sum.plotOn(frame1,Name("Total"));
  
  //data->SetMarkerStyle(21);
  //frame1->SetMarkerStyle(21);
  //-------


  float WH = canv->GetWh();
  float WW = canv->GetWw();
  float l = canv->GetLeftMargin();
  float t = canv->GetTopMargin();
  float r = canv->GetRightMargin();
  float b = canv->GetBottomMargin();

  TString cmsText     = "CMS";
  float cmsTextFont   = 61;  // default is helvetic-bold
  TString extraText   = "Preliminary";
  float extraTextFont = 52;  // default is helvetica-italics
  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float lumiTextSize     = 0.6;
  float lumiTextOffset   = 0.2;
  float cmsTextSize      = 0.75;
  float cmsTextOffset    = 0.1;  // only used in outOfFrame version
  float relPosX    = 0.045;
  float relPosY    = 0.035;
  float relExtraDY = 1.2;
  
  TString lumiText = "46.1 fb^{-1}  (13 TeV)";
  
  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    


  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;
  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  float posX_=0.16;
  //  float posX_=0.22;
  //  float posY_=0.88;
  //  float posY_=0.99;
  float posY_=.977;

  latex.SetTextFont(cmsTextFont);
  latex.SetTextSize(cmsTextSize*t);
  latex.SetTextAlign(23);
  latex.DrawLatex(posX_, posY_, cmsText);

  
  //  posY_ =   1-t+lumiTextOffset*t;
  posY_ = 0.967;
  posX_ = posX_ + .15;

  latex.SetTextFont(extraTextFont);
  latex.SetTextSize(extraTextSize*t);
  latex.SetTextAlign(23);
  latex.DrawLatex(posX_, posY_, extraText);      

  //-------------TLegend WOrk---------------
  TLegend *leg = new TLegend(0.50,0.60,0.89,0.85);//(0.50,0.70,0.89,0.85)
  leg->SetFillColor(0);
  //leg->SetBorderMode(0);
  //leg->SetFrameFillStyle(0);
  //leg->AddEntry(Massframe,"Data","lep");

  TLegendEntry *entry=leg->AddEntry("frame1","Data","lep");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(8);
  entry->SetMarkerSize(1);
  //entry=leg->AddEntry("Total","Background","l");
  
  leg->AddEntry(Total,"Signal+Background","l");
  leg->AddEntry(V,"Model Signal","l");
  leg->AddEntry(Background2,"Background","l");
  leg->SetTextAlign(23);
  leg->SetTextFont(42);
  leg->SetTextSize(extraTextSize*t);
  leg->SetBorderSize(0);
  leg->Draw();
  //cout <<"The value of slope is "<<slope<<endl;
}






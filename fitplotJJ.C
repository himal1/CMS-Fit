{
 //setting the enviroment
  gROOT->SetStyle("Plain");

  gStyle->SetOptStat("kTRUE");
  gStyle->SetOptTitle(kFALSE);

  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  //-------------------------------
  using namespace RooFit;

  Int_t nbins=25;
  Float_t low = 40.;
  Float_t high = 140.;

  Float_t momentum = 0.;
  cout << momentum << endl;

  RooRealVar m4m("m4m","4 muon mass / GeV",low,high);
  RooRealVar pt4m("pt4m","4 muon p_T",momentum,300.);
 
  RooRealVar deltay("deltay","4 muon Delta Y",0.,3.);

  RooRealVar mj1("mj1","2 muon mass high pT",3.0,3.2);
  RooRealVar mj2("mj2","2 muon mass low pT",2.95,3.25);

  RooRealVar tiso("tiso","track isolation",-1.0,0.5);

  RooArgSet args(m4m, pt4m, deltay,mj1,mj2,tiso);

  RooDataSet *data =  new RooDataSet("data","data",args,0);
  //  data = RooDataSet::read("data.dat",args);
  data = RooDataSet::read("blind.dat",args);


  // write the candidate events (186)
data.Print("v") ;
data.write("himal.dat");
// 




  // make the Higgs
  //RooRealVar meanH("meanH","{4 muon}_{mass} / GeV",125.00);
  //RooRealVar sigmaH("sigmaH","sigma 4mu mass",1.5);  
  //RooGaussian gaussH("gaussH","gaussian PDF", m4m, meanH, sigmaH);//for lower level fit
  //---------------higgs new ------------------
  RooRealVar meanH1("meanH1","{4 muon}_{mass} / GeV",1.24824e+02);
  RooRealVar sigmaH1("sigmaH1","sigma 4mu mass",1.76759);
  RooGaussian gaussH1("gaussH1","gaussian PDF", m4m, meanH1, sigmaH1);

  RooRealVar sigmaH2("sigmaH2","sigma 4mu mass",8.20949e-01);
  RooGaussian gaussH2("gaussH2","gaussian PDF", m4m, meanH1, sigmaH2);
  RooRealVar hfrac("hfrac","Higgs fraction",2.58215e-01);
  RooAddPdf gaussH("gaussH","Higgs Sum",RooArgList(gaussH1,gaussH2),hfrac);

  
  // make the Z
  //RooRealVar meanZ("meanZ","{4 muon}_{mass} / GeV",91.00);
  //RooRealVar sigmaZ("sigmaZ","sigma 4mu mass",1.0);  
  //RooGaussian gaussZ("gaussZ","gaussian PDF", m4m, meanZ, sigmaZ);//for lower level fit

  RooRealVar meanZ("meanZ","{4 muon}_{mass} / GeV",91.1876);
  RooRealVar sigmaZ("sigmaZ","sigma 4mu mass",0.5);
  RooRealVar widthZ("widthZ","width 4mu mass",2.4952);
  RooVoigtian gaussZ("gaussZ","Z Peak",m4m,meanZ,widthZ,sigmaZ);
  

  // Background - exponential
  RooRealVar slope("slope","Slope",-0.49,-10,10);
  RooExponential backge("backge","Background", m4m, slope);


  RooRealVar a0("a0","a0",2.) ;
  RooRealVar a1("a1","a1",0.) ;
  RooRealVar a2("a2","a2",0.) ;

  //  RooExponential p2("p2","Exp 2", m4m, a0);

//RooPolynomial p2("p2","p2",m4m,RooArgList(a0,a1,a2),0) ;
//  RooPolynomial p2("p2","p2",m4m,RooArgList(a0,a1,a2),0) ;

  RooUniform p2("p2","p2",m4m);
  
  RooRealVar frac("frac","fraction",0.1,0.,1.);

  RooAddPdf backg("backg","Background Sum",RooArgList(backge,p2),frac);


  // add hypothesis
  RooRealVar nsig("nsig","nsignal",0.1,0.,10.);
  RooRealVar nsigz("nsigz","nsignalz",0.1,0,100);
  RooRealVar nbkg("nbkg","nbackground",500,0,4000);
  RooRealVar nexp("nexp","nbackgroundexp",500,0,4000);
  RooRealVar npol("npol","nbackgroundpol",5,-5,100);

  //  RooAddPdf sum("sum","S+BH",RooArgList(gaussH,backg),RooArgList(nsig,nbkg));
  //  RooAddPdf sum("sum","S+BH+BZ",RooArgList(gaussH,gaussZ,backg),RooArgList(nsig,nsigz,nexp));
  RooAddPdf sum("sum","S+BH+BZ",RooArgList(gaussH,gaussZ,backg),RooArgList(nsig,nsigz,nexp));
  sum.fitTo(*data,Extended(kTRUE),Range(low,high));

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
  canv->SetLogy();

  //-------

  RooPlot* Massframe = m4m.frame(nbins);    //initially nbins

  Massframe->SetAxisRange(low,high,"X");
  

  //  Massframe->GetXaxis()->SetNdivisions(0);

  Massframe->GetXaxis()->SetTitle("4 Muon Mass  (GeV/c^{2})");

  Massframe->GetYaxis()->SetTitle("Events / 4 GeV");

  Massframe->GetYaxis()->SetTitleOffset(1.0);
  Massframe->GetYaxis()->SetTitleSize(0.055);
  Massframe->SetMaximum( 105 );
  //Massframe->SetMinimum(-0.0000001)
  //Massframe->GetYaxis()->SetRangeUser(1,105.);
  
  //Massframe->SetMarkerStyle(21);
  //gr->SetMarkerSize(1.3);


  gStyle->SetEndErrorSize(0);
  Massframe->SetMarkerSize(1.0);
  data->plotOn(Massframe,DataError(RooAbsData::SumW2)) ;
  sum.plotOn(Massframe,Name("Total"));
  //sum.plotOn(Massframe,Name("Background"),Components(backge),LineStyle(kDashed));



  nsigz=2.;
  sum.plotOn(Massframe,Name("Z"),Components(gaussZ),LineStyle(kDashed),LineColor(kRed));
  
  //sum.plotOn(Massframe,Name("Background2"),Components(p2),LineStyle(kDashed));

  nsig=2.;
  sum.plotOn(Massframe,Name("Signal"),Components(gaussH),LineStyle(kDashed),LineColor(kRed)); 
  
  Massframe->SetMarkerStyle(21);
  Massframe->Draw();
  RooPlot* frame1 = m4m.frame(nbins);
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
  
  TString lumiText = "37.5 fb^{-1}  (13 TeV)";
  
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
  TLegend *leg = new TLegend(0.5,0.70,0.89,0.85);
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
  
  leg->AddEntry(Total,"Background","l");
  leg->AddEntry(Z,"Model Boson Signals","l");
  //leg->AddEntry(Signal,"Fixed H-Boson(2-event)","l");
  leg->SetTextAlign(23);
  leg->SetTextFont(42);
  leg->SetTextSize(extraTextSize*t);
  leg->SetBorderSize(0);
  leg->Draw();
  cout <<"The value of slope is "<<slope<<endl;
}






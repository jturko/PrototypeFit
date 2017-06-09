
#include "NeutronFit_BC537.hh"

NeutronFit_BC537::NeutronFit_BC537(int run_num) :
    fRunNum(run_num),
    fSimFile(NULL),
    fExpFile(NULL),
    fSimHist(NULL),
    fExpHist(NULL),
    fSimTree(NULL),
    fEdepBranch(NULL),
    fEkinBranch(NULL),
    fPtypeBranch(NULL),
    fEdepVector(NULL),
    fEkinVector(NULL),
    fPtypeVector(NULL),
    fFitFunc(NULL),
    fRebin(false)
{
 
    double energy_vector[] = 
    {
        0.687649    ,0.651889    ,0.596412    ,0.526792    ,0.488742    ,0.371838    ,0.266086    ,0.18371      ,0.18371     ,0.126838    ,
        0.100888    ,0.075058    ,0.059957    ,0.816573    ,0.695436    ,0.595748    ,0.516107    ,0.45421      ,0.389146    ,1.974584    ,
        1.900621    ,1.784647    ,1.636554    ,1.46825     ,1.292083    ,1.119299    ,0.958821    ,0.816573     ,2.964447    ,2.860877    ,
        2.698133    ,2.48962     ,2.251518    ,2.000698    ,1.752682    ,1.520026    ,1.217669    ,1.052637     ,0.808785    ,0.662011    ,
        3.308838    ,2.957477    ,2.608323    ,2.278822    ,1.981274    ,1.722613    ,1.50509     ,1.03544      ,4.300737    ,4.157133    ,
        3.931182    ,3.641081    ,3.308838    ,7.908938    ,7.643525    ,7.226025    ,6.690194    ,6.076821     ,5.428531    ,4.784732    ,
        4.177603    ,3.629763    ,3.153886    ,17.945689   ,12.526461   ,20.907393   ,20.634645   ,20.196496    ,19.615953   ,18.9222     ,
        18.148086   ,17.327541   ,16.493228   ,15.674682   ,14.897078   ,14.180666   ,13.540783   ,12.988332    ,12.338569   
    };
    fEnergy = energy_vector[fRunNum];

    double cutoff_low_vector[] =
    {
        0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,40  ,
        40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,
        40  ,40  ,40  ,40  ,40  ,40  ,40  ,40  ,210 ,210 ,205 ,205 ,205 ,240 ,250 ,260 ,250 ,280 ,280 ,290 ,
        290 ,260 ,280 ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   
    };
    fCutoffLow = cutoff_low_vector[fRunNum];
    
    double cutoff_high_vector[] =
    {
        190 ,170 ,150 ,130 ,120 ,90  ,60  ,45  ,45  ,35  ,35  ,30  ,30  ,190 ,170 ,130 ,120 ,100 ,85  ,500 ,
        460 ,420 ,360 ,320 ,270 ,215 ,180 ,140 ,850 ,800 ,735 ,650 ,550 ,460 ,380 ,310 ,230 ,190 ,140 ,110 ,
        940 ,810 ,690 ,560 ,460 ,380 ,310 ,185 ,1070,1010,945 ,850 ,750 ,5000,4900,4500,4000,3500,3050,2600,
        2150,1700,1400,4500,3100,3700,3750,5200,3000,3000,3000,4000,4000,3750,3550,3400,3500,3400,3000       
    };
    fCutoffHigh = cutoff_high_vector[fRunNum];
 
    fExpFile = TFile::Open("~/data/prototype2013.root"); 

    std::string hist_name = "WhiteCal" + std::to_string(fRunNum);
    std::string title = std::to_string(fEnergy) + " MeV";
    fExpHist = (TH1F*)(fExpFile->Get(hist_name.c_str())->Clone());
    fExpHist->SetNameTitle(hist_name.c_str(),title.c_str());
    fExpHist->Scale(10000./fExpHist->Integral());
    fExpHist->SetLineColor(kBlack);
    fExpHist->GetXaxis()->UnZoom(); 
    fExpHist->GetYaxis()->UnZoom();
    fExpHist->GetXaxis()->SetTitleOffset(0.75);
    fExpHist->GetXaxis()->SetTitleSize(0.05);
    fExpHist->SetStats(false);
    fExpBinNum = fExpHist->GetNbinsX();
    fExpBinHigh = fExpHist->GetBinLowEdge(fExpBinNum+1);
    fExpBinLow = fExpHist->GetBinLowEdge(1);
    //fCutoffHigh = fExpBinHigh;
    //fCutoffLow = fExpBinLow;
    ApplyCutoffLow(fCutoffLow,"exp");

    std::string name = "~/data/smearing/prototype/G4_RAW/Sim" + std::to_string(fRunNum) + "/g4out.root";
    //std::string name = "/nessa/geant4/joey/data/smearing/prototype/G4_RAW/Sim" + std::to_string(fRunNum) + "/g4out.root";
    fSimFile = TFile::Open(name.c_str());     

    fSimTree = (TTree*)(fSimFile->Get("ntuple/ntuple")); 
    fNumEntries = fSimTree->GetEntries();
    
    fSimTree->SetBranchAddress("eDepVector",&fEdepVector,&fEdepBranch);
    fSimTree->SetBranchAddress("eKinVector",&fEkinVector,&fEkinBranch);
    fSimTree->SetBranchAddress("particleTypeVector",&fPtypeVector,&fPtypeBranch);
    
    fProtonCoeff[0] = 0.74; fProtonCoeff[1] = 3.2; fProtonCoeff[2] = 0.20; fProtonCoeff[3] = 0.97;
    fDeuteronCoeff[0] = 0.75; fDeuteronCoeff[1] = 2.80; fDeuteronCoeff[2] = 0.25; fDeuteronCoeff[3] = 0.93;
    fCarbonCoeff[0] = 0.05; fCarbonCoeff[1] = 0.0; fCarbonCoeff[2] = 0.0;fCarbonCoeff[3] = 0.0;
    fAlphaCoeff[0] = 0.14; fAlphaCoeff[1] = 0.59; fAlphaCoeff[2] = 0.065; fAlphaCoeff[3] = 1.01;
    fBeCoeff[0] = 0.0821; fBeCoeff[1] = 0.0; fBeCoeff[2] = 0.0; fBeCoeff[3] = 0.0;
    fBCoeff[0] = 0.0375; fBCoeff[1] = 0.0; fBCoeff[2] = 0.0; fBCoeff[3] = 0.0;
    
    fSmearingCoeff[0] = 0.123; fSmearingCoeff[1] = 0.125; fSmearingCoeff[2] = 0.0075;

    fParameters[0] = 0.639;
    fParameters[1] = 1.462;
    fParameters[2] = 0.373;
    fParameters[3] = 0.968;
    fParameters[4] = 0.0;
    fParameters[5] = 0.123;
    fParameters[6] = 0.125;
    fParameters[7] = 0.0074;
    SetParameters(fParameters);
  
    fOffset = 0;

    //fSimSortMax = 200000;
    
    //if(fSimTree->GetEntries() > fExpHist->GetEntries()) fSimSortMax = fExpHist->GetEntries();
    //else fSimSortMax = fSimTree->GetEntries();
    
    if(fSimTree->GetEntries() >= 2e5) fSimSortMax = 2e5;
    else fSimSortMax = fSimTree->GetEntries();

    std::cout << "Run# = " << fRunNum << " ; Energy = " << fEnergy << " MeV ; cutoff(low,high) = (" << fCutoffLow << ","; 
    std::cout << fCutoffHigh << ") " << " ; #evts ratio = " << double(fSimSortMax)/double(fExpHist->GetEntries()) << std::endl;
    
    fSimSortMax = fNumEntries;

    TF1 * tmpFit = NULL;
    tmpFit = fExpHist->GetFunction("fit");
    if(tmpFit) tmpFit->SetBit(TF1::kNotDraw);
    tmpFit = fExpHist->GetFunction("");
    if(tmpFit) tmpFit->SetBit(TF1::kNotDraw);

    fExpHist->Rebin(10);
    //if(fExpBinNum == 50100) fExpHist->Rebin(10);

}

NeutronFit_BC537::~NeutronFit_BC537() {}

void NeutronFit_BC537::SetParameters(double * par)
{
    fDeuteronCoeff[0] = par[0];
    fDeuteronCoeff[1] = par[1];
    fDeuteronCoeff[2] = par[2];
    fDeuteronCoeff[3] = par[3];
    fCarbonCoeff[0] = par[4];
    fSmearingCoeff[0] = par[5]; 
    fSmearingCoeff[1] = par[6]; 
    fSmearingCoeff[2] = par[7]; 

    fParameters[0] = par[0];
    fParameters[1] = par[1];
    fParameters[2] = par[2];
    fParameters[3] = par[3];
    fParameters[4] = par[4];
    fParameters[5] = par[5];
    fParameters[6] = par[6];
    fParameters[7] = par[7];
}

void NeutronFit_BC537::Sort(double a1, double a2, double a3, double a4, double carbon, double A, double B, double C)
{
    double par[8];
    par[0] = a1;
    par[1] = a2;
    par[2] = a3;
    par[3] = a4;
    par[4] = carbon;
    par[5] = A;    
    par[6] = B;    
    par[7] = C;    

    Sort(par);
}

void NeutronFit_BC537::Sort(double * par)
{
    fRandom.SetSeed(1);
    gErrorIgnoreLevel = kError;    
    
    SetParameters(par);
    //PrintParameters();
    
    fExpBinNum = fExpHist->GetNbinsX();
    if(fSimHist) { delete fSimHist; fSimHist = NULL; }
    //fSimHist = new TH1F("fSimHist","fSimHist",fExpBinNum,-10,5000); 
    fSimHist = new TH1F("fSimHist","fSimHist",fExpBinNum,fExpBinLow,fExpBinHigh); 
    int nHits = 0;
    double light = 0.;
    double centroidEkin = 0.;    
    double centroidEres = 0.;    
    
    //clock_t overalstart = clock();
    //double resolutiontime = 0.;
    //double runningtime = 0.;
    
    int counter = 0;
    
    for(int i=0; i<fSimSortMax; i++)
    {
        counter++;
        if( counter%50000==0 ) std::cout << "sorting " << fEnergy << " MeV... " << "evt " << counter << "/" << fSimSortMax << "; " << double(counter)/double(fSimSortMax)*100 << "% complete \r"  << std::flush; 
     
        fEdepBranch->GetEntry(i);   
        fEkinBranch->GetEntry(i);   
        fPtypeBranch->GetEntry(i);   
        nHits = fEdepVector->size();
        light = 0.;
        for(int j=0; j<nHits; j++)
        {
            if(fPtypeVector->at(j) == 2 || fPtypeVector->at(j) == 3) {
                centroidEkin = fEkinVector->at(j);
                centroidEres = fEkinVector->at(j)-fEdepVector->at(j);
            }
            else if(fPtypeVector->at(j) == 4) {
                centroidEkin = LightOutput(fEkinVector->at(j), fProtonCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fProtonCoeff);
            }
            else if(fPtypeVector->at(j) == 6) {
                centroidEkin = LightOutput(fEkinVector->at(j), fDeuteronCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fDeuteronCoeff);
            }
            else if(fPtypeVector->at(j) == 7) {
                centroidEkin = LightOutput(fEkinVector->at(j), fCarbonCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fCarbonCoeff);
            }
            else if(fPtypeVector->at(j) == 8) {
                centroidEkin = LightOutput(fEkinVector->at(j), fAlphaCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fAlphaCoeff);
            }
            else if(fPtypeVector->at(j) == 9) {
                centroidEkin = LightOutput(fEkinVector->at(j), fBeCoeff);
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fBeCoeff);
            }
            else if(fPtypeVector->at(j) == 10) {
                centroidEkin = LightOutput(fEkinVector->at(j), fBCoeff );
                centroidEres = LightOutput(fEkinVector->at(j)-fEdepVector->at(j), fBCoeff);
            }
            else { 
                centroidEkin = 0.; 
                centroidEres = 0.; 
            } 
            //clock_t resolutionstart = clock();
            if(centroidEkin>0.){
                light += 1000.*fRandom.Gaus(centroidEkin, Resolution(centroidEkin,fSmearingCoeff));
            }
            if(centroidEres>0.){
                light -= 1000.*fRandom.Gaus(centroidEres, Resolution(centroidEres,fSmearingCoeff));
            } 
            //clock_t resolutionend = clock();
            //resolutiontime += (double)(resolutionend - resolutionstart);
        }//end scatters loop       
        if(light>0.) fSimHist->Fill(light+fOffset);
        //if(light>0.) fSimHist->Fill(light);
    }//end event loop
    
    //clock_t overalend = clock();
    //std::cout << "Overall time: " << (int)(overalend - overalstart) << std::endl;
    //std::cout << "Resolution time: " << (int)resolutiontime << std::endl;

    fExpHist->SetBinContent(fExpBinNum+1,0);

    ApplyCutoffLow(fCutoffLow,"sim");    
    fSimHist->Scale(fExpHist->Integral(fExpHist->FindBin(fCutoffLow),fExpHist->FindBin(fCutoffHigh),"width")/fSimHist->Integral(fSimHist->FindBin(fCutoffLow),fSimHist->FindBin(fCutoffHigh),"width"));
    fSimHist->SetStats(false);
    std::cout << "sorting " << fEnergy << " MeV... done!                                                   " << std::endl;

    std::string title = std::to_string(fEnergy) + " MeV - Run " + std::to_string(fRunNum) + "; #chi^{2} = " + std::to_string(DoChi2());
    fExpHist->SetTitle(title.c_str());
    
    //if(fFitFunc) { delete fFitFunc; fFitFunc = NULL; }
    //fFitFunc = new TF1("fFitFunc",this,&NeutronFit_BC537::HistCompare,fCutoffLow,fCutoffHigh,8);
    //fFitFunc->SetNpx(100);
    //fFitFunc->SetParameters(fParameters[0],fParameters[1],fParameters[2],fParameters[3],fParameters[4],fParameters[5],fParameters[6],fParameters[7]);
    //fFitFunc->SetParLimits(0,0.2,1);
    //fFitFunc->SetParLimits(1,0.5,10);
    //fFitFunc->SetParLimits(2,0.05,0.4);
    //fFitFunc->SetParLimits(3,0.8,1.2);
    //fFitFunc->SetParLimits(4,0,0.1);
    //fFitFunc->SetParLimits(5,0,0.3);
    //fFitFunc->SetParLimits(6,0,0.3);
    //fFitFunc->SetParLimits(7,0,0.05);

}

double NeutronFit_BC537::HistCompare(double * x, double * par) 
{
    if(!fSimHist) Sort();
    if(DidParametersChange(par)) Sort(par);
    double xx = x[0];
    int bin = fSimHist->GetXaxis()->FindBin(xx);
    double content = fSimHist->GetBinContent(bin);
    return content;

}

TF1 * NeutronFit_BC537::Fit()
{
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Combination");
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit","Simplex");
    TVirtualFitter::SetPrecision(1.0e-10);
    TVirtualFitter::SetMaxIterations(10000);
    TFitResultPtr res = fExpHist->Fit("fFitFunc","RSV");
    return fFitFunc;
}

bool NeutronFit_BC537::DidParametersChange(double * par)
{
    for(int i=0; i<8; i++) 
    {
        if(TMath::Abs(par[i] - fParameters[i]) > 1e-4) return true;
    }
    return false;
}




#ifndef NEUTRONFIT_BC537_H
#define NEUTRONFIT_BC537_H
#endif

#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TColor.h"
#include "TApplication.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom3.h"

#include "TMath.h"

#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

class NeutronFit_BC537 
{

public:
    NeutronFit_BC537(int run_num=0);
    ~NeutronFit_BC537();

    void PrintParameters() {
        std::cout << "     a1 = " << fParameters[0] << std::endl;
        std::cout << "     a2 = " << fParameters[1] << std::endl;
        std::cout << "     a3 = " << fParameters[2] << std::endl;
        std::cout << "     a4 = " << fParameters[3] << std::endl;
        std::cout << " carbon = " << fParameters[4] << std::endl;
        std::cout << "      A = " << fParameters[5] << std::endl;
        std::cout << "      B = " << fParameters[6] << std::endl;
        std::cout << "      C = " << fParameters[7] << std::endl;
        std::cout << " offset = " << fOffset << " keVee " << std::endl;
    }

    void Sort(double * par);
    void Sort(double a1=0.639, double a2=1.462, double a3=0.373, double a4=0.968, double carbon=0, double A=0.123, double B=0.125, double C=0.0074);

    void SetSmearingCoeff(double A=0.123, double B=0.125, double C=0.0074) {
        fSmearingCoeff[0] = A;
        fSmearingCoeff[1] = B;
        fSmearingCoeff[2] = C;
    }
    double GetSmearingCoeff(int i) { 
        if(i>=0 && i<=2) return fSmearingCoeff[i]; 
        else { std::cout << "error in GetSmearingCoeff() !" << std::endl; return 0; }
    }    

    void SetOffset(double offset) {
        fOffset = offset;
    }
    double GetOffset() { return fOffset; }

    void SetParameters(double * par);

    double LightOutput(double E, double * par) {
        return ( par[0]*E-par[1]*(1.0-TMath::Exp(-par[2]*TMath::Power(E,par[3]))) );
        //return ( fOffset/1000. + par[0]*E-par[1]*(1.0-TMath::Exp(-par[2]*TMath::Power(E,par[3]))) );
    }
    double Resolution(double E, double * par) {
        return (E*TMath::Sqrt(TMath::Power(par[0],2)+TMath::Power(par[1],2)/E+TMath::Power(par[2]/E,2)))/(2.*TMath::Sqrt(2.*TMath::Log(2.))) ;
    }

    void ApplyCutoffLow(double cutoff, std::string str) {
        if(str == "sim")      for(int i=0; i<fSimHist->FindBin(cutoff); i++) fSimHist->SetBinContent(i,0.);
        else if(str == "exp") for(int i=0; i<fExpHist->FindBin(cutoff); i++) fExpHist->SetBinContent(i,0.);
        else                  std::cout << "--->  Error with ApplyCutoffLow!" << std::endl;
    }
    void ApplyCutoffHigh(double cutoff, std::string str) {
        if(str == "sim")      for(int i=fSimHist->FindBin(cutoff); i<fSimHist->GetNbinsX(); i++) fSimHist->SetBinContent(i,0.);
        else if(str == "exp") for(int i=fExpHist->FindBin(cutoff); i<fExpHist->GetNbinsX(); i++) fExpHist->SetBinContent(i,0.);
        else                  std::cout << "--->  Error with ApplyCutoffHigh!" << std::endl;
    }
    
    double DoChi2() {
        TH1F * h_e = (TH1F*)fExpHist->Clone();
        TH1F * h_s = (TH1F*)fSimHist->Clone();
        for(int i=h_e->FindBin(fCutoffHigh); i<h_e->GetNbinsX(); i++) h_e->SetBinContent(i,0.);
        for(int i=h_s->FindBin(fCutoffHigh); i<h_s->GetNbinsX(); i++) h_s->SetBinContent(i,0.);
        for(int i=0; i<h_e->FindBin(fCutoffLow); i++) h_e->SetBinContent(i,0.);
        for(int i=0; i<h_s->FindBin(fCutoffLow); i++) h_s->SetBinContent(i,0.);
        fChi2 = h_e->Chi2Test(h_s,"CHI2/NDF");
        delete h_e;
        delete h_s;

        return fChi2;
    }

    void Draw() {
        fExpHist->SetLineColor(kBlack);
        if(fSimHist) fSimHist->SetLineColor(kRed);
        else { std::cout << "no sim hist sorted yet!" << std::endl; return; }
        fExpHist->GetXaxis()->SetRangeUser(fCutoffLow-30,fCutoffHigh+150);
        double ymax = 1.25*fExpHist->GetBinContent(fExpHist->GetMaximumBin());
        fExpHist->GetYaxis()->SetRangeUser(0.1,ymax);

        fExpHist->Draw();
        //fFitFunc->Draw("same");
        fSimHist->Draw("same");   
        //fExpHist->Draw("same");   
    }

    double GetEnergy() { return fEnergy; }
    int GetRunNum() { return fRunNum; }

    void SetCutoffHigh(double cut) { fCutoffHigh = cut; }
    void SetCutoffLow(double cut) { fCutoffLow = cut; }
    double GetCutoffHigh() { return fCutoffHigh; }    
    double GetCutoffLow() { return fCutoffLow; }    

    TH1F * GetSimHist() { return fSimHist; }
    TH1F * GetExpHist() { return fExpHist; }
    
    double GetSimEntries() { return fNumEntries; }

    double HistCompare(double * x, double * par);
    TF1 * Fit();
    bool DidParametersChange(double * par);

    TF1 * GetFitFunc() { return fFitFunc; }

    TF1 * fFitFunc;

    double fEnergy;
    int fRunNum;

    double fCutoffHigh;
    double fCutoffLow;
    
    TFile * fSimFile;
    TFile * fExpFile;

    TH1F * fSimHist;
    TH1F * fExpHist;

    TTree * fSimTree;
    TBranch * fEdepBranch;
    TBranch * fEkinBranch;
    TBranch * fPtypeBranch;
    std::vector<double> * fEdepVector;
    std::vector<double> * fEkinVector;
    std::vector<int> * fPtypeVector;

    double fProtonCoeff[4];
    double fDeuteronCoeff[4];
    double fCarbonCoeff[4];
    double fAlphaCoeff[4];
    double fBeCoeff[4];
    double fBCoeff[4];
    double fSmearingCoeff[3];
    
    double fParameters[8];

    int fSimSortMax;
    int fNumEntries;
    TRandom3 fRandom;

    double fOffset;

    double fChi2;
    int fExpBinNum;
    double fExpBinHigh;
    double fExpBinLow;
   
    bool fRebin;
    

};


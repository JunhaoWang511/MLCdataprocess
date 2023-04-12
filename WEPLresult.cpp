#include <iostream>
#include <fstream>
#include <math.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TAxis.h>
using namespace std;

void WEPLres(void)
{
    TCanvas *can = new TCanvas();
    gStyle->SetOptFit(1111);
    double energy[5]={1.004,0.8559,0.6933,0.5108,0.293};
    double energy_error[5]={0.005279,0.007719,0.009154,0.01086,0.01396};
    double energy_error_error[5]={1.2e-5,1.3e-5,1.6e-5,2e-5,3e-5};
    double WEPL[5]={0,50,100,150,200};
    double WEPL_error[5]={0};
    TF1* poly2=new TF1("poly2","[0]*x*x+[1]*x+[2]",0,2);
    poly2->SetParameters(-100,-100,200);
    poly2->SetParLimits(0,-1000,1000);
    poly2->SetParLimits(1,-1000,1000);
    poly2->SetParLimits(2,-1000,1000);
    TGraphErrors* gr1=new TGraphErrors(5,energy,WEPL,energy_error,WEPL_error);
    gr1->SetTitle("Calibration Curve;Energy/AU;WEPL/mm");
    gr1->SetMarkerStyle(8);
    gr1->SetMarkerSize(1);
    gr1->Draw("AP");
    gr1->Fit("poly2","","",0,2);
    can->Print("WEPL-Energy.pdf");
    can->Clear();

    double WEPLres[5]={0};
    double WEPLres_error[5]={0};
    double derivative[5]={0};
    int i;
    for(i=0;i<5;i++){
        derivative[i]=fabs(poly2->Derivative(energy[i]));
        WEPLres[i]=derivative[i]*energy_error[i];
        WEPLres_error[i]=derivative[i]*energy_error_error[i];
    }
    TGraphErrors* gr2=new TGraphErrors(5,WEPL,WEPLres,WEPL_error,WEPLres_error);
    gr2->SetTitle("WEPL-Resolution;WEPL/mm;WEPL-resolution/mm");
    //gr2->GetYaxis()->SetRangeUser(1,4);
    gr2->SetMarkerStyle(7);
    gr2->SetMarkerSize(1);
    gr2->Draw("AP");
    can->Print("WEPL-resolution.pdf");
}
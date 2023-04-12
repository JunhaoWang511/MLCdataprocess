#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

void DrawAverageW() 
{
    TCanvas* can = new TCanvas();
    TF1 f("fun","sqrt(pow(0.011*(260-fmod(260-4-x,260./[0])),2)+pow(1.8*[1]*fmod(260-4-x,260./[0]),2))",0,260);//[0] is nlayer, [1] is deltaE
    ROOT::Math::WrappedTF1 wf(f);
    ROOT::Math::GaussIntegrator ig;
    ig.SetFunction(wf);
    ig.SetRelTolerance(0.0001);

    int i=0;
    double SigW[100]={0};
    double DelE[100]={0};
    for(int i=0; i<100; i++)
    {
        DelE[i]=0.05*i;
    }

    f.SetParameter(0,1);
    for(i=0; i<100; i++)
    {
        f.SetParameter(1,DelE[i]/100);
        SigW[i]=ig.Integral(0,260)/260;
        cout<<SigW[i]<<'\t'<<endl;
    }
    TGraph* gra1 = new TGraph(100,DelE,SigW);
    gra1->SetLineWidth(2);
    gra1->SetLineColor(2);
    gra1->SetLineStyle(1);
    gra1->SetMarkerColor(2);
    gra1->SetTitle("n=1");

    f.SetParameter(0,2);
    for(i=0; i<100; i++)
    {
        f.SetParameter(1,DelE[i]/100);
        SigW[i]=ig.Integral(0,260)/260;
    }
    TGraph* gra2 = new TGraph(100,DelE,SigW);
    gra2->SetLineWidth(2);
    gra2->SetLineColor(3);
    gra2->SetLineStyle(1);
    gra2->SetMarkerColor(3);
    gra2->SetTitle("n=2");

    f.SetParameter(0,3);
    for(i=0; i<100; i++)
    {
        f.SetParameter(1,DelE[i]/100);
        SigW[i]=ig.Integral(0,260)/260;
    }
    TGraph* gra3 = new TGraph(100,DelE,SigW);
    gra3->SetLineWidth(2);
    gra3->SetLineColor(4);
    gra3->SetLineStyle(1);
    gra3->SetMarkerColor(4);
    gra3->SetTitle("n=3");

    f.SetParameter(0,4);
    for(int i=0; i<100; i++)
    {
        f.SetParameter(1,DelE[i]/100);
        SigW[i]=ig.Integral(0,260)/260;
    }
    TGraph* gra4 = new TGraph(100,DelE,SigW);
    gra4->SetLineWidth(2);
    gra4->SetLineColor(5);
    gra4->SetLineStyle(1);
    gra4->SetMarkerColor(5);
    gra4->SetTitle("n=4");

    f.SetParameter(0,5);
    for(int i=0; i<100; i++)
    {
        f.SetParameter(1,DelE[i]/100);
        SigW[i]=ig.Integral(0,260)/260;
    }
    TGraph* gra5 = new TGraph(100,DelE,SigW);
    gra5->SetLineWidth(2);
    gra5->SetLineColor(6);
    gra5->SetLineStyle(1);
    gra5->SetMarkerColor(6);
    gra5->SetTitle("n=5");
    
    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle(";Energy resolution,#deltaE,%;Average WEPL resolution,mm");
    mg->Add(gra1);
    mg->Add(gra2);
    mg->Add(gra3);
    mg->Add(gra4);
    mg->Add(gra5);
    mg->GetYaxis()->SetRangeUser(1,5);
    mg->GetXaxis()->SetRangeUser(0,5);
    mg->Draw("alp");

    TLegend* legend = new TLegend(0.1,0.7,0.3,0.9);
    legend->AddEntry(gra1,"n=1","l");
    legend->AddEntry(gra2,"n=2","l");
    legend->AddEntry(gra3,"n=3","l");
    legend->AddEntry(gra4,"n=4","l");
    legend->AddEntry(gra5,"n=5","l");
    legend->Draw();

    can->Print("Average_WEPL_Resolution.pdf");
}
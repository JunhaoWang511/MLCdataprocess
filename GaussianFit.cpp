//Fit a TH1D from a root file, named as "200MeV_layeri", i is a serial number from 1 to layer number n.
#include<iostream>
#include<fstream>
#include<math.h>
using namespace std;

void GaussianFit(TH1D* his)//
{
    int maxbin,iminus=1,iplus=1;
    double binwidth,maxcont,levelcont,tempcont;
    his->Rebin(10);
    his->Draw();
    maxbin=his->GetMaximumBin();
    maxcont=his->GetBinContent(maxbin);
    binwidth=his->GetBinWidth(0);
    levelcont=0.3*maxcont;
    tempcont=maxcont;
    while((tempcont>levelcont)||(his->GetBinContent(maxbin+iplus)>levelcont))
    {
        tempcont=his->GetBinContent(maxbin+iplus);
        iplus++;
    }
    tempcont=maxcont;
    while((tempcont>levelcont)||(his->GetBinContent(maxbin-iminus)>levelcont))
    {
        tempcont=his->GetBinContent(maxbin-iminus);
        iminus++;
    }
    his->Fit("gaus","","",binwidth*(maxbin-iminus),binwidth*(maxbin+iplus));
}

void LayerFit(TFile* infile,int nlayer)
{
    TCanvas* can=new TCanvas();
    gStyle->SetOptFit(1111);
    TH1D* his; 
    int i=1;
    string name="200MeV_layer";
    string namei;
    string filename=std::to_string(nlayer)+"layer_";
    string filenamei;
    for(i;i<=1;i++)
    {
        namei=name+std::to_string(i)+";1";
        filenamei=filename+name+std::to_string(i)+".pdf";
        his=(TH1D*)infile->Get(namei.c_str());
        his->GetXaxis()->SetRangeUser(0,300000);
        GaussianFit(his);
        can->Print(filenamei.c_str());
    }
}

/*void PositionFit(TFile* infile)
{
    TCanvas* can=new TCanvas();
    gStyle->SetOptFit(1111);
    TH1D* his; 
    int i;
    string name="MeV stop position";
    for(i=0;i<10;i++)
    {
        string namei="";
        string filename="";
        namei=std::to_string((i+1)*20)+name;
        filename=namei+".pdf";
        his=(TH1D*)infile->Get(namei.c_str());
        GaussianFit(his);
        can->Print(filename.c_str());
    }
}*/


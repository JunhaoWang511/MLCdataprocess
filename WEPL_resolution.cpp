#include<iostream>
#include<fstream>
#include<math.h>
using namespace std;

void DrawSigmaW(double* deltaE)
{
    TCanvas* can = new TCanvas();
    TF1* f1 = new TF1("fun1","sqrt(pow(0.011*(260-fmod(260-4-x,260./[0])),2)+pow(1.8*[1]*fmod(260-4-x,260./[0]),2))",0,260);//[0] is nlayer, [1] is deltaE
    f1->SetParameter(0,1);
    f1->SetParameter(1,deltaE[0]);
    f1->SetLineWidth(1);
    f1->SetLineColor(2);
    f1->SetLineStyle(1);

    TF1* f2 = new TF1("fun2","sqrt(pow(0.011*(260-fmod(260-4-x,260./[0])),2)+pow(1.8*[1]*fmod(260-4-x,260./[0]),2))",0,260);
    f2->SetParameter(0,2);
    f2->SetParameter(1,deltaE[1]);
    f2->SetLineWidth(1);
    f2->SetLineColor(3);
    f2->SetLineStyle(2);

    TF1* f3 = new TF1("fun3","sqrt(pow(0.011*(260-fmod(260-4-x,260./[0])),2)+pow(1.8*[1]*fmod(260-4-x,260./[0]),2))",0,260);
    f3->SetParameter(0,5);
    f3->SetParameter(1,deltaE[2]);
    f3->SetLineWidth(1);
    f3->SetLineColor(4);
    f3->SetLineStyle(3);

    /*TF1* f4 = new TF1("fun4","sqrt(pow(0.011*(260-fmod(260-4-x,260./[0])),2)+pow(1.8*[1]*fmod(260-4-x,260./[0]),2))",0,260);
    f1->SetParameter(0,);
    f1->SetParameter(1,deltaE[3]);
    TF1* f5 = new TF1("fun5","sqrt(pow(0.011*(260-fmod(260-4-x,260./[0])),2)+pow(1.8*[1]*fmod(260-4-x,260./[0]),2))",0,260);
    f5->SetParameter(0,);
    f5->SetParameter(1,deltaE[4]);*/

    f1->Draw();
    f2->Draw("same");
    f3->Draw("same");
    f1->SetTitle(";WEPL,mm;WEPL resolution,mm");
    f1->GetYaxis()->SetRangeUser(1,4);

    auto legend = new TLegend(0.6,0.7,0.9,0.9);
    legend->AddEntry("fun1","n=1, #deltaE = 0.57%","l");
    legend->AddEntry("fun2","n=2, #deltaE = 0.67%","l");
    legend->AddEntry("fun3","n=5, #deltaE = 1.05%","l");
    legend->Draw();

    can->Print("Predicted_WEPL_Resolution.pdf");
}

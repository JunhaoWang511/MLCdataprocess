// correction: divide scintillator surface into 100*100 bins, sum and average each bin's light collection(2D), include each bins' distribution(3D).
#include <iostream>
#include <fstream>
#include <math.h>
#include <TROOT.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF2.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TFile.h>
using namespace std;
//add a change
#define NUM1 400 // bin number
double correct_factor[NUM1-76]={33406.7,33353.3,33591.3,33569.1,33451.8,32588.3,33651.6,32113,32825.6,32467.5,33268.5,32931,33555.7,32505.6,32906.9,32877.4,33820.6,31535.4,33378.8,33360.4,33127.3,33006.8,32144.7,32900.1,33792.7,33609.9,32436.7,32960.9,33853.7,33251.4,33228.7,33969.4,33478.7,33401.1,32170.3,33064.3,33436.6,33389.2,33317.1,32505.8,32760.1,32837.8,33052.2,32681,33330.1,33653.7,32260.1,32784.9,32686.9,33515.4,33239.1,33629.2,33127.9,33383.7,33011.5,33212.2,32924.6,33277.1,33689.9,32888.8,32520.4,32782.7,33469.3,32226.6,33826.4,33331.6,33302.8,33065.1,33079.9,32454.7,33028.5,33099.6,33478.5,33171.3,32490.6,33046.1,33076.2,33511,32835.2,33621.9,33993.2,32817,33976.4,33175,34067.5,33075.9,33853,33905,33312.2,32597.8,33095,33140.5,34071.6,33418,33918.9,33673.9,33363.3,33673.5,33024.6,32907.4,34152.9,34489.5,33853.9,33002.7,33514.4,32669.6,33323,34201.3,33266.5,33421.7,34557.8,33613.5,34496.3,32944.4,34051.5,34197,34608.7,33387.4,32679.8,33465.4,32892.8,33668.4,33837.6,34460.5,33252.7,33520.6,34907.8,34466.9,35064.8,33406.8,33666.3,33617.8,34278,33217.3,32536.1,34375,33863.1,34809.2,34114.1,34414.5,34381.4,33164.1,34062.6,34006.2,33854.3,34456.8,33917.7,34290.1,34120.3,34504.2,33716.7,34677.5,34512.9,34825.9,33577.5,34226.9,34871.3,32931.1,33935.8,34330.1,34302.2,35048,34628.1,34510,35190.6,34969,34368.4,33641.2,34748,36014.8,33174.9,35728,34981.1,34376,34645.1,34171.3,35625.3,34807.5,34013.6,34048.6,35825.2,33567,35663.2,34914.6,34750.4,34744.8,36016.6,34396.1,34496.3,35147.7,33820.7,35361.9,35687.2,34166.4,35620.8,36048.6,34712.1,35194.2,35375.8,35503.8,36090.7,34040.7,35636.7,35158.4,35433.2,35262.7,35431.1,34948.8,34777.7,35642.2,35464.6,35116.3,34720.7,35744.8,35272.8,35211.2,36612.7,36438,35650.4,37041.1,34873,35519.5,35588.3,36728.3,35361.2,36237.7,35949.7,36123.2,36707.7,34700.2,36884.4,35856.9,36675.7,36292.6,36891.9,36925.5,36000,36463.4,37251.4,36362,37260.4,35743.1,35611,36713.1,37494.3,37048.4,36586.5,36735.9,36603.7,37057.4,36058.4,36189.7,36392.7,37852.1,37718.9,37123.4,37976.6,37249,36925.8,37057.2,35502.5,37403.4,38456.7,38253.5,36612.3,37284,36913.4,37863.9,37001.5,37046.3,38157.7,38561.4,38667.1,37840,37423.7,38085.6,38094.4,37173.6,38814.3,38819.4,38175.6,37096.8,37459.7,37106.8,37318.3,38070.7,37837.6,38463.8,39840.2,39078.8,40006,39016.2,38758.8,39996.6,39157.9,39390.6,39863.3,38032,39791.3,39166.4,38531.1,39503,39565,38762.4,40001.3,39898.2,40221.3,39734.6,40896.3,39936.5,41352.5,40596,39739.9,41225.4,42032,40579.1,41083.1,40182.7,39759.7,40270.2,40928.4,40490.2,41435.8,40061.1};//use each bins' average photon number as correct factor(using 200 MeV data).
// #define NUM2 50000   //event number
// along Z is x-axis,along X is y-axis
int BinSerial(double x1, double z1, double x2, double z2) // input Pos3 and Pos4, return which bin proton will hit.(1-400)
{
    double posX, posZ;
    int serial;
    posX = 1.35 * x2 - 0.35 * x1;
    posZ = 1.35 * z2 - 0.35 * z1;
    if (fabs(posX) > 5 || fabs(posZ) > 5)
        return (0);
    serial = floor((posZ + 5.) / (10. / sqrt(NUM1))) * sqrt(NUM1) + ceil((posX + 5.) / (10. / sqrt(NUM1)));
    return (serial);
}
int BinSerial1(double x1, double z1, double x2, double z2) // remove bins along edge(1-324)
{
    double posX, posZ;
    int serial;
    posX = 1.35 * x2 - 0.35 * x1;
    posZ = 1.35 * z2 - 0.35 * z1;
    if (fabs(posX) > 4.5 || fabs(posZ) > 4.5)
        return (0);
    serial = floor((posZ + 4.5) / (10. / sqrt(NUM1))) * (sqrt(NUM1)-2) + ceil((posX + 4.5) / (10. / sqrt(NUM1)));
    return (serial);
}

double *AverageSignal(TTree *tr) // return a double array containing average light collection of each bin.
{
    int EventNum = tr->GetEntries();
    cout << EventNum << endl;
    double x1, x2, z1, z2;
    int photon;
    long long sumSignal[NUM1] = {0};
    int hitnumber[NUM1] = {0}, i, bin;
    double *aveSignal = new double[NUM1];
    int hit1 = 0, hit2 = 0, hit3 = 0, hit = 0;

    TBranch *br1, *br2, *br3, *br4, *br5;
    br1 = tr->GetBranch("Pos3X");
    br2 = tr->GetBranch("Pos3Z");
    br3 = tr->GetBranch("Pos4X");
    br4 = tr->GetBranch("Pos4Z");
    br5 = tr->GetBranch("Photon1");

    br1->SetAddress(&x1);
    br2->SetAddress(&z1);
    br3->SetAddress(&x2);
    br4->SetAddress(&z2);
    br5->SetAddress(&photon);

    for (i = 0; i < EventNum; i++)
    {
        br1->GetEntry(i);
        br2->GetEntry(i);
        br3->GetEntry(i);
        br4->GetEntry(i);
        br5->GetEntry(i);

        if (x1 == 0. && x2 == 0. && z1 == 0. && z2 == 0.)
        {
            hit1++;
            // cout<<x1<<endl;
            continue;
        }
        if (sqrt(pow((x2 - x1), 2) + pow((z2 - z1), 2)) / 5.0 > 0.1)
        {
            hit2++;
            continue;
        }
        bin = BinSerial(x1, z1, x2, z2);
        if (bin == 0)
        {
            hit3++;
            continue;
        }
        else
        {
            sumSignal[bin - 1] += photon;
            hitnumber[bin - 1]++;
            hit++;
        }
    }

    for (i = 0; i < NUM1; i++)
    {
        if (hitnumber[i] != 0)
        {
            aveSignal[i] = double(sumSignal[i]) / double(hitnumber[i]);
        }
        else
            aveSignal[i] = 0;
        // cout<<aveSignal[i]<<'\t';
        // cout << hitnumber[i] << '\t';
    }
    // cout << "total hit number is hit1=" << hit1 << " + hit2=" << hit2 << " + hit3=" << hit3 << " + hit=" << hit << " =  "<< hit1 + hit2 + hit3 + hit << endl;

    return aveSignal;
}

double position_rms(TProfile2D *pro)
{
    double photon[NUM1 - 76] = {0};
    int i, j, bin = 0;
    for (i = 1; i < sqrt(NUM1) - 1; i++)
        for (j = 1; j < sqrt(NUM1) - 1; j++)
        {
            photon[bin] = pro->GetBinContent(i, j);
            bin++;
        }
    return TMath::RMS(NUM1 - 76, photon) / TMath::Mean(NUM1 - 76, photon);
}
double position_rms1(TH2D *pro)
{
    double photon[NUM1 - 76] = {0};
    int i, j, bin = 0;
    for (i = 1; i < sqrt(NUM1) - 1; i++)
        for (j = 1; j < sqrt(NUM1) - 1; j++)
        {
            photon[bin] = pro->GetBinContent(i, j);
            bin++;
        }
    return TMath::RMS(NUM1 - 76, photon) / TMath::Mean(NUM1 - 76, photon);
}

void DrawSignal2D(TFile *infile) // draw average signal distribution
{
    TTree *tr = (TTree *)infile->Get("PSDposition_PMTcollection");
    TCanvas *can = new TCanvas();
    TF2 *f = new TF2("fit", "[0]*x^2+[1]*y^2+[2]*x*y+[3]*x+[4]*y+[5]", -5, 5, -5, 5);
    // TF2 *f = new TF2("f2", "[0]*x^2+[1]*x+[2]", -5, 5, -5, 5);
    TF2 *fscale = new TF2("fscale", "4.56e+02*x^2-2.65e+02*y^2+1.45e+01*x*y+4.11e+03*x+4.58e-02*y+2.63e+05", -5, 5, -5, 5);
    double width = 5 - 10. / sqrt(NUM1);
    double error = 0;
    int errorbin = 0;
    TH2D *his = new TH2D("his", "AverageSignal", sqrt(NUM1) - 2, -width, width, sqrt(NUM1) - 2, -width, width);
    double *x = AverageSignal(tr);
    for (int i = 1; i < sqrt(NUM1) - 1; i++)
        for (int j = 1; j < sqrt(NUM1) - 1; j++)
        {
            // his->Fill((-5 + 5. / sqrt(NUM1)) + 10. / sqrt(NUM1) * i, (-5 + 5. / sqrt(NUM1)) + 10. / sqrt(NUM1) * j, x[int(sqrt(NUM1)) * i + j] / fscale->Eval((-5 + 5. / sqrt(NUM1)) + 10. / sqrt(NUM1) * i, (-5 + 5. / sqrt(NUM1)) + 10. / sqrt(NUM1) * j));

            his->Fill((-5 + 5. / sqrt(NUM1)) + 10. / sqrt(NUM1) * i, (-5 + 5. / sqrt(NUM1)) + 10. / sqrt(NUM1) * j, x[int(sqrt(NUM1)) * i + j]);
            // cout<<x[int(sqrt(NUM1)) * i + j]<<'\t';
        }
    // cout << "total error is " << error / errorbin << endl;
    his->Fit("fit");
    //  his->Draw("COLZ");
    his->Draw("LEGO2");
    // cout<<"RMS is "<<position_rms1(his);
    // cout<<"bin number is "<<his->GetNbinsX()<<endl;
    //  can->Print("distribution.pdf");
}
void DrawSignal3D(TFile *infile)
{
    int binnum=0;
    TCanvas *can = new TCanvas();
    TTree *tr = (TTree *)infile->Get("PSDposition_PMTcollection");
    int EventNum = tr->GetEntries();
    cout << EventNum << endl;
    double x1, x2, z1, z2;
    int photon;
    // long long sumSignal[NUM1] = {0};
    int hitnumber[NUM1] = {0}, i, j, bin, num = 0;
    // double *aveSignal = new double[NUM1];
    int hit1 = 0, hit2 = 0, hit3 = 0, hit4 = 0, hit = 0;
    double posX, posZ;

    TBranch *br1, *br2, *br3, *br4, *br5;
    br1 = tr->GetBranch("Pos3X");
    br2 = tr->GetBranch("Pos3Z");
    br3 = tr->GetBranch("Pos4X");
    br4 = tr->GetBranch("Pos4Z");
    br5 = tr->GetBranch("Photon5");

    br1->SetAddress(&x1);
    br2->SetAddress(&z1);
    br3->SetAddress(&x2);
    br4->SetAddress(&z2);
    br5->SetAddress(&photon);

    double width = 5 - 10. / sqrt(NUM1);
    // TH3D *his = new TH3D("his", "SignalDis", sqrt(NUM1) - 2, -width, width, sqrt(NUM1) - 2, -width, width, 1000000, 0, 1000000);
    TProfile2D *prohis = new TProfile2D("his", "SignalDis", sqrt(NUM1) - 2, -width, width, sqrt(NUM1) - 2, -width, width, 0, 1000000);

    //TF2 *fscale = new TF2("fscale", "4.56497e+02*x^2-2.67454e+02*y^2+1.87243e+01*x*y+4.10670e+03*x+-2.05897e+00*y+2.63454e+05", -5, 5, -5, 5);
    TF2 *fscale = new TF2("fscale", "1.37079e+02*x^2-3.83890e+01*y^2+6.23422*x*y+1.10124e+03*x-7.82735*y+4.63912e+04", -5, 5, -5, 5); // Fit result: require photon> 250000
    TH1I *phodis = new TH1I("pho", "pho", 10000, 0, 500000);
    phodis->SetTitle("photon distribution for 200 MeV");
    TH1I *phodisbin = new TH1I("phob", "phob", 5000, 0, 500000);
    phodisbin->SetTitle("photon distribution in a bin");


    for (i = 0; i < EventNum; i++)
    {
        tr->GetEntry(i);
        if (photon < 1000)
        {
            continue;
        }
        phodis->Fill(photon);

        if (x1 == 0. && x2 == 0. && z1 == 0. && z2 == 0.)
        {
            hit1++;
            // cout<<x1<<endl;
            continue;
        }
        if (sqrt(pow((x2 - x1), 2) + pow((z2 - z1), 2)) / 5.0 > 0.1)
        {
            hit2++;
            continue;
        }
        bin = BinSerial1(x1, z1, x2, z2);
        if (bin == 0)
        {
            hit3++;
            continue;
        }
        /*if (photon < 30000)//select events with photon>250000
        {
            hit4++;
            continue;
        }*/
        else
        {
            if(bin==100){phodisbin->Fill(photon);}
            posX = 1.35 * x2 - 0.35 * x1;
            posZ = 1.35 * z2 - 0.35 * z1;
            // sumSignal[bin - 1] += photon;
            hitnumber[bin - 1]++;
            //prohis->Fill(posZ, posX, photon/fscale->Eval(posZ,posX));
            prohis->Fill(posZ, posX, photon);
            //prohis->Fill(posZ, posX, photon/correct_factor[bin-1]);
            hit++;
        }
    }
    //phodisbin->Draw();
    // cout<<hit4<<endl;
    //phodis->Draw();
    prohis->Draw("lego2");

    //TF2 *f = new TF2("fit", "[0]*x^2+[1]*y^2+[2]*x*y+[3]*x+[4]*y+[5]", -5, 5, -5, 5);
    //prohis->Fit("fit");

    /*for (i = 1; i < 19; i++)
    {
        for (j = 1; j < 19; j++)
        {
            correct_factor[binnum]= prohis->GetBinContent(i, j);
            cout<<correct_factor[binnum]<<",";
            binnum++;
        }

    }*/
    //write each bins' signal to correct facor. 

    cout << "RMS is " << position_rms(prohis)<<endl;
}

void cor_dis(TFile *infile) // draw corrected photon distribution, and fit with gaussian function to get mean and sigma.
{
    //TCanvas *can = new TCanvas();
    TTree *tr = (TTree *)infile->Get("PSDposition_PMTcollection");
    int EventNum = tr->GetEntries();
    cout << EventNum << endl;
    double x1, x2, z1, z2;
    int photon;
    // long long sumSignal[NUM1] = {0};
    int hitnumber[NUM1] = {0}, i, j, bin, num = 0;
    // double *aveSignal = new double[NUM1];
    int hit1 = 0, hit2 = 0, hit3 = 0, hit = 0;
    double posX, posZ, cor_signal;
    double corr_factor[NUM1] = {0};

    // TF2 *fscale = new TF2("fscale", "4.56497e+02*x^2-2.67454e+02*y^2+1.87243e+01*x*y+4.10670e+03*x+-2.05897e+00*y+2.63454e+05", -5, 5, -5, 5);
    TF2 *fscale = new TF2("fscale", "1.37079e+02*x^2-3.83890e+01*y^2+6.23422*x*y+1.10124e+03*x-7.82735*y+4.63912e+04", -5, 5, -5, 5);//require photon>250000
    for (int i = 1; i < sqrt(NUM1) - 1; i++)
        for (int j = 1; j < sqrt(NUM1) - 1; j++)
        {
            corr_factor[int(sqrt(NUM1)) * i + j] = fscale->Eval((-5 + 5. / sqrt(NUM1)) + 10. / sqrt(NUM1) * i, (-5 + 5. / sqrt(NUM1)) + 10. / sqrt(NUM1) * j);
        }//this method (corr_factor) is abandoned, I apply correct_factor

    TBranch *br1, *br2, *br3, *br4, *br5;
    br1 = tr->GetBranch("Pos3X");
    br2 = tr->GetBranch("Pos3Z");
    br3 = tr->GetBranch("Pos4X");
    br4 = tr->GetBranch("Pos4Z");
    br5 = tr->GetBranch("Photon1");

    br1->SetAddress(&x1);
    br2->SetAddress(&z1);
    br3->SetAddress(&x2);
    br4->SetAddress(&z2);
    br5->SetAddress(&photon);

    TH1D *his = new TH1D("his", "corrected_distribution", 10000, 0, 2);
    // TH2D *his2D = new TH2D("his2D", "test", 200, 2.5e5, 4e5, 200, -5, 5);

    for (i = 0; i < EventNum; i++)
    {
        tr->GetEntry(i);

        if (x1 == 0. && x2 == 0. && z1 == 0. && z2 == 0.)
        {
            hit1++;
            // cout<<x1<<endl;
            continue;
        }
        if (sqrt(pow((x2 - x1), 2) + pow((z2 - z1), 2)) / 5.0 > 0.1)
        {
            hit2++;
            continue;
        }
        bin = BinSerial1(x1, z1, x2, z2);
        if (bin == 0)
        {
            hit3++;
            continue;
        }
        else
        {
            posX = 1.35 * x2 - 0.35 * x1;
            posZ = 1.35 * z2 - 0.35 * z1;
            hitnumber[bin - 1]++;
            cor_signal = photon /correct_factor[bin-1];

            his->Fill(cor_signal);
            // his2D->Fill(cor_signal, photon);
            //his2D->Fill(photon, posZ);
            hit++;
        }
    }
    his->Draw();
    // his2D->Draw("COLZ");
}

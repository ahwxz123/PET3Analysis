#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
//#include "io.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TDatime.h"
#include "TF1.h"
#include "TMath.h"
using namespace std;
const double PI=3.1415926;

const unsigned int Nslotbyte=64;
const unsigned int NPixels=29184;
const unsigned int Nbin=1024;
const int NCX=608;
const int NZX=48;
const int NZperMod=24;
const int NCperMod=16;
string FILELIST="FILELIST.txt";

TH1D *hist1dD(char *histname,int Xstep,double xlow,double xup,char*xtitle,char *ytitle);
TH1F *hist1d(char *histname,int Xstep,double xlow,double xup,char*xtitle,char *ytitle);
TH1S *hist1dS(char *histname, Int_t Xstep,Double_t xlow,Double_t xup,char*xtitle,char *ytitle);
TH2S *hist2dS(char *hist2dname,int Xstep,double xlow,double xup,int Ystep,double ylow,double yup,char *xtitle,char *ytitle);
TH2F *hist2d(char *hist2dname,int Xstep,double xlow,double xup,int Ystep,double ylow,double yup,char *xtitle,char *ytitle);
void MycanvasSetting(TCanvas *cc,double leftmargin,double rightmargin,double topmargin,double bottommargin); 
void MyPadSetting(TVirtualPad *cc,double leftmargin,double rightmargin,double topmargin,double bottommargin);
void HistPFX1dSetting(  TH1F *hh,  char *hist1dname, double xlow, double xup, double ylow, double yup,  char *xtitle,   char *ytitle );
unsigned int dataconvert(unsigned int  data);
void getFiles(string filepath, vector<string>& files);

int main(int argc,char **argv){

	 string filepath=argv[1];
     gStyle->SetOptStat(0);
     gStyle->SetOptFit(11);
     TStopwatch *ts1=new TStopwatch();
     ts1->Start();
      TDatime *tdate =new TDatime();
      int systime = tdate->GetTime();
      int sysdate = tdate->GetDate();
     
      char buf[100];
      int index1,index2,bankindex1,bankindex2;
 
      ULong64_t data=0;
      ULong64_t dataeff1[2];      
      unsigned int dataeff[2];
      int Nshift=0;
      unsigned int Nslotbyte=64;      
      int Dtflag,multiflag,energy,pos;
      int finetime,coarsetime,reserved;
      int Ovflow;
      int axial,transaxial;
	  int X,Y,E;
	  int coarse,fine;
	  double timetick;
      double coarseLSB=10.;  //ns;
      double fineLSB=0.3125;  //ns;
      int kind,tdiff,Ax1,Tx1,Ax2,Tx2,Adc1,Adc2;
      int index;
      double dtemp,cc;

      char countmapFile[100],energyvsringFile[100],energyvsringpfxFile[100];
      char energyvsangleFile[100],energyvsanglepfxFile[100],ringangleFile[100],timeFile[100],energypixelFile[100];
      char slicesFile[100],t1addt2File[100],t1difft2File[100];
      char peakandresolutionmapFile[100];
      char energybankFile1[100],energybankFile2[100],energybankFile3[100];
      char pixelsenergyFile1[100],pixelsenergyFile2[100],pixelsenergyFile3[100];
      sprintf(countmapFile, "%s/Countsmap_%d.pdf",filepath.c_str(),sysdate);
      sprintf(energyvsringFile, "%s/EnergyvsRing_%d.pdf",filepath.c_str(),sysdate);
      sprintf(energyvsringpfxFile, "%s/EnergyvsRing_PFX_%d.pdf",filepath.c_str(),sysdate);
      sprintf(energyvsangleFile, "%s/EnergyvsAngle_%d.pdf",filepath.c_str(),sysdate);
      sprintf(energyvsanglepfxFile, "%s/EnergyvsAngle_PFX_%d.pdf",filepath.c_str(),sysdate);
      sprintf(ringangleFile, "%s/RingMap&AngleMap_%d.pdf",filepath.c_str(),sysdate);
      sprintf(timeFile, "%s/Time_Spectrum_%d.pdf",filepath.c_str(),sysdate);
      sprintf(slicesFile, "%s/Slice_Spectrum_%d.pdf",filepath.c_str(),sysdate);
      sprintf(t1addt2File, "%s/T1ADDT2_Spectrum_%d.pdf",filepath.c_str(),sysdate);
      sprintf(t1difft2File, "%s/T1DIFFT2_Spectrum_%d.pdf",filepath.c_str(),sysdate);
      sprintf(peakandresolutionmapFile, "%s/PeakMap&ResolutionMap_%d.pdf",filepath.c_str(),sysdate);
      sprintf(energypixelFile, "%s/EnergyvsPixelsNumber_%d.png",filepath.c_str(),sysdate);
      sprintf(energybankFile1, "%s/BankEnergy_%d.pdf[",filepath.c_str(),sysdate);
      sprintf(energybankFile2, "%s/BankEnergy_%d.pdf",filepath.c_str(),sysdate);
      sprintf(energybankFile3, "%s/BankEnergy_%d.pdf]",filepath.c_str(),sysdate);
      sprintf(pixelsenergyFile1, "%s/PixelsEnergy_%d.pdf[",filepath.c_str(),sysdate);
      sprintf(pixelsenergyFile2, "%s/PixelsEnergy_%d.pdf",filepath.c_str(),sysdate);
      sprintf(pixelsenergyFile3, "%s/PixelsEnergy_%d.pdf]",filepath.c_str(),sysdate);


      char countmapPNG[100],energyvsringPNG[100],energyvsringpfxPNG[100];
      char energyvsanglePNG[100],energyvsanglepfxPNG[100],ringanglePNG[100],timePNG[100],energypixelPNG[100];
      char slicesPNG[100],t1addt2PNG[100],t1difft2PNG[100];
      char peakandresolutionmapPNG[100];
	  char energypixelFile1[100],energypixelFile2[100],energypixelFile3[100];
      //char energybankPNG1[100],energybankPNG2[100],energybankPNG3[100];
      sprintf(countmapPNG, "%s/Countsmap_%d.png",filepath.c_str(),sysdate);
      sprintf(energyvsringPNG, "%s/EnergyvsRing_%d.png",filepath.c_str(),sysdate);
      sprintf(energyvsringpfxPNG, "%s/EnergyvsRing_PFX_%d.png",filepath.c_str(),sysdate);
      sprintf(energyvsanglePNG, "%s/EnergyvsAngle_%d.png",filepath.c_str(),sysdate);
      sprintf(energyvsanglepfxPNG, "%s/EnergyvsAngle_PFX_%d.png",filepath.c_str(),sysdate);
      sprintf(ringanglePNG, "%s/RingMap&AngleMap_%d.png",filepath.c_str(),sysdate);
      sprintf(timePNG, "%s/Time_Spectrum_%d.png",filepath.c_str(),sysdate);
      sprintf(slicesPNG, "%s/Slice_Spectrum_%d.png",filepath.c_str(),sysdate);
      sprintf(t1addt2PNG, "%s/T1ADDT2_Spectrum_%d.png",filepath.c_str(),sysdate);
      sprintf(t1difft2PNG, "%s/T1DIFFT2_Spectrum_%d.png",filepath.c_str(),sysdate);
      sprintf(peakandresolutionmapPNG, "%s/PeakMap&ResolutionMap_%d.png",filepath.c_str(),sysdate);
      sprintf(energypixelPNG, "%s/EnergyvsPixelsNumber_%d.png",filepath.c_str(),sysdate);
      //sprintf(energybankPNG1, "BankEnergy.png[");
      //sprintf(energybankPNG2, "BankEnergy.png");
      //sprintf(energybankPNG3, "BankEnergy.png]");
      sprintf(energypixelFile1, "%s/PixelEnergy.pdf[",filepath.c_str());
      sprintf(energypixelFile2, "%s/PixelEnergy.pdf",filepath.c_str());
      sprintf(energypixelFile3, "%s/PixelEnergy.pdf]",filepath.c_str());



      TF1 *f1=new TF1("f1","gaus",-1e5,1e5);

      TH2F *CountsMap;
      TH1F *AngleMap;
      TH1F *RingMap;
      TH2F *EnergyvsRing;
      TH2F *EnergyvsAngle;
      TH1S *Energy[NPixels];
      TH1F *Slices;
	  TH1F *T1ADDT2;
	  TH1F *T1DIFFT2;
	  TH1D *TIME;
	  TH1D *EnergyBank[2][38];

      TH1F *EnergyvsRing_PFX;
      TH1F *EnergyvsAngle_PFX;

      TH2F *EnergyvsPixel;
      TH2F *PeakMap, *ResolutionMap;
      Slices = hist1d("Sensitivity",96,0,96,"silce","Entries");
	  T1ADDT2= hist1d("(T1+T2)/2",1216,0,1216,"(T1+T2)/2.","Entries");
	  T1DIFFT2= hist1d("|T1-T2|",608,0,608,"|T1-T2|","Entries");
	  TIME =hist1dD("TIME",100,-300,300,"TimeDiff (ps)","Entries");
      PeakMap = hist2d("PeakMap",608,0,608,48,0,48,"TranAxial","Axial");
      ResolutionMap = hist2d("ResolutionMap",608,0,608,48,0,48,"TranAxial","Axial");
      
      
     sprintf(buf,"CountsMap");
     CountsMap= hist2d(buf,608,0,608,48,0,48,"TransAxial","Axial");
     sprintf(buf,"RingMap");
     RingMap= hist1d(buf,48,0,48,"RingSector","Entries");
     RingMap->SetLineColor(2);
     sprintf(buf,"AngleMap");
     AngleMap= hist1d(buf,608,0,608,"AngleSector","Entries");
     AngleMap->SetLineColor(2);
     sprintf(buf,"EnergyvsRing");
     EnergyvsRing= hist2d(buf,48,0,48,500,10,1010,"RingSector","Energy");
     sprintf(buf,"EnergyvsAngle");
     EnergyvsAngle= hist2d(buf,608,0,608,500,10,1010,"AngleSector","Energy");
     for(int j=0 ;j <NPixels; j++){
        sprintf(buf,"Energy_A%d_T%d",j/608,j%608);
        Energy[j]= hist1dS(buf,500,10,1010,"Energy","Entries"); 
		Energy[j]->SetTitle(buf);
     }

     EnergyvsPixel = hist2d("EnergyvsPixel",NPixels,0,NPixels,500,10,1010,"PixelNo.","Entries");

     for(int i=0 ;i<2;i++){
        for(int j=0 ;j<38;j++){
            sprintf(buf,"EnergyBank_%d_%d",i,j);
            EnergyBank[i][j] = hist1dD(buf,500,10,1010,"Energy","Entries");
        }
     }


     string txtpath = filepath+"/"+"FILELIST.txt";
     ifstream fintest(txtpath.c_str());
	 if(fintest.fail()) {
        cout<<"  No 'FILELIST.txt' in current directory !!"<<endl;
        return 1;
     }
     fintest.close();

     vector<string> files;
    
     getFiles(txtpath,files);
	 int Nfiles=files.size();
     //for(int i=0;i<Nfiles;i++)cout<< files[i].c_str()<<" ------"<<endl;

     string rootpath = filepath +"/"+"Coincidences_raw.root" ;
     TFile *rfile=new TFile(rootpath.c_str(),"recreate");
     TTree *coin = new TTree("coin","coin");
     coin->Branch("kind"    ,&kind    ,"kind/I"    );
     coin->Branch("tdiff"   ,&tdiff   ,"tdiff/I"   );
     coin->Branch("Adc1"    ,&Adc1    ,"Adc1/I"    );
     coin->Branch("Adc2"    ,&Adc2    ,"Adc2/I"    );
     coin->Branch("Tx1"     ,&Tx1     ,"Tx1/I"     );
     coin->Branch("Ax1"     ,&Ax1     ,"Ax1/I"     );
     coin->Branch("Tx2"     ,&Tx2     ,"Tx2/I"     );
     coin->Branch("Ax2"     ,&Ax2     ,"Ax2/I"     );



     ifstream fin; 
for(int MM=0;MM<Nfiles;MM++){
    //sprintf(buf,"%s/%s%d.bin",filepath.c_str(),filename,MM);
    cout<< files[MM].c_str()<<" ------"<<endl;
	//fin.open( (filepath+"/"+files[MM]).c_str(),ios::binary);
	fin.open( files[MM].c_str(),ios::binary);
	fin.seekg( 0,ios::beg);
    while(fin.good()){
       
         fin.read((char*)(&dataeff),2*sizeof(unsigned int));
//         fin>>dataeff[0]>>dataeff[1];
         
         //cout<<hex<<dataeff[0]<<"   "<<dataeff[1]<<"   "<<endl;
         dataeff1[0] = dataconvert(dataeff[1]);
         dataeff1[1] = dataconvert(dataeff[0]);

         //cout<<"   "<<((data>>32)&0xFFFFFFFF)<<endl;
         //cout<<hex<<dataeff1[1]<<"   "<<dataeff1[0]<<"   "<<endl;
         data = ((dataeff1[1] )<<32) | (dataeff1[0] &0xFFFFFFFF);
         //cout<<data<<endl;
         kind         = (data>>63) & 0x1;
         tdiff        = (data>>52) & 0x7FF;
         Adc1         = (data>>42) & 0x3FF;
         Adc2         = (data>>32) & 0x3FF;
         Tx1          = (data>>22) & 0x3FF;
         Ax1          = (data>>16) & 0x3F;
         Tx2          = (data>>6)  & 0x3FF;
         Ax2          = (data)     & 0x3F;
//         coin->Fill();
			
//         cout<<"  X "<<X<<"   Y "<<Y<<"   E "<<E<<endl;
         if(kind==0)continue;
         if(Adc1<0||Adc1>=Nbin||Adc2<0||Adc2>=Nbin)continue;
		 if(Tx1<0||Tx1>=NCX||Tx2<0||Tx2>=NCX)continue;
		 if(Ax1<0||Ax1>=NZX||Ax2<0||Ax2>=NZX)continue; 
  	
	     CountsMap->Fill(Tx1,Ax1);
	     CountsMap->Fill(Tx2,Ax2);
	     RingMap->Fill(Tx1);
	     RingMap->Fill(Tx2);
	     AngleMap->Fill(Ax1);  
	     AngleMap->Fill(Ax2);  
	     
	     EnergyvsRing->Fill(Ax1,Adc1);
	     EnergyvsRing->Fill(Ax2,Adc2);
	     EnergyvsAngle->Fill(Tx1,Adc1);
	     EnergyvsAngle->Fill(Tx2,Adc2);
	   
	     bankindex1= (Ax1/NZperMod);
	     bankindex2= Tx1/NCperMod;
	     EnergyBank[bankindex1][bankindex2] ->Fill(Adc1);
	     bankindex1= (Ax2/NZperMod);
	     bankindex2= Tx2/NCperMod;
	     EnergyBank[bankindex1][bankindex2] ->Fill(Adc2);
	    }  

	    fin.close();
		cout<<"file--->"<<files[MM].c_str()<<"--has been read done !!"<<endl;
    }

    rfile->cd();
    coin->Write();
    rfile->Close();

    TCanvas *c1=new TCanvas("c1","c1",1200,800);
    TCanvas *c2=new TCanvas("c2","c2",1200,800);
    TCanvas *c3=new TCanvas("c3","c3",1200,800);
    TCanvas *c4=new TCanvas("c4","c4",1200,800);
    TCanvas *c5=new TCanvas("c5","c5",1200,800);
    TCanvas *c6=new TCanvas("c6","c6",1600,800);
    TCanvas *c7=new TCanvas("c7","c7",1000,800);

    MycanvasSetting(c1,0.12,0.14,0.12,0.12);
    MycanvasSetting(c2,0.12,0.14,0.12,0.12);
    MycanvasSetting(c3,0.12,0.14,0.12,0.12);
    MycanvasSetting(c4,0.12,0.14,0.12,0.12);
    MycanvasSetting(c5,0.12,0.14,0.12,0.12);
    MycanvasSetting(c6,0.12,0.14,0.12,0.12);
    MycanvasSetting(c7,0.12,0.14,0.12,0.12);

    TVirtualPad *tvp[3];
    c1->cd();
	CountsMap->Draw("colz");
    c1->Print(countmapFile);
    c1->Print(countmapPNG);
 
    c2->cd();
    EnergyvsRing ->Draw("colz");
    c2->Print(energyvsringFile);
    c2->Print(energyvsringPNG);

    c2->Clear();
	c2->cd();
	EnergyvsRing ->ProfileX();
    sprintf(buf,"EnergyvsRing_pfx");
    EnergyvsRing_PFX= (TH1F*)gDirectory->Get(buf);
    HistPFX1dSetting(EnergyvsRing_PFX,buf,0,48,425,650,"RingSector","Energy");
    EnergyvsRing_PFX ->Draw("");
    c2->Print(energyvsringpfxFile);
    c2->Print(energyvsringpfxPNG);


    c3->cd();
    EnergyvsAngle->Draw("colz");
    c3->Print(energyvsangleFile);
    c3->Print(energyvsanglePNG);
    
    c3->Clear();
	c3->cd();
    EnergyvsAngle ->ProfileX();
    sprintf(buf,"EnergyvsAngle_pfx");
    EnergyvsAngle_PFX= (TH1F*)gDirectory->Get(buf);
    HistPFX1dSetting(EnergyvsAngle_PFX,buf,0,608,425,650,"AngleSector","Energy");
    EnergyvsAngle_PFX->Draw("");
    c3->Print(energyvsanglepfxPNG);


    c4->Divide(2,1);
    TVirtualPad *tv[2];
    for(int i=0;i<2;i++){
        if(i<1) {

            tv[i]=c4->cd(i+1);
            MyPadSetting(tv[i],0.12,0.12,0.12,0.15);
            RingMap->Draw();
        }
        else{
            tv[i]=c4->cd(i+1);
            MyPadSetting(tv[i],0.12,0.12,0.12,0.15);
            AngleMap->Draw();        
        }
    }
    c4->Print(ringangleFile);
    c4->Print(ringanglePNG);

    c5->cd();
    TIME->SetLineColor(4);
    TIME->SetFillColor(4);
    TIME->SetFillStyle(3001);
    TIME->Fit(f1);
    double aa=f1->GetParameter(1);
    double bb=f1->GetParameter(2);
    TIME->Fit(f1,"","",aa-5*bb,aa+5*bb);
    c5->Print(timeFile);
    c5->Print(timePNG);

    c5->Clear();
	c5->cd();
    Slices->SetLineColor(5);
    Slices->SetFillColor(5);
    Slices->SetFillStyle(3002);
    Slices->Draw();
	c5->Print(slicesFile);
    c5->Print(slicesPNG);

    c5->Clear();
	c5->cd();
    T1ADDT2->SetLineColor(5);
    T1ADDT2->SetFillColor(5);
    T1ADDT2->SetFillStyle(3002);
    T1ADDT2->Draw();
	c5->Print(t1addt2File);
    c5->Print(t1addt2PNG);

    c5->Clear();
	c5->cd();
    T1DIFFT2->SetLineColor(5);
    T1DIFFT2->SetFillColor(5);
    T1DIFFT2->SetFillStyle(3002);
    T1DIFFT2->Draw();
	c5->Print(t1difft2File);
    c5->Print(t1difft2PNG);

    TVirtualPad *tvpad[2];

    c6->Divide(1,2);
    tvpad[0]=c6->cd(1);
    MyPadSetting(tvpad[0],0.12,0.12,0.12,0.15);
    EnergyvsPixel->Draw("colz");
    tvpad[1]=c6->cd(2);
    MyPadSetting(tvpad[1],0.12,0.12,0.12,0.15);
    EnergyvsPixel->ProfileX();
    TH1F *EnergyvsPixel_PFX = (TH1F *)gDirectory->Get("EnergyvsPixel_pfx");
    HistPFX1dSetting(EnergyvsPixel_PFX,"EnergyPeak",0,NPixels,425,650,"PixelsNo.","EnergyPeak");
    EnergyvsPixel_PFX->Draw();
    c6->Print(energypixelFile);
    c6->Print(energypixelPNG);



      
 //   for(int i=0 ;i<NPixels;i++){
 //       c7->Clear();
 //       c7->cd();
 //       Energy[i]->Fit(f1,"","",0,1023);
 //       peakV[i] = f1->GetParameter(1);
 //       sigma[i] = f1->GetParameter(2);
 //       Energy[i]->Fit(f1,"","",peakV[i]-1.5*sigma[i],peakV[i]+2.5*sigma[i]);
 //       peakV[i] = f1->GetParameter(1);
 //       sigma[i] = f1->GetParameter(2);
 //     }
//
//
//    for(int i=0 ;i<608; i++){
//        for( int j=0; j<48; j++){
//            PeakMap->SetBinContent(i+1,j+1,peakV[j+48*i]);
//            ResolutionMap->SetBinContent(i+1,j+1,2.354*sigma[j+48*i]/peakV[j+48*i]);
//        }
//    }
    
    c1->Clear();
    c1->Divide(1,2);
   
    tvpad[0]=c1->cd(1);
    MyPadSetting(tvpad[0],0.12,0.12,0.12,0.15);
    PeakMap->GetZaxis()->SetRangeUser(400,600);
    PeakMap->Draw("colz");
    tvpad[1]=c1->cd(2);
    MyPadSetting(tvpad[1],0.12,0.12,0.12,0.15);
    ResolutionMap->GetZaxis()->SetRangeUser(0,0.3);
    ResolutionMap->Draw("colz");
    c1->Print(peakandresolutionmapPNG);
    c1->Print(peakandresolutionmapFile);

    c5->Print(energybankFile1);
    for(int i=0 ;i<38;i++){
        c5->Clear();
        c5->cd();
        EnergyBank[0][i]->SetLineColor(2);    
        EnergyBank[1][i]->SetLineColor(4);    
        EnergyBank[0][i]->SetFillColor(2);    
        EnergyBank[1][i]->SetFillColor(4);    
        EnergyBank[0][i]->SetFillStyle(3004);    
        EnergyBank[1][i]->SetFillStyle(3005);    
        dtemp = TMath::Max( EnergyBank[0][i]->GetMaximum() /EnergyBank[0][i]->GetEntries() ,EnergyBank[1][i]->GetMaximum() /EnergyBank[1][i]->GetEntries() );    
		EnergyBank[0][i]->GetXaxis()->SetRangeUser(0.0, dtemp);
		EnergyBank[0][i]->DrawNormalized();
        EnergyBank[1][i]->DrawNormalized("same");
        c5->Print(energybankFile2);
    }
    c5->Print(energybankFile3);
/* 
	c7->Print(energypixelFile1);
	TVirtualPad *tp[16];
    for(int i=0 ;i<NPixels;i++){
		if(i%16==0){
			c7->Clear();
	         c7->Divide(4,4);
		}
		tp[i%16]=c7->cd(i%16+1);
        cc = Energy[i]->GetBinCenter( Energy[i]->GetMaximumBin());
		Energy[i]->Draw();
		//Energy[i]->Fit(f1,"","",cc-60,cc+100);
        //peakV[i] = f1->GetParameter(1);
        //sigma[i] = f1->GetParameter(2);
        //Energy[i]->Fit(f1,"","",peakV[i]-1.5*sigma[i],peakV[i]+2.5*sigma[i]);
        //peakV[i] = f1->GetParameter(1);
        //sigma[i] = f1->GetParameter(2);
        
		if(i%16==15)c7->Print(energypixelFile2);  		

      }
	c7->Print(energypixelFile3);
*/

//    sprintf(buf,"%s/Peak_value.raw",filepath.c_str());
//    FILE *fileout=fopen(buf,"wb");
//    fwrite(&peakV,sizeof(double),NPixels,fileout);
//    fclose(fileout);
// 
//    sprintf(buf,"%s/Resolution.raw",filepath.c_str());
//    FILE *fileout1=fopen(buf,"wb");
//    fwrite(&sigma,sizeof(double),NPixels,fileout);
//    fclose(fileout1);
   

	sprintf(buf,"%s/Histogram_Results.root",filepath.c_str());
    TFile *rf=new TFile(buf,"recreate");
    rf->cd();
 
    TIME->Write();
	Slices->Write();
	T1ADDT2->Write();
	T1DIFFT2->Write();
	EnergyvsPixel->Write();
    PeakMap->Write();
    ResolutionMap->Write();
    CountsMap->Write();
    RingMap->Write();
    AngleMap->Write();
    EnergyvsRing->Write();
    EnergyvsAngle->Write();
    EnergyvsRing_PFX->Write();
    EnergyvsAngle_PFX->Write();
    for(int i=0 ;i<2;i++){
        for(int j=0 ;j<38;j++){
            EnergyBank[i][j]->Write();
        }
    }

    rf->Close();




return 0;
}

   void MyPadSetting(TVirtualPad *cc,double leftmargin,double rightmargin,double topmargin,double bottommargin){      
        cc->SetLeftMargin(leftmargin);
        cc->SetRightMargin(rightmargin);
        cc->SetTopMargin(topmargin);
        cc->SetBottomMargin(bottommargin);	  
        cc->SetBorderMode(0);
        cc->SetBorderSize(0);
        cc->SetFrameFillColor(10);
        cc->SetFrameBorderMode(0);
        cc->SetFrameBorderSize(0);
        cc->SetFrameLineWidth(2);        
      }

   void MycanvasSetting(TCanvas *cc,double leftmargin,double rightmargin,double topmargin,double bottommargin){      
        cc->SetLeftMargin(leftmargin);
        cc->SetRightMargin(rightmargin);
        cc->SetTopMargin(topmargin);
        cc->SetBottomMargin(bottommargin);	  
        cc->SetBorderMode(0);
        cc->SetBorderSize(0);
        cc->SetFrameFillColor(10);
        cc->SetFrameBorderMode(0);
        cc->SetFrameBorderSize(0);
           
      }
   TH1D *hist1dD(char *histname, Int_t Xstep,Double_t xlow,Double_t xup,char*xtitle,char *ytitle){      
      TH1D *Myhist1d=new TH1D(histname,histname,Xstep,xlow,xup);
      //Myhist1d->SetMinimum(ylow);
      //Myhist1d->SetMaximum(yup);
      Myhist1d->SetLineWidth(2.);
      Myhist1d->GetXaxis()->SetTitle(xtitle);
      Myhist1d->GetYaxis()->SetTitle(ytitle);
      Myhist1d->GetXaxis()->SetTitleSize(0.055);
      Myhist1d->GetXaxis()->SetTitleOffset(1.1);
      Myhist1d->GetXaxis()->SetLabelSize(0.045);
      Myhist1d->GetYaxis()->SetTitleSize(0.055);
      Myhist1d->GetYaxis()->SetTitleOffset(1.1);
      Myhist1d->GetYaxis()->SetLabelSize(0.045);
      Myhist1d->GetXaxis()->SetNdivisions(510);
      Myhist1d->GetXaxis()->CenterTitle();
      Myhist1d->GetYaxis()->CenterTitle();
      return Myhist1d;
     }

   TH1S *hist1dS(char *histname, Int_t Xstep,Double_t xlow,Double_t xup,char*xtitle,char *ytitle){      
      TH1S *Myhist1d=new TH1S(histname,"",Xstep,xlow,xup);
      //Myhist1d->SetMinimum(ylow);
      //Myhist1d->SetMaximum(yup);
      Myhist1d->SetLineWidth(2.);
      Myhist1d->GetXaxis()->SetTitle(xtitle);
      Myhist1d->GetYaxis()->SetTitle(ytitle);
      Myhist1d->GetXaxis()->SetTitleSize(0.055);
      Myhist1d->GetXaxis()->SetTitleOffset(1.1);
      Myhist1d->GetXaxis()->SetLabelSize(0.045);
      Myhist1d->GetYaxis()->SetTitleSize(0.055);
      Myhist1d->GetYaxis()->SetTitleOffset(1.1);
      Myhist1d->GetYaxis()->SetLabelSize(0.045);
      Myhist1d->GetXaxis()->SetNdivisions(510);
      Myhist1d->GetXaxis()->CenterTitle();
      Myhist1d->GetYaxis()->CenterTitle();
      return Myhist1d;
     }

   TH1F *hist1d(char *histname, Int_t Xstep,Double_t xlow,Double_t xup,char*xtitle,char *ytitle){      
      TH1F *Myhist1d=new TH1F(histname,"",Xstep,xlow,xup);
      //Myhist1d->SetMinimum(ylow);
      //Myhist1d->SetMaximum(yup);
      Myhist1d->SetLineWidth(2.);
      Myhist1d->GetXaxis()->SetTitle(xtitle);
      Myhist1d->GetYaxis()->SetTitle(ytitle);
      Myhist1d->GetXaxis()->SetTitleSize(0.055);
      Myhist1d->GetXaxis()->SetTitleOffset(1.1);
      Myhist1d->GetXaxis()->SetLabelSize(0.045);
      Myhist1d->GetYaxis()->SetTitleSize(0.055);
      Myhist1d->GetYaxis()->SetTitleOffset(1.1);
      Myhist1d->GetYaxis()->SetLabelSize(0.045);
      Myhist1d->GetXaxis()->SetNdivisions(510);
      Myhist1d->GetXaxis()->CenterTitle();
      Myhist1d->GetYaxis()->CenterTitle();
      return Myhist1d;
     }


 TH2F *hist2d(char *hist2dname,int Xstep,double xlow,double xup,int Ystep,double ylow,double yup,char *xtitle,char *ytitle){   
      TH2F *Myhist2d = new TH2F(hist2dname,"",Xstep,xlow,xup,Ystep,ylow,yup);
      //Myhist2d->SetMinimum(ylow);
      //Myhist2d->SetMaximum(yup);
      Myhist2d->SetLineWidth(2.);
      Myhist2d->GetXaxis()->SetTitle(xtitle);
      Myhist2d->GetYaxis()->SetTitle(ytitle);      
      Myhist2d->GetXaxis()->SetTitleSize(0.055);
      Myhist2d->GetXaxis()->SetTitleOffset(1.1);
      Myhist2d->GetXaxis()->SetLabelSize(0.055);
      Myhist2d->GetYaxis()->SetTitleSize(0.055);
      Myhist2d->GetYaxis()->SetTitleOffset(1.1);
      Myhist2d->GetYaxis()->SetLabelSize(0.055);
      Myhist2d->GetXaxis()->SetNdivisions(505);
      Myhist2d->GetYaxis()->SetNdivisions(505);
      Myhist2d->GetXaxis()->CenterTitle();
      Myhist2d->GetYaxis()->CenterTitle();
      //Myhist2d->GetXaxis()->SetNdivisions(512);
    return Myhist2d;
  }



 TH2S *hist2dS(char *hist2dname,int Xstep,double xlow,double xup,int Ystep,double ylow,double yup,char *xtitle,char *ytitle){   
      TH2S *Myhist2d = new TH2S(hist2dname,"",Xstep,xlow,xup,Ystep,ylow,yup);
      //Myhist2d->SetMinimum(ylow);
      //Myhist2d->SetMaximum(yup);
      Myhist2d->SetLineWidth(2.);
      Myhist2d->GetXaxis()->SetTitle(xtitle);
      Myhist2d->GetYaxis()->SetTitle(ytitle);      
      Myhist2d->GetXaxis()->SetTitleSize(0.055);
      Myhist2d->GetXaxis()->SetTitleOffset(1.1);
      Myhist2d->GetXaxis()->SetLabelSize(0.055);
      Myhist2d->GetYaxis()->SetTitleSize(0.055);
      Myhist2d->GetYaxis()->SetTitleOffset(1.1);
      Myhist2d->GetYaxis()->SetLabelSize(0.055);
      Myhist2d->GetXaxis()->SetNdivisions(505);
      Myhist2d->GetYaxis()->SetNdivisions(505);
      Myhist2d->GetXaxis()->CenterTitle();
      Myhist2d->GetYaxis()->CenterTitle();
      //Myhist2d->GetXaxis()->SetNdivisions(512);
    return Myhist2d;
  }



  void HistPFX1dSetting(  TH1F *hh,  char *hist1dname, double xlow, double xup, double ylow, double yup,  char *xtitle,   char *ytitle ){   
      hh->SetNameTitle(hist1dname, ""); 
      hh->SetMarkerStyle(20);
      hh->SetMarkerColor(4);
      hh->SetMarkerSize(0.5);
      hh->SetLineColor(4);
      //hh->SetMarkerSize(1.5);
      //hh->SetBins( Xstep, xlow, xup);  
      hh->SetMinimum(ylow);
      hh->SetMaximum(yup);
      hh->GetYaxis()->SetRangeUser(ylow, yup);
      hh->GetXaxis()->SetRangeUser(xlow, xup);
      hh->GetXaxis()->SetTitle(xtitle);
      hh->GetYaxis()->SetTitle(ytitle);      
      hh->GetXaxis()->SetTitleSize(0.055);
      hh->GetXaxis()->SetTitleOffset(1.1);
      hh->GetXaxis()->SetLabelSize(0.055);
      hh->GetYaxis()->SetTitleSize(0.055);
      hh->GetYaxis()->SetTitleOffset(1.1);
      hh->GetYaxis()->SetLabelSize(0.055);
      hh->GetXaxis()->SetNdivisions(510);
      hh->GetYaxis()->SetNdivisions(510);
      hh->GetXaxis()->CenterTitle();
      hh->GetYaxis()->CenterTitle();
  
      return ;
  }
 
   unsigned int dataconvert(unsigned int  data){
     return  (((data&0xFF000000)>>24 ) | ((data&0x00FF0000)>>8 ) | ((data&0x0000FF00)<<8 )| ((data&0x000000FF)<<24 ) );

  };

 void getFiles(string path, vector<string>& files){
      ifstream fin;
      string name;
      fin.open(path.c_str());
      while( fin.good() ){
        getline(fin , name);
        files.push_back(name);
      }
      fin.close();
 }



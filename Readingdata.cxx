#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TPad.h"
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
#include "TText.h"
#include "TSpectrum.h"
#include "TError.h"

using namespace std;
const double PI=3.14159265359;
const double R=415.;
const unsigned int NPixels=29184;
const unsigned int Nbin=1024;
const unsigned int ELUTN=29884416;
const unsigned int TLUTN=58368000;   //NCX*NZX*2000;
const int HalfWindow=1000;
const int NCX=608;
const int NZX=48;
const int Nmod=38;
const int Nbank=76;
const int TxInMod=16;
const int TxInBank=16;
const int AxInBank=24;
const int AxInMod=48;
const int AxInSubmod=8;
const int TxInSubmod=8;
const int halfNmod=19;

double peak0[NPixels];
double peak1[NPixels];
double res[NPixels];
double tmean[NPixels];
const unsigned int FakerFlagThreshold=6400;
string FILELIST="FILELIST.txt";

TH1D *hist1dD(char *histname,int Xstep,double xlow,double xup,char*xtitle,char *ytitle);
TH1F *hist1d(char *histname,int Xstep,double xlow,double xup,char*xtitle,char *ytitle);
TH1S *hist1dS(char *histname, Int_t Xstep,Double_t xlow,Double_t xup,char*xtitle,char *ytitle);
TH2S *hist2dS(char *hist2dname,int Xstep,double xlow,double xup,int Ystep,double ylow,double yup,char *xtitle,char *ytitle);
TH2F *hist2d(char *hist2dname,int Xstep,double xlow,double xup,int Ystep,double ylow,double yup,char *xtitle,char *ytitle);
void MycanvasSetting(TCanvas *cc,double leftmargin,double rightmargin,double topmargin,double bottommargin); 
void MyPadSetting(TVirtualPad *cc,double leftmargin,double rightmargin,double topmargin,double bottommargin);
void HistPFX1dSetting(  TH1F *hh,  char *hist1dname, double xlow, double xup, double ylow, double yup,  char *xtitle,   char *ytitle );
unsigned int PET3RawdataConvert(unsigned int  data);
void getFiles( string path, vector<string>& files );
double TxToAngle(int Tx);
int PET3Txshift(int Tx);
double GetX(int Tx);
double GetY(int Tx);
double GetZ(int Ax);


void ShowVersion(int a_returnCode){
    cout<<endl;
    cout<<"   PET3Calibration version  0.5  "<<endl;
    cout<<"     :add reading depth flag and option.  2018/12/25 "<<endl;
    cout<<"     :Adc range 1-1022.                   2018/12/29 "<<endl;
    cout<<"     :input & output exits ro not.        2019/01/02 "<<endl;
    cout<<"     Use ROOT histogram.                  2019/01/10 "<<endl;
    cout<<endl;
    exit(a_returnCode);
}
void ShowWrongDataType(int a_returnCode){
    cout<<endl;
    cout<<"***************************************************************************************************************"<<endl;
    cout<<"*********** Cannot read the data-type flag ! Please recheck input file type (Coincidence?? Singles??)  ***"<<endl;
    cout<<endl;
    exit(a_returnCode);
}
void ShowHelp(int a_returnCode){
    cout<<endl;
    cout<<endl;
    cout<<"***************************************************************************************************************"<<endl;
    cout<<" It's a application for Minfounders in MINFOUND-PET3(BGO) calibration. "<<endl;
    cout<<" Auther: Xiaozhuang Wang "<<endl;
    cout<<" Email: xiaozhuang.wang@minfound.com "<<endl;
    cout<<" Date: Dec. 20th, 2018 "<<endl; 
    cout<<endl;
    cout<<" PET3Calibration Usage "<<endl;
    cout<<" Options: "<<endl;
    cout<<"          -i/-il  path_to_input_file "<<endl;
    cout<<"          -o      path_of_output "<<endl;
    cout<<"         [-hist   yes/no ] "<<endl;
    cout<<"         [-pdf    yes/no ] "<<endl;
    cout<<"         [-png    yes/no ] "<<endl;
    cout<<"         [-txt    yes/no ] "<<endl;
    cout<<"         [-type   coin/single ] "<<endl;
    cout<<"         [-event  delay/prompt/both ] "<<endl;
    cout<<"         [-calib  energy/time/both ] "<<endl;
    cout<<"         [-iter   xxx ] "<<endl;
    cout<<"         [-dp     number ] "<<endl;
    cout<<"         [-all    ] "<<endl;
    cout<<"         [-v      ] "<<endl;
    cout<<"         [-h      ] "<<endl;
    cout<<endl;
    cout<<endl;
    cout<<" [Main setting]: "<<endl;
    cout<<"         -i  path_to_input_file    : give an input datafile of PET3 datafile(s) Within absolute PATH. "<<endl;
    cout<<"                                     option -in/-input is same as -i "<<endl;
    cout<<"         -il path_to_input_txtfile : give an input txt-like file With absolute PATH, txtfile includes PET3 data format files. "<<endl;
    cout<<"         -o  path_to_output        : give an absolute PATH for output, it must exist before analysis. "<<endl;
    cout<<"                                     option -out/-output is same as -o  "<<endl;
    cout<<"         -hist XXX                 : save histogram in rootfiles or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default positive)"<<endl;
    cout<<"         -pdf  XXX                 : save pdf file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default negetive)"<<endl;
    cout<<"         -png  XXX                 : save png file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default positive)"<<endl;
    cout<<"         -txt  XXX                 : save txt file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default negetive)"<<endl;
    cout<<"         -type  YYY                : data file type. YYY: coin/coins/single/singles. (default coincidence data)"<<endl; 
    cout<<"         -event YYY                : event kind. YYY:delay/prompt/both. (default both) "<<endl;
    cout<<"         -calib YYY                : calibration mode. YYY:energy/time/both. (default energy) "<<endl;
    cout<<"         -iter YYY                 : Iteration time in TIME calibration, YYY should ben >0 int , default( 3 ) "<<endl;
    cout<<"         -dp YYY                   : read part events of inputs file,YYY: int >0, default(0xFFFFFFFF)"  <<endl;
    cout<<"         -all                      : save all outputs. raw,hist,pdf,png,binary. "<<endl; 
    cout<<"         -v                        : print the version of this application "<<endl; 
    cout<<"                                     option -version/--version is same as -v  "<<endl;
    cout<<"         -h                        : print Help information "<<endl; 
    cout<<"                                     option -help/--help is same as -h  "<<endl;
    cout<<"***************************************************************************************************************"<<endl;
    cout<<endl;
    cout<<endl;
    exit(a_returnCode);
}


int main(int argc, char** argv)
{
   string path =(string) argv[1];
   char buf[100],buf2[100];
   int datatest,datatest1,datatest2; 
   ifstream filein;
   ofstream fileout;
   sprintf(buf,"%s/PET3TimeLUT.binary",path.c_str());
   sprintf(buf2,"%s/PET3TimeLUT_Correction.binary",path.c_str());
   fileout.open(buf2,ios::binary);
   filein.open(buf,ios::binary);
        filein.read((char*)(&datatest),sizeof(int));
        datatest1 = datatest;
        fileout.write((char*)(&datatest1),sizeof(int));
        
        filein.read((char*)(&datatest),sizeof(int));
        datatest1 = datatest;
        fileout.write((char*)(&datatest1),sizeof(int));
        
        filein.read((char*)(&datatest),sizeof(int));
        datatest1 = datatest;
        fileout.write((char*)(&datatest1),sizeof(int));
        
        filein.read((char*)(&datatest),sizeof(int));
        datatest1 = datatest;
        fileout.write((char*)(&datatest1),sizeof(int));
         
   for(int i=0;i<NCX*(NZX+16);i++){
        
        filein.read((char*)(&datatest),sizeof(int));
        datatest1 = datatest;
        //if(datatest==0){
        //    datatest1 =0x0A0A;
        //}
        datatest2 = -datatest1;
        datatest1 = (datatest2 &0xFFFF0000);
        fileout.write((char*)(&datatest1),sizeof(int));
   }
   filein.close();
   fileout.close();
/*	
	unsigned int timelutflag=0;
    fwrite(&timelutflag,sizeof(int),1,timelut);
	fwrite(&timelutflag,sizeof(int),1,timelut);
	fwrite(&timelutflag,sizeof(int),1,timelut);
	timelutflag = 0xFFFE;
	timelutflag = PET3RawdataConvert(timelutflag);
	fwrite(&timelutflag,sizeof(int),1,timelut);
    
   for(int i=0;i<Nbank;i++){
   	for(int j=0;j<AxInSubmod;j++){
   		for(int k=0;k<TxInSubmod;k++){
   			zeros = PET3RawdataConvert(zeros);
   			if(EnergyCalibOn)fwrite(&zeros,sizeof(int),1,energylut);
   			if(TimeCalibOn)fwrite(&zeros,sizeof(int),1,timelut);
   		    //filetxt<<"   index"<<counter<<"   value  "<<zeros<<"   aa  "<<0<<"   bb   "<<0<<endl;	
   			//counter++;
   		}		
   	}
   	for(int j=0;j<NZX;j++){
   		for(int k=0;k<TxInSubmod;k++){
   			temp1 = i*TxInSubmod+k;
   			aa=0;
   			bb=0;///////////////////////////initialized;
   		    
            Ax = j; 
            Tx = temp1;
            peakt0 = peak0[j+temp1*NZX]; 
   		    peakt1 = peak1[j+temp1*NZX];
   		    sigmat = res[j+temp1*NZX];
   		    gain =0;
   		    shift=0;
   		    par1=0;
   		    par2=0;
            tdelay=(tmean[j+temp1*NZX]);
            tpar=0;
///////////////////////////////////////////////////////////////////////////////////////////////////////
               if( peakt0>=1. && peakt1<1024. && sigmat>=1 && sigmat<1024.  ){
   				aa = static_cast<unsigned short>(  (511.)/(peakt0) *8192);
   			    bb =static_cast<unsigned short>(0x0);
                   gain = (511.)/(peakt0);
                   shift = 0.;
               }
               par1=aa;
               par2=bb;
   /////////////////////////////////////// //////////////////////////////////////////////////////////////                
               
   			unsigned int wxz = (static_cast<unsigned int>(aa))*65536 + static_cast<unsigned int>(bb&0xFFFF);
   			//filetxt<<"   index"<<counter<<"   value  "<<wxz<<"   aa  "<<static_cast<unsigned short>((wxz>>16)&0xFFFF)<<"   bb   "<<static_cast<short>((wxz)&0xFFFF)<<endl;
   			wxz = PET3RawdataConvert(wxz);
   			if(EnergyCalibOn)fwrite(&wxz,sizeof(unsigned int),1,energylut);
   			
           
            tpar= static_cast<short>(tdelay);
            wxz = static_cast<unsigned int>(tpar&0xFFFF);
            wxz = PET3RawdataConvert(wxz);
            if(TimeCalibOn)fwrite(&wxz,sizeof(unsigned int),1,timelut);
   			//counter+i+;
            tr->Fill();
   		}		
   	}
   
   	for(int j=0;j<AxInSubmod;j++){
   	    	for(int k=0;k<TxInSubmod;k++){
   	    		zeros = PET3RawdataConvert(zeros);
   	    		if(EnergyCalibOn)fwrite(&zeros,sizeof(int),1,energylut);
   	    		if(TimeCalibOn)fwrite(&zeros,sizeof(int),1,timelut);
   	    		//filetxt<<"   index"<<counter<<"   value  "<<zeros<<"   aa  "<<0<<"   bb   "<<0<<endl;
   	    		//counter++;
   	    	}
   	    }
   }
   if(EnergyCalibOn)fclose(energylut);
   if(TimeCalibOn)fclose(timelut);
   //filetxt.close();
    
    
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
      Myhist1d->GetYaxis()->SetTitleOffset(0.9);
      Myhist1d->GetYaxis()->SetLabelSize(0.045);
      Myhist1d->GetXaxis()->SetNdivisions(510);
      Myhist1d->GetXaxis()->CenterTitle();
      Myhist1d->GetYaxis()->CenterTitle();
      return Myhist1d;
     }

   TH1S *hist1dS(char *histname, Int_t Xstep,Double_t xlow,Double_t xup,char*xtitle,char *ytitle){      
      TH1S *Myhist1d=new TH1S(histname,histname,Xstep,xlow,xup);
      //Myhist1d->SetMinimum(ylow);
      //Myhist1d->SetMaximum(yup);
      Myhist1d->SetLineWidth(2.);
      Myhist1d->GetXaxis()->SetTitle(xtitle);
      Myhist1d->GetYaxis()->SetTitle(ytitle);
      Myhist1d->GetXaxis()->SetTitleSize(0.055);
      Myhist1d->GetXaxis()->SetTitleOffset(1.1);
      Myhist1d->GetXaxis()->SetLabelSize(0.045);
      Myhist1d->GetYaxis()->SetTitleSize(0.055);
      Myhist1d->GetYaxis()->SetTitleOffset(0.95);
      Myhist1d->GetYaxis()->SetLabelSize(0.045);
      Myhist1d->GetXaxis()->SetNdivisions(510);
      Myhist1d->GetXaxis()->CenterTitle();
      Myhist1d->GetYaxis()->CenterTitle();
      return Myhist1d;
     }

   TH1F *hist1d(char *histname, Int_t Xstep,Double_t xlow,Double_t xup,char*xtitle,char *ytitle){      
      TH1F *Myhist1d=new TH1F(histname,histname,Xstep,xlow,xup);
      //Myhist1d->SetMinimum(ylow);
      //Myhist1d->SetMaximum(yup);
      Myhist1d->SetLineWidth(2.);
      Myhist1d->GetXaxis()->SetTitle(xtitle);
      Myhist1d->GetYaxis()->SetTitle(ytitle);
      Myhist1d->GetXaxis()->SetTitleSize(0.055);
      Myhist1d->GetXaxis()->SetTitleOffset(1.1);
      Myhist1d->GetXaxis()->SetLabelSize(0.045);
      Myhist1d->GetYaxis()->SetTitleSize(0.055);
      Myhist1d->GetYaxis()->SetTitleOffset(0.95);
      Myhist1d->GetYaxis()->SetLabelSize(0.045);
      Myhist1d->GetXaxis()->SetNdivisions(510);
      Myhist1d->GetXaxis()->CenterTitle();
      Myhist1d->GetYaxis()->CenterTitle();
      return Myhist1d;
     }


 TH2F *hist2d(char *hist2dname,int Xstep,double xlow,double xup,int Ystep,double ylow,double yup,char *xtitle,char *ytitle){   
      TH2F *Myhist2d = new TH2F(hist2dname,hist2dname,Xstep,xlow,xup,Ystep,ylow,yup);
      //Myhist2d->SetMinimum(ylow);
      //Myhist2d->SetMaximum(yup);
      Myhist2d->SetLineWidth(2.);
      Myhist2d->GetXaxis()->SetTitle(xtitle);
      Myhist2d->GetYaxis()->SetTitle(ytitle);      
      Myhist2d->GetXaxis()->SetTitleSize(0.055);
      Myhist2d->GetXaxis()->SetTitleOffset(1.1);
      Myhist2d->GetXaxis()->SetLabelSize(0.055);
      Myhist2d->GetYaxis()->SetTitleSize(0.055);
      Myhist2d->GetYaxis()->SetTitleOffset(0.95);
      Myhist2d->GetYaxis()->SetLabelSize(0.055);
      Myhist2d->GetXaxis()->SetNdivisions(505);
      Myhist2d->GetYaxis()->SetNdivisions(505);
      Myhist2d->GetXaxis()->CenterTitle();
      Myhist2d->GetYaxis()->CenterTitle();
      //Myhist2d->GetXaxis()->SetNdivisions(512);
    return Myhist2d;
  }



 TH2S *hist2dS(char *hist2dname,int Xstep,double xlow,double xup,int Ystep,double ylow,double yup,char *xtitle,char *ytitle){   
      TH2S *Myhist2d = new TH2S(hist2dname,hist2dname,Xstep,xlow,xup,Ystep,ylow,yup);
      //Myhist2d->SetMinimum(ylow);
      //Myhist2d->SetMaximum(yup);
      Myhist2d->SetLineWidth(2.);
      Myhist2d->GetXaxis()->SetTitle(xtitle);
      Myhist2d->GetYaxis()->SetTitle(ytitle);      
      Myhist2d->GetXaxis()->SetTitleSize(0.055);
      Myhist2d->GetXaxis()->SetTitleOffset(1.1);
      Myhist2d->GetXaxis()->SetLabelSize(0.055);
      Myhist2d->GetYaxis()->SetTitleSize(0.055);
      Myhist2d->GetYaxis()->SetTitleOffset(0.95);
      Myhist2d->GetYaxis()->SetLabelSize(0.055);
      Myhist2d->GetXaxis()->SetNdivisions(505);
      Myhist2d->GetYaxis()->SetNdivisions(505);
      Myhist2d->GetXaxis()->CenterTitle();
      Myhist2d->GetYaxis()->CenterTitle();
      //Myhist2d->GetXaxis()->SetNdivisions(512);
    return Myhist2d;
  }



  void HistPFX1dSetting(  TH1F *hh,  char *hist1dname, double xlow, double xup, double ylow, double yup,  char *xtitle,   char *ytitle ){   
      hh->SetNameTitle(hist1dname, hist1dname); 
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
      hh->GetYaxis()->SetTitleOffset(0.95);
      hh->GetYaxis()->SetLabelSize(0.055);
      hh->GetXaxis()->SetNdivisions(510);
      hh->GetYaxis()->SetNdivisions(510);
      hh->GetXaxis()->CenterTitle();
      hh->GetYaxis()->CenterTitle();
  
      return ;
  }
  
 unsigned int PET3RawdataConvert(unsigned int  data){
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
double TxToAngle(int Tx=0){
     if(Tx<0||Tx>=NCX){
        return -1;
     }
     else{
        return  Tx*PI*1./NCX;
     }
}
int PET3Txshift(int Tx=0){
    if(Tx <592)return Tx+TxInMod;
    else  return Tx+TxInMod-NCX;
}

double GetX(int Tx=0){
    if(Tx<0||Tx>=NCX){
        return -999999;
    }
    else{
        double x;
        int nmod=int(Tx/TxInMod);
        int npix=int(Tx%TxInMod);
        double xbase = (7-npix)*4.2+2.1;
        double ybase = R;
        x=xbase*cos( nmod*PI*2./Nmod) - ybase*sin( nmod*PI*2./Nmod) ; 
        return x;
    }
}

double GetY(int Tx=0){
    if(Tx<0||Tx>=NCX){
        return -999999;
    }
    else{
        double y;
        int nmod=int(Tx/TxInMod);
        int npix=int(Tx%TxInMod);
        double xbase = (7-npix)*4.2+2.1;
        double ybase = R;
        y=xbase*sin( nmod*PI*2./Nmod) +ybase*cos( nmod*PI*2./Nmod) ; 
        return y;
    }
}


double GetZ(int Ax=0){
    if(Ax<0||Ax>=NZX){
        return -999999.;
    }
    else{
        double z;
        z = (Ax-AxInMod)*4.2+2.1;
        return z;
    }
}

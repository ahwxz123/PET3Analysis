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
using namespace std;
const double PI=3.14159265359;
const double R=415.;
const unsigned int NPixels=29184;
const unsigned int Nbin=1024;
const unsigned int ELUTN=29884416;
const int NCX=608;
const int NZX=48;
const int NSubAx=6;
const int NSubTx=76;
const int Nmod=38;
const int TxInSubMod=8;
const int AxInSubMod=8;
const int TxInBank=16;
const int AxInBank=24;
const int TxInMod=16;
const int AxInMod=24;
const int halfNmod=19;

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
    cout<<"   PET3Analysis version  1.1  "<<endl;
    cout<<"     :add reading depth flag and option.                       2018/12/25 "<<endl;
    cout<<"     :Adc range 1-1022.                                        2018/12/29 "<<endl;
    cout<<"     :input & output exits or not.                             2019/01/02 "<<endl;
    cout<<"     :Coincidences 0xFFFFFFC Flags;32 events in one flag.      2019/01/11 "<<endl;
    cout<<"     :Big Endian       .                                       2019/01/11 "<<endl;
    cout<<"     :Singles no flag.                                         2019/01/11 "<<endl;
    cout<<"     :--------------------       .                             2019/01/02 "<<endl;
    cout<<"     :Add Submodule Energy;                                    2019/01/14 "<<endl;
    cout<<"     :Singles High32 Low32 exchange.                           2019/01/18 "<<endl;
    cout<<"     :SetStats flag On, 2D floop map flag-off                    2019/01/21 "<<endl;
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
    cout<<" It's a application for Minfounders in MINFOUND-PET3(BGO) files analysis. "<<endl;
    cout<<" Auther: Xiaozhuang Wang "<<endl;
    cout<<" Email: xiaozhuang.wang@minfound.com "<<endl;
    cout<<" Date: Dec. 27th, 2018 "<<endl; 
    cout<<endl;
    cout<<" PET3Analysis Usage "<<endl;
    cout<<" Options: "<<endl;
    cout<<"          -i/-il  path_to_input_file "<<endl;
    cout<<"          -o      path_of_output "<<endl;
    cout<<"         [-raw    yes/no ] "<<endl;
    cout<<"         [-hist   yes/no ] "<<endl;
    cout<<"         [-pdf    yes/no ] "<<endl;
    cout<<"         [-png    yes/no ] "<<endl;
    cout<<"         [-bin    yes/no ] "<<endl;
    cout<<"         [-sino   yes/no ] "<<endl;
    cout<<"         [-type   coin/single ] "<<endl;
    cout<<"         [-event  delay/prompt/both ] "<<endl;
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
    cout<<"         -raw  XXX                 : save raw data in rootfiles or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default negetive)"<<endl;
    cout<<"         -hist XXX                 : save histogram in rootfiles or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default positive)"<<endl;
    cout<<"         -pdf  XXX                 : save pdf file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default negetive)"<<endl;
    cout<<"         -png  XXX                 : save png file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default positive)"<<endl;
    cout<<"         -bin  XXX                 : save binary file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default negetive)"<<endl;
    cout<<"         -sino  XXX                : save sino file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default positive)"<<endl;
    cout<<"         -type  YYY                : data file type. YYY: coin/coins/single/singles. (default coincidence data)"<<endl; 
    cout<<"         -event YYY                : event kind. YYY:delay/prompt/both. (default prompt) "<<endl;
    cout<<"         -dp    YYY                : Read only some part of events,YYY is a int number. (default 0xFFFFFFFF) "<<endl;
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
    
    if(argc==1) ShowHelp(0);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////                          Parameterized variables with default values                                   ////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    string input_file="";
    vector<string> path_to_input_files;
    vector<string> txtstring;
    string path_to_output_file="";
    string data_type="coin";
    const string single_data_type="single";
    const string coin_data_type="coin";
    unsigned int ReadDepth=0xFFFFFFFF;
    int event_type=1;           // 0  delay,  1  prompt,   2  both;
    int Percentage=10;          //10%; 
    bool Save_Raw=false;
    bool Save_Hist=true;
    bool Save_Pdf=false;
    bool Save_Png=true;
    bool Save_Bin=false;
    bool Save_Sino=true;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////                                         Reading Command-Line parameters                                ////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(int i=1;i<argc;i++){
        string option = (string) argv[i];
         
        if(option=="-h" || option=="--help" ||option=="-help") {
            ShowHelp(0);
        }
        else if(option=="-v"||option=="-version"||option=="--version"){
            ShowVersion(0); 
        } 
        else if(option=="-i"||option=="-in" || option=="-input"){
            if(argv[i+1]==NULL) {
                cerr<<"******** argument missing for option:"<<option<<endl;
                exit(1);
            }
            else{
                input_file = argv[i+1];
                path_to_input_files.push_back(input_file);
                     
                i++;
            }

            FILE *filetest;
            filetest=fopen(input_file.c_str(),"r");
            if( filetest==NULL){
                cerr<<"*********  input file does not exist!!!!  Please recheck the path or file name !!!"<<endl;
                exit(1);
            }
            
        }
        else if(option=="-il"){
            if(argv[i+1]==NULL) {
                cerr<<"******** argument missing for option:"<<option<<endl;
                exit(1);
            }
            else{
                input_file =(string) argv[i+1];
                ifstream txtfileread;
                string fname;
                txtfileread.open(input_file.c_str());
                if(! txtfileread){
                    txtfileread.close();
                    cerr<<"*********  input TXT file does not exist!!!!  Please recheck the path or TXT file name !!!"<<endl;
                    exit(1);
                }

                while(txtfileread.good()){
                    getline(txtfileread, fname);
                    path_to_input_files.push_back(fname);
                    FILE  *filetest;
                    filetest=fopen(fname.c_str(),"r");
                    if(filetest=fopen==NULL){
                        cerr<<"*********  input file "<<fname<<"  does not exist!!!!  Please recheck the path or file name !!!"<<endl;
                        exit(1);
                    }
                         
                }
                txtfileread.close();
            
                i++;
            }
        }
        else if(option=="-o"||option=="-out"||option=="-output"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                exit(1);
            }
            else{
                path_to_output_file = (string) argv[i+1];
                if(access(path_to_output_file.c_str(),0)== -1){
                   cerr<<"******** outpath does not exist!!! Pleaser recheck the Path "<<endl;
                   exit(1);
                }
                i++;
            }
        }
        else if(option=="-raw"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara =="yes"||optionpara=="YES"||optionpara=="Yes"||optionpara=="Y"||optionpara=="y"){
                    Save_Raw = true;                 
                }
                else if(optionpara =="no"||optionpara=="NO"||optionpara=="No"||optionpara=="N"||optionpara=="n"){
                    Save_Raw = false;                 
                }
                else{
                    cerr<<"******** invalid paramters for option:  "<<option<<"  !"<<endl;
                    exit(1);
                }
    
                i++;
            } 
        }
        else if(option=="-hist"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara =="yes"||optionpara=="YES"||optionpara=="Yes"||optionpara=="Y"||optionpara=="y"){
                    Save_Hist = true;                 
                }
                else if(optionpara =="no"||optionpara=="NO"||optionpara=="No"||optionpara=="N"||optionpara=="n"){
                    Save_Hist = false;                 
                }
                else{
                    cerr<<"******** invalid paramters for option:  "<<option<<"  !"<<endl;
                    exit(1);
                }
    
                i++;
            } 
        }
        else if(option=="-pdf"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara =="yes"||optionpara=="YES"||optionpara=="Yes"||optionpara=="Y"||optionpara=="y"){
                    Save_Pdf = true;                 
                }
                else if(optionpara =="no"||optionpara=="NO"||optionpara=="No"||optionpara=="N"||optionpara=="n"){
                    Save_Pdf = false;                 
                }
                else{
                    cerr<<"******** invalid paramters for option:  "<<option<<"  !"<<endl;
                    exit(1);
                }
    
                i++;
            } 
        }
        else if(option=="-png"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara =="yes"||optionpara=="YES"||optionpara=="Yes"||optionpara=="Y"||optionpara=="y"){
                    Save_Png = true;                 
                }
                else if(optionpara =="no"||optionpara=="NO"||optionpara=="No"||optionpara=="N"||optionpara=="n"){
                    Save_Png = false;                 
                }
                else{
                    cerr<<"******** invalid paramters for option:  "<<option<<"  !"<<endl;
                    exit(1);
                }
    
                i++;
            } 
        }
        else if(option=="-bin"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara =="yes"||optionpara=="YES"||optionpara=="Yes"||optionpara=="Y"||optionpara=="y"){
                    Save_Bin = true;                 
                }
                else if(optionpara =="no"||optionpara=="NO"||optionpara=="No"||optionpara=="N"||optionpara=="n"){
                    Save_Bin = false;                 
                }
                else{
                    cerr<<"******** invalid paramters for option:  "<<option<<"  !"<<endl;
                    exit(1);
                }
                
                i++;
            } 
        }
        else if(option=="-sino"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara =="yes"||optionpara=="YES"||optionpara=="Yes"||optionpara=="Y"||optionpara=="y"){
                    Save_Sino = true;                 
                }
                else if(optionpara =="no"||optionpara=="NO"||optionpara=="No"||optionpara=="N"||optionpara=="n"){
                    Save_Sino = false;                 
                }
                else{
                    cerr<<"******** invalid paramters for option:  "<<option<<"  !"<<endl;
                    exit(1);
                }
                
                i++;
            } 
        }
        else if(option=="-type"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara == "singles"||optionpara=="single"||optionpara=="Singles"||optionpara=="Single"){
                    data_type = single_data_type;                 
                }
                else if(optionpara == "coins"||optionpara=="coin"||optionpara=="Coins"||optionpara=="Coin"){
                    data_type = coin_data_type;                 
                }
                else{
                    cerr<<"******** invalid paramters for option:"<<option<<endl;
                    exit(1);
                }
            } 
    
            i++;
        }
        
        else if(option=="-dp"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                ReadDepth = atoi(optionpara.c_str());
                if(ReadDepth<0 ){
                    cerr<<"******** invalid paramters for option:"<<option<<endl;
                    cerr<<"******** it should be a >0 integral number     !!!  recommend among 1~8  ..."<<endl;
                    exit(1);
                }
                else{
                    cout<<"******** Only Read  "<<ReadDepth <<" events in data files !!!"<<endl;
                }
            } 
    
            i++;
        }
        else if(option=="-event"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara == "both"){
                    event_type = 2;                 
                }
                else if(optionpara == "prompt"||optionpara == "pp"){
                    event_type = 1;                 
                }
                else if(optionpara == "delay"||optionpara == "de"){
                    event_type = 0;                 
                }
                else{
                    cerr<<"******** invalid paramters for option:"<<option<<endl;
                    exit(1);
                }
            } 
    
            i++;
        }
        else if(option=="-all"||option=="-All"||option=="-ALL"){
            Save_Raw=true;
            Save_Hist=true;
            Save_Pdf=true;
            Save_Png=true;
            Save_Bin=true;
            Save_Sino=true;
            
            event_type=2;
        } 
        else{
            cerr<<endl;
            cerr<<"***************** Unknown options: "<<option<<"   ! Please recheck your spells"<<endl;
            cerr<<endl;
            exit(1);
        }

    }
    
/////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                 Parameters in file decode                                      ////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


      //gStyle->SetOptStat(0);
	  gStyle->SetPadLeftMargin(0.15);
	  gStyle->SetPadRightMargin(0.15);
	  gStyle->SetPadTopMargin(0.12);
	  gStyle->SetPadBottomMargin(0.15);
      gErrorIgnoreLevel=kError,kBreak,kSysError,kFatal;
      TStopwatch *ts1=new TStopwatch();
	  ts1->Start(); 

	  TDatime *tdate =new TDatime();
      int systime = tdate->GetTime();
      int sysdate = tdate->GetDate();
      char buf[100],buf2[100];

      const unsigned int CoincidenceFlag=0xFFFFFFFC;

	  const int energylow=2;
	  const int energyhigh=1020;
	  const int Ntimebin=2048;
      unsigned int fakerflag=0;
      int counts;
      int index;
      int index1,index2,bankindex1,bankindex2;
      unsigned int eindex1, eindex2;
	  int temp1,temp2;
	  double cc;

      int Adc1= -1,Adc2 = -1, tdiff= -1,  DeltaT = -1;
      int Ax1= -1,  Ax2= -1,    Tx1= -1,  Tx2= -1;
      int co=-1,co1=-1,co2=-1,kind;
      int multhit, ovflow;
      int Ax,Tx, Adc, coarse,fine,timezone;
      int Txmax,Txmin;
      double timetick;
      double coarseLSB=10.;       //ns;
      double fineLSB=0.3125;      //ns;
      
      double x1,y1,x2,y2,z1,z2;
      double theta,S,coeff;

      int Nshift=0;
      unsigned int Nslotbyte=64;      
      unsigned int Nbyte=Nslotbyte/2;      
	  int data[Nslotbyte];
      ULong64_t datacon[Nslotbyte] ;
      ULong64_t dataeff[Nbyte] ;
      unsigned int datatest;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                              histgream define                                               ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<"---------------- Start histograms definations---------------"<<endl;



      TH2F *Countmap1 = hist2d("Countmap1",NCX,0,Nmod,NZX,0,2,"Trans_Axial","Axial");
      TH2F *Countmap2 = hist2d("Countmap2",NCX,0,Nmod,NZX,0,2,"Trans_Axial","Axial");
      TH2F *Countmap = hist2d("Countmap",NCX,0,Nmod,NZX,0,2,"Trans_Axial","Axial");
      Countmap1->GetXaxis()->SetNdivisions(505);
      Countmap2->GetXaxis()->SetNdivisions(505);
      Countmap ->GetXaxis()->SetNdivisions(505);

      TH1F *Axial1 = hist1d("Axial_1st",NZX,0,NZX,"Axial","Entries");
      TH1F *Axial2 = hist1d("Axial_2nd",NZX,0,NZX,"Axial","Entries");
      TH1F *Axial = hist1d("Axial",NZX,0,NZX,"Axial","Entries");
      TH1F *TransAxial1 = hist1d("Trans_Axial_1st",NCX,0,NCX,"Axial","Entries");
      TH1F *TransAxial2 = hist1d("Trans_Axial_2nd",NCX,0,NCX,"Axial","Entries");
      TH1F *TransAxial = hist1d("Trans_Axial",NCX,0,NCX,"Axial","Entries");
      TH1F *Energy1= hist1d("Energy1",Nbin,0,Nbin,"Energy first","Entries");
      TH1F *Energy2= hist1d("Energy2",Nbin,0,Nbin,"Energy second","Entries");
      TH1F *Energy = hist1d("Energy",Nbin,0,Nbin,"Energy","Entries");
      TH1F *Time = hist1d("Time",Ntimebin,-Ntimebin/2,Ntimebin/2,"Time","Entries");
      
      TH2S *EnergyvsRing = hist2dS("EnergyvsRing",NZX,0,NZX,Nbin,0,Nbin,"Axial","Energy");
      TH2S *EnergyvsAngle = hist2dS("EnergyvsAngle",NCX,0,NCX,Nbin,0,Nbin,"TransAxial","Energy");

      TH1F *EnergyModule[Nmod];
      TH1S *EnergySubMod[NSubAx][NSubTx];
      TH1F *Timecouple[halfNmod];
      for(int i=0;i<Nmod;i++){
          sprintf(buf,"Energy_mod%d",i);
          EnergyModule[i] =hist1d(buf,Nbin,0,Nbin,"Energy","Entries");
          if(i<halfNmod){
              sprintf(buf,"time_couple_%dvs%d",i,i+halfNmod);
              Timecouple[i] =hist1d(buf,400,-200,200,"Time","Entries");
          }
      }
      
      for(int i=0;i<NSubAx;i++){
        for(int j=0;j<NSubTx;j++){
            sprintf(buf,"Energy_SubMod_Ax%d_Tx%d",i,j);
            EnergySubMod[i][j] = hist1dS(buf,Nbin,0,Nbin,"Energy","Entries");
        }
     }
     TH2S *TimevsAngle =hist2dS("TimevsAngle",NCX,0,NCX,400,-200,200,"Angle","Time");
     TH2S *TimevsRing =hist2dS("TimevsRing",NZX,0,NZX,400,-200,200,"Axial","Time");
     
     TH2S *ThetavsS =hist2dS("ThetavsS",334,-350.7,350.7,304,0,PI,"S","Theta");

     Countmap ->SetStats(kFALSE);
     Countmap1->SetStats(kFALSE);
     Countmap2->SetStats(kFALSE);
     EnergyvsAngle ->SetStats(kFALSE);
     EnergyvsRing  ->SetStats(kFALSE);
     TimevsAngle ->SetStats(kFALSE);
     TimevsRing  ->SetStats(kFALSE);
     ThetavsS    ->SetStats(kFALSE);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                              string for pdf/png                                              ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ////////////// string name for saving pdf/png;
    cout<<"---------------- Start path/pdf/png definations---------------"<<endl;
      string CountMapName             = path_to_output_file+"/"+"CountMap";
      string AxialANDTransName        = path_to_output_file+"/"+"Axial_AND_Trans_Hist";
      string EnergyTotalName          = path_to_output_file+"/"+"Energy_Total_Hist";
      string EnergyModuleName         = path_to_output_file+"/"+"Energy_Module_Hist";
      string EnergyBankName           = path_to_output_file+"/"+"Energy_Bank_Hist";
      string EnergyvsAngleName        = path_to_output_file+"/"+"Energy_vs_angle_Hist";
      string EnergyvsRingName         = path_to_output_file+"/"+"Energy_vs_ring_Hist";
      string Energvs2DimensionsName   = path_to_output_file+"/"+"Energy_vs_ring_Energy_vs_angle_Hist";
      string TimeTotalName            = path_to_output_file+"/"+"Time_Total_Hist";
      string TimeModuleName           = path_to_output_file+"/"+"Time_Module_Hist";
      string TimeBankName             = path_to_output_file+"/"+"Time_Bank_Hist";
      string TimevsAngleName          = path_to_output_file+"/"+"Time_vs_angle_Hist";
      string TimevsRingName           = path_to_output_file+"/"+"Time_vs_ring_Hist";
      string ThetavsSName             = path_to_output_file+"/"+"Theta_vs_S";


      string CountMapPdf             =CountMapName +".pdf";
      string AxialANDTransPdf        =AxialANDTransName+".pdf";
      string EnergyModulePdf         =EnergyModuleName+".pdf";
      string EnergyTotalPdf          =EnergyTotalName+".pdf";
      string EnergyBankPdf           =EnergyBankName+".pdf";
      string EnergyvsAnglePdf        =EnergyvsAngleName+".pdf";
      string EnergyvsRingPdf         =EnergyvsRingName+".pdf";
      string Energvs2DimensionsPdf   =Energvs2DimensionsName+".pdf";
      string TimeTotalPdf            =TimeTotalName+".pdf";
      string TimeModulePdf           =TimeModuleName+".pdf";
      string TimeBankPdf             =TimeBankName+".pdf";
      string TimevsAnglePdf          =TimevsAngleName+".pdf";
      string TimevsRingPdf           =TimevsRingName+".pdf";
      string ThetavsSPdf             =ThetavsSName+".pdf";
      
      string CountMapPng             =CountMapName +".png";
      string AxialANDTransPng        =AxialANDTransName+".png";
      string EnergyModulePng         =EnergyModuleName+".png";
      string EnergyTotalPng          =EnergyTotalName+".png";
      string EnergyBankPng           =EnergyBankName+".png";
      string EnergyvsAnglePng        =EnergyvsAngleName+".png";
      string EnergyvsRingPng         =EnergyvsRingName+".png";
      string Energvs2DimensionsPng   =Energvs2DimensionsName+".png";
      string TimeTotalPng            =TimeTotalName+".png";
      string TimeModulePng           =TimeModuleName+".png";
      string TimeBankPng             =TimeBankName+".png";
      string TimevsAnglePng          =TimevsAngleName+".png";
      string TimevsRingPng           =TimevsRingName+".png";
      string ThetavsSPng             =ThetavsSName+".png";

      string CountMapBin             =CountMapName +".binary";
      string EnergyModuleBin         =EnergyModuleName +".binary";
      string TimeModuleBin           =TimeModuleName+".binary";
      string EnergyvsAngleBin        =EnergyvsAngleName+".binary";
      string EnergyvsRingBin         =EnergyvsRingName+".binary";
      string TimevsAngleBin          =TimevsAngleName+".binary";
      string TimevsRingBin           =TimevsRingName+".binary";
      string ThetavsSBin             =ThetavsSName+".binary";


      FILE *binfile;

      //Generate root file;
      string root_output_file = path_to_output_file +"/" +"Rawdata.root";
      string root_hist_file = path_to_output_file +"/" +"Histogram.root";
      TFile *rfile =new TFile(root_output_file.c_str(),"recreate");
      TTree *raw = new TTree("raw","raw");
      if(Save_Raw){
          if(data_type==coin_data_type){
              raw->Branch("kind", &kind,"kind/I");
              raw->Branch("Ax1", &Ax1,"Ax1/I");
              raw->Branch("Tx1", &Tx1,"Tx1/I");
              raw->Branch("Adc1", &Adc1,"Adc1/I");
              raw->Branch("Ax2", &Ax2,"Ax2/I");
              raw->Branch("Tx2", &Tx2,"Tx2/I");
              raw->Branch("Adc2", &Adc2,"Adc2/I");
              raw->Branch("tdiff", &tdiff,"tdiff/I");
              
              raw->Branch("x1", &x1, "x1/D");
              raw->Branch("x2", &x2, "x2/D");
              raw->Branch("y1", &y1, "y1/D");
              raw->Branch("y2", &y2, "y2/D");
              raw->Branch("theta", &theta, "theta/D");
              raw->Branch("S", &S, "S/D");
          }
          if(data_type==single_data_type){
              raw->Branch("multhit"   ,&multhit       ,"multhit/I"     );
              raw->Branch("ovflow"    ,&ovflow        ,"ovflow/I"      );
              raw->Branch("Ax"        ,&Ax            ,"Ax/I"          );
              raw->Branch("Tx"        ,&Tx            ,"Tx/I"          );
              raw->Branch("Adc"       ,&Adc           ,"Adc/I"         );
              raw->Branch("coarse"    ,&coarse        ,"coarse/I"      );
              raw->Branch("fine"      ,&fine          ,"fine/I"        );
              raw->Branch("timezone"  ,&timezone      ,"timezone/I"    ); //  ns;
              raw->Branch("timetick"  ,&timetick      ,"timetick/D"    ); //  ns;
          }
      }







    cout<<"---------------- Start reading data files---------------"<<endl;

	int size = path_to_input_files.size();
    cout<<"----Reading inputs files :"<<endl;      
    for(int i=0;i<size;i++){
	     cout<<"    "<<path_to_input_files[i].c_str()<<endl;      
	}
    int PercantRead;
    unsigned int Readcounts;
    ULong64_t ReadLimit=64;

	ifstream fin;
    for(int MM=0; MM<size; MM++){
        fin.open(path_to_input_files[MM].c_str(),ios::binary);
        //fin.seekg( 0,ios::end);
        //ReadLimit =  fin.tellg()*Percentage/100./sizeof(unsigned int)  ;
        fin.seekg( Nshift*4,ios::beg);
        Readcounts=0;
        if(data_type==coin_data_type){
            
            while(fin.good() &&Readcounts<= ReadDepth ){
                fin.read((char*)(&datatest), sizeof(unsigned int));
                Readcounts +=1;
                //datatest = PET3RawdataConvert(datatest);
                if(CoincidenceFlag!= datatest){
                    fakerflag++;
                    if(fakerflag>=FakerFlagThreshold)ShowWrongDataType(1);
                   cout<<"-----------------------------------faker flag" <<endl;
                    continue;
                }
                else{
                    fin.read((char*)(&data), Nslotbyte*sizeof(unsigned int));
                    Readcounts +=Nslotbyte;
                    for(int i=0; i<Nbyte; i++){
                        //datacon[2*i]   = PET3RawdataConvert(data[2*i]);
                        //datacon[2*i+1] = PET3RawdataConvert(data[2*i+1]);
                        datacon[2*i]   =data[2*i];
                        datacon[2*i+1] =data[2*i+1];
                        dataeff[i] = ((datacon[2*i] )<<32) | (datacon[2*i+1] &0xFFFFFFFF);
                    }

                    for(int i=0; i<Nbyte; i++){
    	                   
                        kind         = (dataeff[i]>>63) & 0x1;
                        tdiff        = (dataeff[i]>>52) & 0x7FF;
                        Adc1         = (dataeff[i]>>42) & 0x3FF;
                        Adc2         = (dataeff[i]>>32) & 0x3FF;
                        Tx1          = (dataeff[i]>>22) & 0x3FF;
                        Ax1          = (dataeff[i]>>16) & 0x3F;
                        Tx2          = (dataeff[i]>>6)  & 0x3FF;
                        Ax2          = (dataeff[i])     & 0x3F;
               
                        //cout<<"dataeff "<<hex<<dataeff[i]<<endl;
                        x1 = GetX(Tx1);
                        x2 = GetX(Tx2);
                        y1 = GetY(Tx1);
                        y2 = GetY(Tx2);
                        z1 = GetZ(Ax1);
                        z2 = GetZ(Ax2);
                        
                        
                        if(((Tx1)+(Tx2))/2. ==159.5|| ((Tx1)+(Tx2))/2. ==463.5 ){
                            theta = PI/2.;
                            S = x1; 
                        }
                        else{
                            theta = atan((y1-y2)/(x1-x2))+PI/2.;
                            //if(theta>=0)theta= theta +PI/2.;
                            //else        theta=PI/2.+theta;
                            S = sqrt( ((x1+x2)/2.)*((x1+x2)/2.)+  ((y1+y2)/2.)*((y1+y2)/2.)  );
                            coeff = -x1*(y2-y1)/(x2-x1)+y1;
                            if(coeff<=0) S = -1*S;
                        }   
                        
                        
                        if(kind!=event_type &&event_type<2) continue;
    	                if(Save_Raw)raw->Fill();
                        if(Ax1>=NZX || Ax2>=NZX || Tx1>=NCX||Tx2>=NCX)continue;    
                        if(Adc1==0 || Adc1>=Nbin-1 || Adc2==0||Adc2>=Nbin-1 )continue;    
                        
                        if(Save_Hist){
                            Axial1 ->Fill(Ax1); 
                            Axial2 ->Fill(Ax2); 
                            Axial ->Fill(Ax1); 
                            Axial ->Fill(Ax2); 
                            TransAxial1 ->Fill(Tx1); 
                            TransAxial2 ->Fill(Tx2); 
                            TransAxial ->Fill(Tx1); 
                            TransAxial ->Fill(Tx2); 
                            Energy1->Fill(Adc1); 
                            Energy2->Fill(Adc2); 
                            Energy->Fill(Adc1); 
                            Energy->Fill(Adc2); 
                       
                            Countmap1 ->Fill(Tx1*1./TxInBank,Ax1*1./AxInBank);
                            Countmap2 ->Fill(Tx2*1./TxInBank,Ax2*1./AxInBank);
                            Countmap  ->Fill(Tx1*1./TxInBank,Ax1*1./AxInBank);
                            Countmap  ->Fill(Tx2*1./TxInBank,Ax2*1./AxInBank);
                      
                            Energy->Fill(Adc1);
                            Energy->Fill(Adc2);
                            EnergyvsRing->Fill(Ax1,Adc1);
                            EnergyvsRing->Fill(Ax2,Adc2);
                      
                            EnergyvsAngle->Fill(Tx1,Adc1);
                            EnergyvsAngle->Fill(Tx2,Adc2);
                      
                            index1 = Ax1/AxInSubMod;
                            index2 = Tx1/TxInSubMod;
                            EnergySubMod[index1][index2]->Fill(Adc1);
                            index1 = Ax2/AxInSubMod;
                            index2 = Tx2/TxInSubMod;
                            EnergySubMod[index1][index2]->Fill(Adc2);
                             
                            index1 = Tx1/TxInMod;
                            index2 = Tx2/TxInMod;
                            EnergyModule[index1]->Fill(Adc1);
                            EnergyModule[index2]->Fill(Adc2);
                            Txmax=max(Tx1/TxInMod,Tx2/TxInMod);
                            Txmin=min(Tx1/TxInMod,Tx2/TxInMod);
                            if(Txmax-Txmin==halfNmod){
                                if(Tx1/TxInMod == Txmin){
                                    Timecouple[Txmin]->Fill(tdiff);
                                }
                                if(Tx2/TxInMod == Txmin){
                                    Timecouple[Txmin]->Fill(-tdiff);
                                }
                            }
                            Time->Fill(-tdiff);
                            Time->Fill( tdiff);

                            TimevsAngle ->Fill(Tx1,tdiff);
                            TimevsAngle ->Fill(Tx2,-tdiff);
                            
                            TimevsRing ->Fill(Ax1,tdiff);
                            TimevsRing ->Fill(Ax2,-tdiff);
                            if(Save_Sino)ThetavsS->Fill(S,theta);
                        }
                       
                       
                    }
                }
            }
        }
        if(data_type==single_data_type){
            while(fin.good() && Readcounts<= ReadDepth ){
                fin.read((char*)(&data),4*Nslotbyte);
                Readcounts += Nslotbyte;
                
                for(int i=0; i<Nbyte; i++){
                    //datacon[2*i]   = PET3RawdataConvert(data[2*i+1]);
                    //datacon[2*i+1] = PET3RawdataConvert(data[2*i]);
                    datacon[2*i]   = (data[2*i+1]);
                    datacon[2*i+1] = (data[2*i]);
                    dataeff[i] = ((datacon[2*i+1] )<<32) | (datacon[2*i] &0xFFFFFFFF);
                    //cout<<hex<<datacon[2*i+1]<<" "<<datacon[2*i]<<endl;
                }
                for(int i=0;i<Nbyte;i++){
                    if(( (dataeff[i]>>48) & 0xFFFF) == 0x8000){
                        timezone = dataeff[i]&0xFFFFFFFF ;
                    }
                    else if(( (dataeff[i]>>48) & 0xFFFF) == 0){
                        multhit      = (dataeff[i]>>47) & 0x1;
                        ovflow       = (dataeff[i]>>46) & 0x1;
                        coarse       = (dataeff[i]>>34) & 0xFFF;
                        fine         = (dataeff[i]>>26) & 0x1F;
                        Adc          = (dataeff[i]>>16) & 0x3FF;
                        Tx           = (dataeff[i]>>6)  & 0x3FF;
                        Ax           =  dataeff[i]      & 0x3F;
                        timetick = coarse*coarseLSB+fine*fineLSB;  
                        //cout<<"  Adc "<<Adc<<"  Tx "<<Tx<<"  Ax "<<" coarse "<<coarse<<"  fine"<<fine<<endl;
                        if(Save_Raw)raw->Fill();
    	                if(Save_Hist){
                            if(Ax<NZX&&Tx<NCX){ 
                                Axial->Fill(Ax);
                                TransAxial->Fill(Tx);
                                Countmap ->Fill(Tx*1./TxInBank,Ax*1./AxInBank);
                                EnergyvsRing->Fill(Ax,Adc);
                                EnergyvsAngle->Fill(Tx,Adc);
                                Energy ->Fill(Adc);
                        
                                index1 = Ax/AxInSubMod;
                                index2 = Tx/TxInSubMod;
                                EnergySubMod[index1][index2]->Fill(Adc);
                                index = Tx/TxInBank;
                                EnergyModule[index]->Fill(Adc);
                            }
                        }
          
                        
                    }
                    else{
                        fakerflag++;
                        if(fakerflag>=FakerFlagThreshold)ShowWrongDataType(1);
                    }
                }

            }
        }

        fin.close();
	    cout<<" file : "<< path_to_input_files[MM].c_str() <<"   read done !!!"<<endl;
    }      
    
    rfile->cd();
    if(Save_Raw)raw->Write();
    rfile->Close();	



    cout<<"---------------- Start Saving .root files --------------"<<endl;

    TFile *hfile =new TFile(root_hist_file.c_str(),"recreate");
    hfile ->cd();
    if(Save_Hist){
       Axial->Write(); 
       TransAxial->Write(); 
       Countmap->Write();
       EnergyvsRing->Write();
       EnergyvsAngle->Write();
       Energy->Write();
       for(int i=0;i<Nmod;i++){
            EnergyModule[i]->Write();
       }
      for(int i=0 ;i<NSubAx;i++){
        for(int j=0 ;j<NSubTx;j++){
            EnergySubMod[i][j]->Write();
        }
      }
       if(data_type==coin_data_type){
           Countmap1->Write();
           Countmap2->Write();
           Axial1->Write();
           Axial2->Write();
           TransAxial1->Write();
           TransAxial2->Write();
           Energy1->Write();
           Energy2->Write();
           Time->Write();
           for(int i=0;i<halfNmod;i++){
               Timecouple[i]->Write();
           }
           
           TimevsAngle->Write();
           TimevsRing->Write();
           if(Save_Sino)ThetavsS->Write();
       }
    }
    hfile->Close();






    cout<<"---------------- Start saving pdf/png---------------"<<endl;

    TCanvas *c1=new TCanvas("c1","c1",1500,1000);
    TCanvas *c2=new TCanvas("c2","c2",1500,1000);
    MycanvasSetting(c1,0.15,0.14,0.15,0.15);
    MycanvasSetting(c2,0.15,0.14,0.15,0.15);

    //TVirtualPad *tv[6], *tvp[3],*tvpad[2];
     if(data_type==coin_data_type){
         if(Save_Sino){
             c1->Clear();
             c1->cd();
             ThetavsS->Draw("colz");
             if(Save_Pdf)c1->Print(ThetavsSPdf.c_str());
             if(Save_Png)c1->Print(ThetavsSPng.c_str());
         }
         c1->Clear();
         c1->Divide(2,3);
         c1->cd(1);
         Axial1->Draw();
         c1->cd(2);
         TransAxial1->Draw();
         c1->cd(3);
         Axial2->Draw();
         c1->cd(4);
         TransAxial2->Draw();
         c1->cd(5);
         Axial->Draw();
         c1->cd(6);
         TransAxial->Draw();
         if(Save_Pdf)c1->Print(AxialANDTransPdf.c_str());
         if(Save_Png)c1->Print(AxialANDTransPng.c_str());

         c1->Clear();
         c1->Divide(1,3);
         c1->cd(1);
         Energy1->Draw();
         c1->cd(2);
         Energy2->Draw();
         c1->cd(3);
         Energy->Draw();
         if(Save_Pdf)c1->Print(EnergyTotalPdf.c_str());
         if(Save_Png)c1->Print(EnergyTotalPng.c_str());
  
         c2->Clear();
         c2->cd();
         Time->Draw();
         if(Save_Pdf)c2->Print(TimeTotalPdf.c_str());
         if(Save_Png)c2->Print(TimeTotalPng.c_str());
  
         c2->Clear();
         c2->Divide(1,3);
         c2->cd(1);
         Countmap1->Draw("colz");
         c2->cd(2);
         Countmap2->Draw("colz");
         c2->cd(3);
         Countmap->Draw("colz");
         if(Save_Pdf)c2->Print(CountMapPdf.c_str());
         if(Save_Png)c2->Print(CountMapPng.c_str());
  
         c1->Clear();
         c1->Divide(1,2);
         c1->cd(1);
         EnergyvsRing ->Draw("colz");
         c1->cd(2);
         EnergyvsAngle ->Draw("colz");
         if(Save_Pdf)c1->Print(Energvs2DimensionsPdf.c_str());
         if(Save_Png) c1->Print(Energvs2DimensionsPng.c_str());
         
         c2->Clear();
         c2->Divide(5,8);
         for(int i=0;i<Nmod;i++){
            c2->cd(i+1);
            EnergyModule[i]->Draw();
         }
         if(Save_Pdf)c2->Print(EnergyModulePdf.c_str());
         if(Save_Png)c2->Print(EnergyModulePng.c_str());
         
         //c1->Clear();
         //for(int i=0;i<Nmod;i++){
         //   c2->cd(i+1);
         //   EnergyModule[i]->Draw();
         //}
         //if(Save_Pdf)c2->Print(EnergyModulePdf.c_str());
         //if(Save_Png)c2->Print(EnergyModulePng.c_str());
         //
         //c2->Clear();
         //c2->Divide(5,8);
         //for(int i=0;i<Nmod;i++){
         //   c2->cd(i+1);
         //   EnergyModule[i]->Draw();
         //}
         //if(Save_Pdf)c2->Print(EnergyModulePdf.c_str());
         //if(Save_Png)c2->Print(EnergyModulePng.c_str());

         c1->Clear();
         c1->Divide(5,4);
         TF1 *fitf[20];
            TText text1[20];
            TText text2[20];
         for(int i=0;i<20;i++){
            sprintf(buf,"fitf%d",i);
            fitf[i]=new TF1("fitf","gaus",-1e4,1e4);
            text1[i].SetNDC(kTRUE);
            text2[i].SetNDC(kTRUE);
            text1[i].SetTextSize(0.06);
            text2[i].SetTextSize(0.06);
            text1[i].SetTextColor(1);
            text2[i].SetTextColor(1);
         }
         float fmean[20] ,fsigm[20];
         TVirtualPad *vpd[20];
         for(int i=0;i<halfNmod;i++){
            vpd[i]=c1->cd(i+1);
            MyPadSetting(vpd[i],0.15,0.14,0.15,0.15);
            vpd[i]->cd();
            double mean = Timecouple[i]->GetBinCenter(Timecouple[i]->GetMaximumBin());
            double rms  = Timecouple[i]->GetRMS();
            fitf[i]->SetParameters(Timecouple[i]->GetMaximum(),mean,rms); 
            fitf[i]->SetParameters(Timecouple[i]->GetMaximum(),mean,rms); 
            Timecouple[i]->Draw();
            Timecouple[i]->Fit(fitf[i],"Q","",mean-rms,mean+rms);
            fmean[i] = fitf[i]->GetParameter(1);
            fsigm[i] = fitf[i]->GetParameter(2);
            sprintf(buf,"mean : %3.2f", fmean[i]);
            sprintf(buf2,"sigma: %3.2f", fsigm[i]);
            fitf[i]->SetParameters(fitf[i]->GetParameter(0),fmean[i],fsigm[i]); 
            fitf[i]->Draw("same");
            //text1->DrawTextNDC(0.15,0.85,buf );
            //text2->DrawTextNDC(0.75,0.85,buf2);
            text1[i].SetText(0.17,0.3,buf );
            text2[i].SetText(0.55,0.3,buf2);
            text1[i].Draw("same");
            text2[i].Draw("same");

         }
         if(Save_Pdf)c1->Print(TimeModulePdf.c_str());
         if(Save_Png)c1->Print(TimeModulePng.c_str());
    
         c2->Clear();
         c2->cd();
         TimevsAngle->Draw("colz");
         if(Save_Pdf)c2->Print(TimevsAnglePdf.c_str());
         if(Save_Png)c2->Print(TimevsAnglePng.c_str());

         c1->Clear();
         c1->cd();
         TimevsRing->Draw("colz");
         if(Save_Pdf)c1->Print(TimevsRingPdf.c_str());
         if(Save_Png)c1->Print(TimevsRingPng.c_str());
     
     
         if(Save_Bin){
             if(Save_Sino){
                 unsigned int bb[304];
                 binfile = fopen(ThetavsSBin.c_str(),"wb");
                 for(int i=0;i<334;i++){
                    for(int j=0;j<304;j++){
                        bb[j] = ThetavsS ->GetBinContent(i+1,j+1);
                    }
                    fwrite(&bb,sizeof(unsigned int),304,binfile);
                 }
                 fclose(binfile);
             }

             unsigned int timec[400];
             binfile = fopen(TimevsAngleBin.c_str(),"wb");
             for(int i=0;i<NCX;i++){
                for(int j=0;j<400;j++){
                    timec[j] = TimevsAngle ->GetBinContent(i+1,j+1);
                }
                fwrite(&timec,sizeof(unsigned int),400,binfile);
             }
             fclose(binfile);
             
             binfile = fopen(TimevsRingBin.c_str(),"wb");
             for(int i=0;i<NCX;i++){
                for(int j=0;j<400;j++){
                    timec[j] = TimevsRing ->GetBinContent(i+1,j+1);
                }
                fwrite(&timec,sizeof(unsigned int),400,binfile);
             }
             fclose(binfile);
             

         }
     
     }
     if(data_type==single_data_type){

        c1->Clear();
        c1->Divide(2,1);
        c1->cd(1);
        Axial ->Draw();
        c1->cd(2);
        TransAxial->Draw();
        if(Save_Pdf)c1->Print(AxialANDTransPdf.c_str());
        if(Save_Png)c1->Print(AxialANDTransPng.c_str());

        c1->Clear();
        c1->cd();
        Energy->Draw();
        if(Save_Pdf)c1->Print(EnergyTotalPdf.c_str());
        if(Save_Png)c1->Print(EnergyTotalPng.c_str());

        c2->Clear();
        c2->cd();
        Time->Draw();
        if(Save_Pdf)c2->Print(TimeTotalPdf.c_str());
        if(Save_Png)c2->Print(TimeTotalPng.c_str());

        c2->Clear();
        c2->Divide(5,8);
        for(int i=0;i<Nmod;i++){
           c2->cd(i+1);
           EnergyModule[i]->Draw();
        }
        if(Save_Pdf)c2->Print(EnergyModulePdf.c_str());
        if(Save_Png)c2->Print(EnergyModulePng.c_str());
        
        c1->Clear();
        c1->Divide(1,2);
        c1->cd(1);
        EnergyvsRing ->Draw("colz");
        c1->cd(2);
        EnergyvsAngle ->Draw("colz");
        if(Save_Pdf) c1->Print(Energvs2DimensionsPdf.c_str());
        if(Save_Png) c1->Print(Energvs2DimensionsPng.c_str());
         
        c2->Clear();
        c2->cd();
        Countmap->Draw("colz");
        if(Save_Pdf) c2->Print(CountMapPdf.c_str());
        if(Save_Png) c2->Print(CountMapPng.c_str());
        
     }

     if(Save_Bin){
         unsigned int aa[NZX];
         FILE *binfile;
         binfile = fopen(CountMapBin.c_str(),"wb");
         for(int i=0;i<NCX;i++){
            for(int j=0;j<NZX;j++){
                aa[j] = Countmap ->GetBinContent(i+1,j+1);
            }
            fwrite(&aa,sizeof(unsigned int),NZX,binfile);
         }
         fclose(binfile);
         
         unsigned int cc[1024];
         binfile = fopen(EnergyvsRingBin.c_str(),"wb");
         for(int i=0;i<NZX;i++){
            for(int j=0;j<Nbin;j++){
                cc[j] = EnergyvsRing ->GetBinContent(i+1,j+1);
            }
            fwrite(&cc,sizeof(unsigned int),Nbin,binfile);
         }
         fclose(binfile);
 
         binfile = fopen(EnergyvsAngleBin.c_str(),"wb");
         for(int i=0;i<NCX;i++){
            for(int j=0;j<Nbin;j++){
                cc[j] = EnergyvsAngle ->GetBinContent(i+1,j+1);
            }
            fwrite(&cc,sizeof(unsigned int),Nbin,binfile);
         }
         fclose(binfile);
     }
    


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

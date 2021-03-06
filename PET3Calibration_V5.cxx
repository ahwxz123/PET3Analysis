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
    cout<<"   PET3Calibration version  1.2  "<<endl;
    cout<<"     :add reading depth flag and option.                         2018/12/25 "<<endl;
    cout<<"     :Adc range 1-1022.                                          2018/12/29 "<<endl;
    cout<<"     :input & output exits ro not.                               2019/01/02 "<<endl;
    cout<<"     :Use ROOT histogram.                                        2019/01/10 "<<endl;
//    cout<<"     :Cocincidence flag 0xFFFFFFFC; 16vevents in one flag;       2019/01/10 "<<endl;
    cout<<"     :Cocincidence flag 0xFFFFFFFC; 32vevents in one flag;       2019/01/11 "<<endl;
    cout<<"     :Big Endian, No need using PET3RawdataConvert .             2019/01/11 "<<endl;
    cout<<"     :Singles no flag.                                           2019/01/11 "<<endl;
    cout<<"     :Header flag 0xA0A, no zero in lut .                        2019/01/13 "<<endl;
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
    int event_type=1;           // 0  delay,  1  prompt,   2  both;
    int calib_mode=0;           // 0  energy, 1  time  ,   2  both;
    int Niter=3;                // time calibration times;
    int Percentage=0x7FFFFFFF ;          ;
    bool Save_Raw=false;
    bool Save_Hist=true;
    bool Save_Pdf=false;
    bool Save_Png=false;
    bool Save_Txt=false;
   
    bool EnergyCalibOn=true;
    bool TimeCalibOn=false;
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
            if(filetest==NULL){
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
                txtfileread.open(input_file.c_str(),ios::in);
                if(! txtfileread){
                    cerr<<"*********  input TXT file does not exist!!!!  Please recheck the path or TXT file name !!!"<<endl;
                    exit(1);
                }

                while(txtfileread.good()){
                    getline(txtfileread, fname);
                    path_to_input_files.push_back(fname);
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
        else if(option=="-txt"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara =="yes"||optionpara=="YES"||optionpara=="Yes"||optionpara=="Y"||optionpara=="y"){
                    Save_Txt = true;                 
                }
                else if(optionpara =="no"||optionpara=="NO"||optionpara=="No"||optionpara=="N"||optionpara=="n"){
                    Save_Txt = false;                 
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
        else if(option=="-calib"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                if(optionpara == "both"){
                    calib_mode = 2;     
                    EnergyCalibOn=true;
                    TimeCalibOn=true;
                }
                else if(optionpara == "energy"||optionpara == "e"){
                    calib_mode = 0;                 
                    EnergyCalibOn=true;
                    TimeCalibOn=false;
                }
                else if(optionpara == "time"||optionpara == "t"){
                    calib_mode = 1;                 
                    EnergyCalibOn=false;
                    TimeCalibOn=true;
                }
                else{
                    cerr<<"******** invalid paramters for option:"<<option<<endl;
                    exit(1);
                }
            } 
    
            i++;
        }
    
        else if(option=="-iter"){
            if(argv[i+1]==NULL){
                cerr<<"******** argument missing for option:"<<option<<endl;
                cerr<<"******** using default option !"<<endl;
            }
            else{
                string optionpara = (string) argv[i+1];
                Niter = atoi(optionpara.c_str());
                if(Niter<0 ){
                    cerr<<"******** invalid paramters for option:"<<option<<endl;
                    cerr<<"******** it should be a >0 integral number     !!!  recommend among 1~8  ..."<<endl;
                    exit(1);
                }
                else{
                    cout<<"******** TIME calibration iteration times ---> "<<Niter <<"  "<<endl;
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
                Percentage = atoi(optionpara.c_str());
                if(Percentage<0 ){
                    cerr<<"******** invalid paramters for option:"<<option<<endl;
                    cerr<<"******** it should be a >0 integral number     !!!"<<endl;
                    exit(1);
                }
                else{
                    cout<<"******** Reading "<<Percentage <<"\%  events"<<endl;
                }
            } 
    
            i++;
        }
        else if(option=="-all"||option=="-All"||option=="-ALL"){
            Save_Hist=true;
            Save_Pdf=true;
            Save_Png=true;
            Save_Txt=true;
            
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
      cout<<"--------------------   Starting ----------------"<<endl;

      gStyle->SetOptStat(0);
	  gStyle->SetPadLeftMargin(0.15);
	  gStyle->SetPadRightMargin(0.15);
	  gStyle->SetPadTopMargin(0.12);
	  gStyle->SetPadBottomMargin(0.15);
      gErrorIgnoreLevel = kError, kBreak, kSysError, kFatal;
      TStopwatch *ts1=new TStopwatch();
	  ts1->Start(); 

	  TDatime *tdate =new TDatime();
      int systime = tdate->GetTime();
      int sysdate = tdate->GetDate();
      char buf[100];

      const unsigned int CoincidenceFlag=0xFFFFFFFC;
	  const int energylow=2;
	  const int energyhigh=1020;
	  const int Ntimebin=2048;
      unsigned int fakerflag=0;
      int counts;
      int index;
      int index1=0;
      int index2=0;
      int bankindex1=0;
      int bankindex2=0;
      
      unsigned int eindex1, eindex2;
	  //int bankpeak1[NPixels],bankpeak2[NPixels];
	  int temp1,temp2;
	  double cc;

      int Adc1= -1,Adc2 = -1, Tdiff= -1,  DeltaT = -1;
      int Ax1= -1,  Ax2= -1,    Tx1= -1,  Tx2= -1;
      int co=-1,co1=-1,co2=-1,kind;
	  int tdiff;
      double tdiff_cor1,tdiff_cor2;
      int multhit, ovflow;
      int Ax,Tx, Adc, coarse,fine,timezone;
      int Txmax,Txmin;
      double timetick;
      double coarseLSB=10.;       //ns;
      double fineLSB=0.3125;      //ns;
     
      double peakt0,peakt1,sigmat;


      double x1,y1,x2,y2,z1,z2;
      double theta,S,coeff;
      
      
      int npeak;


      int Nshift=0;
      unsigned int Nslotbyte=64;      
      unsigned int Nbyte=Nslotbyte/2;      
	  int data[Nslotbyte];
      Long64_t datacon[Nslotbyte] ;
      Long64_t dataeff[Nbyte] ;
      Long64_t datatest;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                              histgream define                                               ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      cout<<"--------------------   Defination ----------------"<<endl;
      TF1 *fitfun[NPixels];

      TH1S *Energy[NPixels];
      TH1S *Time[NPixels];
 
      for(int i=0; i<NPixels; i++){
          sprintf(buf,"fun%d",i);
          fitfun[i] =new TF1(buf,"gaus",-1e3,1e3);
          if(EnergyCalibOn){
              sprintf(buf,"Energy_Ax%d_Tx%d",i%NZX,i/NZX);
              Energy[i] = hist1dS(buf,Nbin,0,Nbin,"Energy (ch)","Entries");
          }
          if(TimeCalibOn){
              sprintf(buf,"Time_Ax%d_Tx%d", i%NZX,i/NZX);
              Time[i] = hist1dS(buf,2*HalfWindow,-HalfWindow,HalfWindow,"Time (ch)","Entries");
          }   
      }  
    
      TH1F *EnergyTotal   = hist1d("EnergyTotal",Nbin,0,Nbin,"Energy (ch)","Entries"); 
      TH2F *EnergyvsAngle = hist2d("EnergyvsAngle",NCX,0,NCX,Nbin,0,Nbin," Angle ","Energy (ch)");
      TH2F *EnergyvsRing  = hist2d("EnergyvsRing" ,NZX,0,NZX,Nbin,0,Nbin," Ring ","Energy (ch)");

      TH1F *TimeRes     = hist1d("Timeres",2*HalfWindow,-HalfWindow,HalfWindow,"Time (ch)","Entries"); 
      TH2F *TimevsAngle = hist2d("TimevsAngle",NCX,0,NCX,2*HalfWindow,-HalfWindow,HalfWindow," Angle ","Time (ch)");
      TH2F *TimevsRing  = hist2d("TimevsRing" ,NZX,0,NZX,2*HalfWindow,-HalfWindow,HalfWindow," Ring ","Time (ch)");

      TH1F *TimeRes_cor     = hist1d("Timeres_cor",2*HalfWindow,-HalfWindow,HalfWindow,"Time (ch)","Entries"); 
      TH2F *TimevsAngle_cor = hist2d("TimevsAngle_cor",NCX,0,NCX,2*HalfWindow,-HalfWindow,HalfWindow," Angle ","Time (ch)");
      TH2F *TimevsRing_cor  = hist2d("TimevsRing_cor" ,NZX,0,NZX,2*HalfWindow,-HalfWindow,HalfWindow," Ring ","Time (ch)");
      TH2F *EnergyLutMap    = hist2d("EnergyLutMap" ,NCX,0,NCX,NZX,0,NZX," Transxial ","Axial");
      TH2F *TimeLutMap      = hist2d("TimeLutMap" ,NCX,0,NCX,NZX,0,NZX," Transxial ","Axial");

      TimeRes->SetLineColor(1);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                              string for pdf/png                                              ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ////////////// string name for saving pdf/png;
      string CalibResultsName        = path_to_output_file+"/"+"PET3Calibration_lut_Iteration"+to_string(Niter)+".root";
      string CalibHistName           = path_to_output_file+"/"+"PET3Calibration_histogram_Iteration"+to_string(Niter)+".root";
      string EnergyLUTFileName       = path_to_output_file+"/"+"PET3EnergyLUT.binary";
      string TimeLUTFileName         = path_to_output_file+"/"+"PET3TimeLUT.binary";
      string EnergySpectrumName      = path_to_output_file+"/"+"Energy_pixel_spectrum_Hist";
      string TimeSpectrumName        = path_to_output_file+"/"+"Time_pixel_spectrum_Hist";
      string EnergyTotalName         = path_to_output_file+"/"+"Energy_Total_Hist";
      string TimeTotalName           = path_to_output_file+"/"+"Time_Total_Hist";
      string Energy2DName            = path_to_output_file+"/"+"Energy_versusAngleRing";
      string Time2DName              = path_to_output_file+"/"+"Time_versusAngleRing";
      
      string EnergySpectrumPdf1      = EnergySpectrumName+".pdf[";
      string EnergySpectrumPdf2      = EnergySpectrumName+".pdf";
      string EnergySpectrumPdf3      = EnergySpectrumName+".pdf]";
      string TimeSpectrumPdf1        = TimeSpectrumName+".pdf[";
      string TimeSpectrumPdf2        = TimeSpectrumName+".pdf";
      string TimeSpectrumPdf3        = TimeSpectrumName+".pdf]";
      
      string Energy2DPng             = Energy2DName+".png";
      string Time2DPng               = Time2DName+".png";
      string EnergyTotalPng          = EnergyTotalName+".png";
      string TimeTotalPng            = TimeTotalName+".png";
      
      string EnergyLUTTxt            = EnergySpectrumName+".txt";
      string TimeLUTTxt              = TimeSpectrumName+".txt";
     
/////////////////////////////////////////////////////////////////////////////////////////////     
      cout<<"--------------------   Reading data files ----------------"<<endl;

	int size = path_to_input_files.size();
    cout<<"----Reading inputs files :"<<endl;      
    for(int i=0;i<size;i++){
	     cout<<"    "<<path_to_input_files[i].c_str()<<endl;      
	}
    int PercantRead;
    ULong64_t Readcounts;
    ULong64_t ReadLimit= Percentage;
	  
	ifstream fin;
    for(int MM=0; MM<size; MM++){
        fin.open(path_to_input_files[MM].c_str(),ios::binary);
        fin.seekg( 0,ios::end);
        fin.seekg( Nshift*4,ios::beg);
        Readcounts=0;
        if(data_type==coin_data_type){
            while(fin.good()  && Readcounts<= ReadLimit ){
                fin.read((char*)(&datatest), sizeof(unsigned int));
                Readcounts +=1;
                //datatest = PET3RawdataConvert(datatest);
                if(CoincidenceFlag!= datatest){
                    fakerflag++;
                    if(fakerflag>=FakerFlagThreshold)ShowWrongDataType(1);
                    continue;
                }
                else{ 
                    fin.read((char*)(&data),Nslotbyte*sizeof(unsigned int));
                    Readcounts +=Nslotbyte;
                    for(int i=0; i<Nbyte; i++){
                        //datacon[2*i]   = PET3RawdataConvert(data[2*i]);
                        //datacon[2*i+1] = PET3RawdataConvert(data[2*i+1]);
                        datacon[2*i]   = data[2*i];
                        datacon[2*i+1] = data[2*i+1];
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
                    
                        if(kind!=event_type &&event_type<2) continue;
                        if(Ax1>=NZX||Ax2>=NZX||Tx1>=NCX||Tx2>=NCX) continue;
                        if(Adc1<=1 ||Adc1>=Nbin-1 ||Adc2<=1 ||Adc2>=Nbin-1)continue;
                        if(tdiff>=HalfWindow)continue;   // time is too large;
                        index1= Tx1*NZX +Ax1;
                        index2= Tx2*NZX +Ax2;
                        if(EnergyCalibOn){
                            Energy[index1]->Fill(Adc1);
                            Energy[index2]->Fill(Adc2);
                            
                            EnergyTotal->Fill(Adc1);
                            EnergyTotal->Fill(Adc2);
                    
                            EnergyvsAngle->Fill(Tx1,Adc1);
                            EnergyvsAngle->Fill(Tx2,Adc2);
                            EnergyvsRing->Fill(Ax1,Adc1);
                            EnergyvsRing->Fill(Ax2,Adc2);
                    }

                        if(TimeCalibOn){
                            Time[index1]->Fill( -tdiff); 
                            Time[index2]->Fill(  tdiff); 

                            TimevsAngle->Fill(Tx1,-tdiff);
                            TimevsAngle->Fill(Tx2, tdiff);
                            TimevsRing ->Fill(Ax1,-tdiff);
                            TimevsRing ->Fill(Ax2, tdiff);
                            
                            TimeRes ->Fill( tdiff);
                            TimeRes ->Fill(-tdiff);
                       }
                    }
                }
            }
        }
        if(data_type==single_data_type){
            while(fin.good() && Readcounts<= ReadLimit ){
                fin.read((char*)(&data),4*Nslotbyte);
                Readcounts +=Nslotbyte;
                 
                for(int i=0; i<Nbyte; i++){
                    //datacon[2*i]   = PET3RawdataConvert(data[2*i+1]);
                    //datacon[2*i+1] = PET3RawdataConvert(data[2*i]);
                    datacon[2*i]   = data[2*i+1];
                    datacon[2*i+1] = data[2*i];
                    dataeff[i] = ((datacon[2*i+1] )<<32) | (datacon[2*i] &0xFFFFFFFF);
                    //cout<<hex<< data[2*i]<<"   "<<data[2*i+1] <<endl;
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
                        if(Ax>=NZX||Tx>=NCX) continue;
                       
                        index = Tx*NZX*Nbin+Ax*Nbin+Adc;
                        Energy[index]->Fill(Adc);
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
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                 Fill histogram                                                      ////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(int i=0;i<NPixels;i++){
        if(TimeCalibOn){
            double coeff0 = Time[i]->GetMaximum();
            //double coeff1 = Time[i]->GetBinCenter(Time[i]->GetMaximumBin());
            double coeff1 = Time[i]->GetMean();
            double coeff2 = 20.;
            fitfun[i]->SetParameter(0,coeff0);
            fitfun[i]->SetParameter(1,coeff1);
            Time[i]->Fit(fitfun[i],"QNC","",coeff1-5*coeff2,coeff1+5*coeff2);
            coeff1 = fitfun[i]->GetParameter(1);
            coeff2 = fitfun[i]->GetParameter(2);
            Time[i]->GetXaxis()->SetRangeUser(coeff1-5*coeff2,coeff1+5*coeff2);
            Time[i]->Fit(fitfun[i],"QNC","",coeff1-5*coeff2,coeff1+5*coeff2);
            tmean[i] =tmean[i]+ fitfun[i]->GetParameter(1);
            //if(i==15000)cout<<i<<"  "<< Time[i]->Integral()<<"  tmean "<<tmean[i]<<" TIME RMS "<<Time[i]->GetRMS()<<endl;;
        }
    }
    



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                 Calibration                                                         ////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
      cout<<"--------------------   Starting Calibration ----------------"<<endl;
    double *xpos;
    if(EnergyCalibOn){
        TSpectrum *sp=new TSpectrum(3,25);
        
        for(int i=0;i<NPixels;i++){
            npeak  =  sp-> Search(Energy[i],25.,"",0.6);
            xpos   =  sp-> GetPositionX();
            sort(xpos+0,xpos+npeak);
            peak0[i] = xpos[npeak -1];
			fitfun[i]->SetParameters(Energy[i]->GetBinContent(xpos[npeak-1]), peak0[i], 25.);
            Energy[i]->Fit(fitfun[i],"QN","",peak0[i]-30,peak0[i]+30 );
            peak0[i] = fitfun[i]->GetParameter(1);
            res[i] = fitfun[i]->GetParameter(2);
          
        }
    }
    
    if(TimeCalibOn){

        for(int NN=0;NN<Niter;NN++){

	    //cout<<" -----   Time calibration Iteration -->"<<NN+1<<" begins ....."<<endl;
            TimeRes_cor ->Reset();
            for(int i=0; i<NPixels; i++){
		        Time[i]->Reset();
            } 
            for(int MM=0; MM<size; MM++){
                fin.open(path_to_input_files[MM].c_str(),ios::binary);
                fin.seekg( 0,ios::end);
                fin.seekg( Nshift*4,ios::beg);
                Readcounts=0;
                if(data_type==coin_data_type){
                    while(fin.good()  && Readcounts<= ReadLimit ){
                        fin.read((char*)(&datatest),sizeof(unsigned int));
                        Readcounts +=1;
                        //datatest = PET3RawdataConvert(datatest);
                        if(CoincidenceFlag!=datatest){
                            fakerflag++;
                            if(fakerflag>=FakerFlagThreshold)ShowWrongDataType(1);
                            continue;
                        }
                        else{
                            fin.read((char*)(&data),Nslotbyte*sizeof(unsigned int));
                            Readcounts +=Nslotbyte;
                            for(int i=0; i<Nbyte; i++){
                                //datacon[2*i]   = PET3RawdataConvert(data[2*i]);
                                //datacon[2*i+1] = PET3RawdataConvert(data[2*i+1]);
                                datacon[2*i]   = data[2*i];
                                datacon[2*i+1] = data[2*i+1];
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
                            
                                if(kind!=event_type &&event_type<2) continue;
                                if(Ax1>=NZX||Ax2>=NZX||Tx1>=NCX||Tx2>=NCX) continue;
                                if(Adc1<=1||Adc1>=Nbin-1||Adc2<=1||Adc2>=Nbin-1) continue;
                              
                                if(tdiff>=HalfWindow)continue;   // time is too large;
                                index1 = Tx1*NZX +Ax1;
                                index2 = Tx2*NZX +Ax2;
                                tdiff_cor1 = static_cast<double>(-tdiff) -tmean[index1];
                                tdiff_cor2 = static_cast<double>( tdiff) -tmean[index2];
                                
                                if(abs(tdiff_cor1)>=HalfWindow)continue;   // time is too large;
                                if(abs(tdiff_cor2)>=HalfWindow)continue;   // time is too large;
                                //cout<<"  tdiff "<<tdiff<<"  tdiff_cor1 "<<tdiff_cor1+HalfWindow<< "  tdiff_cor2 "<<tdiff_cor2+HalfWindow<<endl;
                                
                                Time[index1]->Fill( tdiff_cor1);
                                Time[index2]->Fill( tdiff_cor2);

                                if(NN==Niter-1){
                                    TimevsAngle_cor->Fill(Tx1,tdiff_cor1);
                                    TimevsAngle_cor->Fill(Tx2,tdiff_cor2);
                                    TimevsRing_cor ->Fill(Ax1,tdiff_cor1);
                                    TimevsRing_cor ->Fill(Ax2,tdiff_cor2);
                                    TimeRes_cor ->Fill(tdiff_cor1);
                                    TimeRes_cor ->Fill(tdiff_cor2);
                                }
                            }
                        }
                }
            }      
            fin.close();
	        //cout<<"       reread file : "<< path_to_input_files[MM].c_str() <<" "<<endl;
        }

        //ofstream fout;
        //sprintf(buf,"%s/lut%d.txt",path_to_output_file.c_str(),NN);
        //fout.open(buf);
        for(int i=0;i<NPixels;i++){
            double coeff0 = Time[i]->GetMaximum();
            //double coeff1 = Time[i]->GetBinCenter(Time[i]->GetMaximumBin());
            double coeff1 = Time[i]->GetMean();
            double coeff2 = 20.;
            if(abs(coeff1)<200){
            //fitfun[i]->SetParameters(coeff0,coeff1,coeff2);
            fitfun[i]->SetParameter(0,coeff0);
            fitfun[i]->SetParameter(1,coeff1);
            Time[i]->Fit(fitfun[i],"QNC","",coeff1-5*coeff2,coeff1+5*coeff2);
            coeff1 = fitfun[i]->GetParameter(1);
            coeff2 = fitfun[i]->GetParameter(2);
            Time[i]->GetXaxis()->SetRangeUser(coeff1-5*coeff2,coeff1+5*coeff2);
            Time[i]->Fit(fitfun[i],"QNC","",coeff1-5*coeff2,coeff1+5*coeff2);
            tmean[i] =tmean[i]+ fitfun[i]->GetParameter(1);
            }
        }
      
    /* 
       string savetest                = path_to_output_file+"/"+"PET3Calibration_COR2times_histogram_iter"+to_string(NN)+".root";
       TFile *stfile = new TFile(savetest.c_str(), "recreate");
       stfile->cd();
       for(int i=0;i<NPixels;i++){
           Time[i] ->Write();
       }
       TimeRes_cor->Write();
       stfile->Close();
*/
    cout<<"------ time Calibration iteration "<< NN+1<<" has been done ..."<<endl;
      }


    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                              Save LUT                                                  ////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      cout<<"--------------------   Saving results ----------------"<<endl;

    const unsigned int HeaderFlag=0xA0A;
    double gain,shift;
    double tdelay;
    unsigned short par1;
    short par2,tpar;
	unsigned short aa;
	short bb;
	unsigned int zeros=0;
	FILE *energylut, *timelut;
	energylut=fopen(EnergyLUTFileName.c_str(),"wb");
	timelut  =fopen(TimeLUTFileName.c_str(),"wb");

	unsigned int energylutflag=0;
	fwrite(&energylutflag,sizeof(int),1,energylut);
	fwrite(&energylutflag,sizeof(int),1,energylut);
	fwrite(&energylutflag,sizeof(int),1,energylut);
	energylutflag = 0xFFFF;
	energylutflag = PET3RawdataConvert(energylutflag);
	fwrite(&energylutflag,sizeof(int),1,energylut);
	
	unsigned int timelutflag=0;
    fwrite(&timelutflag,sizeof(int),1,timelut);
	fwrite(&timelutflag,sizeof(int),1,timelut);
	fwrite(&timelutflag,sizeof(int),1,timelut);
	timelutflag = 0xFFFE;
	timelutflag = PET3RawdataConvert(timelutflag);
	fwrite(&timelutflag,sizeof(int),1,timelut);

    TFile *rf=new TFile( CalibResultsName.c_str(),"recreate");
    TTree *tr=new TTree("calib","calib");
    if(Save_Hist){
        tr->Branch("Ax"   , &Ax, "Ax/I");
        tr->Branch("Tx"   , &Tx, "Tx/I");
        tr->Branch("peak0", &peakt0, "peakt0/D"  );
        tr->Branch("peak1", &peakt1, "peakt1/D"  );
        tr->Branch("sigma", &sigmat, "sigmat/D"  );
        tr->Branch("gain",  &gain,  "gain/D"   );
        tr->Branch("shift", &shift, "shift/D"  );
        tr->Branch("par1",  &par1,  "par1/s"   );
        tr->Branch("par2",  &par2,  "par2/S"   );
        tr->Branch("tdelay",&tdelay,  "tdelay/D"   );   //LSB;
        tr->Branch("tpar",  &tpar,  "tpar/S"   );
    }
    
    
    
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
   			aa=HeaderFlag;
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
            tpar=HeaderFlag;
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
   			if(EnergyCalibOn){
                fwrite(&wxz,sizeof(unsigned int),1,energylut);
   		        EnergyLutMap ->SetBinContent(Tx+1,Ax+1,gain);	
            }
            tpar= static_cast<short>(tdelay);
            wxz = static_cast<unsigned int>(tpar&0xFFFF);
            wxz = PET3RawdataConvert(wxz);
            if(TimeCalibOn){
                fwrite(&wxz,sizeof(unsigned int),1,timelut);
                TimeLutMap ->SetBinContent(Tx+1,Ax+1,tdelay);	
   			}
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
    
    
    rf->cd();
    
    rf = tr->GetCurrentFile();
    tr->Write();
    rf->Close();	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if(Save_Hist){
    	TFile *rfile =new TFile(CalibHistName.c_str(),"recreate");
        

	rfile->cd();
        TimeRes->Write(); 
        TimevsAngle->Write();
        TimevsRing ->Write();
        
	    TimeRes_cor->Write(); 
        TimevsAngle_cor->Write();
        TimevsRing_cor ->Write();
        if(EnergyCalibOn){
            for(int i=0; i<NPixels;i++){
                Energy[i] ->Write();
            }
        }
        if(TimeCalibOn){
            for(int i=0; i<NPixels;i++){
                Time[i] ->Write();
            }
        }
        EnergyLutMap ->Write();
        TimeLutMap->Write();
	rfile->Close();	
    }
   
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c1=new TCanvas("c1","c1",1500,1000);
    TCanvas *c2=new TCanvas("c2","c2",1500,1000);
    MycanvasSetting(c1,0.15,0.14,0.15,0.15);
    MycanvasSetting(c2,0.15,0.14,0.15,0.15);

    TVirtualPad *tv[16];

    if(Save_Png){
        if(EnergyCalibOn){
            c1->Clear();
            c1->cd();
          
            EnergyTotal->Draw();
            c1->Print(EnergyTotalPng.c_str());

            c2->Clear();
            c2->Divide(1,2);
            c2->cd(1);
            EnergyvsAngle->Draw("colz");
            c2->cd(2);
            EnergyvsRing ->Draw("colz");
            c2->Print(Energy2DPng.c_str());
        }
        if(TimeCalibOn){
            c1->Clear();
            c1->cd();
            TimeRes_cor->GetXaxis()->SetRangeUser(-200,200);
            TimeRes_cor->Draw();
            TimeRes->Draw("same");
            c1->Print(TimeTotalPng.c_str());

            c2->Clear();
            c2->Divide(1,2);
            c2->cd(1);
            TimevsAngle->Draw("colz");
            c2->cd(2);
            TimevsRing ->Draw("colz");
            c2->Print(Time2DPng.c_str());
        }
    }

    if(Save_Pdf){
        if(EnergyCalibOn){
            
            c1->Print(EnergySpectrumPdf1.c_str());
            c1->Clear();
            c1->cd();
            for(int i=0;i<NPixels;i++){
                if(i%16==0){
                    c1->Clear();
                    c1->Divide(4,4);
                }
                tv[i] = c1->cd(i%16+1);
                Energy[i]->GetXaxis()->SetRangeUser(10,1020);
                Energy[i]->Draw();
            
                if(i%16==15){
                   c1->Print(EnergySpectrumPdf2.c_str());
                }
            }
            c1->Print(EnergySpectrumPdf3.c_str());
        }  
  
        if(TimeCalibOn){
            c2->Print(TimeSpectrumPdf1.c_str());
            c2->Clear();
            c2->cd();
            for(int i=0;i<NPixels;i++){
                if(i%16==0){
                    c2->Clear();
                    c2->Divide(4,4);
                }
                tv[i] = c2->cd(i%16+1);
                Time[i]->GetXaxis()->SetRangeUser(-200, 200);
                Time[i]->Draw();
            
                if(i%16==15){
                   c2->Print(TimeSpectrumPdf2.c_str());
                }
cout<<" i" <<i<<"----------------------------"<<endl;
            }
            c2->Print(TimeSpectrumPdf3.c_str());
        }  
    }



     if(Save_Txt){
        ofstream energytxtfile;
        ofstream timetxtfile;
  
        if(EnergyCalibOn)energytxtfile.open(EnergyLUTTxt.c_str());
        if(TimeCalibOn)  timetxtfile.open(TimeLUTTxt.c_str());

	    if(EnergyCalibOn){
             energytxtfile<<setw(10)<<0<<endl;;
             energytxtfile<<setw(10)<<0<<endl;;
	         energytxtfile<<setw(10)<<0xFFFF<<endl;;
             energytxtfile<<setw(10)<<0<<endl;;
	    }
	    if(TimeCalibOn){
             timetxtfile<<setw(10)<<0<<endl;;
             timetxtfile<<setw(10)<<0<<endl;;
	         timetxtfile<<setw(10)<<0xFFFE<<endl;;
             timetxtfile<<setw(10)<<0<<endl;;
	    }

       for(int i=0;i<Nbank;i++){
       	for(int j=0;j<AxInSubmod;j++){
       		for(int k=0;k<TxInSubmod;k++){
       			zeros = PET3RawdataConvert(zeros);
       			if(EnergyCalibOn)energytxtfile<<setw(10)<<zeros<<endl;;
       			if(TimeCalibOn)timetxtfile<<setw(10)<<zeros<<endl;;
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
                tdelay =tmean[j+temp1*NZX]; 
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
                   if( peakt0>=1. && peakt1<1024. && sigmat>=1 && sigmat<1024.  ){
       				aa = static_cast<unsigned short>(  (511.)/(peakt0) *8192);
       			    bb =static_cast<unsigned short>(0x0);
                   }
       /////////////////////////////////////// //////////////////////////////////////////////////////////////                
                   
       			if(EnergyCalibOn)energytxtfile<<"Ax "<<Ax<<" Tx "<<Tx<<"  aa "<<aa<<"  bb"<<bb<<endl;;
                if(TimeCalibOn)timetxtfile<<"Ax "<<Ax<<" Tx "<<Tx<<"  tdelay "<<tdelay<<endl;
       		}		
       	}
       
       	for(int j=0;j<AxInSubmod;j++){
       	    	for(int k=0;k<TxInSubmod;k++){
       	    		zeros = PET3RawdataConvert(zeros);
       			    if(EnergyCalibOn)energytxtfile<<setw(10)<<zeros<<endl;;
       			    if(TimeCalibOn)timetxtfile<<setw(10)<<zeros<<endl;;
       	    	}
       	    }
       }
     


     if(EnergyCalibOn)energytxtfile.close();
     if(TimeCalibOn)timetxtfile.close();
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

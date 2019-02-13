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
    cout<<"   PET3Calibration version  0.1  "<<endl;
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
    cout<<"         [-dp     percentage ] "<<endl;
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
    cout<<"         -dp YYY                   : read part percentage(%) of inputs file,YYY: 0-100, default(10)"  <<endl;
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
    string data_type="single";
    const string single_data_type="single";
    const string coin_data_type="coin";
    int event_type=1;           // 0  delay,  1  prompt,   2  both;
    int calib_mode=0;           // 0  energy, 1  time  ,   2  both;
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


      gStyle->SetOptStat(0);
	  gStyle->SetPadLeftMargin(0.15);
	  gStyle->SetPadRightMargin(0.15);
	  gStyle->SetPadTopMargin(0.12);
	  gStyle->SetPadBottomMargin(0.15);
      TStopwatch *ts1=new TStopwatch();
	  ts1->Start(); 

	  TDatime *tdate =new TDatime();
      int systime = tdate->GetTime();
      int sysdate = tdate->GetDate();
      char buf[100];

      int Ehist[Nbin],temp;

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
      int multhit, ovflow;
      int Ax,Tx, Adc, coarse,fine,timezone;
      int Txmax,Txmin;
      double timetick;
      double coarseLSB=10.;       //ns;
      double fineLSB=0.3125;      //ns;
     
      double peakt0,peakt1,sigmat;


      double x1,y1,x2,y2,z1,z2;
      double theta,S,coeff;
      
      
      unsigned int *Edata=new unsigned int[ELUTN];
      memset(Edata,0,sizeof(unsigned int)*ELUTN);
      unsigned int *Tdata=new unsigned int[TLUTN];
      memset(Tdata,0,sizeof(unsigned int)*TLUTN);
      int npeak;



      int Nshift=0;
      unsigned int Nslotbyte=128;      
      unsigned int Nbyte=64;      
	  int data[Nslotbyte];
      ULong64_t datacon[Nslotbyte] ;
      ULong64_t dataeff[Nbyte] ;
      ULong64_t datatest;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                              histgream define                                               ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      TF1 *fitfun[NPixels];

      TH1S *Energy[NPixels];
      TH1S *Time[NPixels];
 
      for(int i=0; i<NPixels; i++){
          sprintf(buf,"fun%d",i);
          fitfun[i] =new TF1(buf,"gaus",-1e3,1e3);
          if(EnergyCalibOn){
              sprintf(buf,"Energy_Ax%d_Tx%d",i%NZX,i/NZX);
              Energy[i] = hist1dS(buf,Nbin,0,Nbin,"Energy","Entries");
          }
          if(TimeCalibOn){
              sprintf(buf,"Time_Ax%d_Tx%d", i%NZX,i/NZX);
              Time[i] = hist1dS(buf,2*HalfWindow,-HalfWindow,HalfWindow,"Time","Entries");
          }   
      }  

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                              string for pdf/png                                              ////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ////////////// string name for saving pdf/png;
      string CalibResultsName        = path_to_output_file+"/"+"PET3Calibration_lut.root";
      string CalibHistName           = path_to_output_file+"/"+"PET3Calibration_histogram.root";
      string EnergyLUTFileName       = path_to_output_file+"/"+"PET3EnergyLUT.binary";
      string TimeLUTFileName         = path_to_output_file+"/"+"PET3TimeLUT.binary";
      string EnergySpectrumName      = path_to_output_file+"/"+"Energy_spectrum_Hist";
      string TimeSpectrumName        = path_to_output_file+"/"+"Time_spectrum_Hist";
      
      string EnergySpectrumPdf1      = EnergySpectrumName+".pdf[";
      string EnergySpectrumPdf2      = EnergySpectrumName+".pdf";
      string EnergySpectrumPdf3      = EnergySpectrumName+".pdf]";
      string TimeSpectrumPdf1        = TimeSpectrumName+".pdf[";
      string TimeSpectrumPdf2        = TimeSpectrumName+".pdf";
      string TimeSpectrumPdf3        = TimeSpectrumName+".pdf]";
      
      string EnergySpectrumPng       = EnergySpectrumName+".png";
      string TimeSpectrumPng         = TimeSpectrumName+".png";
      
      string EnergyLUTTxt            = EnergySpectrumName+".txt";
      string TimeLUTTxt              = TimeSpectrumName+".txt";
      

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
        //ReadLimit =  fin.tellg()*Percentage/100./sizeof(unsigned int)  ;
        //cout<<"----------------------------------------"<<endl;
        //cout<<"--------------------Only"<<ReadLimit<<"  events will read !! "<<endl;
        fin.seekg( Nshift*4,ios::beg);
        Readcounts=0;
        if(data_type==coin_data_type){
           //for(int KK=0;KK<1e5;KK++){
            while(fin.good()  && Readcounts<= ReadLimit ){
                fin.read((char*)(&data),Nslotbyte*sizeof(unsigned int));
                //for(int i=0; i<Nslotbyte; i++){
                //    fin>>data[i];
                //}
                Readcounts +=Nslotbyte;
                for(int i=0; i<Nbyte; i++){
                    datacon[2*i]   = PET3RawdataConvert(data[2*i+1]);
                    datacon[2*i+1] = PET3RawdataConvert(data[2*i]);
                    dataeff[i] = ((datacon[2*i+1] )<<32) | (datacon[2*i] &0xFFFFFFFF);
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
                    index1= Tx1*NZX*Nbin +Ax1*Nbin+Adc1;
                    index2= Tx2*NZX*Nbin +Ax2*Nbin+Adc2;
                    Edata[index1]++;
                    Edata[index2]++;
                    
                    if(tdiff>=HalfWindow)continue;   // time is too large;
                    if(TimeCalibOn){
                        index1 = Tx1*NZX*2*HalfWindow+Ax1*2*HalfWindow+tdiff+HalfWindow;
                        Tdata[index1] ++;
                        index2 = Tx2*NZX*2*HalfWindow+Ax2*2*HalfWindow-tdiff+HalfWindow;
                        Tdata[index2] ++;
                   }
                }
                
            }
        }
        if(data_type==single_data_type){
            while(fin.good() && Readcounts<= ReadLimit ){
                fin.read((char*)(&data),4*Nslotbyte);
                Readcounts +=Nslotbyte;
                 
                for(int i=0; i<Nbyte; i++){
                    datacon[2*i]   = PET3RawdataConvert(data[2*i+1]);
                    datacon[2*i+1] = PET3RawdataConvert(data[2*i]);
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
                       
                        index = Tx*NZX*Nbin+Ax*Nbin+Adc;
                        Edata[index]++;
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
        if(EnergyCalibOn){
            for(int j=0;j<Nbin;j++){
                Energy[i]->SetBinContent(j+1,Edata[i*Nbin+j]);
            }
        }
        if(TimeCalibOn){
            for(int j=0;j<2*HalfWindow;j++){
                Time[i]->SetBinContent(j+1,Tdata[i*2*HalfWindow+j]);
            }
        }
    }




/////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                 Calibration                                                         ////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double *xpos;
    if(EnergyCalibOn){
        TSpectrum *sp=new TSpectrum(3,25);
        
        for(int i=0;i<NPixels;i++){
            npeak  =  sp-> Search(Energy[i],25.,"",0.6);
            xpos   =  sp-> GetPositionX();
            sort(xpos+0,xpos+npeak);
            peak0[i] = xpos[npeak -1];
			fitfun[i]->SetParameters(Energy[i]->GetBinContent(xpos[npeak-1]), peak0[i], 25.);
            Energy[i]->Fit(fitfun[i],"Q","",peak0[i]-30,peak0[i]+30 );
            peak0[i] = fitfun[i]->GetParameter(1);
            res[i] = fitfun[i]->GetParameter(2);
          
        }
    }
    
    if(TimeCalibOn){
        for(int i=0;i<NPixels;i++){
            //double tx = Time[i]->GetBinCenter(Time[i]->GetMaximumBin());        
            //fitfun[i]->SetParameters(Time[i]->GetMaximum(), tx,20);
            //Time[i] ->Fit( fitfun[i],"Q","",tx-50,tx+50  );
            //tmean[i] =  fitfun[i] ->GetParameter(1);
            
            tmean[i] = Time[i]->GetMean(); 
        }
    }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                              Save LUT                                                  ////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
    
    rf->cd();
    
    
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
    
    
    
    tr->Write();
    rf->Close();	

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
    TFile *rfile =new TFile(CalibHistName.c_str(),"recreate");
    if(Save_Hist){
        rfile->cd();
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
    }
    rfile->Close();	
   
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *c1=new TCanvas("c1","c1",1500,1000);
    TCanvas *c2=new TCanvas("c2","c2",1500,1000);
    MycanvasSetting(c1,0.15,0.14,0.15,0.15);
    MycanvasSetting(c2,0.15,0.14,0.15,0.15);

    TVirtualPad *tv[16];
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
                Time[i]->Draw();
            
                if(i%16==15){
                   c2->Print(TimeSpectrumPdf2.c_str());
                }
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
    



delete [] Edata;
delete [] Tdata;
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

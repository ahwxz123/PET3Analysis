

***************************************************************************************************************
 It's a application for Minfounders in MINFOUND-PET3(BGO) files analysis. 
 Auther: Xiaozhuang Wang 
 Email: xiaozhuang.wang@minfound.com 
 Date: Dec. 20th, 2018 

 PET3Analysis Usage 
 Options: 
          -i/-il  path_to_input_file 
          -o      path_of_output 
         [-raw    yes/no ] 
         [-hist   yes/no ] 
         [-pdf    yes/no ] 
         [-png    yes/no ] 
         [-bin    yes/no ] 
         [-sino   yes/no ] 
         [-type   coin/single ] 
         [-event  delay/prompt/both ] 
         [-all    ] 
         [-v      ] 
         [-h      ] 


 [Main setting]: 
         -i  path_to_input_file    : give an input datafile of PET3 datafile(s) Within absolute PATH. 
                                     option -in/-input is same as -i 
         -il path_to_input_txtfile : give an input txt-like file With absolute PATH, txtfile includes PET3 data format files. 
         -o  path_to_output        : give an absolute PATH for output, it must exist before analysis. 
                                     option -out/-output is same as -o  
         -raw  XXX                 : save raw data in rootfiles or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default negetive)
         -hist XXX                 : save histogram in rootfiles or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default positive)
         -pdf  XXX                 : save pdf file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default negetive)
         -png  XXX                 : save png file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default positive)
         -bin  XXX                 : save binary file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default negetive)
         -sino  XXX                : save sino file or not. XXX: y/yes/Yes/YES/Y/n/no/No/NO/N. (default positive)
         -type  YYY                : data file type. YYY: coin/coins/single/singles. (default coincidence data)
         -event YYY                : event kind. YYY:delay/prompt/both. (default both) 
         -all                      : save all outputs. raw,hist,pdf,png,binary. 
         -v                        : print the version of this application 
                                     option -version/--version is same as -v  
         -h                        : print Help information 
                                     option -help/--help is same as -h  
***************************************************************************************************************



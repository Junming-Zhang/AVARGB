Name of the program: AVARGB.m
Title of manuscript:Automatic velocity analysis of seismic data containing multiples based on RGB space mapping 
Author name: Junming Zhang
Affiliation:  Jilin University, China
Address:  Jilin University, People’s Republic of China, China.
Email address: zhangjm19@mails.jlu.edu.cn

The MATLAB package contains one main script, four subfunctions, one packet of seismic data and a library function package. The three subfunctions and the library functions are required for the main script.
The following are the main subroutines and describe the purpose and input/output data of each.
1.	MultiplePrediction.m This subfunction is used to predict multiples, we provide two methods for predicting free surface multiples and internal multiples respectively. The input data of this code is full wavefield data, and the output is one type of predicted multiple. The type of multiple can be selected by the user.
2.	Localsimilarity.m This subfunction is used to get predictability spectrum by calculate the local similarity of the full wavefield data velocity spectrum and predicted multiple velocity spectrum. The input data of this subfunction are two velocity spectrum and the output data are predictability spectrum.
3.	RGBmapping.m  This subfunction is used to the steps of RGB space mapping. The input data are the full wavefield data velocity spectrum and predictability spectrum. The output results are the RGB space and the RGB values of each energy group.
4.	Pick.m    This subfunction is used to use RGB space mapping results to obtain the final velocity model. The input data is the RGB space and the RGB value of each energy group. The output is the final velocity model.
5.	Library function package.  This function package includes some essential library functions in the main scripts and subfunction, such as the smooth function and the clip function.
The Matlab package contains synthetic seismic data used for testing. In addition, the data packet also contains the reference velocity model of the data for comparison with analysis results. 
Program language: Matlab
Software required: Matlab 2014b and above
Program size: 1.65 MB 

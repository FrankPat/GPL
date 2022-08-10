# GPL
General Power Law

![image](https://user-images.githubusercontent.com/62480664/183958477-8181f74c-fb8b-465c-bf91-f1f13de2ef8d.png)

GPL is a computer program written in C++ that calculates the shape of a valley cross-section based on a list of x,y values describing the valley cross section. Since its introduction by Svensson (1959), the power law curve $y = a x^b$ (with $x$ and $y$ horizontal and vertical direction respectively) has been widely used in morphological analysis of glacial trough cross profiles. The numerical constants a and b are obtained by a linear regression analysis of the logarithmic form of the power law curve:

$\log y = \log a + b \log x$

The value $b$ then gives a measure for the form of the cross section. However, over the years this form of the power law has endured a lot of criticism. This criticism is well founded, since this particular form of the power law is not suitable for curve fitting in morphological analyses. Therefore, a general power law is proposed, of the form

$y - y_0 = a | x - x_0 |^b$

with $x_0$, $y_0$ the coordinates of the origin of the cross profile. A unique and unbiased solution for this equation is obtained with a general least squares method, thereby minimizing the error between the cross profile data and the curve, and not between the logarithmic transform of the data and its regression line. This provides a robust way to characterize trough cross-section forms.

Synopsis of GPL

GPL is a program that estimates the shape of a valley cross-section according to the general power law: $y - y_0 = a | x - x_0 |^b$
It was written by Frank Pattyn and Wim Van Huele, and described in "POWER LAW OR POWER FLAW?", by Pattyn and Van Huele. 1998. Earth Surface Processes and Landforms, 23: 761-767. GPL inputs an ascii list of $x-y$ coordinates, one point couple per line, and outputs the parameters $x_0$, $y_0$, $a$ and $b$ corresponding to the best fit curve.

Compilation

The source code gpl.c can be compiled with following command:

cc -o gpl gpl.c -lm

This creates the executable version gpl on a UNIX machine.

Input data file
The program input consist of an ASCII data file with x,y coordinates of the cross-sectional profile. On the first line of this file, the number of data point couples (n) is given. The n consecutive lines are x and y coordinates of each point separated by blanks or tabs. As an example we list here the contents of the file exam1.dat:

16
0 1160
100 954
300 945
500 814
700 648
900 620
1100 573
1300 471
1500 420
1700 485
1900 500
2100 650
2300 754
2500 919
2700 1008
2950 1182

![image](https://user-images.githubusercontent.com/62480664/183962122-e7b5c2f2-ca62-4203-a958-ed864f443f1e.png)

In this example, all measurements of the cross-sectional profile are of equal weight. This means that the covariance matrix in the least-squares is the unity matrix. However, it is also possible to introduce an error measure (or weight factor) in the curve fitting process. This is illustrated in the file exam3.dat:

19
0 1400 1
200 1124 1
400 887 1
600 791 1
800 765 1
1000 652 1
1200 614 1
1400 579 1
1600 562 1
2000 679 100
2200 648 100
2400 485 0.01
3600 526 1
3800 559 1
4200 652 1
4400 772 1
4600 878 1
4800 998 1
5300 1420 1

The cross sectional profile looks like this:

![image](https://user-images.githubusercontent.com/62480664/183962178-1cf4eca3-8c9c-4c84-8d01-0e1ed13df461.png)

The file exam2.dat represents the similar cross-sectional profile as exam3.dat, but without the error measures. In the exam3.dat example, all errors are set to 1, except the central point with a lower error (and thus higher accuracy and weight) of 0.01, and the two points before (the small bump in the profile) with an error of 100 (low weight). The resulting best fit curve will be more or less forced to go through the central point, and the points representing the bump will have less weight in the calculation. If measurements of the cross-sectional profile are made in the field, it is also possible to use those accuracy measurements.

Running the program
The program runs in two ways on a UNIX machine, i.e. interactively so that the program prompts for input and output file names and options, or command line driven. The Macintosh version works only interactively.

COMMAND LINE

gpl inputfile outputfile [-Aaccuracy] [-Ox0/y0] [-C]

The options between square brackets are optional

examples:

gpl exam1.dat exam1.rep -A0.001

This runs the program with exam1.dat as input file. The output file is exam1.rep. The parameter b in the power law equation is estimated with an accuracy of 0.001. Default accuracy is 0.01.

gpl exam1.dat exam1.rep -A0.001 -O1500.0/400.0

This runs the program with exam1.dat as input file. The output file is exam1.rep. The parameter b in the power law equation is estimated with an accuracy of 0.001. The fitted curve is FORCED to go through the inflection point (or origin of the curve) at x=1500.0, y=400.0.

gpl exam3.dat exam3.rep -C

This runs the program with exam3.dat as input file. The output file is exam3.rep. The parameter b in the power law equation is estimated with a default accuracy of 0.01. The input file exam3.dat has a third column with error estimates for each measurement along the cross-sectional profile.

INTERACTIVE

$ gpl

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*             G P L :  G e n e r a l   P o w e r   L a w
*
* Copyright (c) 1998, Frank Pattyn and Wim Van Huele
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

-> Input file name  : exam1.dat
-> Output file name  : exam1.rep
-> Accuracy level power estimate  : 0.001
-> Fix origin x0 / y0 ? (1 = yes; 0 = no)  : 0
-> Use covariance matrix ? (1 = yes; 0 = no)  : 0

In both cases (interative or command-line driven) the program outputs on screen the iterative result of the fitting process, that works its way progressively to the b-value that corresponds with the lowest RMS error between the fitted curve and the observed values. Higher accuracy level estimates will increase the number of iterations needed to obtain this best fit.

Iteration   Power    estimate RMSE
-------------------------------------
  1        2.00000     47.67525
  2        2.10000     50.43143
  3        2.00000     47.67524
  4        1.90000     45.13112
  5        1.80000     42.88868
  6        1.70000     41.06033
  7        1.60000     39.77859

...

 39        1.47344     39.16108
 40        1.47500     39.16086
 41        1.47656     39.16087
 42        1.47500     39.16086
 
 This figure shows the iterative process that seeks the b-value corresponding to the lowest RMS error

![image](https://user-images.githubusercontent.com/62480664/183962300-9ec28403-5e00-4dee-be12-6eaf21334ab0.png)

If we plot the RMS error as a function of the b-value, we can clearly observe that a unique optimal b-value is found corresponding to the lowest RMS error

![image](https://user-images.githubusercontent.com/62480664/183962363-ea1d1452-502c-451a-81c2-8c8417347a78.png)

The figure below displays a detail of the figure above.

![image](https://user-images.githubusercontent.com/62480664/183962418-38425b8a-1cc2-45d9-a7f2-838c688b3ccf.png)

Program output
The program creates one output file, which is listed hereunder for the experiment described above. The contents of the file exam1.rep are:

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*             G P L :  G e n e r a l   P o w e r   L a w
*
* Copyright (c) 1998, Frank Pattyn and Wim Van Huele
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

Power Estimate : y - y0 = a * |x - x0|^b
----------------------------------------
correlation coefficient : 0.973359
RMSE                    : 39.160858
coefficient a           : 1.512631e-02
coefficient b           : 1.475000
coefficient x0          : 1.439432e+03
coefficient y0          : 4.337806e+02

    orig-x      orig-y       new-y      error-y   log|x-x0| log(y-y0) log(newy-y0)

   0.000000  1160.000000  1122.534424  -37.465576  7.272004  6.587852  6.534884
 100.000000   954.000000  1053.136108   99.136108  7.200001  6.254251  6.428679
 300.000000   945.000000   921.699341  -23.300659  7.038285  6.236799  6.190149
 500.000000   814.000000   800.816467  -13.183533  6.845276  5.940748  5.905460
 700.000000   648.000000   691.624634   43.624634  6.605882  5.367001  5.552355
 900.000000   620.000000   595.714966  -24.285034  6.290517  5.226925  5.087191
1100.000000   573.000000   515.550415  -57.449585  5.827274  4.936051  4.403908
1300.000000   471.000000   455.793030  -15.206970  4.937578  3.616830  3.091607
1500.000000   420.000000   440.215485   20.215485  4.103765  2.623262  1.861732
1700.000000   485.000000   489.143097    4.143097  5.562863  3.936118  4.013902
1900.000000   500.000000   562.040344   62.040344  6.132460  4.192973  4.854057
2100.000000   650.000000   652.109070    2.109070  6.493100  5.376294  5.386001
2300.000000   754.000000   756.288330    2.288330  6.757592  5.769006  5.776127
2500.000000   919.000000   872.716675  -46.283325  6.966560  6.184601  6.084354
2700.000000  1008.000000  1000.107544   -7.892456  7.139318  6.353011  6.339171
2950.000000  1182.000000  1173.324341   -8.675659  7.320241  6.617696  6.606033

This output file is a tabulated file, which can easily be read by a spreadsheet program, such as Excel. In the output file the coefficients a, b, x0, y0 are given together with the correlation coefficient and the RMS error corresponding to the best fit. The tabulated list gives the original x and y coordinates, the new y-coordinates (corresponding to the best fit, the error between the fit and the original values, and their logarithmic transforms. The next graph shows the output of orig-y and new-y versus orig-x of the file exam1.rep:

![image](https://user-images.githubusercontent.com/62480664/183962500-565c6658-c911-4922-bf62-01bac65188ad.png)

The following figure displays the results of both the exam2 (red line) and the exam3 (black line) experiment. Using the above described weight factors shows that the bump in the profile is neglected and the profile is fitted through the central point (black line).

![image](https://user-images.githubusercontent.com/62480664/183962564-1f63422f-aac6-4b76-81cb-6148d347c926.png)

Besides the sample packages exam1.dat, exam2.dat and exam3.dat, two other files are available, one with a perfect V-shaped valley (linear.dat) and one containing data of a perfect parabolic valley. Running the gpl program with these data files should give you b-values of 1.0 and 2.0 respectively.

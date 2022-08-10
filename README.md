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

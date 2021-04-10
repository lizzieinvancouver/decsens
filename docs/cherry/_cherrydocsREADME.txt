Started 29 March 2021


<><><><><><><><><><><><>
<> From 25 March 2021 <>
<> slope_interpretation <>
<><><><><><><><><><><><>

From Jonathan:

It was nice talking to you two on Monday. I attached some thoughts on how to interpret the regression slope (i.e. temperature sensitivity)—or more accurately I attached some thoughts that may help you decide how to interpret the regression slope. One section makes no assumption about the distribution of X, other than X > 0. The other section assumes normality. In that case, if the variance of X is much smaller than its mean, the slope of the line approximating Y = s/X is the derivative evaluated at the mean, - s/E[X]^2, which I guess isn't very surprising. But there's a simple formula for the error.


<><><><><><><><><><><><>
<> 31 March 2021
<> slope_interpretation_cjc.pdf <>
<><><><><><><><><><><><>

From Jonathan:

To answer Cat's questions:

1. Average percent change calculates the percent change from x to mu for each x, and then averages over x (i.e. E[ (x - mu) / x ] ). Average absolute percent change calculates the percent change from x to mu for each x, takes the absolute value, and then averages over x (i.e. E[ |( x - mu) / x | ]). The only difference between the two measures is the absolute values—the first can be negative, the second can't. You can think of the first as a measure of bias, except while E[x - mu] = 0,  E[(x-mu)/x] < 0

2. Yes. Sorry for the typo, Y_1 should be X_2.

3. O(sigma^4/mu^4) represents a constant that is small if sigma^4/mu^4 is small. I think you can ignore the O(sigma^4/mu^4) term because, if mean temperatures are ten times the variance, sigma/mu = 1/10 and sigma^4/mu^4 = 1/10000.
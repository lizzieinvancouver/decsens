Started 9 April 2021


<><><><><><><>
5 (or around then) April 2021

Attached:
leafout.R
leafout.stan

Great! Here's my first attempt at modeling a moving temperature window:

Let n denote the day of leafout. Suppose the mean temperature is calculated over day 0 to day b, but the plant doesn't start accumulating temperature until day a, where 0 < a < n < b. That is, the plant accumulates temperature from day a until leafout (day n).

If "a" is unknown, we can model it as a random variable. For example, a ~ Geometric(p). In that case, the assumption is that every day the tree flips a coin with weight p. If the coin is heads, the tree starts accumulating temperature. We could make p proportional to alpha_0, the expected temperature on day 0; that would allow warmer temperatures to delay the start of accumulation. 

See simulation and Stan code attached. Note that Stan can only handle continuous parameters so we have to "average out" a to do inference on p, which makes Stan much slower. (A chain takes 40 minutes on my computer.) If we want to go in this direction, I can look into using reduce_sum to speed up the model. 


<><><><><><><>
8 April 2021

Attached:
leafout_sin.R
leafout_sin_a.stan
leafout_sin_p.stan


I took the liberty of rewriting the model for sinusoidal temperatures instead of linear temperatures, see attached. leafout_sin_a.stan assumes the date plants start accumulating temperatures is known and takes a few minutes to fit. leafout_sin_p.stan assumes plants start accumulating temperatures randomly according to a geometric distribution with probability p and takes a few hours to fit.

Sinusoidal temperatures add an additional parameter (assuming the frequency is known to be annual). Instead of intercept and slope, we have annual average temperature, amplitude, and phase. According to the article, it sounds like we should make p proportional to the amplitude—the higher the amplitude, the cooler the winters, and the faster the tree should be accumulating temperature.



<><><><><><><><><><><>
10 June 2021 (by Cat)


Here is information on the three data sources for cherry blossom data. 

(1) PEP725 data from Germany and Norway 
	
	Data is from: http://www.pep725.eu/data_download/data_selection.php
	Species: Prunus avium 

(2) NOAA data (original publication: https://www.researchgate.net/profile/This-Rutishauser/publication/238552037_A_280Year_Long_Series_of_Phenological_Observations_of_Cherry_Tree_Blossoming_Dates_for_Switzerland/links/5f7b20f3299bf1b53e0e51ec/A-280Year-Long-Series-of-Phenological-Observations-of-Cherry-Tree-Blossoming-Dates-for-Switzerland.pdf) 
	
	Data is from: https://www.ncei.noaa.gov/pub/data/paleo/historical/europe/switzerland/swiss-cherry-phenology2003.txt
	Species: Prunus avium

(3) Meteo Swiss from Liestal site (47.481389, 7.730556; 350m)
	
	Data is from: https://gate.meteoswiss.ch/idaweb/system/stationList.do
		(ask Cat for login information if needed)
	Species: Prunus avium




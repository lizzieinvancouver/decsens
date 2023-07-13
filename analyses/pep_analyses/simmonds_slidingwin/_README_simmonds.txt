Started 13 July 2023
By Lizzie
So written just a few (3-4) years after last working on this.

<><><><><><><><><><><><><><><>
<> Overview <>
<><><><><><><><><><><><><><><>

Code in this folder explores the climwin package from R including implementations of it from Simmonds et al. 2019 "Cue identification in phenology: A case study of the predictive performance of current statistical tools." It was written by Lizzie and Cat Chamberlain to try to figure how well it actually worked.

We dropped this work in late 2020 or early 2021 I believe as no one had time to lead it.

<><><><><><><><><><><><><><><>
<> Quick review of files <>
<><><><><><><><><><><><><><><>
2023May_movingwindow_auerbach.pdf -- email from Auerbach about his student working on this; they took a robust statistics approach as opposed to a more biological approach. [May only be on Lizzie's computer.]


* All of the below are best guesses by Lizzie in 2023; we should update if we figure out that they are wrong.* 

betpen_climate_slidingwin.R and 
fagsyl_climate_slidingwin.R 
	code to make mean climate data for PEP data (by Cat)

bp_sw_simmonds.R and 
fs_sw_simmonds.R 
	uses code from Emily Simmond's sliding window approach: https://github.com/emilygsimmonds/Cue_Identification (by Cat)

fsbetpen.sh and
fsfagsyl.sh
	Code to run moving window code on server (bash script)

Params_SW.R 
 	Function to extract parameters from sliding window analysis (sourced in bp_sw_simmonds. R I believe)

Run_SW.R 
	Function to run sliding window analyses (sourced in bp_sw_simmonds. R I believe

sandboxy_sw.R (Lizzie's code) 
	This file runs  the sliding window from climwin for one site for Betpen (PEP) and makes some notes on what happens in Run_SW f(x)

sims_runsw.R (Lizzie's code)
	Reads in simulated budburst and climate data that I made, then runs the sliding windows on it, then plots the output. (See Figures folder) 

sims_sw.R (Lizzie's code) 
	Here I made the simulated climate data (and maybe budburst). 


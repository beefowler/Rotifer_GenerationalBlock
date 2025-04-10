# Rotifer_GenerationalBlock
Code associated with manuscript, "Timing the Initiation of Sex" 



Authors: Bethany L.F. Stevens, Silke van Daalen, Tirzah J. Blomquist, Kristin E. Gribble, Michael G. Neubert
contact bstevens@ucsb.edu for questions. 

Summary: We develop a model of cyclically parthenogenic rotifer populations with a novel formulation of a ``mictic block'' that prevents sexual reproduction by females that are not sufficiently distant, genealogically, from a stem ancestor that was produced sexually. We consider the success of strains with different mixis phenotypes in different environments, including consecutive seasons of stochastic length as in Serra et al. 2005. We also consider invasion success of different strategies and explore the potential for a stable polymorphism. 

All scripts in this repository were written by BLFS, SVD, and MGN in MatLAB. 

delay_example.m 
    This script generates Figure 2 in the manuscript, demonstrating the difference in egg production between strains with temporal and generational mictic blocks. The first two panels arise from single season simulations of two scenarios described in the text, for which we plot egg production over time. The third panel allows for consecutive seasons and plots the ratio of the number of eggs in each scenario for each consecutive season. 

gen_rotifer_onestrain.m 
    This is the function that contains the main system of differential equations described in the manuscript in equations 8-18 for a single monomorphic population of rotifers. It calculates the derivative with respect to time for each of the groups of interest in the population at time, t. The number of these groups and hence the number of state variables that need to be included depends on the strain's phenotype value of G_i. The state variable vector, x, is organized as (Resting Eggs, Mictic Adults, Juveniles of Generation 1, Adults of Generation 1, Juveniles of Generation 2, Adult sof Generation 2, ..., Juveniles of Genreation G_i, Adults of Genreation G_i). Input parameters are described at the top of the script. Juveniles and Adults of Generation 0 have determinstic beavior and can be calculated just from t, phi, tau, q in functions at bottom of the script. These functions are included in the calculation of the total population size at time, t, which affects desnity dependent birth rate and mixis ratios. 
    This function is called by many subsequent scripts. 
   
gen_rotifer_twostrain.m 
    This is similar to the function above, but modified to allow for two competing strains of rotifers. Input phenotype parameters, m_i, T_i, and G_i can now be size 2 to have different phenotypes competing. State variable, x, now includes that which would be used above for the first strain followed by the second. Both strains contribute to the total population size to control density dependent birth rate and mixis ratios.
    This function is called by many subsequent scripts 

Figure4_WithinSeason.m 
    This script simulates the within-season dynamics of a monomorphic population by calling gen_rotifer_onestrain within dde23 delay differential equation solver. The top of the script defines the parameter values, then simulates and plots for two scenarios with different mixis values. 




Data_for_Figur5.mat
    This script has results of stochastic simulations included in Figure 5. 


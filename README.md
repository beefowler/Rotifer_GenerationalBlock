# Rotifer_GenerationalBlock
Code associated with manuscript, "Timing the Initiation of Sex" 

Authors: [ Redacted for Review ]

Summary: We develop a model of cyclically parthenogenic rotifer populations with a novel formulation of a ``mictic block'' that prevents sexual reproduction by females that are not sufficiently distant, genealogically, from a stem ancestor that was produced sexually. We consider the success of strains with different mixis phenotypes in different environments, including consecutive seasons of stochastic length as in Serra et al. 2005. We also consider invasion success of different strategies and explore the potential for a stable polymorphism. 

All scripts in this repository were written in MatLAB. 

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

Figure5_part1.m   
    This script simulates 40 independent experiments of 40 consecutive seasons of stachastic season length for an assortment of mixis phenotypes as plotted in Figure 5 of the manuscript. It does this by simulating the model with dde23 and adding the number of eggs produced in each season to create a record of cumulative_egg_production for each phenotye and each experimemt. 

Data_for_Figure5.mat  
    This script has results of stochastic simulations included in Figure 5.  
    mixis = values of m_i used in experiment (x axis)  
    G_vals = values of G_i used in experiments  
    threshold = values of T_i used in experiment   
    mixis_and_G = all combined pairs of phenotypes for the first set of experiments where T_i  = 0  
    mixis_and_thresh = all combined pairs of phenotype values for the second set of experiments where G_i - 0   
    num_seasons = number of consecutive seasons per experiment  
    num_experiments = number of independent experiments per phenotype  
    cumulative_egg_production = cumulative eggs produced in each experiment with T_i = 0 for each phenotype organized according to "mixis_and_G"   
    cumulative_egg_production_2 = cumulative eggs produced in each experiment with G_i = 0 for each phenotype organized according to "mixis_and_thresh"   
    
Figure5_part2.m     
    This script takes the outputs of Figure5_part1.m (Data_for_Figure5.mat) and generates the desired figure. 
    
Figure6_fixedseasons.m   
    This script explores how monomorphic populations of different phenotypes behave in short, medium, and long seasons of fixed length. Again we use dde23 applied to gen_rotifer_onestrain to run these simulations. Lines 1-186 carry this out for G_i = 0 for each of three season lengths (Could have used a for loop, but we didn't), then the same steps are carried out for G_i = 8 to create the second column of the figure.   

Plot_PairwiseInvasionMatrices.m   
    Plots the results of the invasion experiments. Pulls from folders with outputs for individual invasion experiments, looks for the mean invasion proportion at the end of the experiments, and plots success if greater than or equal to 0.05 and failure if less than 0.05. One section of the code for eachof the Pairwise Invasion Plots in Figure S5, two of which make up Figure 7 in the main text.   

Figure8_NoBlockResident.m   
    This script generates figure 8 from the results of the No Block Resident Invasion Experiments. 

Figure9_TwoStrainTrajectories.m   
    This script simulates two strains simultaneously for a single season. We use it to explore the dynamics of two strains that can mutualy invade eachother in our invasion experiments and to see their relative success over the course of a single season. We simulate this scenario using dde23 and gen_rotifer_twostrains.m This script produces Figure 9. 

    

    

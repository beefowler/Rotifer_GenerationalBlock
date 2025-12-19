# Rotifer_GenerationalBlock
Code associated with manuscript, "Timing the Initiation of Sex: delay mechanisms alter fitness outcomes in a rotifer population model" published in Journal of Theoretical Biology 

Authors: Bethany Stevens, Silke van Daalen, Tirzah Blomquist, Kristin Gribble, and Michael Neubert 

https://doi.org/10.1016/j.jtbi.2025.112333

Summary: We develop a model of cyclically parthenogenic rotifer populations with a novel formulation of a ``mictic block'' that prevents sexual reproduction by females that are not sufficiently distant, generationally, from a stem ancestor that was produced sexually. We consider the success of strains with different mixis phenotypes in different environments, including consecutive seasons of stochastic length as in Serra et al. 2005. We also consider invasion success of different strategies and explore the potential for coexistence between distinct phenotypes. 

All scripts in this repository were written in MatLAB. 

gen_rotifer_onestrain.m   
    This is the function that contains the main system of differential equations described in the manuscript in equations 1-10 for a single monomorphic population of rotifers. It calculates the derivative with respect to time for each of the groups of interest in the population at time, t. The number of these groups and hence the number of state variables that need to be included depends on the strain's phenotype value of G_i. The state variable vector, x, is organized as: (Resting Eggs, Mictic Adults, Juveniles of Generation 1, Adults of Generation 1, Juveniles of Generation 2, Adult sof Generation 2, ..., Juveniles of Generation G_i, Adults of Generation G_i). Input parameters are described at the top of the script. Juveniles and Adults of Generation 0 have determinstic beavior and can be calculated just from t, phi, tau, q in functions at bottom of the script. These functions are included in the calculation of the total population size at time, t, which affects density dependent birth rate and mixis ratios. 
    This function is called by many subsequent scripts. 
   
gen_rotifer_twostrain.m   
    This is similar to the function above, but modified to allow for two competing strains of rotifers. Input phenotype parameters, m_i, T_i, and G_i can now be size 2 to have different phenotypes competing. State variable, x, now includes that which would be used above for the first strain followed by the second. Both strains contribute to the total population size to control density dependent birth rate and mixis ratios.
    This function is called by many subsequent scripts 

Figure3_WithinSeason.m   
    This script simulates the within-season dynamics of a monomorphic population by calling gen_rotifer_onestrain within dde23 delay differential equation solver. The top of the script defines the parameter values, then simulates and plots for two scenarios with different mixis values. 

Figure4_part1_Fitness.m   
    This script simulates 40 independent experiments of 40 consecutive seasons of stachastic season length for an assortment of mixis phenotypes as plotted in Figure 4 of the manuscript. It does this by simulating the model with dde23 and calculating the average fitness across seasons for each phenotye and each experiment. These results are saved to Data_for_Figure4_Fitness.mat.  

Data_for_Figure4_Fitness.mat  
    This file has results of stochastic simulations included in Figure 5.  
    mixis = values of m_i used in experiment (x axis)  
    G_vals = values of G_i used in experiments  
    threshold = values of T_i used in experiment   
    mixis_and_G = all combined pairs of phenotypes for the first set of experiments where T_i  = 0  
    mixis_and_thresh = all combined pairs of phenotype values for the second set of experiments where G_i - 0   
    num_seasons = number of consecutive seasons per experiment  
    num_experiments = number of independent experiments per phenotype  
    cumulative_egg_production = cumulative eggs produced in each experiment with T_i = 0 for each phenotype organized according to "mixis_and_G"   
    cumulative_egg_production_2 = cumulative eggs produced in each experiment with G_i = 0 for each phenotype organized according to "mixis_and_thresh"   
    fintess_all = fitness for each experiment with Ti = 0, organized as above
    fitness_all_2 = fitness for each experiment with G_i = 0, organized as above
    
Figure4_part2_Fitness.m     
    This script takes the outputs of Figure4_part1.m (Data_for_Figure4_Fitness.mat) and generates the desired figure. 
    
Figure5_fixedseasons.m   
    This script explores how monomorphic populations of different phenotypes behave in short, medium, and long seasons of fixed length. Again we use dde23 applied to gen_rotifer_onestrain to run these simulations. Lines 1-186 carry this out for G_i = 0 for each of three season lengths, then the same steps are carried out for G_i = 8 to create the second column of the figure.   

Invasion_Baseline_part1
    This script conducts the invasion experiments needed to create the baseline pairwise invasibility plots for variable m when G_i = 0 and T_i = 0.   
    For each pair of phenotypes, the script simulates 20 independent experiments with 40 consecutive seasons each using gen_rotifer_twostrains and dde23. The invader is introduced at the start of each experiement at frequency 0.05. Seasons have randomly chosen length between 10 and 51 days. If either population goes extinct before end of experiment, experiment is ended. Final invader proprotion is saved to invasion_prop variable within a matlab file named for the combination of phenotypes in that experiment. 

PairwiseInvasion_M_Gis3.m 
    Almost identical to Invasion_Baseline_part1.m, this script conducts the same set of experiments for G_i = 3 and saves results to PairwiseInvasion_M_Gis3. Mostly only made it a new script so we could run simultaneously on the same computer. 
    
InvasionofDelayTypes.m 
    Following the same structure as Invasion_Baseline_part1.m, this script conducts invasion experiemnts for strains with mi = 0.11 and threshold =0. Strains with generational blocks (G_i = 0:8) are  to a resident without a generational block and outputs are saved to Invasion_NoDelay_byG. 

PairwiseInvasion_G.m 
    Following the same structure as InvasionofDelayTypes.m, this script conducts pariwise invasion experiments for strains with variable G values. In this case all phenotypes can be either resident or invader in the experiments. Results are saved to PairwiseInvasion_G directory. We ran this script both for mixis = 0.17 and mixis = 0.11 and results are labelled accordingly. 

Plot_PairwiseInvasionMatrices.m   
    Plots the results of the invasion experiments. Pulls from folders with outputs for individual invasion experiments, looks for the mean invasion proportion at the end of the experiments, and plots success if greater than or equal to 0.05 and failure if less than 0.05. One section of the code for each of the Pairwise Invasion Plots in Figure S11, two of which make up Figure 6 in the main text.   

Figure7_NoBlockResident.m   
    This script generates boxchart portion of figure 7 from the results of the No Block Resident Invasion Experiments. 

CompareFitness_ManySeasons.m 
    This script calculates relative ftness of invader strategies compared to No Block Resident to add on top of Figure 7. The dynamics of a resident are simulated for 1000 seasons and fitness is calculuated. Then we simulate the introduction of an invader at low frequences at the start of each season. Invader dynamics do not impact resident dynamics by calling gen_rotifer_invader_negligible.m. Relative fitness is calculated as the ratio of the mean fitness of the invader to the resident across all seasons. 



The following scripts contribute to the paper Supplement 

delay_example.m   
    This script generates Figure S1 in the manuscript, demonstrating the difference in egg production between strains with temporal and generational mictic blocks. The first two panels arise from single season simulations of two scenarios described in the text, for which we plot egg production over time. The third panel allows for consecutive seasons and plots the ratio of the number of eggs in each scenario for each consecutive season.

FigureS_birthdaycloseup.m
    This script generates Figures S4

CheckForBestT_Fitness.m 
    This script generates Figure S5 and corresponding data is saved to Data_CheckForBestT.mat

FigureS6_part1.m 
    Generates data for Figure S6
FigureS6_part2.m
    Generates Figure S6 from Data_for_FigureS6.mat

Sensitivity_Analysis.m and Sensitivity_Analysis_Bmax.m 
    Generate Figures S7 and S8

Invasion_G_and_m.m
    Runs experiments for Figure S12. Outputs are summarized in Invasion_with_2_params_outputmatrix.mat

FigureS_TwoStrainTrajectories.m   
    This script simulates two strains simultaneously for a single season. We use it to explore the dynamics of two strains that can mutualy invade eachother in our invasion experiments and to see their relative success over the course of a single season. We simulate this scenario using dde23 and gen_rotifer_twostrains.m This script produces Figure S13B. 

NoLateClones_SI.m 
    Generates Figure S14
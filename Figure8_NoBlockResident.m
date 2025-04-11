% This script produces a boxchard figure for invasion results of No Delay
% Resident by invader pehnotypes with Generational Block 

outputlist = dir('Invasion_NoDelay_byG/*.mat'); 

G_vals = [0:8]; 
N = length(G_vals);
g_invader = nan(1,N);
all_invasion_prop = nan(20,N); 


for i = 1:length(G_vals) 

     load(['Invasion_NoDelay_byG/outputs_nodelay_resident_by_G' num2str(G_vals(i)) '.mat'])
     g_invader(i) = G_i(2); 
     all_invasion_prop(:,i) = invasion_prop; 

end

figure
b = boxchart(all_invasion_prop, 'BoxFaceColor', 'k', 'MarkerColor', 'k');
b.MarkerStyle = 'o';
xticklabels(G_vals)
 
ylabel('Final invader proportion')
xlabel('Invader generational block, G_2')
title(['Mixis ratio: m_i = ' num2str(pheno_1(1)) ' Mixis threshold: T_i =' num2str(pheno_1(2))])
set(findall(gcf,'-property','FontSize'),'FontSize',14)
fontname('Arial')
box on 
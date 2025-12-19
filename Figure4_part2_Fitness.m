clear all

% Create Figure 4 from results of stochastic experiments run in
% Figure4_part1.m 

load('Data_for_Figure4_Fitness.mat')

% Get colors
bluemono = [8, 47, 97;
            8, 69, 148;
            33, 113, 181;  
            66, 146, 198; 
            107, 174, 214; 
            158, 202, 225; ]./255; 
colorvals =  bluemono(end:-1:1,:);


fitness_nonans = cellfun(@replacenans, fitness_all, 'UniformOutput', false);
fitness_all = fitness_nonans; 


%%
Average_fitness1 = cellfun(@geomean, fitness_all)
Average_fitness1 = mean(Average_fitness1')


% now make a figure
G = length(G_vals); 
m = length(mixis); 

figure(1)
subplot(1,2,2)
%expected_egg_production = mean(cumulative_egg_production')./num_seasons;
for i = 1:length(G_vals) %want to have lines connecting dots, stupidly hard now 
    hold on 
    plot(mixis_and_G(1,i:G:end),Average_fitness1(i:G:end), 'color', colorvals(i,:), 'linewidth', 3)
    scatter(mixis_and_thresh(1,i:G:end), Average_fitness1(i:G:end), 80, colorvals(i,:), 'filled')

end
%scatter(mixis_and_G(1,:), Average_fitness1, 80, mixis_and_G(2,:), 'filled')
colormap(colorvals)
title(['T_1 = 0'])
ylabel({'Geometric mean fitness'})
xlabel('Mixis ratio, m_1')
h = colorbar ;
h.Label.String = 'Generational block, G_1';
h.Ticks = G_vals
box on 
%ylim([0 108])

%%
%second panel

fitness_nonans_2 = cellfun(@replacenans, fitness_all_2, 'UniformOutput', false);
fitness_all_2 = fitness_nonans_2; 
Average_fitness2 = cellfun(@geomean, fitness_all_2)
Average_fitness2 = mean(Average_fitness2')


subplot(1,2,1)
t = length(threshold); 
m = length(mixis); 
colorvals =  bluemono(end:-1:1,:);
%expected_egg_production = mean(cumulative_egg_production_2')./num_seasons;
for i = 1:6 %want to have lines connecting dots, stupidly hard now 
    hold on 
    plot(mixis_and_thresh(1,i:t:end),Average_fitness2(i:t:end), 'color', colorvals(i,:), 'linewidth', 3)
    scatter(mixis_and_thresh(1,i:t:end), Average_fitness2(i:t:end), 80, colorvals(i,:), 'filled')

end




%scatter(mixis_and_thresh(1,:), Average_fitness2, 80, mixis_and_thresh(2,:), 'filled')
colormap(colorvals)
title(['G_1 = 0'])
ylabel({'Geometric mean fitness'})
xlabel('Mixis ratio, m_1')
h = colorbar ;
h.Label.String = {'Mixis threshold density, T_1'; '(Individuals L^{-1})'};
h.Ticks = threshold
%ylim([0 108])

box on 

set(findall(gcf,'-property','FontSize'),'FontSize',15)


function C = replacenans(C)
    C(isnan(C)) = 0;
    C(C<0) = 0; 
end

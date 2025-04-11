clear all

% Create Figure 5 from results of stochastic experiments run in
% Figure5_part1.m 

load('Data_for_Figure5.mat')

% Get colors
bluemono = [8, 69, 148;
            33, 113, 181;  
            66, 146, 198; 
            107, 174, 214; 
            158, 202, 225; 
            198, 219, 239]./255; 
colorvals =  bluemono(end:-1:1,:);



% now make a figure
G = length(G_vals); 
m = length(mixis); 

figure(1)
subplot(1,2,2)
expected_egg_production = mean(cumulative_egg_production')./num_seasons;
for i = 1:length(G_vals) %want to have lines connecting dots, stupidly hard now 
    hold on 
    plot(mixis_and_G(1,i:G:end),expected_egg_production(i:G:end), 'color', colorvals(i,:), 'linewidth', 3)
end
scatter(mixis_and_G(1,:), expected_egg_production, 80, mixis_and_G(2,:), 'filled')
colormap(colorvals)
title(['T_1 = 0'])
ylabel({'Expected egg production'; 'per season (L^{-1})'})
xlabel('Mixis ratio, m_1')
h = colorbar ;
h.Label.String = 'Generational block, G_1';
h.Ticks = G_vals
box on 
ylim([0 108])


%second panel
subplot(1,2,1)
t = length(threshold); 
m = length(mixis); 
colorvals =  bluemono(end:-1:1,:);
expected_egg_production = mean(cumulative_egg_production_2')./num_seasons;
for i = 1:t %want to have lines connecting dots, stupidly hard now 
    hold on 
    plot(mixis_and_thresh(1,i:t:end),expected_egg_production(i:t:end), 'color', colorvals(i,:), 'linewidth', 3)
end
scatter(mixis_and_thresh(1,:), expected_egg_production, 80, mixis_and_thresh(2,:), 'filled')
colormap(colorvals)
title(['G_1 = 0'])
ylabel({'Expected egg production'; 'per season (L^{-1})'})
xlabel('Mixis ratio, m_1')
h = colorbar ;
h.Label.String = 'Mixis threshold, T_1';
h.Ticks = threshold
ylim([0 108])

box on 

set(findall(gcf,'-property','FontSize'),'FontSize',15)
clear all

% Explore relative fitness of invader vs resident in negligible invader
% scenarios 

% Define the parameters
bmax=2.1; q=0.5; K=500; c=0.2; H=31; tau = 1;
% H is the final day of the season that hatching can happening
B_max = 1; %maximum number of resting eggs that survive to next season. 
%Hatching schedule
global sk x_width phi
t0=0; %0 is season begins on first hatch day
sk = t0:1:H; %days that hatching will occur

%phenotypes - Use best phenotype in absense of a block 
m_i = .11;
T_i = 0;

G_vals = [0:8]; %generational delays to try for invader 
G_num = length(G_vals);

num_seasons = 500; 
% initialize outputs 
R_out_resident_all = nan(length(G_vals),num_seasons); 
R_out_invader_all = nan(length(G_vals),num_seasons); 
realized_seasons = nan(1, num_seasons); 


pheno_1 = [m_i, T_i];
pheno_2 = [m_i, T_i];



% Start Simulations with season 1 

R_y = B_max; %number of resting eggs that we start with.

%Assume all eggs are resident
phi_1 = (R_y)./length(sk); % figure out how many resting eggs will emerge each hatch day

% And also consider an imaginary population of 5% invader
phi_2 = (B_max* .05)./length(sk); % figure out how many resting eggs will emerge each hatch day
phi = [phi_1 phi_2];

for season = 1:num_seasons

    x_hists = zeros(x_width-2,1); %past state space, not including generation 1. Resets every season

    Tn = randi([10 51]); %end of season, randomly selected from range
    realized_seasons(season) = Tn;

    %First timestep need history function.
    timestep = 1;
    tspan = 0:timestep:Tn;

    m_i = [pheno_1(1) pheno_2(1)];
    T_i = [pheno_1(2) pheno_2(2)];

    for g = 1:length(G_vals)

        %Resident generational delay, Second value will be changed below
        G_i = [0 NaN];
        G_i(2) = G_vals(g);  %this one cycles through the values for invader generational delay

        %adjust xwidth to generation delay for both phenotypes
        if G_i(1) ~= 0 && G_i(2) ~=0
            x_width =  2+2*(G_i(1)+1) +  2+2*(G_i(2)+1);   % (R, M, J, A, J, A, ...) 0 to G_i
            half_x_width = 2*(G_i(1)+1);
        elseif G_i(1) == 0 && G_i(2) ~=0
            x_width =  6 + 2+2*(G_i(2)+1);
            half_x_width = 4;
        elseif G_i(1) ~=0 && G_i(2) ==0
            x_width =  2+2*(G_i(1)+1) + 6;
            half_x_width = 2*(G_i(1)+1);
        else
            x_width = 12;
            half_x_width = 4;
        end

        %run simulation
        sol=dde23(@(t,x,x_hists) gen_rotifer_invader_negligible(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @(t) history_fun(t, sk, x_width), tspan);

        R_out_invader = sol.y(half_x_width+1,end);
        if R_out_invader <=0
            R_out_invader = realmin 
        end
        R_out_invader_all(g, season) = R_out_invader;

        R_out_resident = sol.y(1, end);
        if R_out_resident <= 0 
            keyboard
        end
        R_out_resident_all(g, season) = R_out_resident;
        
    end
    R_y = min(sol.y(1,end), B_max);
    % resident number of eggs changes, invader is always the same;
    phi = [(R_y)./length(sk) (B_max* .05)./length(sk)];

    disp([season phi])
end


save('CompareFitness_ManySeasons2.mat')


map = [238, 204, 102;
    238, 153, 170;
    102, 153, 204; 
    153, 119. 0;
    153, 68, 85; 
    0, 68, 136
]./255;

hold on 
yyaxis left
R0_invader = R_out_invader_all(:, 2:end)./(B_max*.05);
R0_resident = R_out_resident_all(:, 2:end)./min(R_out_resident_all(:, 1:end-1), B_max);


plot(G_vals+1, geomean(R0_invader')./geomean(R0_resident'), '--', 'linewidth', 2, 'color', map(6,:))
ylabel('Relative Fitness of Invader')
xlabel('Invader generational block, G_2')
title(['Mixis ratio: m_i = ' num2str(pheno_1(1)) ' Mixis threshold: T_i =' num2str(pheno_1(2))])

set(findall(gcf,'-property','FontSize'),'FontSize',14)
fontname('Arial')

ax = gca;
ax.YColor = map(6,:);

xlim([0 8])
ylim([0 1.4])


%% Combine this with Figure 8


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

hold on 
yyaxis right
b = boxchart(all_invasion_prop, 'BoxFaceColor', 'k', 'MarkerColor', 'k');
b.MarkerStyle = 'o';
xticklabels(G_vals)
 
ylabel('Final invader proportion')
xlabel('Invader generational block, G_2')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
fontname('Arial')
box on 

ax = gca;
ax.YColor = 'k';

%% history functions 

function v = history_fun(t, sk, x_width) %not really necessary now that we are not including generation 0 in diffe solver
clear zeros
  if t< sk(1)
      v =  zeros(x_width-4,1); 
  elseif t == sk(1)
      v = zeros(x_width-4,1); 
  else
      v = NaN; 
  end
end


function J_0 = how_many_J(t, phi, sk, tau, q)
  [~,last_birthday_ind] = find(sk<=t, 1, 'last'); 
  time_between_clutches = sk(2); 
  if tau <= time_between_clutches %maturation time is less than time between birthdays
    %all individuals from previous clutch will be gone
    time_since_clutch = t-sk(last_birthday_ind);
    %either they have all matured, or not
    if time_since_clutch >= tau
        J_0 = 0; 
    else
        J_0 = phi.*exp(-q.*time_since_clutch); 
    end

  else  %we will have some overlapping generations
     [~, clutches_to_consider] = find(sk > t-tau & sk <= t);
     J_0 = 0; 
     for c = clutches_to_consider
         time_since_clutch = t - sk(c); 
         J_0 = J_0 + phi.*exp(-q.*time_since_clutch); %add clutch members that haven't died
     end
  end

  if t<0 %careful with history
      J_0 = 0; 
  end

end


function A_0 = how_many_A(t, phi, sk, tau, q)

    [~, clutches_to_consider] = find(sk <= t-tau);
    maturation_days = sk(clutches_to_consider) + tau; 
    A_0 = sum(phi.*exp(-q.*tau).*exp(-q.*(t-maturation_days)));  %# of juvs who matured decaying exponential for amount of time since maturation

end


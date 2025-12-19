clear all


%%   keep G_i = 0 and vary T_i and m_i 


% Define the parameters 
bmax= 2.1; q=0.5; K=500; c=0.2; H=31; tau = 1; 
% H is the final day of the season that hatching can happening 
B_max = 1; %maximum number of resting eggs that survive to next season. In paper it is 1????
%Hatching schedule 
global sk x_width phi
t0=0; %0 is season begins on first hatch day 
sk = t0:1:H; %days that hatching will occur


num_experiments = 40; 
num_seasons = 40; 


%choose generation delay =0 and resize state variable
G_i = 0;
x_width = max(2+2*(G_i+1), 6);   % (R, M, J, A, J, A, ...) 0 to G_i
J_ind = [3:2:x_width];
A_ind = [4:2:x_width];


%phenotypes to try 
mixis = [0:0.05:.5]; %mixis values along x axis 
threshold = [0:.1:1]; % values of T_i for each of the lines when Gi=0 

[M,T] = meshgrid(mixis, threshold);
n = length(mixis)*length(threshold);
MM = reshape(M, [1, n]);
TT = reshape(T, [1, n]);
mixis_and_thresh = [MM; TT];


% create output 
a_num = length(mixis_and_thresh); 
cumulative_egg_production_2 = zeros(a_num, num_experiments); 
fitness_all_2 = cell(a_num, num_experiments); 

for x = 1:num_experiments
    x
    for a=1:a_num

        m_i=mixis_and_thresh(1,a);
        T_i=mixis_and_thresh(2,a);

        R_y = B_max; %number of resting eggs that we start with. Can't tell what this should be?
        phi = R_y ./length(sk); % figure out how many resting eggs will emerge each hatch day
        season_fitness = nan(1, num_seasons); 

        for season = 1:num_seasons
            if season>1 
                this_seasons_Ry = min(sol.y(1,end), B_max);
            else
                this_seasons_Ry = R_y; 
            end
            
            x_hists = zeros(x_width-2,1); %past state space, not including generation 1. Resets every season

            Tn = randi([10 51]); %end of season

            timestep = 1;
            tspan = 0:timestep:Tn;

            %run simulation
            sol=dde23(@(t,x,x_hists) gen_rotifer_onestrain(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @history_fun, tspan);

            phi = min(sol.y(1,end), B_max) ./length(sk); %number of eggs to start next season
            
            season_fitness(season) = sol.y(1,end)./this_seasons_Ry;

            if sol.y(1,end) == 0 %if no eggs survive, no need to continue to next season
                break 
            end

            cumulative_egg_production_2(a, x) = cumulative_egg_production_2(a, x) + sol.y(1,end);

        end
       fitness_all_2{a, x} = season_fitness;

    end
end


 save('Data_CheckForBestT.mat', 'mixis', 'threshold', 'mixis_and_thresh', 'num_seasons', 'num_experiments', 'cumulative_egg_production_2', 'fitness_all_2')


%%


load('Data_CheckForBestT.mat')

% Get colors
bluemono = [8, 69, 148;
            33, 113, 181;  
            66, 146, 198; 
            107, 174, 214; 
            158, 202, 225; 
            198, 219, 239]./255; 
colorvals =  bluemono(end:-1:1,:);



fitness_nonans_2 = cellfun(@replacenans, fitness_all_2, 'UniformOutput', false);
fitness_all_2 = fitness_nonans_2; 
%%
Average_fitness2 = cellfun(@geomean, fitness_all_2)
Average_fitness2 = mean(Average_fitness2')

%second panel
subplot(1,2,1)
t = length(threshold); 
m = length(mixis); 
colorvals =  bluemono(end:-1:1,:);
expected_egg_production = mean(cumulative_egg_production_2')./num_seasons;
%for i = 1:t %want to have lines connecting dots, stupidly hard now 
 %   hold on 
  %  plot(mixis_and_thresh(1,i:t:end),Average_fitness2(i:t:end), 'color', colorvals(i,:), 'linewidth', 3)
%end
scatter(mixis_and_thresh(1,:), Average_fitness2, 80, mixis_and_thresh(2,:), 'filled')
colormap(colorvals)
title(['G_1 = 0'])
ylabel({'Geometric mean fitness'})
xlabel('Mixis ratio, m_1')
h = colorbar ;
h.Label.String = 'Mixis threshold, T_1';
h.Ticks = threshold
%ylim([0 108])

box on 

%% history functions 


function v = history_fun(t) %not really necessary now that we are not including generation 0 in diffe solver
global sk x_width phi
  if t< sk(1)
      v =  zeros(x_width-2, 1); 
  elseif t == sk(1)
      v = zeros(x_width-2, 1); 
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

function C = replacenans(C)
    C(isnan(C)) = 0;
    C(C<0) = 0; 
end

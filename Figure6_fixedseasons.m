clear all
figure

%Let's see how different strains perform in fixed seasons. 


% Define the parameters 
bmax=2.1; q=0.5; K=500; c=0.2; H=31; tau = 1; 
% H is the final day of the season that hatching can happening 
B_max = 1; %maximum number of resting eggs that survive to next season
%Hatching schedule 
global sk x_width phi
t0=0; %0 is season begins on first hatch day 
sk = t0:1:H; %days that hatching will occur
 
%choose generation delay and presize things
G_i = 0;
x_width = max(2+2*(G_i+1), 6);   % (R, M, J, A, J, A, ...) 0 to G_i
J_ind = [3:2:x_width];
A_ind = [4:2:x_width];


%These are the three season lengths to explore 
Tn1 = 15;
Tn2 = 30;
Tn3 = 45; 

%phenotypes to try 
mixis = [0:0.05:1];
threshold = 0:20:100 
[M,T] = meshgrid(mixis, threshold);
n = length(mixis)*length(threshold);
MM = reshape(M, [1, n]);
TT = reshape(T, [1, n]);
mixis_and_thresh = [MM; TT];


%save each season results
season = 1; %no consecutive seasons here
a_num = length(mixis_and_thresh); 

realized_season_length = nan(a_num, 1, 1); 
realized_egg_production = nan(a_num, 1, 1); 


x_width = max(2+2*(G_i+1), 6);   % (R, M, J, A, J, A, ...) 0 to G_i
J_ind = [3:2:x_width];
A_ind = [4:2:x_width]; 


%% Panel A - G_i = 0, Short Season 

% we could have put it all together into a for loop, but we didn't, so that
% we could adjust plot scales and titles as we went etc.  

for a=1:a_num %for each phenotype

    m_i=mixis_and_thresh(1,a); 
    T_i=mixis_and_thresh(2,a);

    R_y = B_max; %number of resting eggs that we start with. 
    phi = R_y ./length(sk); % figure out how many resting eggs will emerge each hatch day

        x_hists = zeros(x_width-2,1); %past state space, not including generation 1. Resets every season

        Tn = Tn1; %end of season, only expected value for now

        %First timestep need history function.
        timestep = 1;
        tspan = 0:timestep:Tn;

        %run simulation
        sol=dde23(@(t,x,x_hists) gen_rotifer_onestrain(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @history_fun, tspan);
        
        realized_egg_production(a, season) = sol.y(1,end);
        phi = min(sol.y(1,end), B_max) ./length(sk); %number of eggs to start next season 

end

% now make a figure
subplot(3,2,1)
bluemono = [8, 69, 148;
            33, 113, 181;  
            66, 146, 198; 
            107, 174, 214; 
            158, 202, 225; 
            198, 219, 239]./255; 
t = length(threshold); 
colorvals = bluemono(end:-1:1,:);
for i = 1:t %want to have lines connecting dots  
    hold on 
    plot(mixis_and_thresh(1,i:t:end),realized_egg_production(i:t:end)', 'color', colorvals(i,:))
end
scatter(mixis_and_thresh(1,:), realized_egg_production, 40, mixis_and_thresh(2,:), 'filled')
colormap(colorvals)
ylabel('Cumulative egg production')
xlabel('Mixis ratio')
title('G_i = 0')
ylim([0 1])

%% Panel C - Gi = 0, Medium Season 

realized_egg_production = nan(a_num, 1, 1); 

for a=1:a_num

    m_i=mixis_and_thresh(1,a); %mixis_and_thresh(1,a);
    T_i=mixis_and_thresh(2,a);

    R_y = B_max; %number of resting eggs that we start with. 
    phi = R_y ./length(sk); % figure out how many resting eggs will emerge each hatch day

        x_hists = zeros(x_width-2,1); %past state space, not including generation 1. Resets every season

        Tn = Tn2; %end of season, only expected value for now

        %First timestep need history function.
        timestep = 1;
        tspan = 0:timestep:Tn;

        %run simulation
        sol=dde23(@(t,x,x_hists)  gen_rotifer_onestrain(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @history_fun, tspan);
        
        realized_egg_production(a, season) = sol.y(1,end);
        phi = min(sol.y(1,end), B_max) ./length(sk); %number of eggs to start next season 


end


subplot(3,2,3)
t = length(threshold); 
m = length(mixis); 
for i = 1:t %want to have lines connecting dots
    hold on 
    plot(mixis_and_thresh(1,i:t:end),realized_egg_production(i:t:end)', 'color', colorvals(i,:))
end
scatter(mixis_and_thresh(1,:), realized_egg_production, 40, mixis_and_thresh(2,:), 'filled')
ylabel('Cumulative egg production')
xlabel('Mixis ratio')
ylim([0 250])

%% Panel E - Gi = 0, Long Season 

realized_egg_production = nan(a_num, 1, 1); 

for a=1:a_num

    m_i=mixis_and_thresh(1,a); %mixis_and_thresh(1,a);
    T_i=mixis_and_thresh(2,a);

    R_y = B_max; %number of resting eggs that we start with. 
    phi = R_y ./length(sk); % figure out how many resting eggs will emerge each hatch day

        x_hists = zeros(x_width-2,1); %past state space, not including generation 1. Resets every season

        Tn = Tn3; %end of season, only expected value for now

        %First timestep need history function.
        timestep = 1;
        tspan = 0:timestep:Tn;

        %run simulation
        sol=dde23(@(t,x,x_hists)  gen_rotifer_onestrain(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @history_fun, tspan);
        
        realized_egg_production(a, season) = sol.y(1,end);
        phi = min(sol.y(1,end), B_max) ./length(sk); %number of eggs to start next season 


end

subplot(3,2,5)
t = length(threshold); 
m = length(mixis); 
for i = 1:t %want to have lines connecting dots  
    hold on 
    plot(mixis_and_thresh(1,i:t:end),realized_egg_production(i:t:end)', 'color', colorvals(i,:))
end
scatter(mixis_and_thresh(1,:), realized_egg_production, 40, mixis_and_thresh(2,:), 'filled')
ylabel('Cumulative egg production')
xlabel('Mixis ratio')
h = colorbar ;
h.Label.String = 'Threshold value';
ylim([0 400])

%% One more column with different block phenotype 
% Copied from above but now G_i is set to 8

%choose generation delay and presize things
G_i = 8;
x_width = 2+2*(G_i+1);   % (R, M, J, A, J, A, ...) 0 to G_i
J_ind = [3:2:x_width];
A_ind = [4:2:x_width];

%% Panel B - Gi = 8, Short Season 

%reset output values
realized_season_length = nan(a_num, 1, 1); 
realized_egg_production = nan(a_num, 1, 1); 


for a=1:a_num %for each phenotype 

    m_i=mixis_and_thresh(1,a); %mixis_and_thresh(1,a);
    T_i=mixis_and_thresh(2,a);

    R_y = B_max; %number of resting eggs that we start with. 
    phi = R_y ./length(sk); % figure out how many resting eggs will emerge each hatch day

        x_hists = zeros(x_width-2,1); %past state space, not including generation 1. Resets every season

        Tn = Tn1; %end of season, only expected value for now

        %First timestep need history function.
        timestep = 1;
        tspan = 0:timestep:Tn;

        %run simulation
        sol=dde23(@(t,x,x_hists)  gen_rotifer_onestrain(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @history_fun, tspan);
        
        realized_egg_production(a, season) = sol.y(1,end);
        phi = min(sol.y(1,end), B_max) ./length(sk); %number of eggs to start next season 

end

%plot
subplot(3,2,2)
t = length(threshold); 
m = length(mixis); 
for i = 1:t %want to have lines connecting dots  
    hold on 
    plot(mixis_and_thresh(1,i:t:end),realized_egg_production(i:t:end)', 'color', colorvals(i,:))
end
scatter(mixis_and_thresh(1,:), realized_egg_production, 40, mixis_and_thresh(2,:), 'filled')
colormap(colorvals)
ylabel('Cumulative egg production')
xlabel('Mixis ratio')
h = colorbar ;
h.Label.String = 'Threshold value';
title('G_i = 8')
ylim([0 1])


%% Panel D - G_i = 8 Medium Season 
realized_egg_production = nan(a_num, 1, 1); 

for a=1:a_num %for each phenotype 

    m_i=mixis_and_thresh(1,a); %mixis_and_thresh(1,a);
    T_i=mixis_and_thresh(2,a);

    R_y = B_max; %number of resting eggs that we start with. 
    phi = R_y ./length(sk); % figure out how many resting eggs will emerge each hatch day

        x_hists = zeros(x_width-2,1); %past state space, not including generation 1. Resets every season

        Tn = Tn2; %end of season, only expected value for now

        %First timestep need history function.
        timestep = 1;
        tspan = 0:timestep:Tn;

        %run simulation
        sol=dde23(@(t,x,x_hists)  gen_rotifer_onestrain(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @history_fun, tspan);
        
        realized_egg_production(a, season) = sol.y(1,end);
        phi = min(sol.y(1,end), B_max) ./length(sk); %number of eggs to start next season 


end

%plot 
subplot(3,2,4)
t = length(threshold); 
m = length(mixis); 
for i = 1:t %want to have lines connecting dots  
    hold on 
    plot(mixis_and_thresh(1,i:t:end),realized_egg_production(i:t:end)', 'color', colorvals(i,:))
end
scatter(mixis_and_thresh(1,:), realized_egg_production, 40, mixis_and_thresh(2,:), 'filled')
colormap(colorvals)
ylabel('Cumulative egg production')
xlabel('Mixis ratio')
h = colorbar ;
h.Label.String = 'Threshold value';
ylim([0 250])



%% Panel F - G_i = 8 Long Season 

realized_egg_production = nan(a_num, 1, 1); 

for a=1:a_num %for each phenotype 

    m_i=mixis_and_thresh(1,a); 
    T_i=mixis_and_thresh(2,a);

    R_y = B_max; %number of resting eggs that we start with.
    phi = R_y ./length(sk); % figure out how many resting eggs will emerge each hatch day

        x_hists = zeros(x_width-2,1); %past state space, not including generation 1. Resets every season

        Tn = Tn3; %end of season, only expected value for now

        %First timestep need history function.
        timestep = 1;
        tspan = 0:timestep:Tn;

        %run simulation
        sol=dde23(@(t,x,x_hists)  gen_rotifer_onestrain(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @history_fun, tspan);
        
        realized_egg_production(a, season) = sol.y(1,end);
        phi = min(sol.y(1,end), B_max) ./length(sk); %number of eggs to start next season 

end


%plot 
subplot(3,2,6)
t = length(threshold); 
m = length(mixis); 
for i = 1:t %want to have lines connecting dots  
    hold on 
    plot(mixis_and_thresh(1,i:t:end),realized_egg_production(i:t:end)', 'color', colorvals(i,:))
end
scatter(mixis_and_thresh(1,:), realized_egg_production, 40, mixis_and_thresh(2,:), 'filled')
colormap(colorvals)
ylabel('Cumulative egg production')
xlabel('Mixis ratio')
h = colorbar ;
h.Label.String = 'Threshold value';
ylim([0 400])





%% history functions 


function v = history_fun(t) %not really necessary now that we are not including generation 0 in diffe solver
global sk x_width phi
  if t< sk(1)
      v =  zeros(x_width-2, 1); 
  elseif t == sk(1)
      v = zeros(x_width-2, 1); 
      %v(3) = phi; 
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


clear all

% I want one plot that shows the tradeoffs between two strategies

% Define the parameters 
bmax=2.1; q=0.5; K=500; c=0.2; H=31; tau = 1; 
% H is the final day of the season that hatching can happening 
B_max = 1; %maximum number of resting eggs that survive to next season. 
%Hatching schedule 
t0=0; %0 is season begins on first hatch day 
sk = t0:1:H; %days that hatching will occur


%choose mixis phenotypes for both strains
pheno_1 = [0.17 0];
pheno_2 = [0.23 0];
G_i = [3 4];

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

num_seasons = 1;

R_y = B_max; %number of resting eggs that we start with.

%split eggs equally
phi_1 = (R_y*.5)./length(sk); % figure out how many resting eggs will emerge each hatch day
phi_2 = (R_y* .5)./length(sk); % figure out how many resting eggs will emerge each hatch day
phi = [phi_1 phi_2];

J_ind_1 = [3:2:half_x_width];
A_ind_1 = [4:2:half_x_width];
J_ind_2 = [half_x_width+3:2:x_width];
A_ind_2 = [half_x_width+4:2:x_width];


m_i = [pheno_1(1) pheno_2(1)];
T_i = [pheno_1(2) pheno_2(2)];

for season = 1:num_seasons

    x_hists = zeros(x_width-4,1); %past state space, not including generation 0 for either strain. Resets every season

    Tn = 51 %end of season
    timestep = 1;
    tspan = 0:timestep:Tn;

    m_i = [pheno_1(1) pheno_2(1)];
    T_i = [pheno_1(2) pheno_2(2)];


    %run simulation
    sol=dde23(@(t,x,x_hists) gen_rotifer_twostrains(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @(t) history_fun(t,sk,x_width), tspan);


    R_out = sol.y(1,end) + sol.y(half_x_width+1, end);
    R_y = min(R_out, B_max); %number of eggs to start next season

    if R_out <= 0
        prop_2 = NaN;
        break
    end
    prop_1 = sol.y(1,end)./R_out;
    prop_2 = sol.y(half_x_width+1, end)./R_out;

    phi = [R_y*prop_1./length(sk) R_y*prop_2./length(sk)];

    if phi ==0 %if no eggs survive, no need to continue to next season
        phi
        break
    elseif prop_2 == 0 || prop_1 == 0 %if either prop is 0, no need to continue. invasion unsucessful (or super success)
        prop_2
        break
    end

end %seasons

%% plotting  within season dynamics
t = sol.x; 
%need to add stem juveniles and adults back into state variable to plot
%total pop size 
As = nan(2,length(t));
Js = nan(2, length(t));
for j = 1:length(t) %could make this more efficient
    t_i = t(j);
    Js(1,j) =  how_many_J(t_i, phi(1), sk, tau, q);
    As(1,j) = how_many_A(t_i, phi(1), sk, tau, q);
    Js(2,j) =  how_many_J(t_i, phi(2), sk, tau, q);
    As(2,j) = how_many_A(t_i, phi(2), sk, tau, q);
end

y = sol.y; 
y1 = [sol.y(1:2, :); Js(1,:); As(1,:); sol.y(3:half_x_width,:)];
y2 = [sol.y(half_x_width+1:half_x_width+2, :); Js(2,:); As(2,:); sol.y(half_x_width+3:end,:)]; 
hold on 

gold = [245 173 45]./255;

figure
hold on
plot(t, y1(1,:), 'color', gold, 'linewidth', 2) %Resting Eggs
hold on 
plot(t, sum(y1(2:end, :)), 'linewidth', 2, 'color', 'k')
hold on
plot(t, y2(1,:), '--', 'color', gold, 'linewidth', 2) %Resting Eggs
hold on 
plot(t, sum(y2(2:end, :)), '--', 'linewidth', 2, 'color', 'k')
ylabel('Individuals L^{-1}')
xlabel('Time (d)')
ylim([0 250])
xlim([0 51]); 

h = legend({'Resting eggs'; 'Adults'; ''; '';});
h.Location = 'NorthWest';

box on 

set(findall(gcf,'-property','FontSize'),'FontSize',14)
fontname('Arial')


%% history functions 


function v = history_fun(t, sk, x_width) %not really necessary now that we are not including generation 0 in diffe solver
  if t< sk(1)
      v =  zeros(x_width-4, 1); 
  elseif t == sk(1)
      v = zeros(x_width-4, 1); 
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

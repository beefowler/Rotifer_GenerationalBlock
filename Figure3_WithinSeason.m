% Create Figure 3 of the manuscript showing within season dynamics of
% resting eggs, mictic females, amictic females, and juveniles. 


clear all
figure


%choosing colors    %paul tol map 
map = [238, 204, 102;
    238, 153, 170;
    102, 153, 204; 
    153, 119. 0;
    153, 68, 85; 
    0, 68, 136
]./255;
gold = [245 173 45]./255;


% Define the parameters
bmax=2.1; q=0.5; K=500; c=0.2; H=31; tau = 1;
% H is the final day of the season that hatching can happening
B_max = 1; %maximum number of resting eggs that survive to next season. 
%Hatching schedule
global sk x_width phi
t0=0; %0 is season begins on first hatch day
sk = t0:1:H; %days that hatching will occur


%phenotypes
m_i_vals = [.3 .7];
T_i = 100;
Gvals = [0 0];  


for a = 1:2 %cycle through two mixis options
    G_i = Gvals(a);
    m_i = m_i_vals(a)

    x_width = max(2+2*(G_i+1), 6);   % (R, M, J, A, J, A, ...) 0 to G_i
    J_ind = [3:2:x_width];
    A_ind = [4:2:x_width];


    R_y = B_max; %number of resting eggs that we start with. Can't tell what this should be?
    phi = R_y ./length(sk); % figure out how many resting eggs will emerge each hatch day

    x_hists = zeros(x_width-2,1); %past state space, not including generation 1. Resets every season

    Tn = 51; %end of season, randomly selected from range

    %First timestep need history function.
    timestep = 1;
    tspan = 0:timestep:Tn;

    %run simulation
    sol=dde23(@(t,x,x_hists) gen_rotifer_onestrain(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @history_fun, tspan);

    % some plotting for within season dynamics
    t = sol.x;
    %need to add stem juveniles and adults back into state variable to plot
    As = nan(1,length(t));
    Js = nan(1, length(t));
    for j = 1:length(t) %could make this more efficient
        t_i = t(j);
        Js(j) =  how_many_J(t_i, phi, sk, tau, q);
        As(j) = how_many_A(t_i, phi, sk, tau, q);
    end
    y = sol.y;
    y = [sol.y(1:2, :); Js; As; sol.y(3:end,:)]; 

    subplot(2,1,a)
    hold on
    plot(t, y(1,:), 'color', gold, 'linewidth', 3) %Resting Eggs
    hold on
    plot(t, y(2,:), '-.', 'color', map(5,:),'LineWidth', 3) %magenta for mictic
    %plot(t, sum(y(J_ind,:)), '-.', 'color', map(3,:), 'LineWidth', 1.5) %"g for juveniles"
    plot(t, sum(y(A_ind,:)), '-.', 'color', map(3,:),'LineWidth', 3) %blue for adults amictic
    plot(t, sum(y(2:end, :)), 'linewidth', 2, 'color', 'k')
    ylabel('Individuals L^{-1}')
    xlabel('Time (d)')

    ylim([ 0 340])
    xlim([0 42])

    title(['m_1= ' num2str(m_i)])
    h = legend({'Resting eggs'; 'Mictic adults'; 'Amictic adults'; 'Total population'; 'Density threshhold'});
    h.Location = 'NorthWest';

    hold on
    plot([t(1) t(end)], [T_i T_i], '--', 'Color', [.3 .3 .3])
    box on
end

set(findall(gcf,'-property','FontSize'),'FontSize',14)
fontname('Arial')

ylim([ 0 340])
xlim([0 42])

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


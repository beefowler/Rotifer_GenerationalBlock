clear all

%How does our change to the mictic block affect the otpimal invadable
%phenotypes for the same environmental paramters as Serra et al?

% Define the parameters 
bmax=2.1; q=0.5; K=500; c=0.2; H=31; tau = 1; 
% H is the final day of the season that hatching can happening 
B_max = 1; %maximum number of resting eggs that survive to next season. 
%Hatching schedule 
global phi half_x_width
t0=0; %0 is season begins on first hatch day 
sk = t0:1:H; %days that hatching will occur

num_experiments = 20; 

%choose generation delay for both phenotypes
G_i = [0 0];
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

%phenotypes to try
mixis = [0:.01:0.2];

threshold = [0];
[M,T] = meshgrid(mixis, threshold);
n = length(mixis)*length(threshold);
MM = reshape(M, [1, n]);
TT = reshape(T, [1, n]);
mixis_and_thresh = [MM; TT];


num_seasons = 40;
a_num = length(mixis_and_thresh); %need this defined for parfor loop
for a=1:a_num % for each phenotype as resident
    parfor b = 1:a_num %for each phenotype as invader 
        if a == b
            continue %no need to invade self 
        end
    

        pheno_1 = mixis_and_thresh(:,a);
        pheno_2 = mixis_and_thresh(:,b);

        %check if we've already done (so we can restart if script crashes) 
        if exist(['PairwiseInvasion_M/outputs_' num2str(pheno_1(1)) '_' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_' num2str(pheno_2(2)) '.mat'], 'file')
            continue
        end

        % initialize outputs 
        invasion_prop = nan(1,num_experiments);
        total_pop = nan(1,num_experiments); 

        for x = 1:num_experiments

            R_y = B_max; %number of resting eggs that we start with.

            %split eggs
            phi_1 = (R_y*.95)./length(sk); % figure out how many resting eggs will emerge each hatch day
            phi_2 = (R_y* .05)./length(sk); % figure out how many resting eggs will emerge each hatch day
            phi = [phi_1 phi_2];

            for season = 1:num_seasons

                x_hists = zeros(x_width-4,1); %past state space, not including generation 0 for either strain. Resets every season

                Tn = randi([10 51]); %end of season, randomly selected from range

                %First timestep need history function.
                timestep = 1;
                tspan = 0:timestep:Tn;
                %%

                m_i = [pheno_1(1) pheno_2(1)];
                T_i = [pheno_1(2) pheno_2(2)];

                %run simulation
                sol=dde23(@(t,x,x_hists) gen_rotifer_twostrains(t, x, x_hists,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk), tau, @(t) history_fun(t, sk, x_width), tspan);

                %%
                R_out = sol.y(1,end) + sol.y(half_x_width+1, end);

                R_y = min(R_out, B_max); %number of eggs to start next season

                prop_1 = sol.y(1,end)./R_out;
                prop_2 = sol.y(half_x_width+1, end)./R_out;

                phi = [R_y*prop_1./length(sk) R_y*prop_2./length(sk)];

                if phi ==0 %if no eggs survive, no need to continue to next season
                    break
                elseif prop_2 == 0 || prop_1 == 0 %if either prop is 0, no need to continue. invasion unsucessful (or super success)
                    break
                end

            end %seasons

            invasion_prop(x) = prop_2; 
            total_pop(x) = R_y
            %disp([num2str(x) ' out of ' num2str(num_experiments)])
        end
            parsave(pheno_1, pheno_2, invasion_prop, total_pop, num_seasons)

    end
end




%% saving function
function parsave(pheno_1, pheno_2, invasion_prop, total_pop, num_seasons)
    save(['PairwiseInvasion_M/outputs_' num2str(pheno_1(1)) '_' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_' num2str(pheno_2(2)) '.mat'], 'pheno_1', 'pheno_2', 'invasion_prop', 'total_pop', 'num_seasons')
end


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


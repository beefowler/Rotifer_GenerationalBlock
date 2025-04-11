%modifying to allow G_i = 0 to be possible

%this is copy of cp_fun_monomorphic to be modified for generation model
%simulates system of ODEs for given parameters

%Inputs:   (** indicate variables that are dependent on phenotype)
    % x, initial state space of system
        % x is size 2+2*(G_i) (Resting eggs + Myctic Adults, Juveniles, Amictic Adults)
        % times G_i or the number of generations we are keeping track of. 
        % (we won't need to keep counting to infinity, but it may be a lot)
        % In order R, M, J_1, A_1, J_2, A_2, etc. 
        %(output of ode will be size t_span .* width(x))
    % x_hist = history of state space, tau units ago until "now". 
    % tau = time, in hours?, until maturation from juveniles. 
    % bmax = birthrate in absense of density dependence
    % q = death rate of adults
    % K = carying capacity
    % c = cost of making a rest egg. 

    % G_i = number of generations until mixis, first generation that will produce eggs **
    % m_i = mixis ratio for offspring ** 
    % T_i = threshhold of onset of mixis **     
    % tunit = we need to know how many units back in X-hist is tau time ago
    % phi = size of hatching 
    % sk = birthdays 

%Output: dxdt = derivatives of state variables at time t during growing season. 


function dxdt =gen_rotifer_onestrain(t, x, x_hist,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk)

    %initialize output
    x_width = max(2+2*(G_i+1), 6);   % (R, M, J, A, J, A, ...) generation 1 to G_i, and if G_i = 0, still 6
    dxdt=zeros(x_width, 1); 

    %make sure inputs are the right orientation 
    if size(x_hist, 1) ~= x_width-2 
        x_hist = x_hist'; 
    end
    

    %go see how many generation zero there are and add those to the state
    %variables
    J_0 = how_many_J(t, phi, sk, tau, q);
    A_0 = how_many_A(t, phi, sk, tau, q); 

    X_all = [x(1:2); J_0; A_0; x(3:end)]; 
    X_hist_all = [x_hist(1:2); how_many_J(t-tau, phi, sk, tau, q); how_many_A(t-tau, phi, sk, tau, q); x_hist(3:end)];

    %resting eggs for next year
    dxdt(1)=c*(bmax-(bmax-q)*(sum(X_all(2:end)))/K)*x(2); % Resting eggs are collected in the first position in x-variable, depend on total pop and M or x(2)
    %previously this had been sum(x(2:end)), July 14, 2022, but Gen 0 needs
    %to count towards sum. 

    %name indeces for different classes, from j =0 to inf 
    J_ind = [3:2:x_width];
    A_ind = [4:2:x_width]; 
        % note x(end) or x(G_ind) is the only class that will make mictic females

        %check if conditions for mixis were met at time (t - tau)
        if G_i(1) == 0 %if no generational delay, both A bins might produce Mictic offspring at same rate m_i.
            if sum(X_hist_all(2:end))>=T_i  % if total water column population is >= threshhold
                %then mixis will happen for generations greater than G_i & asexual
                %reproduction will happen for "closer to stem" generations.
                m_ij_hist = [m_i m_i]; %both A can produce mixis offspring at mixis ratio
            else
                m_ij_hist = [0 0]; %all bins are zero.
            end
        else %Gi > 0 
            m_ij_hist = zeros(1,G_i); %different m value for each generation, is zero until j = G_i
            if sum(X_hist_all(2:end))>=T_i  % if total water column population is >= threshhold
                %then mixis will happen for generations greater than G_i & asexual
                %reproduction will happen for "closer to stem" generations.
                m_ij_hist = [m_ij_hist m_i]; %final bin of m_ij is mixing ratio
            else
                m_ij_hist = [m_ij_hist 0]; %all bins are zero.
            end
        end

    % let some mictic adults mature (m_ij_hist(end)*x_hist(A_ind(end)) is only nonzero entry)  & let some die 
        b_hist = max(bmax-(bmax-q)*(sum(X_hist_all(2:end))/K), 0); %birth rate tau units ago 
        dxdt(2) = b_hist*(m_ij_hist*X_hist_all(A_ind))*exp(-q*tau) - q*x(2); %first term sums automatically because of matrix multiplication 
    
    % produce juveniles, let some mature & die
        b_t = bmax-(bmax-q)*(sum(X_all(2:end))/K); %birth rate now
        b_t = max(b_t, 0); %just make sure birth rate doesn't go negative 

        dxdt(J_ind(2:end)) = b_t.*X_all(A_ind(1:end-1)) - b_hist.*X_hist_all(A_ind(1:end-1))*exp(-q.*tau) - q.*X_all(J_ind(2:end)); 
   
            %these three terms are: 
            % birth of juveniles now, from previous generation of adults. 
            % maturation of juveniles who were born tau units ago, from previous
            %       generation of adults
            % and death rate of juveniles now.

        %don't forget to pool the Amyctic adults & juveniles that are higher generations
        dxdt(J_ind(end)) = dxdt(J_ind(end)) + b_t.*X_all(A_ind(end)) - b_hist.*X_hist_all(A_ind(end))*exp(-q.*tau); 
            %add inviduals born from "furthest from stem" generations, and let
            %those individuals mature. 

    % let amyctic adults mature & some die      
    
        dxdt(A_ind(2:end)) = b_hist.*(X_hist_all(A_ind(1:end-1)).*(1-m_ij_hist(1:end-1))').*exp(-q*tau) - q.*X_all(A_ind(2:end)); 


        %I wrote this out formally, but 1-m_ij_hist(1:end-1) is always 1. 
        %those individuals will only make amictic offspring. 
    
        %Now add in the maturation of adults that were born tau units ago from the
        % "furthest from stem" generations or highest j bin adults. 
            dxdt(A_ind(end)) = dxdt(A_ind(end)) + b_hist*X_hist_all(A_ind(end)).*(1-m_ij_hist(end))*exp(-q*tau); 
   
        %I'm not completely sure this works mathematically: 
        %change in adults old enough to produce mictic offspring is equal to
        %production of amyctic adults by A(G_i -1) + production of amyctic adults
        % by A(>=G_i) - death. 
   
       %generally the code doesn't work if G_i = 0
        dxdt(3:4) = []; %we aren't calculating these anymore. 
        
end
%%


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


%this is model of rotifer population with generational block and density
%threshold, script modified such that we can have two
%phenotypes simulated simultaneously. 

%simulates system of ODEs for given parameters
%This function only works if there are two populations. 

%Inputs:   (** indicate variables that are dependent on phenotype)
    % t, time in days since start of growing season
    % x, initial state space of system. one phenotype then the other. 
        % x is size 2+2*(G_i(1)) + 2+2*(G_i(2)) (Resting eggs + Myctic Adults, Juveniles, Amictic Adults)
        % times G_i or the number of generations we are keeping track of. 
        % (we won't need to keep counting to infinity, but it may be a lot)
        % In order R, M, J_1, A_1, J_2, A_2, etc. then again R, M, J_1, etc
        %(output of ode will be size t_span .* width(x))
    % x_hist = history of state space, tau units ago until "now". 
    % tau = time in days until maturation from juveniles. 
    % bmax = birthrate in absense of density dependence
    % q = death rate of adults
    % K = carying capacity
    % c = cost of making a rest egg. 

    % G_i = number of generations until mixis, first generation that will
    % produce eggs **   size 2
    % m_i = mixis ratio for offspring ** size 2
    % T_i = threshhold of onset of mixis **     size 2
    % phi = size of hatching, size 2
    % sk = birthdays, %assumed to be the same for both phenotypes 

%Output: dxdt = derivatives of state variables at time t during growing season. 


function dxdt = gen_rotifer_twostrains(t, x, x_hist,tau,bmax,q,K,c,G_i,m_i,T_i, phi, sk)

    %first figure out if there are two populations. 
    num_pops = 2;
        % %If not enough distinct paramters were included, assume same for both populatiosn 
        ind = find([length(G_i) length(m_i) length(T_i) length(phi)] ~= num_pops);
        for i = ind
            if i == 1
                G_i = repmat(G_i, 1, num_pops);
            elseif i == 2
                m_i = repmat(m_i, 1, num_pops);
            elseif i == 3
                T_i = repmat(T_i, 1, num_pops);
            elseif i == 4
                phi = repmat(phi, 1, num_pops);
            end
        end

    x_width_1 = max(2+2*(G_i(1)+1), 6);   % (R, M, J, A, J, A, ...) generation 0 to G_i for each group , if G_i = 0, still need 6 
    x_width_2 = max(2+2*(G_i(2)+1), 6); %if G_i = 0, still need 6 state variables 
    dxdt_1=zeros(x_width_1, 1); 
    dxdt_2=zeros(x_width_2, 1); 
    
    clear dxdt 

    x_1 = x(1:x_width_1-2);   %assume state variable has first phenotype followed by second phenotype 
    x_2 = x(x_width_1-1:end); 

   
    %make sure inputs are the right orientation 
    if size(x_hist, 1) ~= x_width_1+x_width_2 - 4 
        x_hist = x_hist'; 
    end
    
    x_hist_1 = x_hist(1:x_width_1-2);
    x_hist_2 = x_hist(x_width_1-1:end);

    clear x_hist %just so we don't accidentally use this  


    %go see how many generation zero there are and add those to the state
    %variables
    J_0_1 = how_many_J(t, phi(1), sk, tau, q);
    J_0_2 = how_many_J(t, phi(2), sk, tau, q);

    A_0_1 = how_many_A(t, phi(1), sk, tau, q); 
    A_0_2 = how_many_A(t, phi(2), sk, tau, q); 

    X_all_1 = [x_1(1:2); J_0_1; A_0_1; x_1(3:end)]; 
    X_all_2 = [x_2(1:2); J_0_2; A_0_2; x_2(3:end)];

    X_hist_all_1 = [x_hist_1(1:2); how_many_J(t-tau, phi(1), sk, tau, q); how_many_A(t-tau, phi(1), sk, tau, q); x_hist_1(3:end)];
    X_hist_all_2 = [x_hist_2(1:2); how_many_J(t-tau, phi(2), sk, tau, q); how_many_A(t-tau, phi(2), sk, tau, q); x_hist_2(3:end)];

    %resting eggs for next year
    total_pop_size = sum(X_all_1(2:end)) + sum(X_all_2(2:end)); %sum all individuals (not eggs) of both strains
    dxdt_1(1)=c*(bmax-(bmax-q)*(total_pop_size/K))*x_1(2); % Resting eggs are collected in the first position in x-variable, depend on total pop and M or x(2)
    dxdt_2(1)=c*(bmax-(bmax-q)*(total_pop_size/K))*x_2(2); % Density dependence includes ALL individuals, of both phenoytpes 

    %name indeces for different classes, from j =0 to inf 
    J_ind_1 = 3:2:x_width_1;
    A_ind_1 = 4:2:x_width_1; 
    J_ind_2 = 3:2:x_width_2;
    A_ind_2 = 4:2:x_width_2;
        % note x(end) or x(G_ind) is the only class that will make mictic females

    %check if conditions for mixis were met at time (t - tau)
    %Threshhold is based on total population size!
    total_pop_size_hist = sum(X_hist_all_1(2:end)) + sum(X_hist_all_2(2:end)); %sum all individuals (not eggs) of both strains
    

    %First determine mixis rates in strain 1
     if G_i(1) == 0 %if no generational delay, both A bins might produce Mictic offspring at same rate m_i. 
        if sum(total_pop_size_hist)>=T_i(1)  % if total water column population is >= threshhold 
        %then mixis will happen for generations greater than G_i & asexual
        %reproduction will happen for "closer to stem" generations. 
        m_ij_hist_1 = [m_i(1) m_i(1)]; %both A can produce mixis offspring at mixis ratio
        else 
        m_ij_hist_1 = [0 0]; %all bins are zero. 
         end
    else %If Gi > 0, then first value will always be 0. 
    m_ij_hist_1 = zeros(1,G_i(1)); %different m value for each generation, is zero until j = G_i
    if sum(total_pop_size_hist)>=T_i(1)  % if total water column population is >= threshhold 
        %then mixis will happen for generations greater than G_i & asexual
        %reproduction will happen for "closer to stem" generations. 
        m_ij_hist_1 = [m_ij_hist_1 m_i(1)]; %final bin of m_ij is mixing ratio
    else 
        m_ij_hist_1 = [m_ij_hist_1 0]; %all bins are zero. 
    end
     end

     %now determine mixis rates in strain 2
    if G_i(2) == 0 %if no generational delay, both A bins might produce Mictic offspring at same rate m_i. 
        if sum(total_pop_size_hist)>=T_i(2)  % if total water column population is >= threshhold 
        %then mixis will happen for generations greater than G_i & asexual
        %reproduction will happen for "closer to stem" generations. 
        m_ij_hist_2 = [m_i(2) m_i(2)]; %both A can produce mixis offspring at mixis ratio
        else 
        m_ij_hist_2 = [0 0]; %all bins are zero. 
         end
    else %If Gi > 0, then first value will always be 0. 
    m_ij_hist_2 = zeros(1,G_i(2)); %different m value for each generation, is zero until j = G_i
    if sum(total_pop_size_hist)>=T_i(2)  % if total water column population is >= threshhold 
        %then mixis will happen for generations greater than G_i & asexual
        %reproduction will happen for "closer to stem" generations. 
        m_ij_hist_2 = [m_ij_hist_2 m_i(2)]; %final bin of m_ij is mixing ratio
    else 
        m_ij_hist_2 = [m_ij_hist_2 0]; %all bins are zero. 
    end
    end


    % let some mictic adults mature (m_ij_hist(end)*x_hist(A_ind(end)) is only nonzero entry)  & let some die 
        b_hist = bmax-(bmax-q)*(total_pop_size_hist/K); %birth rate tau units ago 
        b_hist = max(b_hist, 0); %check that birth doesn't go negative 
        dxdt_1(2) = b_hist*(m_ij_hist_1*X_hist_all_1(A_ind_1))*exp(-q*tau) - q*x_1(2); %first term sums atuomatically because of matrix multiplication 
        dxdt_2(2) = b_hist*(m_ij_hist_2*X_hist_all_2(A_ind_2))*exp(-q*tau) - q*x_2(2); %first term sums atuomatically because of matrix multiplication 

    % produce juveniles, let some mature & die
        b_t = bmax-(bmax-q)*(total_pop_size/K); %birth rate now
        b_t = max(b_t, 0); %check that b_t doesn't go negative 

        dxdt_1(J_ind_1(2:end)) = b_t.*X_all_1(A_ind_1(1:end-1)) - b_hist.*X_hist_all_1(A_ind_1(1:end-1))*exp(-q.*tau) - q.*X_all_1(J_ind_1(2:end)); 
        dxdt_2(J_ind_2(2:end)) = b_t.*X_all_2(A_ind_2(1:end-1)) - b_hist.*X_hist_all_2(A_ind_2(1:end-1))*exp(-q.*tau) - q.*X_all_2(J_ind_2(2:end)); 

            %these three terms are: 
            % birth of juveniles now, from previous generation of adults. 
            % maturation of juveniles who were born tau units ago, from previous
            %       generation of adults
            % and death rate of juveniles now.

        %don't forget to pool the Amyctic adults & juveniles that are higher generations
        dxdt_1(J_ind_1(end)) = dxdt_1(J_ind_1(end)) + b_t.*X_all_1(A_ind_1(end)) - b_hist.*X_hist_all_1(A_ind_1(end))*exp(-q.*tau); 
        dxdt_2(J_ind_2(end)) = dxdt_2(J_ind_2(end)) + b_t.*X_all_2(A_ind_2(end)) - b_hist.*X_hist_all_2(A_ind_2(end))*exp(-q.*tau); 

            %add inviduals born from "furthest from stem" generations, and let
            %those individuals mature. 

    % let amyctic adults mature & some die      
    
        dxdt_1(A_ind_1(2:end)) = b_hist.*(X_hist_all_1(A_ind_1(1:end-1)).*(1-m_ij_hist_1(1:end-1))').*exp(-q*tau) - q.*X_all_1(A_ind_1(2:end)); 
        dxdt_2(A_ind_2(2:end)) = b_hist.*(X_hist_all_2(A_ind_2(1:end-1)).*(1-m_ij_hist_2(1:end-1))').*exp(-q*tau) - q.*X_all_2(A_ind_2(2:end)); 

        %In this version, 1-m_ij_hist(1:end-1) is not necessarily always 1.
        %Not if G_i = 0, then stem individuals produce mictic offspring if
        %threshold is surpassed etc. 
    
        %Now add in the maturation of adults that were born tau units ago from the
        % "furthest from stem" generations or highest j bin adults. 
            dxdt_1(A_ind_1(end)) = dxdt_1(A_ind_1(end)) + b_hist*X_hist_all_1(A_ind_1(end)).*(1-m_ij_hist_1(end))*exp(-q*tau); 
            dxdt_2(A_ind_2(end)) = dxdt_2(A_ind_2(end)) + b_hist*X_hist_all_2(A_ind_2(end)).*(1-m_ij_hist_2(end))*exp(-q*tau); 

        %change in adults old enough to produce mictic offspring is equal to
        %production of amyctic adults by A(G_i -1) + production of amyctic adults
        % by A(>=G_i) - death. 
   
        dxdt_1(3:4) = []; %we aren't calculating these anymore. %see functions below 
        dxdt_2(3:4) = []; %we aren't calculating these anymore. 

        dxdt = [dxdt_1; dxdt_2]; 
      
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


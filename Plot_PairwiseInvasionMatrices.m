
% Gather outputs together to form Matricees for PIP 


% First do Baseline: G = 0, Variable m

outputpath = 'PairwiseInvasion_M/';

mixis = [0:.01:0.21];
threshold = 0;
[M,T] = meshgrid(mixis, threshold);
n = length(mixis)*length(threshold);
MM = reshape(M, [1, n]);
TT = reshape(T, [1, n]);
mixis_and_thresh = [MM; TT];

mean_output = nan(n); 
std_output = nan(n); 
mean_total_pop = nan(n); 
r_string = nan(1,n*n); 
i_string = nan(1,n*n); 
success_string = nan(n*n);
count = 1; 
for r=1:n
    for i = 1:n

        pheno_1 = mixis_and_thresh(:,r);
        pheno_2 = mixis_and_thresh(:,i);
        
        if exist([outputpath 'outputs_' num2str(pheno_1(1)) '_' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_' num2str(pheno_2(2)) '.mat'])
        load([outputpath 'outputs_' num2str(pheno_1(1)) '_' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_' num2str(pheno_2(2)) '.mat'])

        mean_output(r,i) = mean(invasion_prop, 'omitnan'); 
        std_output(r,i) = std(invasion_prop, 'omitnan'); 

        r_string(count) = mixis(r); 
        i_string(count) = mixis(i); 
        success_string(count) = mean(invasion_prop, 'omitnan')>=0.05;
        end
    count = count+1; 

    end
    disp(r)
end


success_output = mean_output; 
success_output(success_output<0.05) = NaN;
success_output(~isnan(success_output)) = 1;
success_output(isnan(success_output)) = 0;

% Plot
subplot(2,2,1)

%manually choose values to draw curve between x and o on PIP 
Xvalues = [-0.0:0.01:0.2]; 
Yvalues3 = [0.205 0.205 0.205 0.205  0.205 0.205 0.205 ...
    0.195  0.175  0.155 0.125 0.11 ...
    0.105 0.095 0.075 0.065 0.055 0.055 0.055 0.045 0.045];
Yvalues = [Xvalues' Yvalues3'-Xvalues']; 
a = area(Xvalues, Yvalues);
a(2).FaceColor = [0 0 0];
a(1).FaceColor = [1 1 1]; 


ylim([0 0.2])
xlim([0 0.2])
hold on 

scatter(r_string(success_string ==1), i_string(success_string==1), 10, [.5 .5 .5], 'filled')
hold on 
scatter(r_string(success_string ==0), i_string(success_string==0), 25, [.5 .5 .5], 'x')
ylabel('Invader mixis ratio, m_2')
xlabel('Resident mixis ratio, m_1')
title({'No generational block: G_i = 0,'; 'No mixis threshold: T_i = 0'})

%%  Panel for m = 0.11, Variable G 

outputpath = 'PairwiseInvasion_G/';

mixis = 0.11;  %This panel is for m = 0.11
G_vals = [0:7];
[M,G] = meshgrid(mixis, G_vals);
n = length(mixis)*length(G_vals);
MM = reshape(M, [1, n]);
GG = reshape(G, [1, n]);
mixis_and_G = [MM; GG];

mean_output = nan(n); 
std_output = nan(n); 
mean_total_pop = nan(n); 

r_string = nan(1,n*n); 
i_string = nan(1,n*n); 
success_string = nan(n*n);
count = 1; 
for r=1:n
    for i = 1:n

        pheno_1 = mixis_and_G(:,r);
        pheno_2 = mixis_and_G(:,i);
        
        if exist([outputpath 'outputs_' num2str(pheno_1(1)) '_G' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_G' num2str(pheno_2(2)) '.mat'])
        load([outputpath 'outputs_' num2str(pheno_1(1)) '_G' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_G' num2str(pheno_2(2)) '.mat'])

        mean_output(r,i) = mean(invasion_prop, 'omitnan'); 
        std_output(r,i) = std(invasion_prop, 'omitnan'); 

        r_string(count) = G_vals(r);
        i_string(count) = G_vals(i);
        success_string(count) = mean(invasion_prop, 'omitnan')>=0.05;

        end
            count = count+1; 
    end
    disp(r)
end

% Plot
subplot(2,2,2)
cla

%manually choose values to draw curve between x and o on PIP 
Xvalues = [0:7]; 
Yvalues3 = [6.5 5.5 4.5 1.5 -.5 0 0 0 ];
Yvalues = [Xvalues' Yvalues3'-Xvalues']; 
a = area(Xvalues, Yvalues)
a(2).FaceColor = [0 0 0];
a(1).FaceColor = [1 1 1]; 

ylim([0 7])
xlim([0 7])
hold on 

scatter(r_string(success_string ==1), i_string(success_string==1), 10, [.5 .5 .5], 'filled')
hold on 
scatter(r_string(success_string ==0), i_string(success_string==0), 25, [.5 .5 .5], 'x')
ylabel('Invader generational block, G_2')
xlabel('Resident generational block, G_1')
title({'Mixis ratio: m_i = 0.11,'; 'No mixis threshold: T_i = 0'})



%%  Panel for G_i = 3, Variable m 


outputpath = 'PairwiseInvasion_M_Gis3/';

mixis = [0.09:.01:0.21];
threshold = 0;
[M,T] = meshgrid(mixis, threshold);
n = length(mixis)*length(threshold);
MM = reshape(M, [1, n]);
TT = reshape(T, [1, n]);
mixis_and_thresh = [MM; TT];

mean_output = nan(n); 
std_output = nan(n); 
mean_total_pop = nan(n); 
r_string = nan(1,n*n); 
i_string = nan(1,n*n); 
success_string = nan(n*n);
count = 1; 
for r=1:n
    for i = 1:n

        pheno_1 = mixis_and_thresh(:,r);
        pheno_2 = mixis_and_thresh(:,i);

        if exist([outputpath 'outputs_' num2str(pheno_1(1)) '_' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_' num2str(pheno_2(2)) '.mat'])
        load([outputpath 'outputs_' num2str(pheno_1(1)) '_' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_' num2str(pheno_2(2)) '.mat'])

        %mean_total_pop(r,i) = mean(total_pop, 'omitnan');
        mean_output(r,i) = mean(invasion_prop, 'omitnan'); 
        std_output(r,i) = std(invasion_prop, 'omitnan'); 

        r_string(count) = mixis(r); 
        i_string(count) = mixis(i); 
        success_string(count) = mean(invasion_prop, 'omitnan')>=0.05;

        end
    count = count+1; 

    end
    disp(r)
end

% Plot
subplot(2,2,3)

%manually choose values to draw curve between x and o on PIP 
Xvalues = [0 0.09:0.01:0.21]; 
Yvalues3 = [.215 0.215 0.215 0.215 0.215 0.215 0.215 0.205 0.175 0.17 0.145 0.115 0.115 0.115];
Yvalues = [Xvalues' Yvalues3'-Xvalues']; 
a = area(Xvalues, Yvalues);
a(2).FaceColor = [0 0 0];
a(1).FaceColor = [1 1 1]; 


ylim([0.09 0.2])
xlim([0.09 0.2])
hold on 

scatter(r_string(success_string ==1), i_string(success_string==1), 10, [.5 .5 .5], 'filled')
hold on 
scatter(r_string(success_string ==0), i_string(success_string==0), 25, [.5 .5 .5], 'x')
ylabel('Invader mixis ratio, m_2')
xlabel('Resident mixis ratio, m_1')
title({'No generational block: G_i = 3,'; 'No mixis threshold: T_i = 0'})


%%  Panel for m = 0.17, Variable G 

outputpath = 'PairwiseInvasion_G/';
%phenotypes to try, for now using same values as in Figure 4
mixis = 0.17;
G_vals = [0:7];
[M,G] = meshgrid(mixis, G_vals);
n = length(mixis)*length(G_vals);
MM = reshape(M, [1, n]);
GG = reshape(G, [1, n]);
mixis_and_G = [MM; GG];

mean_output = nan(n); 
std_output = nan(n); 
mean_total_pop = nan(n); 

r_string = nan(1,n*n); 
i_string = nan(1,n*n); 
success_string = nan(n*n);
count = 1; 
for r=1:n
    for i = 1:n

        pheno_1 = mixis_and_G(:,r);
        pheno_2 = mixis_and_G(:,i);
        
        if exist([outputpath 'outputs_' num2str(pheno_1(1)) '_G' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_G' num2str(pheno_2(2)) '.mat'])
        load([outputpath 'outputs_' num2str(pheno_1(1)) '_G' num2str(pheno_1(2)) 'by' num2str(pheno_2(1)) '_G' num2str(pheno_2(2)) '.mat'])

        %mean_total_pop(r,i) = mean(total_pop, 'omitnan');
        mean_output(r,i) = mean(invasion_prop, 'omitnan'); 
        std_output(r,i) = std(invasion_prop, 'omitnan'); 

        r_string(count) = G_vals(r);
        i_string(count) = G_vals(i);
        success_string(count) = mean(invasion_prop, 'omitnan')>=0.05;

        end
            count = count+1; 
    end
    disp(r)
 end


% Plot
subplot(2,2,4)
cla

%manually choose values to draw curve between x and o on PIP 
Xvalues = [0:7]; 
Yvalues3 = [7 7.5 6.5 5.5 1.5 -.5 0 0 ];
Yvalues = [Xvalues' Yvalues3'-Xvalues']; 
a = area(Xvalues, Yvalues);
a(2).FaceColor = [0 0 0];
a(1).FaceColor = [1 1 1]; 

ylim([0 7])
xlim([0 7])
hold on 

scatter(r_string(success_string ==1), i_string(success_string==1), 10, [.5 .5 .5], 'filled')
hold on 
scatter(r_string(success_string ==0), i_string(success_string==0), 25, [.5 .5 .5], 'x')
ylabel('Invader generational block, G_2')
xlabel('Resident generational block, G_1')
title({'Mixis ratio: m_i = 0.17,'; 'No mixis threshold: T_i = 0'})
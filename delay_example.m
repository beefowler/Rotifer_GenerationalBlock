%Script to compare "simple" example of difference in delays from two block
%types: temporal block vs generational block. 

%% Up here is where we run the models and plot the results %% 
%    Model functions are in the following section   %

%paul tol color map 
map = [238, 204, 102;
    238, 153, 170;
    102, 153, 204; 
    153, 119. 0;
    153, 68, 85; 
    0, 68, 136
]./255;
gold = [245 173 45]./255; % want gold to be darker 

%choose parameter values 
m = .4; %mortality rate
tau = 1; %maturation delay, tau has to be 1 

tspan = 0:1:40; %time period to run simulation 

figure
subplot(3,1,1)
hold on 

%simulate time block population 
sol_tb = dde23(@(t,x, x_hist) time_block(t, x, x_hist, m), tau, @tb_history, tspan);
plot(sol_tb.x, sol_tb.y, '--','color', 'k', 'linewidth',3.5)

%plot stem adults exponential decay
plot(0:10, exp(-m*(0:10)), 'color', map(3,:),'linewidth', 3.5)

%simulate generational block population 
sol_gb = dde23(@(t,x, x_hist) generation_block(t, x, x_hist, m), tau, @gb_history, tspan);
hold on 
plot(sol_gb.x, sol_gb.y, 'color', map(6,:), 'linewidth',3.5)

%plot time block one again so its on top, don't repeat in legend
plot(sol_tb.x, sol_tb.y, '--', 'color', 'k', 'linewidth', 3.5)


%figure formatting
ylabel('Adults (L^{-1})')
lgd = legend({'All adults, $$x$$'; 'Stem adults, $$x_{0}$$'; 'Non-stem adults, $$x_1$$'}, 'Interpreter', 'LaTeX');
fontsize(lgd,14,'points')
set(lgd, 'FontName', 'Arial')
box on 
xlabel('Time (d)')

ylim([0 4])
xlim([0 10])

%ok and now we want to integrate egg production and get a vector over time
%time block (tb) first, integral from 1 to t
x_vals = sol_tb.x(sol_tb.x>=1);
y_vals = sol_tb.y(sol_tb.x>=1); 
tb_y = zeros(1,length(x_vals)); 
for i = 2:length(x_vals) %start at 2, since we need two x vals to make a trap
    tb_y(i) = trapz(x_vals(1:i), y_vals(1:i));
end
%generational block (gb): 
gb_y = zeros(1,length(sol_gb.x)); 
for i = 2:length(sol_gb.x) %start at 2, since we need two x vals to make a trap
    gb_y(i) = trapz(sol_gb.x(1:i), sol_gb.y(1:i));
end

%plot em up 
subplot(3,1,2)
hold on 
plot(x_vals, tb_y, '--', 'color', map(4,:), 'linewidth', 3.5)
plot(sol_gb.x, gb_y, 'color', gold, 'linewidth', 3.5) 

xlabel('Time (d)')
ylabel('Eggs (L^{-1})')

lgd = legend({'Temporal block, $$y$$'; 'Generational block, $$\tilde{y}$$'}, 'Interpreter', 'LaTeX');
fontsize(lgd,14,'points')
fontname(lgd, 'Arial') %this line doesn't work 
%title(lgd, 'Resting eggs')

set(findall(gcf,'-property','FontSize'),'FontSize',14)
fontname("Arial")
box on

ylim([0 17])
xlim([0 10])


% now we want consecutive seasons for third panel 
%% Mikes code below
r = 1;
T_sp = 10;
n = 40;
s = 0.5;

subplot(3,1,3)
cla

    x0_gb = 1;
    x0_tb = 1;

    y1 = zeros(1,n);
    y2 = y1;

    T = T_sp(j);
    for i = 1:n
        sol = nonstem(x0_gb,m,tau,T);
        y1(i) = genblocked_eggs(r,sol);
        x0_gb = s*y1(i);

        sol = nonstem(x0_tb,m,tau,T);
        y2(i) = timeblocked_eggs(genblocked_eggs(r,sol),x0_tb,r,q,tau,T);
        x0_tb = s*y2(i);
    end

    subplot(3,1,3)
    semilogy(1:n,y2./y1,'.','MarkerSize',20, 'color', gold)
    xlabel('Number of seasons'); ylabel('Fitness ratio')

lgd = legend({'Ratio: $$\frac{y}{\tilde{y}}$$'}, 'Interpreter', 'LaTeX');
fontsize(lgd,14,'points')

ylim([1 100])
xlim([0 40])
yticks('auto')


set(findall(gcf,'-property','FontSize'),'FontSize',14)


%% Below are the models themselves %% 


%one differential equation for fixed-time-block
function dxdt = time_block(t, x, x_hist, m)
    dxdt = exp(-m)*x_hist - m*x; 
end

%one differential equation for generational-block
function dxdt = generation_block(t, x, x_hist, m)
    if t < 1
        dxdt = exp(-m)*x_hist - m*x; 
    elseif t >= 1
        dxdt = exp(-m)*x_hist - m*x + exp(-m*t); 
    end

end

%time block history
function h = tb_history(t)
    if t<0 
        h = 0; 
    elseif t == 0
        h = 1; 
    end
end

%generational block history (could also do this with constant history
%format, but let's keep things parallel)
function h = gb_history(t)
    if t<=0 
        h = 0; 
    end
end



% Mike's version to get panel 3

%From Equation 7b in text
function y = timeblocked_eggs(z,x0,r,q,tau,T)
  y = z + x0*(r/q)*(exp(-q*tau) - exp(-q*T));
end

function y = genblocked_eggs(r,sol)
  y = trapz(sol.x,r*sol.y);
end

function x = stem(x0,q,t)
%stem number of stem adults decreases exponentially 
x = x0*exp(-q*t);
end

% number of nonstem adults
function sol = nonstem(x0,q,tau,T)
  sol = dde23(@ddefun,tau,0,[0 T]);
    function dxdt = ddefun(t,x,xdelay)
       dxdt = x0*exp(-q*t) + exp(-q*tau)*xdelay - q*x;
    end
end



close all
clear all
clc

%Import data and define closed-form solution
%Data for ordered spins (up = 1) and T = 1.0
filename = 'Project5_random_mass_check.xlsx';
sheet = 4;
xlRange = 'A2:B101';

[v,T,vT] = xlsread(filename, sheet, xlRange);
data_r=v(:,1); 
data_final=v(:,2); 
%data_phi=v(:,5);
%data_theta=v(:,7);

%Plot histograms
figure

xbins_r = 0:0.5:20;

[f_r,x_r] = hist(data_r,xbins_r)
dx_r = diff(x_r(1:2));
bar(x_r,f_r/sum(f_r*dx_r),'r')

hold on

xbins_final = 0:0.5:20;

[f_final,x_final] = hist(data_final,xbins_final)
dx_final = diff(x_final(1:2));
bar(x_final,f_final/sum(f_final*dx_final),'b')
legend('initial position','final position')

xlabel('Distance from cluster center (lightyears)', 'fontsize',14) % x-axis label
ylabel('probability','fontsize',14) % y-axis label





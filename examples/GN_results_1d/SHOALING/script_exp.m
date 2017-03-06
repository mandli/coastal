load #bea01gh.asc
load #bea02gh.asc
load ('SWE2.mat')
load ('t135_p1_01')


Time_exp = X_bea01gh(:,1) ./ sqrt(2.0) ;
length(Time_exp);
St_1 = X_bea02gh(1:1000,2) ./ 2.0 ;   % x = 4      % 2
St_2 = X_bea01gh(1:1000,4) ./ 2.0 ;   % x = 10.5   % 4
St_3 = X_bea02gh(1:1000,4) ./ 2.0 ;   % x = 12.5   % 5
St_4 = X_bea01gh(1:1000,5) ./ 2.0 ;   % x = 13.5    % 6
St_5 = X_bea01gh(1:1000,6) ./ 2.0 ;   % x = 15.7    % 7
St_6 = X_bea02gh(1:1000,6) ./ 2.0 ;   % x = 17.3    % 8
St_7 = X_bea01gh(1:1000,7) ./ 2.0 ;   % x = 19.0     % 9

time_plot_exp1 = Time_exp(1:1000,1) -23 - 0.7570;
time_plot_exp2 = Time_exp(1:1000,1) -23  + 0.5080;

num_val = 1100;
time = time(1:num_val,1) ./sqrt(9.81) -23;
time_swe = time_s(1:145,1) + 5 - 23;
% station 1 
stname = 'x = 4 m';
eta8_gn = eta_stat(1:num_val,8);
eta8_s = 30 *eta_stat_swe(1:145,8);
plot(time_swe,eta8_s,time,eta8_gn,time_plot_exp2,St_6);
%shoal_plot(time_swe,eta4,time_plot_exp2,St_1,stname);

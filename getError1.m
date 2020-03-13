function error = getError1(time,F_total,i,Data)

% forcename = {'force_05F.mat','force_1F.mat','force_15F.mat','force_2F.mat','force_25F.mat','force_3F.mat','force_085L.mat','force_090L.mat','force_095L.mat'};
% 
% 
% load(forcename{i});
if i == 1 
      F_exp = Data.F_exp1;
elseif i == 2
      F_exp = Data.F_exp2;

elseif i == 3
      F_exp = Data.F_exp3;

elseif i == 4
      F_exp = Data.F_exp4;

elseif i == 5
      F_exp = Data.F_exp5;

elseif i == 6
      F_exp = Data.F_exp6;
elseif i == 7
      F_exp = Data.F_exp7;

elseif i == 8
      F_exp = Data.F_exp8;      
end

time_exp_data = 0 : 1 : length(F_exp)-1;
time_exp_end = time_exp_data(end);
time_sim_end = time(end);
min_time = min (time_exp_end,time_sim_end);


F_exp_spline = spline(time_exp_data,F_exp,[0:1:min_time]);
F_sim_spline = spline(time,F_total,[0:1:min_time]);

error = 100 * sqrt(sumsqr((F_exp_spline - F_sim_spline)))/(min_time);

    

function sum_error =  objective_fun_rat_data(var,Data,flag_plot)

% %% Data from Janssen at 37.5 degree
    delay_TTP = 4.0; 
    Fmax_data = [8.46 16.3 33 54.2];
    TTP_data  = [38.8 40.5 45.4 54.2] - delay_TTP; 
    RT50_data = [19.2 22.7 25.7 29.8];

    Tdev_freq = [37.2 46.1 55.8 62.3 71.1 73.4 73.9 72.9 66.8];
   eTdev_freq = [6.1  7.5  7.5  7.8  8.2  8.5  8.4  8.3  9.5];
    TTP_freq  = [55.1 52.0 49.3 47.2 45.7 42.8 40.8 40.1 38.1] - delay_TTP; 
    RT50_freq = [35.7 33.2 31.0 28.9 28.2 27.3 27.1 25.4 25.2];
   eRT50_freq = [2.6  2.2  1.9  1.5  1.9  1.3  1.3  1.3  1.2];

%% Data from XB parameters used for the mean sham rat

%     Fmax_data = [31.97 53.06 73.53 92.25];
%     TTP_data  = [44 43 43 43] ; 
%     RT50_data = [25.36 26.42 27.47 28.89];
% 
%     Tdev_freq = [54.49 76.43 92.25 103.9 103.6 107.8 106.1 102.5 97.55];
%    eTdev_freq = [6.1  7.5  7.5  7.8  8.2  8.5  8.4  8.3  9.5];
%     TTP_freq  = [53 48 43 41 40 36 35 34 32] ; 
%     RT50_freq = [35.95 31.35 28.85 25.9 25.48 24.54 22.52 21.64 20.9];
%    eRT50_freq = [2.6  2.2  1.9  1.5  1.9  1.3  1.3  1.3  1.2];

lstyle = {'-g','-b','-r','-k','-c','-y'}; k = 0;

freq = Data.freq;
% Set metabolite concentrations,
% 
% MgATP = 8.0494; % mM
% MgADP = 17.7e-3; % mM
% Pi = 0.59287; % mM

%mean SHAM rat
MgATP = 7.93; % mM
MgADP = 43.17e-3; % mM
Pi = 1; % mM

% TAC rat #1
% MgATP = 4.9551; % mM
% MgADP = 20.9e-3; % mM
% Pi = 2.003; % mM

% Set sarcomere lengths, Units: um
SL0 = 2.2;
% 
% kstiff1 = 2000*var(1); % unit (kPa/um) % manually tuned by the factor of * 0.6
% kstiff2 = 25000*var(1);  % unit (kPa/um) % manually tuned by the factor of * 0.6
kstiff1 = var(1)*20000/7.5; % kPa/um (9/5 BM)
kstiff2 = var(1)*375000/7.5; % kPa/um (9/5 BM)
% kstiff1 = var(1)*5.1535e+03; %kPa/um (Tewari et al)
% kstiff2 = var(1)*1.0974e+05;% kPa/um (Tewari et al)
% 
dr = 0.01; % Power-stroke Size; Units: um
k_passive = var(2)* 10; % mN / mm^2 / micron
Kse = 1400;
SLset = 2.2;
% loading the Ca data file, the variable Freq should also change according to the freq Ca.
for i = 1:length(freq)
%         load(filename{i}); % unit (mM)
%  eval(['Ca'  '= parameters.Ca' num2str(i)]);
%  eval(['T'  '= parameters.T' num2str(i)]);

if i == 10
    SL0 = 1.9;
    SLset = 1.9; 
elseif i == 11
  
    SL0 = 2.0;
    SLset = 2.0; 

elseif i == 12

    SL0 = 2.1;
    SLset = 2.1; 

end
%% Extract Ca coeficient based on the Ca data for diferent simulation 
para_fitted_Ca = [2	3	4	5	6	7	8	9	10;
0.0838	0.1306	0.1802	0.2557	0.3099	0.3613	0.408	0.4539	0.4602;
0.7513	0.8756	1.0274	1.4988	1.6107	1.6741	1.7902	2.1398	1.9832;
2.237	2.0486	1.948	1.852	1.6932	1.6773	1.5988	1.4952	1.4524;
0.1865	0.1815	0.1709	0.1693	0.161	0.1661	0.1425	0.1458	0.1222];
freq_all = para_fitted_Ca(1,:);
A_HR_pchip = pchip(freq_all,para_fitted_Ca(2,:));
A_HR = ppval(A_HR_pchip,freq(i));
B_HR_pchip = pchip(freq_all,para_fitted_Ca(3,:));
B_HR = ppval(B_HR_pchip,freq(i));
C_HR_pchip = pchip(freq_all,para_fitted_Ca(4,:));
C_HR = ppval(C_HR_pchip,freq(i));
Ca0_HR_pchip = pchip(freq_all,para_fitted_Ca(5,:));
Ca0_HR = ppval(Ca0_HR_pchip,freq(i));

stim_f = 1/freq(i);
tspan = 0:0.001:stim_f;

para = [Data.TmpC, MgATP, MgADP, Pi, freq(i), kstiff1, kstiff2, k_passive SLset A_HR B_HR C_HR Ca0_HR Kse];

  init = [zeros(1,10),SL0,0.2]; % Initial conditions for the model
  init(10) = 1;% setting the initial value for nonpermissible state equal to 1
  % Solving the system of diffrential equations
 
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-1);
    [~,ys] = ode15s(@Model_XB_Ca_activation,0:0.001:stim_f*5,init,options,para,var );
    init = ys(end,:);
    [t, Y] = ode15s(@Model_XB_Ca_activation,tspan,init,options,para,var );
       p1_0 = Y(:,1);
       p2_0 = Y(:,4);
       p2_1 = Y(:,5);
       p3_0 = Y(:,7);
       p3_1 = Y(:,8);
       SL = Y(:,11);
       N = Y(:,10);
       P = 1 - N - p3_0 - p2_0 - p1_0;
       UNR = Y(:,12);

% plot(
% Ca_i = spline(T'/max(T),Ca',t)
F_SE = Kse*(SLset - SL);
[Tdev_mod(i,1),TTP_mod(i,1),RT50_mod(i,1)] = get_TTP_RT50(SL,t,Kse,SLset);

% Tdev_mod(i,1)= Tdev_mod(i,1)/1.51;

Lthin = 1200; % nm
Lthick = 1670; % nm
Lbare = 100; % nm Thick filament bare zone
OV_Zaxis = min(Lthick/2,1000*SL/2);
OV_Mline = max(1000*SL/2-(1000*SL-Lthin),Lbare/2);
LOV = OV_Zaxis - OV_Mline;
% N_overlap = LOV/Lthin;
N_overlap_thick = LOV*2/(Lthick - Lbare);

B_process = kstiff2 * dr * p3_0;   % Force due to XB ratcheting
C_process = kstiff1 * ( p2_1 + p3_1 );% Force due to stretching of XBs
F_XB = N_overlap_thick.*( B_process + C_process ); % Active Force

Lsref = 1.9;

% Collagen force
SLcollagen = 2.25; % threshold for collagen activation, microns
PConcollagen = 0.01*7.5; % contriubtion of collagen (??)
PExpcollagen = 70; % expresion of collagen (??), unitless
sigma_collagen  = PConcollagen*(exp(PExpcollagen*(SL - SLcollagen)) - 1).*(SL > SLcollagen);

F_passive  = k_passive*(SL/2-Lsref/2) + sigma_collagen;

% Total passive forces
phi = mod(t+0.0001,stim_f)/stim_f;
Ca = (A_HR./phi).*exp(-B_HR.*(log(phi)+C_HR).^2) + Ca0_HR;

if (flag_plot == 1) && (i == 1 || i == 3 || i == 6 || i == 9 || i == 10 || i == 11)
k = k + 1;

figure(2)
hold on
plot(t, F_XB + F_passive, lstyle{k},'linewidth',2)
legend('2Hz','4Hz','7Hz','10Hz','SL = 1.9','SL = 2');
xlim([0 .250])
figure(3)
hold on
title('SL')
plot(t,SL)
legend('2Hz','4Hz','7Hz','10Hz','SL = 1.9','SL = 2');
figure(222)
hold on
plot(t, F_SE, lstyle{k},'linewidth',2)
legend('2Hz','4Hz','7Hz','10Hz','SL = 1.9','SL = 2');
xlim([0 .250])
figure(3)
hold on
title('SL')
plot(t,SL)
legend('2Hz','4Hz','7Hz','10Hz','SL = 1.9','SL = 2');

figure(44)
hold on
plot(t,N)
plot(t,P)

legend('N','P')

figure(55)
hold on
plot(t,UNR)
legend('2Hz','4Hz','7Hz','10Hz','SL = 1.9','SL = 2');

title('UNR')
ylim([0 1]);
figure(66)
hold on
plot(t,p2_1)
title('p2_1')

ylim([0 1]);
figure(77)
hold on
plot(t,p3_0)
title('p3_0')

ylim([0 1]);
figure(88)
hold on
plot(t,p3_1)
title('p3_1')
ylim([0 1]);
figure(99)
plot(phi,Ca)
end

    
clear Ca T Ca_i

end
%considering the lag in Ca data 
lag_Force = 7;

TTP_mod = TTP_mod*1000 - lag_Force;
RT50_mod = RT50_mod*1000;
% lag_Ca = 15;
% TTP_mod = TTP_mod - lag_Ca;
if flag_plot == 1
    h = figure(4); 
%     set(h,'Units','inches','Position',[0.5 0.5 11 6]);
    subplot(231),
    plot([1.9 2.0 2.1 2.2],RT50_data,'ok','linewidth',2,'Markersize',8,'Markerfacecolor','k'); hold on; 
    plot([1.9 2.0 2.1 2.2],[RT50_mod(10,1) RT50_mod(11,1) RT50_mod(12,1) RT50_mod(3,1)],'-ko','linewidth',2,'Markersize',6,'Markerfacecolor','w')
    l = legend('Data','Model');
    set(l,'interpreter','latex','Position', [.16,.73,.15,.45],'Box','off','fontsize',18,'Orientation','horizontal');
    text(1.65,40,'A','interpreter','latex','fontsize',18)
    xlim([1.8 2.3]); ylim([0 40])
%     set(l,'fontsize',18,'interpreter','latex','location','northwest')
    ylabel('RT$_{50}$ (msec)','fontsize',18,'interpreter','latex');
    xlabel('SL ($\mu$m)','fontsize',18,'interpreter','latex');
    set(gca, 'LineWidth',1.5, 'FontSize',14);
    
%     figure(4); clf; 
    subplot(232),
    plot([1.9 2.0 2.1 2.2],Fmax_data,'ok','linewidth',2,'Markersize',8,'Markerfacecolor','k'); hold on; 
    plot([1.9 2.0 2.1 2.2],[Tdev_mod(10,1) Tdev_mod(11,1) Tdev_mod(12,1) Tdev_mod(3,1)],'-ko','linewidth',2,'Markersize',6,'Markerfacecolor','w')
    xlim([1.8 2.3]); ylim([0 80])
    text(1.65,80,'B','interpreter','latex','fontsize',18)
    l = legend('Data','Model');
    set(l,'interpreter','latex','Position', [.44,.73,.15,.45],'Box','off','fontsize',18,'Orientation','horizontal');
    ylabel('T$_{dev}$ (kPa)','fontsize',18,'interpreter','latex');
    xlabel('SL ($\mu$m)','fontsize',18,'interpreter','latex');
    set(gca, 'LineWidth',1.5, 'FontSize',14);
    
%     figure(5); clf;
    subplot(233),
    plot([1.9 2.0 2.1 2.2],TTP_data,'ok','linewidth',2,'Markersize',8,'Markerfacecolor','k'); hold on; 
    plot([1.9 2.0 2.1 2.2],[TTP_mod(10,1) TTP_mod(11,1) TTP_mod(12,1) TTP_mod(3,1)],'-ko','linewidth',2,'Markersize',6,'Markerfacecolor','w')
    xlim([1.8 2.3]); ylim([0 60])
    l = legend('Data','Model');
    set(l,'interpreter','latex','Position', [.72,.73,.15,.45],'Box','off','fontsize',18,'Orientation','horizontal');
    text(1.65,60,'C','interpreter','latex','fontsize',18)
    ylabel('TTP (msec)','fontsize',18,'interpreter','latex');
    xlabel('SL ($\mu$m)','fontsize',18,'interpreter','latex');
    set(gca, 'LineWidth',1.5, 'FontSize',14);

%     figure(6); clf; 
    subplot(234),
    errorbar([2 3 4 5 6 7 8 9 10],RT50_freq,eRT50_freq,'ok','linewidth',2,'Markersize',8,'MarkerFaceColor','k'); hold on; 
    plot([2 3 4 5 6 7 8 9 10],[RT50_mod(1,1) RT50_mod(2,1) RT50_mod(3,1) RT50_mod(4,1) RT50_mod(5,1) RT50_mod(6,1) RT50_mod(7,1) RT50_mod(8,1) RT50_mod(9,1)],'-ko','linewidth',2,'Markersize',6,'Markerfacecolor','w')
    ylabel('RT$_{50}$ (msec)','fontsize',18,'interpreter','latex'); xlim([1.5 10.5])
    xlabel('Freq (Hz)','fontsize',18,'interpreter','latex'), ylim([0 40])
%     l = legend('Data','Model');
%     set(l,'interpreter','latex','Position', [.67,.72,.15,.45],'Box','off','fontsize',18,'Orientation','horizontal');
    text(-1.2,40,'D','interpreter','latex','fontsize',18)

%     set(l,'interpreter','latex','fontsize',18,'location','southeast')
    set(gca, 'LineWidth',1.5, 'FontSize',14);
    
%     figure(7); clf; 
    subplot(235),
    errorbar([2 3 4 5 6 7 8 9 10],Tdev_freq,eTdev_freq,'ok','linewidth',2,'Markersize',8,'MarkerFaceColor','k'); hold on; 
    plot([2 3 4 5 6 7 8 9 10],[Tdev_mod(1,1) Tdev_mod(2,1) Tdev_mod(3,1) Tdev_mod(4,1) Tdev_mod(5,1) Tdev_mod(6,1) Tdev_mod(7,1) Tdev_mod(8,1) Tdev_mod(9,1)],'-ko','linewidth',2,'Markersize',6,'Markerfacecolor','w')
    ylabel('T$_{dev}$ (kPa)','fontsize',18,'interpreter','latex'); xlim([1.5 10.5])
    xlabel('Freq (Hz)','fontsize',18,'interpreter','latex'), ylim([0 100])
%     l = legend('Data','Model');
%     set(l,'interpreter','latex','Position', [.67,.72,.15,.45],'Box','off','fontsize',18,'Orientation','horizontal');
    text(-1.2,100,'E','interpreter','latex','fontsize',18)
    set(gca, 'LineWidth',1.5, 'FontSize',14);
    
%     figure(8); clf; 
    subplot(236),
    plot([2 3 4 5 6 7 8 9 10],TTP_freq,'ok','linewidth',2,'Markersize',8,'MarkerFaceColor','k'); hold on; 
    plot([2 3 4 5 6 7 8 9 10],[TTP_mod(1,1) TTP_mod(2,1) TTP_mod(3,1) TTP_mod(4,1) TTP_mod(5,1) TTP_mod(6,1) TTP_mod(7,1) TTP_mod(8,1) TTP_mod(9,1)],'-ko','linewidth',2,'Markersize',6,'Markerfacecolor','w')
    ylabel('TTP (msec)','fontsize',18,'interpreter','latex');
    xlabel('Freq (Hz)','fontsize',18,'interpreter','latex'), xlim([1.5 10.5]),ylim([0 60])
%     l = legend('Data','Model');
%     set(l,'interpreter','latex','Position', [.67,.72,.15,.45],'Box','off','fontsize',18,'Orientation','horizontal');
    text(-1.2,60,'F','interpreter','latex','fontsize',18)
    set(gca, 'LineWidth',1.5, 'FontSize',14);
end
% calculate the fitting aeror error 
   error(1) = sum((RT50_data - [RT50_mod(10,1) RT50_mod(11,1) RT50_mod(12,1) RT50_mod(3,1)]).^2);
   error(2) =  sum((Fmax_data - [Tdev_mod(10,1) Tdev_mod(11,1) Tdev_mod(12,1) Tdev_mod(3,1)]).^2);
   error(3) =  sum((TTP_data - [TTP_mod(10,1) TTP_mod(11,1) TTP_mod(12,1) TTP_mod(3,1)]).^2);
   error(4) = sum((RT50_freq - [RT50_mod(1,1) RT50_mod(2,1) RT50_mod(3,1) RT50_mod(4,1) RT50_mod(5,1) RT50_mod(6,1) RT50_mod(7,1) RT50_mod(8,1) RT50_mod(9,1)]).^2);
   error(5) = sum((Tdev_freq - [Tdev_mod(1,1) Tdev_mod(2,1) Tdev_mod(3,1) Tdev_mod(4,1) Tdev_mod(5,1) Tdev_mod(6,1) Tdev_mod(7,1) Tdev_mod(8,1) Tdev_mod(9,1)]).^2);
   error(6) = sum((TTP_freq - [TTP_mod(1,1) TTP_mod(2,1) TTP_mod(3,1) TTP_mod(4,1) TTP_mod(5,1) TTP_mod(6,1) TTP_mod(7,1) TTP_mod(8,1) TTP_mod(9,1)]).^2);
   sqrt_error = sqrt(error);
%    sum_error = -sum(sqrt(error)) ;
   sum_error = -sqrt(error(2)+ error(4) + error(5) + error(6));
   Max_Force= getMaxForce (var); 
   if Max_Force >160
            sum_error = sum_error - 5000;
   end
   if sum_error > - 60
        [hill_19 Ca50_19 hill_23 Ca50_23]= getpCa50 (var);  
        error_Ca50 = Ca50_23 - Ca50_19;
       if (hill_19 > 7) || (hill_23 > 7) || (error_Ca50 <0.09) || (hill_19 < 4.0) || (hill_23 < 4.0) 
            sum_error = sum_error - 5000;
       end    
   end

end


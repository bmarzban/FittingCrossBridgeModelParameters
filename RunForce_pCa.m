% RunForcepCa
clear; close all;
var = [ 1.7386*1.7    1.5022    1.4040    1.8637    1.4793    1.4693*1.05    0.6281*0.2    1.9245    0.5000*0.2   2.3293*1.35]

% set the axis font size for the figure plots
flag_plot_state_var = 1; % if "para.flag_plot = 1" plots every state variable else does not plot.
flag_plot_force = 1;
flag_plot_fig3 = 1;
flag_plot_fig3_USR_UNR = 1;

AxisFontSize = 18; LabelFontSize = 18;

%% Defining the parameters into the parameter (para)  

% Set temperature and initial SL (Sarcomere length)
% From Tewari etal JMCC paper one

TmpC =37.5; % centigrade 

lstyle = {'-g','-b','-r','-k'}; k = 0;

% filename = {'Ca_05.mat','Ca_1.mat','Ca_15.mat','Ca_2.mat','Ca_25.mat','Ca_3.mat'};

% Set metabolite concentrations,


%mean SHAM rat
MgATP = 7.93; % mM
MgADP = 43.17e-3; % mM
Pi = 1; % mM

kstiff1 = var(1)*5.1535e+03; %kPa/um (Tewari et al)
kstiff2 = var(1)*1.0974e+05; % kPa/um (Tewari et al)
% kstiff1 = var(1)* 2827.1 * Q10s(4)^((TmpC-17)/10) % kPa/um (9/5 BM)
% kstiff2 = var(1)* 51871 * Q10s(6)^((TmpC-17)/10) % kPa/um (9/5 BM)
dr = 0.01; % Power-stroke Size; Units: um
k_passive = var(2)* 40/7.5; % mN / mm^2 / micron
Kse = 1400;

SLset = [1.90,2.1,2.3];

% Defining the time vector
tspan = 0:0.0001:0.3; nn = length(tspan);

%%Defining the Ca range for the study 
%  Ca_frction = 1; %Unit uM
Ca_fraction = [0.005:0.01:1.5,1.6:0.1:20,30:5:500]; %Unit uM
% Ca_fraction = [0.2:0.1:50]; %Unit uM

t_T = tspan;

tic;
SL0 = [1.90,2.3]; % Set sarcomere length, Units: um

for j = 1: length (SL0)
    
    init = [zeros(1,10),SL0(j),0.2]; % Initial conditions for the model
    init(10) = 1;% setting the initial value for nonpermissible state equal to 1

    for k = 1: length (Ca_fraction)
    para = [TmpC, MgATP, MgADP, Pi, 1, kstiff1, kstiff2, k_passive SLset(j) 1 1 1 1 Kse];
    % Solving the system of diffrential equations
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
    [Time, Y] = ode15s(@Model_XB_Ca_Sensitivity,tspan,init,options,para,var,Ca_fraction(k));
    p2_1 = Y(:,5);
    p3_0 = Y(:,7);
    p3_1 = Y(:,8);
    N = Y(:,10);
    SL = Y(:,11);
    UNR = Y(:,12);
%% Overlap function (Tewari et al.)
    Lthin = 1200; % nm
    Lthick = 1670; % nm
    Lbare = 100; % nm Thick filament bare zone
    OV_Zaxis = min(Lthick/2,1000*SL/2);
    OV_Mline = max(1000*SL/2-(1000*SL-Lthin),Lbare/2);
    LOV = OV_Zaxis - OV_Mline;
    % N_overlap = LOV/Lthin;
    N_overlap_thick = LOV*2/(Lthick - Lbare);

% Active Force
    dr = 0.01; % Power-stroke Size; Units: um
    B_process = kstiff2 * dr * p3_0;   % Force due to XB ratcheting
    C_process = kstiff1 * ( p2_1 + p3_1 );% Force due to stretching of XBs
    F_XB = N_overlap_thick.*( B_process + C_process ); % Active Force

% (linear) passive force model

%% Titin force
    Lsref = 1.9;

% Collagen force
    SLcollagen = 2.25; % threshold for collagen activation, microns
    PConcollagen = 0.01*7.5; % contriubtion of collagen (??)
    PExpcollagen = 70; % expresion of collagen (??), unitless
    sigma_collagen  = PConcollagen*(exp(PExpcollagen*(SL - SLcollagen)) - 1).*(SL > SLcollagen);



    F_passive  = k_passive*(SL/2-Lsref/2)   + sigma_collagen;

    Ftotal = F_XB + F_passive;
     
             SS_N_on (k,j)= 1 - N(end);
             SS_USR (k,j) =  1 - UNR(end);
             SS_UNR (k,j) =  UNR(end);

             SS_F_XB (k,j)= F_XB(end);
             SS_Ftotal (k,j)= Ftotal(end);

%             Max_F_XB (k,j)= max(F_XB);

    end
toc;

% plots the state variables if the flag_plot == 1 
      
    if flag_plot_state_var ==1 
            
        figure(2)
        plot(t_T  ,p3_0/max(p3_0))
        title('Normailized P30')
        hold on
%         figure(3)
%         plot(t_T  ,P)
%         title('Pu')
%         hold on
        figure(4)
        plot(t_T  ,p3_0)
        title('P30')
        hold on
%         figure(5)
%         plot(t_T  ,p1_0)
%         title('P1o')
%         hold on
%         figure(6)
%         plot(t_T  ,p2_0)
%         title('P2o')
%         hold on
        figure(7)
        plot(t_T  ,N)
        title('Np')
        hold on
        figure(8)
        plot(t_T  ,UNR)
        title('UNR')
        hold on
        figure(9)
        plot(t_T  ,SL)
        title('SL')
        hold on
    end 
end

if flag_plot_force == 1

        figure(10)
        plot(-log10(Ca_fraction*(10^(-6))),SS_F_XB(:,1),'-ro','linewidth',1.5)

        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_F_XB(:,2),'linewidth',1.5)
        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_F_XB(:,3),'linewidth',1.5)

        set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);
         title(sprintf('F active'),'fontsize',14); 
         set ( gca, 'xdir', 'reverse' )

         xlabel('pCa = - Log[Ca]','fontsize',LabelFontSize);
         ylabel('Force Active(mN mm^{-2})','fontsize',LabelFontSize);
        legend('SL = 1.9 um','SL = 2.1 um','SL = 2.3 um');
 figure(120)
        plot(-log10(Ca_fraction*(10^(-6))),SS_Ftotal(:,1),'-ro','linewidth',1.5)

        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_Ftotal(:,2),'linewidth',1.5)
        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_Ftotal(:,3),'linewidth',1.5)

        set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);
         title(sprintf('F total'),'fontsize',14); 
        set ( gca, 'xdir', 'reverse' )

         xlabel('pCa = - Log[Ca]','fontsize',LabelFontSize);
         ylabel('Force Total (kPa)','fontsize',LabelFontSize);
        legend('SL = 1.9 um','SL = 2.1 um','SL = 2.3 um');
end


if flag_plot_fig3 == 1


        figure(12)
        plot(-log10(Ca_fraction*(10^(-6))),SS_N_on(:,1),'-ro','linewidth',1.5)

        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_N_on(:,2),'linewidth',1.5)
        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_N_on(:,3),'linewidth',1.5)

        set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);
         title(sprintf('N on'),'fontsize',14); 
         set ( gca, 'xdir', 'reverse' )

         xlabel('pCa = - Log[Ca]','fontsize',LabelFontSize);
         ylabel('N_on','fontsize',LabelFontSize);
        legend('SL = 1.9 um','SL = 2.1 um','SL = 2.3 um');
end


if flag_plot_fig3_USR_UNR == 1

       
        figure(16)
        plot(-log10(Ca_fraction*(10^(-6))),SS_USR(:,1),'-ro','linewidth',1.5)

        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_USR(:,2),'linewidth',1.5)
        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_USR(:,3),'linewidth',1.5)

        set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);
         title(sprintf('USR'),'fontsize',14); 
         set ( gca, 'xdir', 'reverse' )

         xlabel('pCa = - Log[Ca]','fontsize',LabelFontSize);
         ylabel('USR','fontsize',LabelFontSize);
        legend('SL = 1.9 um','SL = 2.1 um','SL = 2.3 um');
        
       
        figure(18)
        plot(-log10(Ca_fraction*(10^(-6))),SS_UNR(:,1),'-ro','linewidth',1.5)

        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_UNR(:,2),'linewidth',1.5)
        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_UNR(:,3),'linewidth',1.5)

        set(gca, 'LineWidth',1.5, 'FontSize',AxisFontSize);
         title(sprintf('UNR'),'fontsize',14); 
         set ( gca, 'xdir', 'reverse' )

         xlabel('pCa = - Log[Ca]','fontsize',LabelFontSize);
         ylabel('UNR','fontsize',LabelFontSize);
        legend('SL = 1.9 um','SL = 2.1 um','SL = 2.3 um');
end

[hill_19 ec50_19]=doseResponse( (Ca_fraction')*10^(-6) ,SS_Ftotal(:,1))
Ca50_19=-log10(ec50_19)

[hill_23 ec50_23]=doseResponse( (Ca_fraction')*10^(-6),SS_Ftotal(:,2))
Ca50_23=-log10(ec50_23)

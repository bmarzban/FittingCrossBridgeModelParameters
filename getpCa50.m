% RunForcepCa
function [hill_19 Ca50_19 hill_23 Ca50_23]= getpCa50 (var) 

% var = [ 1.7386*1.7    1.5022    1.4040    1.8637    1.4793    1.4693*1.05    0.6281*0.2    1.9245    0.5000*0.2   2.3293*1.35]
TmpC =37.5; % centigrade 
% Set metabolite concentrations,
%mean SHAM rat
MgATP = 7.93; % mM
MgADP = 43.17e-3; % mM
Pi = 1; % mM

% MgATP = 8.0494; % mM
% MgADP = 17.7e-3; % mM
% Pi = 0.59287; % mM

% kstiff1 = var(1)*5.1535e+03; %kPa/um (Tewari et al)
% kstiff2 = var(1)*1.0974e+05; % kPa/um (Tewari et al)

kstiff1 = var(1)*20000/7.5; % kPa/um (9/5 BM)
kstiff2 = var(1)*375000/7.5; % kPa/um (9/5 BM)
dr = 0.01; % Power-stroke Size; Units: um
k_passive = var(2)* 40/7.5; % mN / mm^2 / micron
Kse = 1400;
SLset = [1.90,2.3];
% Defining the time vector
tspan = 0:0.0001:0.3; nn = length(tspan);

%%Defining the Ca range for the study 
% Ca_fraction = [0.005:0.01:0.99,0.1:0.01:1.5,1.6:0.1:20,30:5:1000]; %Unit uM
% Ca_fraction = [0.0005:0.005:0.05,0.06:0.01:3.5,3.6:0.1:20,30:20:1000]; %Unit uM
Ca_fraction = [0.2:0.1:20]; %Unit uM

% Ca_fraction = [0.2:0.1:50]; %Unit uM

SL0 = [1.90,2.3]; % Set sarcomere length, Units: um

for j = 1: length (SL0)
    init = [zeros(1,10),SL0(j),0.2]; % Initial conditions for the model
    init(10) = 1;% setting the initial value for nonpermissible state equal to 1

    for k = 1: length (Ca_fraction)
    para = [TmpC, MgATP, MgADP, Pi, 1, kstiff1, kstiff2, k_passive SLset(j) 1 1 1 1 Kse];
    % Solving the system of diffrential equations
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
    [t, Y] = ode15s(@Model_XB_Ca_Sensitivity,tspan,init,options,para,var,Ca_fraction(k));
    p2_1 = Y(:,5);
    p3_0 = Y(:,7);
    p3_1 = Y(:,8);
    SL = Y(:,11);
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
     
%              SS_N_on (k,j)= 1 - N(end);
%              SS_USR (k,j) =  1 - UNR(end);
%              SS_UNR (k,j) =  UNR(end);
%              SS_F_XB (k,j)= F_XB(end);
    SS_Ftotal (k,j)= Ftotal(end);
%             Max_F_XB (k,j)= max(F_XB);
    end
end

 

[hill_19 ec50_19]=pCa_calculate( (Ca_fraction')*10^(-6),SS_Ftotal(:,1));
Ca50_19=-log10(ec50_19);

[hill_23 ec50_23]=pCa_calculate( (Ca_fraction')*10^(-6),SS_Ftotal(:,2));

Ca50_23=-log10(ec50_23);

figure(1)
        plot(-log10(Ca_fraction*(10^(-6))),SS_Ftotal(:,1),'ro','linewidth',1.5)
        hold on 
        plot(-log10(Ca_fraction*(10^(-6))),SS_Ftotal(:,2),'bo','linewidth',1.5)
   
        title(sprintf('F total'),'fontsize',14); 
        set ( gca, 'xdir', 'reverse' )

         xlabel('pCa = - Log[Ca]');
         ylabel('Force Total (kPa)');
A_19 = min(SS_Ftotal(:,1));
D_19 = max(SS_Ftotal(:,1));
A_23 = min(SS_Ftotal(:,2));
D_23 = max(SS_Ftotal(:,2));
x = Ca_fraction*(10^(-6));
cf_19 = D_19+(A_19-D_19)./(1+(x/ec50_19).^hill_19);
cf_23 = D_23+(A_23-D_23)./(1+(x/ec50_23).^hill_23);
% figure(2)
% hold on
% semilogx(x,cf_19)
% semilogx(x,cf_23)
plot(-log10(Ca_fraction*(10^(-6))),cf_19)
plot(-log10(Ca_fraction*(10^(-6))),cf_23)
legend('SL = 1.9 um','SL = 2.3 um','Fitted cureve-SL = 1.9 um','Fitted Curve- SL = 2.3 um');




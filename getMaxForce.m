% RunForcepCa
function [Max_Force]= getMaxForce (var) 
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
SLset = 2.3;
% Defining the time vector
tspan = 0:0.0001:0.3; nn = length(tspan);

%% the Ca range for the study 
Ca_fraction = 20; %Unit uM
% Ca_fraction = [0.2:0.1:50]; %Unit uM
SL0 = 2.3; % Set sarcomere length, Units: um
    init = [zeros(1,10),SL0,0.2]; % Initial conditions for the model
    init(10) = 1;% setting the initial value for nonpermissible state equal to 1
    para = [TmpC, MgATP, MgADP, Pi, 1, kstiff1, kstiff2, k_passive SLset 1 1 1 1 Kse];
    % Solving the system of diffrential equations
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
    [t, Y] = ode15s(@Model_XB_Ca_Sensitivity,tspan,init,options,para,var,Ca_fraction);
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
    Max_Force = max(Ftotal);
 
end

 




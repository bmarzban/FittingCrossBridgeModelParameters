function [ dYdT, dSL, F_XB, F_passive ] = Model_XB_Ca_Sensitivity( t,y,para, var,Ca_fraction)
MgATP = para(2);
MgADP = para(3);
Pi = para(4);
kstiff1 = para (6);
kstiff2 = para (7);
k_passive = para(8);
SLset = para(9);
% fiting the calcium data file
Kse = para(14);
Ca_i = Ca_fraction;
%% Constants and parameters for the cross bridge model
% All the parameters here are from table 1 of  Tewari et al, 2015(Dynamic cross bridge...) 
% for the rat

K_Pi = 4.00; % unit (mM) K_Pi [Pi] dissociation constant in the paper 
K_T = 0.4897 ; % unit (mM) %K_MgATP [MgATP] dissociation constant in the paper 
K_D = 0.194;  % unit (mM) MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).

% Defining the metabolite dependent coeficient, rapid equilibrium of the
% cross bridge sub-states
g1 = (MgADP/K_D)/(1 + MgADP/K_D + MgATP/K_T);
g2 = (MgATP/K_T)/(1 + MgATP/K_T + MgADP/K_D);
f1 = (Pi/K_Pi)/(1 + Pi/K_Pi);
f2 = 1/(1 + Pi/K_Pi);

K_coop = var(3)*6; % Campbell et al, Biophysical Journal 2018,
k_on = var(4)*50; % Campbell et al, Biophysical Journal 2018
k_off =var(5)*180; % manually tuned parameter!

% transitions between super relaxed state and non relaxed state
k1_SR = var(6)*23.1;
kforce = var(7)*1.2825; 
k2_SR= var(8)*150; % 

% XB parameters
ka      = var(9)*608.28; % myosin-actin attach rate constant, 1/sec
kd      = var(10)*183.88;% myosin-actin detach rate constant, 1/sec
k1      = var(11)*21.19; % transition A1 to A2 rate constant, 1/sec
k_1     = var(12)*21.296; % transition A2 to A1 rate constant, 1/sec
k2      = var(13)*811.72; % transition A2 to A3 rate constant, 1/sec
k_2     = var(14)*43.25;% transition A3 to A2 rate constant, 1/sec
k3      = var(15)*73.698; % transition A3 to P rate constant, 1/sec
alpha1  = 10.0; % Stretch sensing parameter for k1 and k?1, 1/um
alpha2  = 9.1; % Stretch sensing parameter for k2 and k?2, 1/um
alpha3  = 0.1*59.3; % Stretch sensing parameter for k3, 1/um
s3      = 9.9e-3;  % Strain at which k3 is minimum, um

kd = kd * f1; % unit (1/s)
k1 = k1* f2; % % unit (1/s)
k_2 = k_2 * g1;  % unit (1/s)
k3 = k3 * g2;% % unit (1/s)

% visc = 2.4*Q10s(3)^((TmpC-17)/10);
visc =0.1;

%% State Variables,
% Moment of state variables
p1_0 = y(1);
p1_1 = y(2);
p1_2 = y(3);
p2_0 = y(4);
p2_1 = y(5);
p2_2 = y(6);
p3_0 = y(7);
p3_1 = y(8);
p3_2 = y(9);
% N represents the Non permissible state in cross bridge model
N = y(10);
% P represents the permssible state in cross bridge model
P = 1.0 - N - p1_0 - p2_0 - p3_0;
% sarcomere length
SL = y(11);
% U_NR represents the Non relaxed state
U_NR = y(12);
% U_SR is the super relaxed state
U_SR = 1.0 - U_NR;

%% Stretch-sensitive rates
f_alpha1o = (p1_0 - alpha1*p1_1 + 0.5*(alpha1*alpha1)*p1_2);
f_alpha1i = (p1_1 - alpha1*p1_2);

f_alpha0o = (p2_0 + alpha1*p2_1 + 0.5*alpha1*alpha1*p2_2);
f_alpha0i = (p2_1 + alpha1*p2_2);

f_alpha2o = (p2_0 - alpha2*p2_1 + 0.5*(alpha2*alpha2)*p2_2);
f_alpha2i = (p2_1 - alpha2*p2_2);

f_alpha3o = (p3_0 + alpha3*(s3*s3*p3_0 + 2.0*s3*p3_1 + p3_2)); 
f_alpha3i = (p3_1 + alpha3*(s3*s3*p3_1 + 2.0*s3*p3_2));


%% Compute Active & Passive Force
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
F_XB = N_overlap_thick*( B_process + C_process ); % Active Force

% (linear) passive force model

%% Titin force
Lsref = 1.9;

% Collagen force
SLcollagen = 2.25; % threshold for collagen activation, microns
PConcollagen = 0.01*7.5; % contriubtion of collagen (??)
PExpcollagen = 70; % expresion of collagen (??), unitless
sigma_collagen  = PConcollagen*(exp(PExpcollagen*(SL - SLcollagen)) - 1)*(SL > SLcollagen);

F_passive  = k_passive*(SL/2-Lsref/2)   + sigma_collagen;

Ftotal = F_XB + F_passive;
intf = (- Ftotal  + Kse*(SLset - SL));

dSL = intf/visc;

dp1_0 = ka*P*N_overlap_thick*U_NR   - kd*p1_0 - k1*f_alpha1o + k_1*f_alpha0o;
dp1_1 = 1*dSL*p1_0 - kd*p1_1 - k1*f_alpha1i + k_1*f_alpha0i;
dp1_2 = 2*dSL*p1_1 - kd*p1_2 - k1*p1_2 + k_1*p2_2;
dp2_0 =         - k_1*f_alpha0o - k2*f_alpha2o + k_2*p3_0 + k1*f_alpha1o;
dp2_1 = 1*dSL*p2_0 - k_1*f_alpha0i - k2*f_alpha2i + k_2*p3_1 + k1*f_alpha1i;
dp2_2 = 2*dSL*p2_1 - k_1*p2_2       - k2*p2_2 + k_2*p3_2 + k1*p1_2;
dp3_0 =         + k2*f_alpha2o - k_2*p3_0 - k3*f_alpha3o;
dp3_1 = 1*dSL*p3_0 + k2*f_alpha2i - k_2*p3_1 - k3*f_alpha3i;
dp3_2 = 2*dSL*p3_1 + k2*p2_2       - k_2*p3_2 - k3*p3_2;


%% Adding Campbell Ca activatin model (Campbell et al, Biophysical Journal 2018,)

% Calculating the forward and backward flux for the transition between the
% permissible and the non permissible state (Ca activation).
% (These transition rates DO NOT depend on SR versus NR state.)
 
Jon = k_on*Ca_i*N*(1 + K_coop*(1 - N));
Joff = k_off*P*(1 + K_coop*N);
dNp = - Jon + Joff;   

%% transitions between super relaxed state and non relaxed state
% 
% dU_NR = k1_SR * ( 1 + kforce * (F_XB + F_passive)) * U_SR - k2_SR * U_NR;
dU_NR = k1_SR * ( 1 + kforce * F_XB ) * U_SR - k2_SR * U_NR;

dYdT = [dp1_0; dp1_1; dp1_2; dp2_0; dp2_1; dp2_2; dp3_0; dp3_1; dp3_2; dNp; dSL; dU_NR];

end
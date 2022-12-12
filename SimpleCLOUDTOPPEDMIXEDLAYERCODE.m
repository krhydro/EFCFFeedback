%% README

% This MATLAB version code is to modify the original CLASS Python code in order to simplify 
% the atmospheric surface layer for a research purpose, the sensitivity of boundary
% layer clouds to Evaporative Fraction across the Contiminous United
% States.
% This research was submitted on Geophysical Research Letter on December
% 2022.

% Kyoungho Ryu and Guido D Salvucci have written this MATLAB code.
% The Department of Earth and Environment, Boston University, Boston, MA.


% Reference
% The Chemistry Land-surface Atmosphere Soil Slab model (CLASS)
% was introduced in the textbook: 
% Atmospheric boundary layer: Integrating air chemistry and land
% interactions (Vilà-Guerau de Arellano et al., 2015).
 
% The authors freely share the Python code at
% github.com/classmodel/modelpy. 

% Vilà-Guerau de Arellano, J., Chiel C. van Heerwaarden, van Stratum, B. J. H., & van den Dries, K. (2015). 
% Atmospheric boundary layer: Intergrating air chemistry and land interactions. Cambridge University Press.

%% Basic setting
clear
clc

% File Save Location

datapath = strcat('/Users/kyounghoryu/Desktop/CBL/UserFriendly-100622/');
cd(datapath)

pathSave   = strcat('/Users/kyounghoryu/Desktop/CBL/UserFriendly-100622/Test-avgEF.mat');

% Free parameters
beta                     = 0.2;  % Entrainment ratio for virtual heat [-]

pa_CHISTAR   = 1.0; % The parameter for moisture being trasponrted from the sub-cloud to cloud layers
                              % 0.84 used in van Stratum et al., (2014).
% van Stratum, B. J. H., Vilá-Guerau de Arellano, J., van Heerwaarden, C. C., & Ouwersloot, H. G. (2014). 
% Subcloud-layer feedbacks driven by the mass flux of shallow cumulus convection over land. Journal of the Atmospheric Sciences, 71(3), 881–895. https://doi.org/10.1175/JAS-D-13-0192.1

% Surface Resistance terms are located in the loop for the ABL dirunal
% evolution.


% Turn on & off

pa_QB       =1;  % if QB=0, no water vapor bouyancy impact, 1 includes water vapor bouyancy
pa_dhdt    = 1; % if 1, include mass flux correction, if zero, use that in the original CLASS model
pa_Rsfclcl = 1; % if 1, use a surface resistance to estimate LCL
                             % 1: use q_aS and ThetaS for ac

% Set constants

qscl=2500;  % exponential scale for initial humidity profile, used in gammaq
MIX=0.8;    %  the mix between surface parcels at surf q and entrainment layer....
% for MIX=0.8, air is an 80% entrainment layer 20% surface mixture.
% matches LES data well

pa_cefjq         = 0; % The coefficient for jump q, -0.125 was used in the CLASS

pa_DZMIN     = 300; % keep dz_h being a positive value.
advtheta             = 0;    % advection of heat [K s-1]
advq                     = 0;    % advection of moisture [kg kg-1 s-1]


% Esential Constants

julian = [1, 31   ; 32, 59  ; 60, 90  ; 91, 120 ; 121, 151; 152, 181; ...
           182, 212; 213, 243; 244, 273; 274, 304; 305, 334; 335, 365]; 
    
% Constant Input Variables
         Lv          = 2.5e6;             % heat of vaporization [J kg-1]
         k           = 0.4;                   % Von Karman constant [-]
         g           = 9.81;                  % gravity acceleration [m s-2]

% Define Constants
         BOLTZMAN = 1.380658e-23;
         AVOGADRO = .602214199e24;
         MD   = 28.9644e-3;
         MV   = 18.0153e-3;
         RV    = (AVOGADRO)*(BOLTZMAN) / (MV);            % constant f or water vapor (J/(kg-K))
         RD    = (AVOGADRO)*(BOLTZMAN) / (MD);            % constant for dry air (J/(kg-K))
         Cp     = 7./2*(RD);                                                                % specific heat of air (J/(kg-K))
         CPD = Cp;
         
% Lifiting Condensation Level
        EPSILON  = RD/RV; % Unitless

% Time setting
         dt   = 1*60;                % seconds
         t     = 24;                    % for Hours
          
         runtime  = t * 3600;                    % Total run time [s]. It needs to change for runtime
         tsteps     = round(runtime ./ dt);  % Seconds
        

 %%  load Inputs
 % SmapleBoston.mat contains data for the NARR's grid closest to Boston 
 % on July 2003.

 % the NARR stands for the North American Regional Reanalysis
 % Readers can find more details on 
 % Mesinger, F., DiMego, G., Kalnay, E., Mitchell, K., Shafran, P. C., Ebisuzaki, W., Jović, D., Woollen, J., Rogers, E., Berbery, E. H., Ek, M. B., Fan, Y., Grumbine, R., Higgins, W., Li, H., Lin, Y., Manikin, G., Parrish, D., & Shi, W. (2006). North American regional reanalysis. 
 % Bulletin of the American Meteorological Society, 87(3), 343–360. https://doi.org/10.1175/BAMS-87-3-343
 
 
 % The structure of each variable is a grid * # of measurements in a day * # of days in July.
 % The net radiation has half-hour measurements => the structure is 1*48* 30
 % The rest of the inputs is the measurement in the early morning near
 % sunrise. => 1 * 30. It will transform as 1 * 1* 30 while running.
% If users want to run your data, please make the data accordingly.

load('SampleBoston.mat')
 

%% Rn interpolation for a second measurement.
inputRn = avgRn;
 
inputRn_perm=permute(inputRn,[2,1,3]); 
inputRnB=interp1([0.25:0.5:23.75]/24, inputRn_perm, ([1:1:tsteps]-0.5)/tsteps,'linear','extrap');

clear inputRn_perm

inputRn=permute(inputRnB,[2,1,3]);
% Daytime is defined when the net radiation is greater than zero.
inputRn(inputRn<0)=0; 
 
clear inputRnB

'Initial Forcings'

[allgrid, ~, ndays] = size(inputRn);
inputPsfc          = zeros(allgrid, 1, ndays); 
inputRH           = zeros(allgrid, 1, ndays);
inputTheta       = zeros(allgrid, 1, ndays);
inputGamma   = zeros(allgrid, 1, ndays);   
inputGammaq = zeros(allgrid, 1, ndays);

inputRHsfc       = zeros(allgrid, 1, ndays);
inputThetasfc  = zeros(allgrid, 1, ndays);

% Setting EF for a sensitivity test
EF = zeros(allgrid, 1, ndays, 3);
EF(:, 1, :, 2) = avgEF;
EF(:, 1, :, 1) = avgEF - 0.05;
EF(:, 1, :, 3)  = avgEF + 0.05;
EF ( EF <= 0 ) = 0.01;

% Inputs
inputPsfc(:, 1, :)    = avgPsfc; % [Pa]
inputRH(:, 1, :)     = avgRH; % on scale 0 to 1.
inputTheta(:, 1, :) = avgT; % [K]
inputGamma(:, 1, :)  = avgGam; % [K/m]  
inputGammaq(:, 1, :) = avgGamq;

clear  avgGam avgPsfc avgRH avgRn avgT avgGamq 


% Preallocation
[~, ~, ~, iEF] = size(EF);

% save 1/2 hourly values

tstepsSave=48;

all_ac      = zeros(allgrid, tstepsSave, ndays, iEF);
all_Q0    = zeros(allgrid, tstepsSave, ndays, iEF );
all_tau     = zeros(allgrid, tstepsSave, ndays, iEF );
all_DQ     = zeros(allgrid, tstepsSave, ndays, iEF );
all_wq      = zeros(allgrid, tstepsSave, ndays, iEF );
all_rho     = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_M        = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_we       = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_q2h     = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_wstar  = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_dzh      = zeros(allgrid, tstepsSave, ndays, iEF ); 

all_wqe                  = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_wqe_cor         = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_wthetae           = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_wthetae_cor  = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_jumpq              = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_jumptheta       = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_qsat_BLtop    = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_Rn                     = zeros(allgrid, tstepsSave, ndays, iEF ); 
all_def                    = zeros(allgrid, tstepsSave, ndays, iEF ); 

all_h                       = zeros(allgrid, tstepsSave, ndays, iEF );
all_lcl                     = zeros(allgrid, tstepsSave, ndays, iEF );
all_Theta              = zeros(allgrid, tstepsSave, ndays, iEF );
all_qa                    = zeros(allgrid, tstepsSave, ndays, iEF );
all_qaS                 = zeros(allgrid, tstepsSave, ndays, iEF );
all_RHsfc             = zeros(allgrid, tstepsSave, ndays, iEF );
all_wqM               = zeros(allgrid, tstepsSave, ndays, iEF );
all_SHmodeled  = zeros(allgrid, tstepsSave, ndays, iEF );
all_LEmodeled  = zeros(allgrid, tstepsSave, ndays, iEF );

%% MIXXED LAYER

% Daily BL Growth with all EF cases.

for idx_EF = 1 : iEF
        
tstepsRun=1;

% Intergrate mixed-layer equations
h                  = zeros(allgrid, tstepsRun, ndays);
Theta          = zeros(allgrid, tstepsRun, ndays);
jumptheta = zeros(allgrid, tstepsRun, ndays);
qsat_a        = zeros(allgrid, tstepsRun, ndays);
q_a             = zeros(allgrid, tstepsRun, ndays);
jumpq        = zeros(allgrid, tstepsRun, ndays);
rho              = zeros(allgrid, tstepsRun, ndays); 

% Surface properties
qsat_sfc   = zeros(allgrid, tstepsRun, ndays);
RH_sfc     = zeros(allgrid, tstepsRun, ndays);

% Surface Fluxes
LEmodeled  = zeros(allgrid, tstepsRun, ndays);
SHmodeled  = zeros(allgrid, tstepsRun, ndays);
wtheta            = zeros(allgrid, tstepsRun, ndays);
wq                   = zeros(allgrid, tstepsRun, ndays);
thetav            = zeros(allgrid, tstepsRun, ndays);
wthetav         = zeros(allgrid, tstepsRun, ndays);
jumpthetav = zeros(allgrid, tstepsRun, ndays);

% Mixed-layer top properties
P_BLtop        = zeros(allgrid, tstepsRun, ndays);
T_BLtop        = zeros(allgrid, tstepsRun, ndays);
qsat_BLtop  = zeros(allgrid, tstepsRun, ndays);
RH_BLtop    = zeros(allgrid, tstepsRun, ndays);

% LCL
T_lcl       = zeros(allgrid, tstepsRun, ndays);
P_lcl       = zeros(allgrid, tstepsRun, ndays);
lcl             = zeros(allgrid, tstepsRun, ndays);     
qsat_lcl  = zeros(allgrid, tstepsRun, ndays);
RH_lcl   = zeros(allgrid, tstepsRun, ndays);

% Calculate convective velocity scale w*
wstar        = zeros(allgrid, tstepsRun, ndays);
% Virtual heat entrainment flux
wthetave = zeros(allgrid, tstepsRun, ndays);
% Compute mixed-layer tendencies, no shear BL growth
we             = zeros(allgrid, tstepsRun, ndays); 

wthetae   = zeros(allgrid, tstepsRun, ndays);
wqe           = zeros(allgrid, tstepsRun, ndays);

% cumulus parameterization
ac            = zeros(allgrid, tstepsRun, ndays);
M            = zeros(allgrid, tstepsRun, ndays);
wqM      = zeros(allgrid, tstepsRun, ndays);
q2_h     = zeros(allgrid, tstepsRun, ndays);
dz_h      = zeros(allgrid, tstepsRun, ndays);
        
% Calculate entrainment fluxes
htend                    = zeros(allgrid, tstepsRun, ndays);
thetatend             = zeros(allgrid, tstepsRun, ndays);
qtend                    = zeros(allgrid, tstepsRun, ndays);
jumpthetatend  = zeros(allgrid, tstepsRun, ndays);
jumpqtend         = zeros(allgrid, tstepsRun, ndays);
dztend                 = zeros(allgrid, tstepsRun, ndays);

% Modification
wthetae_cor    = zeros(allgrid, tstepsRun, ndays);
wqe_cor           = zeros(allgrid, tstepsRun, ndays);

ThetaS  = zeros(allgrid, tstepsRun, ndays);
q_aS    = zeros(allgrid, tstepsRun, ndays);
        
        
       'Every EF'
        [idx_EF, iEF]      

        
% Assign initial values
h   (:, 1, :)  = 200;  % initial ABL height [m]

% Temperature  
Theta (:, 1, :)         = inputTheta;
jumptheta(:, 1, :) = (inputGamma .* h(:, 1, :) .* beta) ./ (1 + 2 * beta);  % initial temperature jump at h= the entrainment zone. the default = 1 [K]
dz_h(:, 1, :)           = 150;  % Transition layer thickness  from dry to moist convection [m]

% Humidity
qsat_a(:, 1, :)  = calcqsat( Theta(:, 1 , :), inputPsfc );
q_a   (:, 1, :)    =  inputRH .* qsat_a(:, 1, :); % Initial mixed-layer specific humidity [kg kg-1]


% for jumpq
jumpq  (:, 1, :)  = pa_cefjq .* q_a(:, 1, :);    % initial specific humidity jump at h [kg kg-1], -0.125 .* q(:,1)


% Initial Mass Flux Change between BL and the free atmosphere. 
ac  (:, 1, :)  = 0.05 ; % Cloud core fraction
q2_h(:, 1, :)  = 0;

LE_modeled = inputRn(:, :, :) .*  EF(:, :, :, idx_EF) ;
SH_modeled = inputRn(:, :, :) .*  ( 1- EF(:, :, :, idx_EF) ) ; 


'Time Evolution'


% These are used later for calculating variance of q from last time
% step.  here they are initialized, later they are saved at end of time
% step loop.

   wqe_lasttimestep          = wqe.* 0;
   wqM_lasttimestep        = wqM .* 0;
   wqe_cor_lasttimestep = wqe_cor .* 0;

tind=0; % This index (tind..time index) is used later to save 1/2 hour averages

for i1 = 1: tsteps

            i1
i=1;  % i is always 1 to make this code faster.

          [idx_EF, iEF,i1,tsteps] 

          
% making rho vary and needed to move wtheta and wq down
rho(:,i,:) = ( (inputPsfc(:, 1, :) ) ./ (RD .* Theta (:, i, :)));

% surface fluxes chaning time to time via i1.
 wtheta   = SH_modeled(:,i1,:)  ./ (rho(:, i, :) .* Cp);  % Kinematic Heat flux [K * m/s]
 wq         = LE_modeled(:,i1,:)  ./ (rho(:, i, :) .* Lv);  % Kinematic Moisture flux [kg/kg * m/s], the default 0.0001;

% Calculate Virtual Temperatures
thetav (:, i, :)          = Theta (:, i, :)  + pa_QB.*0.61 .* Theta(:, i , :) .* q_a(:, i, :); % [K]
wthetav(:, i, :)        = wtheta(:, i, :)  + pa_QB.*0.61 .* Theta(:, i , :) .* wq (:, i, :); % [K * m/s]
jumpthetav(:, i, :) = ( (Theta(:, i, :) + jumptheta(:, i, :) ) .* (1 + pa_QB.*0.61 .* ( q_a(:, i, :) + jumpq (:, i, :) ) ) ) ...
                                - Theta(:, i , :) .* ( 1+ pa_QB.*0.61 .* q_a(:, i, :) ); % [K]

 
% Surface Properties
qsat_sfc(:, i, :) = calcqsat( Theta(:, i , :), inputPsfc ); % [on scale 0 to 1]
RH_sfc  (:, i, :) = q_a(:, i, :) ./ qsat_sfc(:, i, :);   % [on scale 0 to 1]


%% Lifiting Condensation Level

% bulk aero values range from about 10 to 200, so this is
% reasonable


SURFRA  = 38;  % seconds/meter
SURFRAq = 38; % seconds/meter

ThetaS(:, i , :) = Theta(:, i , :) + pa_Rsfclcl .* wtheta(:, i , :) .* SURFRA; % [K]
q_aS(:, i , :)   = q_a(:, i , :) + pa_Rsfclcl .* wq(:, i , :) .* SURFRAq; %[kg/kg]

T_lcl(:, i ,:)  = ...
              2840 ./ (3.5 .* log( ThetaS(:, i , :) ) - log( (inputPsfc )  ./ 1000 .* q_aS(:, i , :) ./ (EPSILON + q_aS(:, i , :))) - 7.108) + 55;
P_lcl(:, i , :)  = inputPsfc .* ( T_lcl(:, i , :) ./ ThetaS(:, i , :) ).^3.5;

lcl  (:, i , :)  = CPD  .* ThetaS(:, i , :) ./ g .* (1 - (P_lcl(:, i , :) ./ inputPsfc ).^(RD ./ CPD)); % [m]
lcl  (:, i , :)  = max(lcl  (:, i , :), 200); % to prevent imaginary value

qsat_lcl(:, i , :) = calcqsat(T_lcl(:, i , :), P_lcl(:, i , :) );
RH_lcl  (:, i , :) = q_a(:, i , :) ./ qsat_lcl(:, i , :);


%% Mixed Layer

% Calculate convective velocity scale w*,    [m/s], on the order of 1 m/s.
%  ( wthetav(:, i, :) ./ thetav(:, i, :) ) .* g => bouyanct force

    wstar(:, i , :) = ( (g ./ thetav(:, i , :) ) .* wthetav(:, i, :) .* h(:, i, :)  ).^(1/3); % Deardorff, James W 1970
    wstar( wthetav(:, i, :) <= 0 ) = 10.^(-6);
       
    wstar( wstar < 1e-6 )=1e-6;  % added to not divide by zero later which causes NAN and crash
   % it works only one of the dimension is one, other than that, please use
   % different ways.
   
% Mixed-layer top properties, qsat_BLtop is used for the estimation of ac.

P_BLtop(:, i, :)      = inputPsfc - (rho(:, i, :) .* g .* h(:, i, :));   
T_BLtop(:, i, :)      = Theta(:, i , :) - ((g/Cp) .* h(:, i, :));  % [K]
qsat_BLtop(:, i, :) = calcqsat( T_BLtop(:, i, :), P_BLtop(:, i, :) );
RH_BLtop  (:, i, :) = q_a(:, i, :) ./ qsat_BLtop(:, i, :);

% Virtual heat entrainment flux
wthetave(:, i , :) = -beta .* wthetav(:, i , :); % [K* m/s], toward the ABL

% Compute mixed-layer tendencies, no shear BL growth, The entrainment
% velocity, typical values are often in the range of 0.01 m/s to 0.20 m/s.
we(:, i , :) = -wthetave(:, i , :) ./ jumpthetav(:, i , :); % [m/s], toward the FA

% No BL sinking if wtheta <0
we(we(:, i, :) < 0 ) = 0;

% Calculate entrainment fluxes

wthetae(:, i , :)   = -we(:, i , :) .* jumptheta(:, i , :); % [ K * m/s]
wqe    (:, i , :)      = -we(:, i , :) .* jumpq    (:, i , :);  % [ kg/kg * m/s]


%%  The Final version of ac

            % gamma q calcs
            
            
            if i1==1
                Q0=q_a   (:, 1, :);
                Theta0=Theta(:, 1, :);                
                gammatheta=inputGamma;                
            end
            
            
            % note, if q/qo=exp(-z/zo) then dqdz at surface=qo/z
            % this is confirmed in gfdl model with mean scale of 2500
            % and very small spatial standard deviation of about 100 M
            
            if i1==1
                gamsurf=-Q0(:,1,:)  ./ qscl;
            end
            
            HO=200; % initial height
            inputGammaq=gamsurf .* exp(-(h-HO)./qscl);
   
            
            % note that even though altering jumpq to account for moisture
            %above h, if wqM linear from value at h to zero at ZD, then
            % dqdt is constant with height in this zone, and thus the gammaq
            % does not actually change...stays simply dependent on z through above
            
            tau=( h(:, i, :) ./ wstar(:, i, :) ); 
            
            %%%%%%%%%%%%%%%%%%%%
             % variance calcs
            
            % straight from jordi
            % using fluxes from last time step, saved below
            DQ=-(jumpq(:, i , :)); % 
                    
            q2_h(:, i , :) = (wqe_lasttimestep(:, i , :) + wqM_lasttimestep(:, i , :)+wqe_cor_lasttimestep(:, i , :)) .* ...
                    ( DQ(:, i , :) ) .* (1./ wstar(:, i, :) ); % [Unitless]
                
                        
            %%%%%%%
            % ac calcs
            
            
            q2_h( q2_h <=  1e-10 )= 1e-10;
            
            
            sigq=(q2_h(:, i , :).^0.5);

            
            qmix2 = (q_aS(:, i , :));
                    
            qmix1 = ( (q_a(:, i , :))- (1-MIX).*qmix2 ) ./ MIX;
            
            qtest =  MIX.*qmix1+ (1-MIX).*qmix2; 


            % MIX of 0.8 and SFAC of 3 are very consistent with double-hump
            % pdf of q in LES data.... this is a mix of two gaussian
            % distributions, one centered at qmix1 with stand dev of sigq,
            % and one centered at qmis2 with standard dev of sigq/SFAC
            % note from above "qtest" that the average of the two mixtures
            % yields the bulk qa...all of this is validated in LES data
            
            SFAC=3;
            
            ac1=normcdf(qmix1, qsat_BLtop(:, i , :), sigq);
            ac2=normcdf(qmix2, qsat_BLtop(:, i , :), sigq/SFAC);
       
            
            ac= MIX.*ac1 + (1-MIX) .* ac2;
                      
            %%%%%%%%%%%%%%%%%%
            
            
            
            
            M  (:, i, :)    = ac(:, i , :) .* wstar(:, i , :);               % [m/s]
            
            %
            
            M(M>we)=we(M>we); % added to make dhdt pos to h doesnt shrink

            
            
            
            DRIVER=sqrt(q2_h); % 
            
            
            
            
            wqM(:, i, :) =  M(:, i , :) .* DRIVER.*pa_CHISTAR ;  % [kg/kg * m/s]
            
            



%% Terms Additions, It stems from dh/dt = we-M
wthetae_cor(:, i , :) = pa_dhdt .* M(:, i , :) .* jumptheta(:, i , :); % [K * m/s]
wqe_cor(:, i , :)        = pa_dhdt .* M(:, i , :) .* jumpq(:, i , :);        % [ kg/kg * m/s]
  

% adding...needed for calculating variance of q from last time step
   wqe_lasttimestep          = wqe;
   wqM_lasttimestep        = wqM;
   wqe_cor_lasttimestep = wqe_cor;

%% Calculate  tendency of key variables, time defferential eqns

htend     (:, i , :)     = we(:, i , :) - M(:, i , :); % [m/s]

thetatend (:, i , :)            = ( wtheta(:, i , :) - wthetae(:, i , :) - wthetae_cor(:, i , :) ) ./ h(:, i , :) + advtheta;  % [ K * 1/s]
qtend     (:, i , :)               = ( wq(:, i , :)     - wqe(:, i , :)  - wqe_cor(:, i , :)    - wqM(:, i , :) )     ./ h(:, i , :) + advq;  % [kg/kg * 1/s]

jumpthetatend(:, i , :)  = inputGamma(:, i, :) .* ( we(:, i , :) - M(:, i , :) ) - thetatend(:, i , :);  % [ K * 1/s]
jumpqtend    (:, i , :)      = inputGammaq(:, i, :) .* ( we(:, i , :) - M(:, i , :) ) - qtend(:, i , :) ;      % [ kg/kg * 1/s]

%% Tendency of the transition layer thickness

dztend(:, i , :) = ( (lcl(:, i , :) - h(:, i , :) ) - dz_h(:, i , :) ) ./ 7200; % [m]

for idx_days = 1  : ndays
    
dztend( ac(:, i , idx_days) <= 0,  i , idx_days) = 0;

dztend( ( lcl(:, i , idx_days) - h(:, i , idx_days) ) >= 300, i , idx_days) = 0;

end

% Omit windtend

%% intergrate mixed layer

% Intergrate mixed-layer equations

h     (:, i , :)        = real( h     (:, i , :) + dt .* htend     (:, i , :) ); % [m]
h     (:, i , :)        = max( h     (:, i , :), 200); % to prevent h from being a weired value.
Theta (:, i , :)   = real( Theta (:, i , :) + dt .* thetatend (:, i , :) ); % [K]
Theta (:, i , :)   = max( Theta (:, i , :), 100); % to prevent h from being a weired value.

jumptheta(:, i , :) = real( jumptheta(:, i , :) + dt .* jumpthetatend(:, i , :) ); % [K]

q_a   (:, i , :)   = real( q_a   (:, i , :) + dt .* qtend     (:, i , :) ); %[ kg/kg ]
q_a   (:, i , :)   = max ( q_a   (:, i , :), 0.000001); % to prevent h from being a weired value.

jumpq    (:, i , :)   = real( jumpq    (:, i , :) + dt .* jumpqtend    (:, i , :) ); % [kg/kg]  % gds
 
dz_h  (:, i , :)   = real( dz_h  (:, i , :) + dt .* dztend    (:, i , :) ); % [m]

dz_h  (:, i , :)   = min(pa_DZMIN, dz_h(:, i , :)  ); % to prevent dz_h from being a weired value. 

% The code below sums the variables at each time step and then when the
% total time reaches 1800 seconds, starts summing them again for the next
% 1/2 hour.  the code is much much faster this way...

incs=round(1800/dt);

if mod(i1,incs)==1

    tind=tind+1;

end


all_h     (:,tind,:,idx_EF)      = all_h   (:,tind,:,idx_EF)       + h/incs;
all_Theta (:,tind,:,idx_EF)  = all_Theta(:,tind,:,idx_EF) + Theta/incs;
all_qa    (:,tind,:,idx_EF)      = all_qa (:,tind,:,idx_EF)      + q_a/incs;
all_qaS    (:,tind,:,idx_EF)    = all_qaS (:,tind,:,idx_EF)   + q_aS/incs;

all_Q0   (:,tind, :, idx_EF) = all_Q0(:,tind,:,idx_EF) + Q0/incs;

all_tau   (:,tind,:,idx_EF)   = all_tau(:,tind,:,idx_EF) + tau/incs ;
all_DQ    (:,tind,:,idx_EF)  = all_DQ(:,tind,:,idx_EF) + DQ/incs ;
all_wq    (:,tind,:,idx_EF)   =  all_wq(:,tind,:,idx_EF)+wq/incs;
all_wstar (:,tind,:,idx_EF) = all_wstar(:,tind,:,idx_EF) + wstar/incs;

all_q2h   (:,tind,:,idx_EF)  = all_q2h(:,tind,:,idx_EF) + q2_h/incs;
all_ac    (:,tind,:,idx_EF)      = all_ac  (:,tind,:,idx_EF) + ac/incs;
all_wqM   (:,tind,:,idx_EF)  = all_wqM  (:,tind,:,idx_EF)  +wqM/incs;

all_lcl   (:,tind,:,idx_EF)       = all_lcl  (:,tind,:,idx_EF)  +lcl/incs;
all_RHsfc (:,tind,:,idx_EF)  = all_RHsfc (:,tind,:,idx_EF) + RH_sfc/incs;
all_wqe   (:,tind,:,idx_EF)    = all_wqe (:,tind,:,idx_EF)  +wqe/incs ;
all_wqe_cor   (:,tind,:,idx_EF)   = all_wqe_cor   (:,tind,:,idx_EF)  + wqe_cor/incs;
all_SHmodeled (:,tind,:,idx_EF) = all_SHmodeled(:,tind,:,idx_EF) +SH_modeled(:,i1,:)/incs;
all_LEmodeled (:,tind,:,idx_EF) = all_LEmodeled(:,tind,:,idx_EF)  +LE_modeled(:,i1,:) /incs;

all_rho   (:,tind,:,idx_EF)      = all_rho      (:,tind,:,idx_EF) + rho/incs ;
all_dzh   (:,tind,:,idx_EF)     = all_dzh      (:,tind,:,idx_EF) +dz_h/incs;
all_jumpq (:,tind,:,idx_EF)  = all_jumpq (:,tind,:,idx_EF) + jumpq/incs ;

all_jumptheta (:,tind,:,idx_EF)  = all_jumptheta  (:,tind,:,idx_EF) + jumptheta/incs ;
all_M         (:,tind,:,idx_EF)  = all_M    (:,tind,:,idx_EF) + M/incs ;
all_we        (:,tind,:,idx_EF)  =  all_we    (:,tind,:,idx_EF) +we/incs ;

all_qsat_BLtop(:,tind,:,idx_EF)  = all_qsat_BLtop(:,tind,:,idx_EF) +qsat_BLtop/incs;

all_def       (:,tind,:,idx_EF)  = all_def (:,tind,:,idx_EF) +(q_a-qsat_BLtop)/incs ;

all_Rn        (:,tind,:,idx_EF)  = all_Rn (:,tind,:,idx_EF) +inputRn(:,i1,:)/incs; 


end % Time evolution
       
% keep the model settings while running.
clearvars -except beta pa_CHISTAR EF qscl MIX  ...
   pa_QB pa_dhdt pa_Rsfclcl pa_cefjq pa_DZMIN advtheta advq ...
   julian Lv k g RD RV Cp CPD EPSILON dt t runtime tsteps pathSave ...
   inputRn inputPsfc inputRH inputTheta inputGamma inputGammaq ...
   inputThetasfc inputRHsfc LE_modeled SH_modeled allgrid ndays ...
    all_h all_Theta all_qa all_Q0 all_tau all_DQ  all_wq  all_wstar  ...
    all_q2h  all_ac  all_lcl all_RHsfc all_wqM all_wqe  ... 
    all_wqe_cor all_LEmodeled LE_modeled all_SHmodeled SH_modeled ... 
    all_rho rho all_dzh all_jumpq all_jumptheta all_M all_we all_qsat_BLtop ...
    all_def all_Rn all_qaS avgEF lat lon iEF
    
end % idx EF 

% The following variables 
all_EFBLtop= mean( (all_wqM .* Lv .* all_rho), 2, 'omitnan' ) ./ (repmat( mean(inputRn, 2, 'omitnan'), [1,1,1, iEF ]));
all_wqMday = mean (all_wqM, 2, 'omitnan' ) ; 

%% save necessary variables.
% Users can save variables based on their interest.       

'saving'

save(pathSave, 'EF', 'lat', 'lon', ...
'pa_dhdt', 'pa_Rsfclcl', 'pa_QB', 'pa_cefjq', 'pa_CHISTAR', 'pa_DZMIN', ...
'all_EFBLtop','avgEF','all_Rn', 'all_h', 'all_qa', 'all_qaS', 'all_Theta', 'all_dzh', ... , 
'all_q2h', 'all_wqMday', 'all_ac', 'all_qsat_BLtop', 'all_wstar', 'all_lcl', ...
'-v7.3')  



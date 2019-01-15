function biochem_out = biochemical_new(biochem_in,Ci_input)
% This is the Updated Version of the SCOPE Biochem Module
% This function calculates:
%    - stomatal resistance of a leaf or needle (s m-1)
%    - photosynthesis of a leaf or needle (umol m-2 s-1)
%    - fluorescence of a leaf or needle (fraction of fluor. in the dark)
%
% Usage:
% biochem_out = biochemical(biochem_in)
% the function was tested for Matlab R2017b
%
% Calculates net assimilation rate A, fluorescence F and stomatal conductance gs using biochemical model
%
% *****
% Original Contributors
% Authors: 	Joe Berry and Christiaan van der Tol, Ari Kornfeld, contributions of others.
% Sources: 	Farquhar et al. 1980, Collatz et al (1991, 1992).
%
% Update: Dec 22 2017 - Debsunder Dutta
% Update: Mar 23 2018 - Debsunder Dutta
%
%
% Photosynthesis Implementation consistent with the CLM Version 4.5
% A-Ci-gs iterations are implemented
% as per the simple Newton Raphson Scheme as per Sun et al
% Temperature Dependence of the Photosynthetic parameters are implemented
% also as per CLM 4.5 (Bonan et. al 2011 paper)
% The Temperature Dependence of C3 Species is consistent with Leuning et.
% al 2001 paper
% The parameters of the C3 Model are derived from the Bonan Paper but one can
% also obtain the parameters from Leuning paper as well
% Currently the parameters are read from the input excel sheet and passed
% into the function as an excel file
% The parameters for the Temperature dependence function of C3 species includes the
% activation energy, deactivation energy and the entropy terms for Vcmax,
% Jmax, TPU, Rd and only Activation energy for Kc, Ko and Gamma_star
% There is the temperature dependence function and high temperature
% inhibition
% For the C4 species the temperature dependence parameters include Q10 and 
% constants s1, s2, s3, s4, s5, s6
% For the C4 species the corrections are for the Vcmax, Rd, and Ke
%
% We have also Incorporated the Computation of Photorespiration (rate of
% oxygenation vo. The formulation has been adapted from Sharkey 1988 paper
%
%
% Input (units are important):
% structure 'biochem_in' with the following elements:
% Knparams   % [], [], []
% Parameters for empirical Kn (NPQ) model:
% Kn = Kno * (1+beta).*x.^alpha./(beta + x.^alpha);
%              [Kno, Kn_alpha, Kn_beta]
%   or, better, as individual fields:
%   Kno                             Kno - the maximum Kn value ("high light")
%   Kn_alpha, Kn_beta               alpha, beta: curvature parameters
%
% Cs                % [ppmV or umol/mol]        initial estimate of conc. of CO2 in the
%                                               ...boundary layer of the leaf
% Q                 % [umol photons m-2 s-1]    net radiation, PAR
% fPAR              % [0-1]                     fraction of incident light that is absorbed by the leaf (default = 1, for compatibility)
% T                 % [oC or K]                 leaf temperature
% eb                % [hPa = mbar]              intial estimate of the vapour pressure in leaf boundary layer
% O                 % [mmol/mol]                concentration of O2 (in the boundary
%                                               ...layer, but no problem to use ambient)
% p                 % [hPa]                     air pressure
% Vcmax25(Vcmo)     % [umol/m2/s]               maximum carboxylation capacity @ 25 degC
% BallBerrySlope(m) % []                        Ball-Berry coefficient 'm' for stomatal regulation
% BallBerry0        % []                        (OPTIONAL) Ball-Berry intercept term 'b' (if present, an iterative solution is used)
%                                               setting this to zeo disables iteration. Default = 0.01
%
% Type              % ['C3', 'C4']              text parameter, either 'C3' for C3 or any
%                                               ...other text for C4
% tempcor           % [0, 1]                    boolean (0 or 1) whether or not
%                                               ...temperature correction to Vcmax has to be applied.

% effcon            % [mol CO2/mol e-]          number of CO2 per electrons - typically 1/5 for C3 and 1/6 for C4

% RdPerVcmax25 (Rdparam)  % []         respiration as fraction of Vcmax25
% stressfactor [0-1]                   stress factor to reduce Vcmax (for
%                                      example soil moisture, leaf age). Use 1 to "disable" (1 = no stress)
% OPTIONAL PARAMETERS
% Kpep25 (kp)    % [umol/m2/s]         PEPcase activity at 25 deg C (defaults to Vcmax/56
% atheta         % [0-1]               smoothing parameter for transition between Vc and Ve (light- and carboxylation-limited photosynthesis)
% useTLforC3     % boolean             whether to enable low-temperature attenuation of Vcmax in C3 plants (its always on for C4 plants)
% po0            %  double             Kp,0 (Kp,max) = Fv/Fm (for curve fitting)
% g_m            % mol/m2/s/bar        Mesophyll conductance (default: Infinity, i.e. no effect of g_m)
% Note: always use the prescribed units. Temperature can be either oC or K
% Note: input can be single numbers, vectors, or n-dimensional
% matrices
%
% Output:
% structure 'biochem_out' with the following elements:
% A         % [umol/m2/s]           net assimilation rate of the leaves
% Cs        % [umol/m3]             CO2 concentration in the boundary layer
% eta0      % []                    fluorescence as fraction of dark
%                                   ...adapted (fs/fo)
% rcw       % [s m-1]               stomatal resistance
% qE        % []                    non photochemical quenching
% fs        % []                    fluorescence as fraction of PAR
% Ci        % [umol/m3]             internal CO2 concentration
% Kn        % []                    rate constant for excess heat
% fo        % []                    dark adapted fluorescence (fraction of aPAR)
% fm        % []                    light saturated fluorescence (fraction of aPAR)
% qQ        % []                    photochemical quenching
% Vcmax     % [umol/m2/s]           carboxylation capacity after
%                                   ... temperature correction

if nargin < 2 % The internal leaf CO2 concentrations are not provided
    Ci_input = [];
end
%% INPUT VARIBLE PROCESSING and PARSING

global constants

% Environmental or Meteorological Variable Checking
if isfield(biochem_in, 'Cs')
    %assert(all(biochem_in.Cs(:) >=0), 'Negative CO2 (Cs) is not allowed!');
    Cs         = max(0, biochem_in.Cs); % just make sure we don't deal with illegal values
else
    % if Cs is missing, Ci must have been supplied. Forcing Cs = NaN invalidates rcw & gs.
    Cs         = NaN; %biochem_in.Ci;
end

Q             = biochem_in.Q;   % APAR Input Value in umol/m2/s so no need of conversion

assert(all(Q(:) >=0), 'Negative light is not allowed!');

T             = biochem_in.T + 273.15*(biochem_in.T<200); % convert temperatures to K if not already
eb            = biochem_in.eb;    % Read in the other met variables from the Input Structure
O             = biochem_in.O;
p             = biochem_in.p;

% Plant Physiological Inputs (Read in or assume default values)

Type                   = biochem_in.Type;
if isfield(biochem_in, 'Vcmax25')       % new field names Vcmax25, BallBerrySlope, RdperVcmax25
    Vcmax25            = biochem_in.Vcmax25;
    
    BallBerrySlope     = 9; % Default assumption of BB Slope
    if isfield(biochem_in, 'BallBerrySlope')  % with g_m and Ci specified, we don't pass BBslope
        BallBerrySlope = biochem_in.BallBerrySlope;
    end
    RdPerVcmax25       = biochem_in.RdPerVcmax25; % This Value of Rdparam25 is for C3 we specify C4 elsewhere
else
    
    % old field names: Vcmo, m, Rdparam
    Vcmax25           = biochem_in.Vcmo;
    BallBerrySlope    = biochem_in.m;
    RdPerVcmax25      = biochem_in.Rdparam;
end


BallBerry0 = 0.01;           % default valueassumed is 0.01 if not provided
if isfield(biochem_in, 'BallBerry0')
    BallBerry0 = biochem_in.BallBerry0;
end


if isfield(biochem_in, 'effcon') % Assigning Effective conductivity value if not provided
    effcon        = biochem_in.effcon;
elseif strcmpi('C3', Type)
    effcon =  1/5;
else
    effcon = 1/6; % C4
end


% Mesophyll conductance: by default we ignore its effect
%  so Cc = Ci - A/gm = Cin
g_m = Inf; % Default Value

if isfield(biochem_in, 'g_m') % If present Read from structure
    g_m = biochem_in.g_m * 1e6; % convert from mol to umol
end

% SCOPE provides PAR as APAR, so default (for SCOPE) = 1
%    The curve-fitting GUI may not be providing APAR and should therefore explicitly set fPAR
fPAR = 1;                           % fraction of incident light that is absorbed by the leaf
if isfield(biochem_in, 'fPAR')      % Or read from the structure if present
    fPAR = biochem_in.fPAR;
end

% Physiological computation options
tempcor           = biochem_in.tempcor;             % Vcmax Temp correction yes/no [1 or 0]
stressfactor      = biochem_in.stressfactor;        % Stress factor apply [1-no stress or 0- totally stressed]
%model_choice  = biochem_in.Fluorescence_model;

%  NOTE: kpep (kp), atheta parameters in next section

%% Fluorescence Parameters (FLUORESCENCE PARAMETERS NEEDS CLEANING AND UPDATING)
if isfield(biochem_in, 'Knparams')
    Knparams      = biochem_in.Knparams;
elseif isfield( biochem_in, 'Kn0')
    Knparams = [biochem_in.Kn0, biochem_in.Kn_alpha, biochem_in.Kn_beta];
elseif isfield(biochem_in, 'Fluorescence_model') && biochem_in.Fluorescence_model==0
    % default drought values:
    Knparams = [5.01, 1.93, 10];
else
    % default general values (cotton dataset)
    Knparams = [2.48, 2.83, 0.114];
end

if isfield(biochem_in, 'po0')
    po0 = biochem_in.po0;
else
    po0 = [];
end

% Electron transport and fluorescence
Kf          = 0.05;                                    % []  rate constant for fluorescence
Kd          = max(0.8738,  0.0301*(T-273.15)+ 0.0773); % []  rate constant for thermal deactivation at Fm other value Kd = 0.95
Kp          = 4.0;                                     % []  rate constant for photochemisty

%% Some other required Leaf Photosynthesis Parameter definitions (at optimum Temperature)
Tref        = 25+273.15;        % [K]           absolute temperature at 25 oC

Jmax25      = 1.97*Vcmax25;
TPU25       = 0.06*Jmax25;              % Triose Phosphate utilization rate

% All of the Following Values are Adopted from Bonan et. al 2011 paper and
% pertaings to C3 photosythesis

Kc25         = 404.9;                   % [umol mol-1]
Ko25         = 278.4;                   % [mmol mol-1]
spfy25       = 2444.4;                  % specificity (Computed from Bernacchhi et al 2001 paper)
Gamma_star25 = 0.5 .*O.*1E3 ./spfy25;   % [ppm] compensation point in absence of Rd


% note:  rhoa/Mair = L/mol (with the current units) = (28.96/1.2047) = 24.039 L/mol
% and  V/n = RT/P ==>  T = 292.95 K @ 1 atm (using R_hPa = 83.144621; 1 atm = 1013.25 hPa)
% Which is the molar volume (V/n) at STP conditions computes to be 24.038

%   These values are used only for computing rcw
rhoa        = 1.2047;           % [kg m-3]       specific mass of air
Mair        = 28.96;            % [g mol-1]      molecular mass of dry air

%% Temperature Correction Functions
% The following two functions pertains to C3 photosynthesis
    function [fTv] = temperature_functionC3(Tref,Rgas,T,deltaHa)
        % Temperature function
        tempfunc1 = (1 - Tref./T);
        fTv = exp(deltaHa/(Tref*Rgas).*tempfunc1);
    end

    function [fHTv] = high_temp_inhibtionC3(Tref,Rgas,T,deltaS,deltaHd)
        % High Temperature Inhibition Function
        hightempfunc_num = (1+exp((Tref*deltaS-deltaHd)/(Tref*Rgas)));
        hightempfunc_deno = (1+exp((deltaS.*T - deltaHd)./(Rgas.*T)));
        fHTv = hightempfunc_num ./ hightempfunc_deno;
    end

Rgas = constants.R;      % Unit is [J K^-1 mol^-1]

% Apply temperature Corrections
if (tempcor == 1)
    if strcmpi('C3', Type)
        % C3 species
        
        Rd25 = RdPerVcmax25 * Vcmax25;
        
        % Constant parameters for temperature correction of Vcmax
        deltaHa  = biochem_in.TDP.delHaV;                % Unit is  [J K^-1]
        %Rgas = 8.314;      % Unit is [J K^-1 mol^-1]
        deltaS   = biochem_in.TDP.delSV;                 % unit is [J mol^-1 K^-1]
        deltaHd  = biochem_in.TDP.delHdV;                % unit is [J mol^-1]
        fTv = temperature_functionC3(Tref,Rgas,T,deltaHa);
        fHTv = high_temp_inhibtionC3(Tref,Rgas,T,deltaS,deltaHd);
        
        Vcmax = Vcmax25 .* fTv .* fHTv * stressfactor;   % Temperature Corrected Vcmax
        
        % Constant parameters for temperature correction of Jmax
        deltaHa =  biochem_in.TDP.delHaJ;                % Unit is  [J K^-1]
        deltaS  =  biochem_in.TDP.delSJ;                 % unit is [J mol^-1 K^-1]
        deltaHd =  biochem_in.TDP.delHdJ;                % unit is [J mol^-1]
        
        fTv = temperature_functionC3(Tref,Rgas,T,deltaHa);
        fHTv = high_temp_inhibtionC3(Tref,Rgas,T,deltaS,deltaHd);
        
        Jmax = Jmax25 .* fTv .* fHTv;                   % Temperature Corrected Vcmax
        
        % Constant parameters for TPU correction
        deltaHa = biochem_in.TDP.delHaP;                % Unit is  [J K^-1]
        deltaS  = biochem_in.TDP.delSP;                 % unit is [J mol^-1 K^-1]
        deltaHd = biochem_in.TDP.delHdP;                % unit is [J mol^-1]
        
        fTv = temperature_functionC3(Tref,Rgas,T,deltaHa);
        fHTv = high_temp_inhibtionC3(Tref,Rgas,T,deltaS,deltaHd);
        
        TPU = TPU25 .* fTv .* fHTv;                     % Temperature Corrected TPU
        
        
        % Constant parameters for Rd param Correction
        deltaHa   = biochem_in.TDP.delHaR;               % Unit is  [J K^-1]
        deltaS    = biochem_in.TDP.delSR;                % unit is [J mol^-1 K^-1]
        deltaHd   = biochem_in.TDP.delHdR;               % unit is [J mol^-1]
        
        fTv = temperature_functionC3(Tref,Rgas,T,deltaHa);
        fHTv = high_temp_inhibtionC3(Tref,Rgas,T,deltaS,deltaHd);
        
        Rd = Rd25 .* fTv .* fHTv * stressfactor; % Temperature Corrected Rd
        
        % Constant parameters for Kc Correction
        deltaHa = biochem_in.TDP.delHaKc;               % Unit is  [J K^-1]
        
        fTv = temperature_functionC3(Tref,Rgas,T,deltaHa);
        
        Kc = Kc25 .* fTv; % Temperature Corrected Kc
        
        % Constant parameters for Ko Correction
        deltaHa = biochem_in.TDP.delHaKo;               % Unit is  [J K^-1]
        
        fTv = temperature_functionC3(Tref,Rgas,T,deltaHa);
        
        Ko = Ko25 .* fTv;                       % Temperature Corrected Ko
        
        % Constant parameters for Gamma_star Correction
        deltaHa = biochem_in.TDP.delHaT;               % Unit is  [J K^-1]
        
        fTv = temperature_functionC3(Tref,Rgas,T,deltaHa);
        
        Gamma_star = Gamma_star25 .* fTv;       % Temperature Corrected Gamma_star
        
        MM_const = Kc .* (1 + O./Ko);           % These are all temperature adjusted vlues [ppm]
    else % C4 Species
        
        RdPerVcmax25 = 0.025;  % Rd25 for C4 is different than C3
        Rd25 = RdPerVcmax25 * Vcmax25;
        % Constant parameters for temperature correction of Vcmax
        Q10 = biochem_in.TDP.Q10;                           % Unit is  []
        s1  = biochem_in.TDP.s1;                            % Unit is [K]
        s2  = biochem_in.TDP.s2;                            % Unit is [K^-1]
        s3  = biochem_in.TDP.s3;                            % Unit is [K]
        s4  = biochem_in.TDP.s4;                            % Unit is [K^-1]
        
        % Constant parameters for temperature correction of Rd
        
        s5  = biochem_in.TDP.s5;                            % Unit is [K]
        s6  = biochem_in.TDP.s6;                            % Unit is [K^-1]
        
        fHTv = 1 + exp(s1.*(T - s2));
        fLTv = 1 + exp(s3.*(s4 - T));
        Vcmax = (Vcmax25 .* Q10.^(0.1.*(T-Tref)))./(fHTv .* fLTv); % Temp Corrected Vcmax
        
        % Temperature correction of Rd
        
        fHTv = 1 + exp(s5.*(T - s6));
        Rd = (Rd25 .* Q10.^(0.1.*(T-Tref)))./fHTv; % Temp Corrected Rd
        
        % Temperature correction of Ke
        Ke25 = 20000 .* Vcmax25 ;               % Unit is  []
        
        Ke = (Ke25 .* Q10.^(0.1.*(T-Tref)));    % Temp Corrected Ke
        
    end
end

if (tempcor == 0) % If no Temperature Corrections are to be applied 
    if strcmpi('C3', Type)
        Rd25 = RdPerVcmax25 * Vcmax25;
        Vcmax = Vcmax25;
        Jmax = Jmax25;
        TPU = TPU25;
        Rd = Rd25;
        Kc = Kc25;
        Ko = Ko25;
        Gamma_star = Gamma_star25;
        MM_const = Kc .* (1 + O./Ko);
    else % C4 Species
        RdPerVcmax25 = 0.025;  % Rd25 for C4 is different than C3
        Rd25 = RdPerVcmax25 * Vcmax25;
        Rd = Rd25;
        Vcmax = Vcmax25;
        Ke25 = 20000 .* Vcmax25 ;
        Ke = Ke25;
    end   
end

%% PHOTOSYNTHESIS AND ASSIMILATION

% C3 Species Photosynthesis

% Calculation of Initial Value of Ci (internal CO2 concentration)
RH = min(1, eb./satvap(T-273.15) ); % jak: don't allow "supersaturated" air! (esp. on T curves)
warnings = [];

if strcmpi('C3', Type)
    minCi = 0.3;                    % Initial Assumption for internal Leaf Ci concentration relative to Cs
else
    minCi = 0.1;
end
    
    fcount = 0;                         % the number of times we called computeA()
    if  ~isempty(Ci_input)
        Ci = Ci_input;                  % in units of bar.
        if any(Ci_input > 1)
            Ci = Ci_input;              % Assuming Ci input is in ppm
        end
        
    else
        % if b = 0: no need to iterate:
        % Get an intial estimate of Ci
        Ci = BallBerry(Cs, RH, [], BallBerrySlope, BallBerry0, minCi);
        
    end
    
    % Compute the net Photosynthesis using the intial Ci estimate
    A =  computeA(Ci);
    % Recompute the Ci and gs using computed A
    [Ci, gs] = BallBerry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci); % Initial gs computation
       
    %% Incorporate Iterative solution of A-Ci-gs
    
    % Iteration loop for Ci computation
    fCi = 100;
    counter = 0;
    while abs(fCi)>0.0001
        
        [tmpA1, A_out] = computeA(Ci);
        A_final = tmpA1;
        A_out_final = A_out;
        counter = counter+1;
        
        [Ci, gs] = BallBerry(Cs, RH, tmpA1, BallBerrySlope, BallBerry0, minCi, Ci); % Iterative gs computation
        Ci_final = Ci;
        gs_final = gs;
        
        fCi   = newtonfunc_ci(Cs,tmpA1,gs,Ci);

        tmpA2 = computeA(fCi + Ci);
        dCi = (newtonfunc_ci(Cs,tmpA2,gs,fCi + Ci) - fCi)./fCi;
        Ci = Ci - fCi./dCi;
        
    end
    
    A                = A_final; % Final A after the iterative estimation
    Ag               = A_out_final.Ag;
    CO2_per_electron = A_out_final.CO2_per_electron;
    
    Ci = Ci_final;
    %gs = gs_final;
    
    
  %% Compute Assimilation.
    function [A, A_out] = computeA(Ci)
    % global: Type, Vcmax, Gamma_star, MM_consts, Vs_C3, effcon, Je, atheta, Rd
    
    if strcmpi('C3', Type)
        wc = (Vcmax .* (Ci - Gamma_star))./(Ci + MM_const); % This is Rubisco Limiting Step
        
        % Compute the Potential Electron Transport Rate
        thetaPSII = 0.7;            % Curvature parameter
        f = 0.15;
        phiPSII = (1-f);            % Quantum Yeild of PSII
        IPSII = 0.5*phiPSII.*Q;     % light utilized in electron transport by PSII
        
        % Je is the potential electron transport rate computed as the
        % smaller root of the quadratic
        Je = sel_root(thetaPSII, -(IPSII+Jmax), IPSII.*Jmax, -1);
        
        wj = Je.*(Ci-Gamma_star)./(4.*Ci + 8.*Gamma_star);      % This is light limiting step
        
        ws = 3*TPU; % This is the product limiting step (doesn't change with iteration of A-Ci)
        
        thetacj = 0.98;
        thetais = 0.95;
        
        CO2_per_electron = (Ci-Gamma_star)./(Ci + 2.*Gamma_star) .* effcon; % same as previous version
        
        
    else % Compute A for C4
        
        wc = Vcmax;                 % Rubisco Limiting Step
        alpha = 0.05;               % Quantum Yeild [mol mol^-1]
        wj = alpha .* Q;            % Light limiting rate Q is already in mumols m-2 s-1
        ws = Ke .* Ci * 1E-6;       % Product limitation
        
        thetacj = 0.80;
        thetais = 0.95;
        CO2_per_electron = effcon; % same as previous version
        
        
        % Compute Je for C4
        f = 0.15;
        phiPSII = (1-f);            % Quantum Yeild of PSII
        IPSII = 0.5*phiPSII.*Q;     % light utilized in electron transport by PSII
        Je=IPSII;
        
    end
    
    
    % find the smoothed minimum of we, wc and then ws

    wi      = sel_root(thetacj,-(wc+wj),wc.*wj, -1 ); % i.e. sign(Gamma_star - Ci)
    Ag      = sel_root(thetais,-(wi+ws),wi.*ws, -1);
    A       = Ag - Rd;
    
    
    % Compute the Photorespiration (rate of oxygenation vo) from net
    % assimilation
    % The formulations are provided in Sharkey 1988 and Bernachhi et al
    % 2001 Papers
    if strcmpi('C3', Type)
        phi = 2.*Gamma_star./Ci;
        vo = (A + Rd)./(1./phi - 0.5);
    end
    
    
    if nargout > 1
        A_out.A = A;
        A_out.Ag = Ag;
        A_out.Vc = wc;
        A_out.Vs = ws;
        A_out.Ve = wj;
        A_out.Je = Je;
        A_out.CO2_per_electron = CO2_per_electron;
        if strcmpi('C3', Type)
            A_out.vo = vo;
        end
    end
    fcount = fcount + 1; % # of times we called computeA
    
    end
    
    %% Newton Raphson Method
    function fCi = newtonfunc_ci(Cs,An,gs,Ci)
    
    fCi = Cs - 1.6 * An./gs - Ci;
    end
    
    %% Calculation of potential electron transport rate
    if isempty(po0)                      % JAK 2015-12: User can specify po0 from measured data
        po0      = Kp./(Kf+Kd+Kp);       % maximum dark photochemistry fraction, i.e. Kn = 0 (Genty et al., 1989)
    end
    %Je          = 0.5*po0 .* Q .* fPAR;  % potential electron transport rate (JAK: add fPAR)
    Je   = A_out_final.Je;        % Je from the Final A Out after A-Ci-gs Iterations
    
    %% Compute Final Value of gs, rcw and Ja (Actual Electron Transport)
    gs = (1.6 .* A )./ (Cs-Ci);
    
    Ja = Ag ./ CO2_per_electron;        % actual electron transport rate
    
    % stomatal resistance
    rcw      =  (rhoa./(Mair*1E-3))./gs;

    
    %% fluorescence (Replace this part by Magnani or other model if needed)
    ps          = po0.*Ja./Je;               % this is the photochemical yield
    nanPs       = isnan(ps);
    if any(nanPs)
        if numel(po0) == 1
            ps(nanPs) = po0;
        else
            ps(nanPs) = po0(nanPs);  % happens when Q = 0, so ps = po0 (other cases of NaN have been resolved)
        end
    end
    ps_rel   = max(0,  1-ps./po0);       % degree of light saturation: 'x' (van der Tol e.a. 2014)
    
    [eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn]    = Fluorescencemodel(ps, ps_rel, Kp,Kf,Kd,Knparams);
    Kpa         = ps./fs*Kf;
    
    %% CO2 concentration based on the mesophyll conductance if anything other than Inf
    Cc = [];
    if ~isempty(g_m)
        Cc    = (Ci - A./g_m);
    end
    
    %% Collect outputs
    
    biochem_out.A        = A;
    biochem_out.Ag       = A_out_final.Ag;
    biochem_out.Ci       = Ci;
    if ~isempty(Cc)
        biochem_out.Cc   = Cc;
    end
    biochem_out.rcw      = rcw;
    biochem_out.gs       =  gs;
    biochem_out.RH       =  RH;
    biochem_out.warnings = warnings;
    biochem_out.fcount   = fcount;  % the number of times we called computeA()  
    biochem_out.Vcmax    = Vcmax;
    biochem_out.Vc       = A_out_final.Vc;  % export the components of A for diagnostic charts
    biochem_out.Ve       = A_out_final.Ve;
    biochem_out.Vs       = A_out_final.Vs;
    if strcmpi('C3', Type)
        biochem_out.vo       = A_out_final.vo;
    end
    biochem_out.Rd       = Rd;
    biochem_out.CO2_pere = CO2_per_electron;
    
    biochem_out.Ja       = Ja; % Actual Electransport rate
    biochem_out.ps       = ps; % photochemical yield
    biochem_out.ps_rel   = ps_rel;   % degree of ETR saturation 'x' (van der Tol e.a. 2014)
    
    % fluoresence outputs:
    % note on Kn: technically, NPQ = (Fm - Fm')/Fm' = Kn/(Kf + Kd);
    %     In this model Kf + Kd is close to but NOT equal to 1 @ 25C Kf + Kd = 0.8798
    %     vdT 2013 fitted Kn assuming NPQ = Kn, but maybe we shouldn't?
    biochem_out.Kd      = Kd;  % K_dark(T)
    biochem_out.Kn      = Kn;  % K_n(x);  x = 1 - ps/p00 == 1 - Ja/Je
    biochem_out.NPQ     = Kn ./ (Kf + Kd); % why not be honest!
    biochem_out.Kf      = Kf;  % Kf = 0.05 (const)
    biochem_out.Kp0     = Kp;  % Kp = 4.0 (const): Kp, max
    biochem_out.Kp      = Kpa; % Kp,actual
    biochem_out.eta     = eta;
    biochem_out.qE      = qE;
    biochem_out.fs      = fs;  % keep this for compatibility with SCOPE
    biochem_out.ft      = fs;  % keep this for the GUI ft is a synonym for what we're calling fs
    biochem_out.SIF     = fs .* Q;
    biochem_out.fo0     = fo0;
    biochem_out.fm0     = fm0;
    biochem_out.fo      = fo;
    biochem_out.fm      = fm;
    biochem_out.Fm_Fo    = fm ./ fo;  % parameters used for curve fitting
    biochem_out.Ft_Fo    = fs ./ fo;  % parameters used for curve fitting
    biochem_out.qQ      = qQ;
    return;
    
end  % end of function biochemical


%% quadratic formula, root of least magnitude
    function x = sel_root(a,b,c, dsign)
        %  sel_root - select a root based on the fourth arg (dsign = discriminant sign)
        %    for the eqn ax^2 + bx + c,
        %    if dsign is:
        %       -1, 0: choose the smaller root
        %       +1: choose the larger root
        %  NOTE: technically, we should check a, but in biochemical, a is always > 0
        if a == 0  % note: this works because 'a' is a scalar parameter!
            x      = -c./b;
        else
            if any(dsign == 0)
                dsign(dsign == 0) = -1; % technically, dsign==0 iff b = c = 0, so this isn't strictly necessary except, possibly for ill-formed cases)
            end
            %disc_root = sqrt(b.^2 - 4.*a.*c); % square root of the discriminant (doesn't need a separate line anymore)
            %  in MATLAB (2013b) assigning the intermediate variable actually slows down the code! (~25%)
            x = (-b + dsign.* sqrt(b.^2 - 4.*a.*c))./(2.*a);
        end
    end %of min_root of quadratic formula


%% Ball Berry Model
    function [Ci, gs] = BallBerry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input)
        %  Cs  : CO2 at leaf surface
        %  RH  : relative humidity
        %  A   : Net assimilation in 'same units of CO2 as Cs'/m2/s
        %  BallBerrySlope, BallBerry0,
        %  minCi : minimum Ci as a fraction of Cs (in case RH is very low?)
        %  Ci_input : will calculate gs if A is specified.
        if nargin > 6 && ~isempty(Ci_input)
            % Ci is given: try and compute gs
            Ci = Ci_input;
            gs = [];
            if ~isempty(A) && nargout > 1
                gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
            end
        elseif all(BallBerry0 == 0) || isempty(A)
            % EXPLANATION:   *at equilibrium* CO2_in = CO2_out => A = gs(Cs - Ci) [1]
            %  so Ci = Cs - A/gs (at equilibrium)                                 [2]
            %  Ball-Berry suggest: gs = m (A RH)/Cs + b   (also at equilib., see Leuning 1990)
            %  if b = 0 we can rearrange B-B for the second term in [2]:  A/gs = Cs/(m RH)
            %  Substituting into [2]
            %  Ci = Cs - Cs/(m RH) = Cs ( 1- 1/(m RH)  [ the 1.6 converts from CO2- to H2O-diffusion ]
            Ci      = max(minCi .* Cs,  Cs.*(1-1.6./(BallBerrySlope .* RH)));
            gs = [];
        else
            %  if b > 0  Ci = Cs( 1 - 1/(m RH + b Cs/A) )
            % if we use Leuning 1990, Ci = Cs - (Cs - Gamma)/(m RH + b(Cs - Gamma)/A)  [see def of Gamma, above]
            % note: the original B-B units are A: umol/m2/s, ci ppm (umol/mol), RH (unitless)
            %   Cs input was ppm but was multiplied by ppm2bar above, so multiply A by ppm2bar to put them on the same scale.
            %  don't let gs go below its minimum value (i.e. when A goes negative)
            gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0);
            Ci = max(minCi .* Cs,  Cs - 1.6 * A./gs) ;
        end
        
    end % function

    function gs = gsFun(Cs, RH, A, BallBerrySlope, BallBerry0)
        % add in a bit just to avoid div zero. 1 ppm = 1e-6 (note since A < 0 if Cs ==0, it gives a small gs rather than maximal gs
        gs = max(BallBerry0,  BallBerrySlope.* A .* RH ./ (Cs+1e-9)  + BallBerry0);
        gs( isnan(Cs) ) = NaN; 
    end


%% Fluorescence model
% NEEDS Little Cleaning
    function [eta,qE,qQ,fs,fo,fm,fo0,fm0,Kn] = Fluorescencemodel(ps,x, Kp,Kf,Kd,Knparams)
        % note: x isn't strictly needed as an input parameter but it avoids code-duplication (of po0) and it's inherent risks.
        
        Kno = Knparams(1);
        alpha = Knparams(2);
        beta = Knparams(3);
        
        % switch model_choice
        %     case 0, % drought
        %         Kno = 5.01;
        %         alpha = 1.93;
        %         beta = 10;
        %         %Kn          = (6.2473 * x - 0.5944).*x; % empirical fit to Flexas' data
        %         %Kn          = (3.9867 * x - 1.0589).*x;  % empirical fit to Flexas, Daumard, Rascher, Berry data
        %     case 1, healthy (cotton)
        %         Kno = 2.48;
        %         alpha = 2.83;
        %         beta = 0.114;
        %         %p = [4.5531;8.5595;1.8510];
        %         %Kn   = p(1)./(p(3)+exp(-p(2)*(x-.5)));
        % end
        
        % using exp(-beta) expands the interesting region between 0-1
        %beta = exp(-beta);
        x_alpha = exp(log(x).*alpha); % this is the most expensive operation in this fn; doing it twice almost doubles the time spent here (MATLAB 2013b doesn't optimize the duplicate code)
        Kn = Kno * (1+beta).* x_alpha./(beta + x_alpha);
        
        %Kn          = Kn .* Kd/0.8738;          % temperature correction of Kn similar to that of Kd
        
        fo0         = Kf./(Kf+Kp+Kd);        % dark-adapted fluorescence yield Fo,0
        fo          = Kf./(Kf+Kp+Kd+Kn);     % light-adapted fluorescence yield in the dark Fo
        fm          = Kf./(Kf   +Kd+Kn);     % light-adapted fluorescence yield Fm
        fm0         = Kf./(Kf   +Kd);        % dark-adapted fluorescence yield Fm
        fs          = fm.*(1-ps);            % steady-state (light-adapted) yield Ft (aka Fs)
        eta         = fs./fo0;
        qQ          = 1-(fs-fo)./(fm-fo);    % photochemical quenching
        qE          = 1-(fm-fo)./(fm0-fo0);  % non-photochemical quenching
        
        
        %eta         = eta*(1+5)/5 - 1/5;     % this corrects for 29% PSI contribution in PAM data, but it is quick and dirty correction that needs to be improved in the next
        
    end



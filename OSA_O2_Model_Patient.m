% OSA model code for patient case
clear; clc; close all;

% Physiological parameters
HRi = 73; %beats/min (heart rate)
SV = 70; %mL/beat (stroke volume)
CO = HRi*SV/60; %mL/s (cardiac output)
Qpi = CO; %mL/min (pulmonary flow)
Qsi = CO; %mL/min (systemic flow)
R = 62360; %mmHg mL/mol K (ideal gas constant)
T = 310; %K (body temperature)
VT = 500; %mL (tidal volume)
br = 12/60; % breaths/s (breathing rate)
yO2 = 0.21; % (mol fraction of O2 in inspired air)
PB = 760; % mmHg (barometric pressure)
PH2O = 47; % mmHg (vapor pressure at 37C, saturation)
VD = 150; %mL (dead space volume)
Beta_p = 1.4e-9; % mol O2/mL blood mmHg (solubility of O2 in plasma) (Frank et al. A finite element method of oxygen diffusion in the pulmonary capillaries)
CHb = 2.3e-6; %mol/mL (concentration of hemoglobin in blood) (Source: Erratum to: Blood HbO2 and HbCO2 dissociation curves at varied O2, CO2, pH, 2,3-DPG and temperature levels)
nV = 1/22400; %mol/mL (Conversion factor at Body temperature)
PAO2i = 104; %mmHg (initial pulmonary venous oxygen partial pressure)
Vsys = 4600; %mL (volume of blood in systemic circulation)
Vsysa = 770; %mL (volume of blood in systemic arteries + aorta) 
Vsyscap = 330; %mL (volume of blood in systemic capillaries)
tcyclei = Vsys/Qsi; %s (time for blood to pass through systemic circulation)
tsysai = Vsysa/Qsi; %s (time for blood to pass from exit of lungs to systemic capillaries)
tsyscapi = Vsyscap/Qsi; %s (time for blood to pass through systemic capillaries)
DLO2 = 21/60; %mLO2/s mmHg
Aeff = 3000; %cm/3 (total pulmonary capillary cross-sectional area)
Vpc = 130; %mL (volume of blood in pulmonary capillaries)
vpci = Qpi/Aeff;
Lc = Vpc/Aeff; %cm (length of pulmonary capillary compartment)
A = 122*100^2; %cm^2 (total surface area for gas exchange)
l_mem = 2e-4; % cm (alveolar-capillary membrane thickness)
MRO2 = 240*nV/60; % molO2/s (basal metabolic rate for all tissues)
DLO2_b = (DLO2*nV)/Beta_p;
D_O2 = DLO2_b*l_mem/A;
k_l = D_O2*A/l_mem;
pars = [yO2, PB, PH2O, R, T, Beta_p, k_l];


% Lookup table approximation for oxygenation
Cd_T = (linspace(0,2e-6,8000))';
S_T = SatO2_2(Cd_T,Beta_p);
CT_T = 4*CHb.*S_T + Cd_T;

% Load breathing pattern file
prompt = "Enter file name for breathing pattern";
file = input(prompt,'s');
Input = load(file,'t_vec','L_Vol','HR_vec');
C = struct2cell(Input);
t_Data = C{1};
V_Data = C{2};
HR_Data = C{3};

% Solving parameters 
dti = 0.007; % Time step for each iteration of for loop
dzi = vpci*dti; % Differential distance for capillary profile IC
k = floor(Lc/dzi);
n_dti = floor(tcyclei/dti);
a_dti = floor(tsysai/dti);
c_dti = floor(tsyscapi/dti);
prompt2 = "Enter total time for data in minutes";
total_t = input(prompt2); % Total time in seconds
s = floor(total_t*60/dti);

% Initialize variable vectors
t = zeros(1,s);
PAO2 = zeros(1,s);  CAO2 = zeros(1,s);
PpaO2 = zeros(1,s); CpaO2d = zeros(1,s); CpaO2T = zeros(1,s); SpaO2 = zeros(1,s);
PpvO2 = zeros(1,s); CpvO2d = zeros(1,s); CpvO2T = zeros(1,s); SpvO2 = zeros(1,s);
PsaO2 = zeros(1,s); CsaO2d = zeros(1,s); CsaO2T = zeros(1,s); SsaO2 = zeros(1,s);
PsvO2 = zeros(1,s); CsvO2d = zeros(1,s); CsvO2T = zeros(1,s); SsvO2 = zeros(1,s);
PcO2 = zeros(1,s); CcO2 = zeros(1,s);
HR = zeros(1,s); Qp = zeros(1,s); Qs = zeros(1,s);
CpcO2T = zeros(k+1,s); CpcO2d = zeros(k+1,s); SpcO2 = zeros(k+1,s); PpcO2 = zeros(k+1,s);
z1 = zeros(k+1,1); vpc = zeros(1,s);
dt = zeros(1,s); n_dt = zeros(1,s); a_dt = zeros(1,s); c_dt = zeros(1,s); 
tcycle = zeros(1,s); tsyscap = zeros(1,s); tsysa = zeros(1,s);
CpcO2T1 = zeros(k+1,1); CpcO2d1 = zeros(k+1,1);
V = zeros(1,s);

% Initialize variables (t=0)
PAO2(1) = PAO2i; % mmHg (oxygen partial pressure in alveoli)
CAO2(1) = PAO2(1)*Beta_p; % mol O2/mL (oxygen concentration in alveoli membrane)

PpvO2(1) = PAO2i; %mmHg (oxygen partial pressure in pulmonary vein)
CpvO2d(1) = PpvO2(1)*Beta_p; % mol O2/mL blood (dissolved oxygen in pulmonary vein)
SpvO2(1) = SatO2_2(CpvO2d(1),Beta_p); % (saturation of pulmonary venous blood)
CpvO2T(1) = 4*CHb*SpvO2(1) + CpvO2d(1);

CpaO2T(1:n_dti) = CpvO2T(1)-MRO2/Qsi; % (total oxygen concentration in pulmonary artery)
CpaO2d(1:n_dti) = interp1(CT_T,Cd_T,CpaO2T(1)); % mol O2/mL blood (dissolved oxygen in pulmonary artery);
PpaO2(1:n_dti) = CpaO2d(1)/Beta_p; %mmHg (oxygen partial pressure in pulmonary artery)
SpaO2(1:n_dti) = SatO2_2(CpaO2d(1),Beta_p); % (saturation of pulmonary arterial blood)

PsaO2(1:a_dti+1) = PpvO2(1);
CsaO2T(1:a_dti+1) = CpvO2T(1);
CsaO2d(1:a_dti+1) = CpvO2d(1);
SsaO2(1:a_dti+1) = SpvO2(1);

PsvO2(1:c_dti+1) = PpaO2(1);
CsvO2T(1:c_dti+1) = CpaO2T(1);
CsvO2d(1:c_dti+1) = CpaO2d(1);
SsvO2(1:c_dti+1) = SpaO2(1);

for n=2:k+1
    CpcO2T1(1) = CpaO2T(1);
    CpcO2d1(1) = CpaO2d(1);
    z1(n) = z1(n-1)+dzi;
    CpcO2T1(n) = CpcO2T1(n-1)+(k_l/Qpi)*(dzi/Lc)*(CAO2(1)-CpcO2d1(n-1));
    CpcO2d1(n) = interp1(CT_T,Cd_T,CpcO2T1(n));
end

CpcO2T(:,1) = CpcO2T1; 
CpcO2d(:,1) = CpcO2d1;
PpcO2(:,1) = CpcO2d(:,1)/Beta_p;
SpcO2(:,1) = SatO2_2(CpcO2d(:,1),Beta_p);


CcO2(1) = trapz(dzi,CpcO2d(:,1))./Lc; 
PcO2(1) = CcO2(1)/Beta_p; 

V(1) = 3400; 
t(1) = 0;
HR(1) = HRi; Qp(1) = Qpi; Qs(1) = Qsi; vpc(1) = vpci;
dt(1) = dti;
n_dt(1) = n_dti; c_dt(1) = c_dti; a_dt(1) = a_dti;
tcycle(1) = tcyclei; tsysa(1) = tsysai; tsyscap(1) = tsyscapi;

for i=2:s
            %Time
            t(i) = t(i-1)+dt(i-1);
            
            %Alveolar Volume
            V(i) = interp1(t_Data,V_Data,t(i));

            %Heart Rate and Flow
            HR(i) = interp1(t_Data,HR_Data,t(i));
            Qp(i) = HR(i)*SV/60;
            Qs(i) = HR(i)*SV/60;

            %Times
            tcycle(i) = Vsys/Qs(i);
            tsyscap(i) = Vsyscap/Qs(i);
            tsysa(i) = Vsysa/Qs(i);
            dt(i) = tcycle(i)*dti/tcyclei;
            n_dt(i) = floor(tcycle(i)/dt(i));
            c_dt(i) = floor(tsyscap(i)/dt(i));
            a_dt(i) = floor(tsysa(i)/dt(i));
            vpc(i) = Qp(i)/Aeff;
           
            %Calculate arterial oxygenation
            CpaO2T(i+n_dt(i)-1) = CpvO2T(i-1)-MRO2/Qs(i-1);
            CpaO2d(i) = interp1(CT_T,Cd_T,CpaO2T(i));
            PpaO2(i) = CpaO2d(i)/Beta_p;
            SpaO2(i) = SatO2_2(CpaO2d(i),Beta_p);
            
            % Pulmonary capillary oxygenation
            CpcO2T(1,i) = CpaO2T(i);
            CpcO2T(2:end,i) = CpcO2T(1:(end-1),i-1) + (k_l/Qp(i-1))*(dzi/Lc)*(CAO2(i-1)-CpcO2d(1:(end-1),i-1));
            CpcO2d(:,i) = interp1(CT_T,Cd_T,CpcO2T(:,i));
            PpcO2(:,i) = CpcO2d(:,i)/Beta_p;
            SpcO2(:,i) = SatO2_2(CpcO2d(:,i),Beta_p);

            % Venous oxygenation
            CpvO2T(i) = CpcO2T(end,i);
            CpvO2d(i) = interp1(CT_T,Cd_T,CpvO2T(i));
            PpvO2(i) = CpvO2d(i)/Beta_p;
            SpvO2(i) = SatO2_2(CpvO2d(i),Beta_p);
          
            % Systemic capillaries input
            CsaO2T(i+a_dt(i)) = CpvO2T(i);
            CsaO2d(i) = interp1(CT_T,Cd_T,CsaO2T(i));
            PsaO2(i) = CsaO2d(i)/Beta_p;
            SsaO2(i) = SatO2_2(CsaO2d(i),Beta_p);
            
            % Systemic capillaries output
            CsvO2T(i+c_dt(i)) = CsaO2T(i)-MRO2/Qs(i);
            CsvO2d(i) = interp1(CT_T,Cd_T,CsvO2T(i));
            PsvO2(i) = CsvO2d(i)/Beta_p;
            SsvO2(i) = SatO2_2(CsvO2d(i),Beta_p);
    
            % Alveolar oxygen partial pressure
            CcO2(i) = trapz(dzi,CpcO2d(:,i))./Lc;
            PcO2(i) = CcO2(i)/Beta_p;
            PAO2(i) = Alv_MassBalance_P(V(i),V(i-1),t(i),t(i-1),PcO2(i-1),PAO2(i-1),pars);
            CAO2(i) = PAO2(i)*Beta_p; 
end

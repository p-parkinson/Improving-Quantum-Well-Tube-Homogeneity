function [parameters,avspec,curve] = PL_LSWfit(settings,watershed,regions_found,regions,flatdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to fit PL spectra with Lasher-Stern-Wurfe Model (LSW) for a
% band-to-band PL peak and output fit parameters
%%% ____________________________ (Inputs) _____________________________ %%%
% settings = PL map measurement settings (for wavength x-axis in PL Spectrum)
% watershed = watershed transform from "PL_watershed.m"
% regions_found = no. regions found from watershed transform of PL map.
% regions = numbered regions in PL map
% flatdata = PL Spectrum data (reshaped to correct dimension)
%%% ____________________________ (Output) _____________________________ %%%
% parameters = parameters from LSW fitting
%     o QW_Width* = QW Width (to be overwritten with interpolation from
%                   sophisticated nextnano simulations in further analysis)
%     o Transition_Energy* = Peak energy of LSW fit
%     o Temperature* = Carrier emission temperature
%     o Gamma* = Urbach energy
%     o Sigma* = FWHM of Gaussian component
%     o Fermi_Energy* = Fermi energy
%     o Intensity = Peak intensity (counts, arb. u)
%     o cow* = (c)entre (o)f (w)avelength
%     o rms_r = residuals of LSW fitting
%         *Includes +/- error from parameter Jacobian 
% avspec = averaged PL spectrum per region
% curve = LSW fitting curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function
opts=optimset('Display','off'); %suppress fitting messages
% Import dependencies
x = settings.wl; 
w = watershed;

% Constants
h = 6.626*10^(-34);
c = 2.998*10^(8);
eV = 1.602*10^(-19);
k_B=8.6173324*1e-5; %eV/K

% LSW Constants
hbar=1.054571*1e-34;
m0=9.10938356*1e-31;
me = 0.063; % i.e. factor of m0
mh = 0.34; % i.e. factor of m0
E_GaAs = 1.42; %(Initial Fitting Parameter)

% Convert wavelength (nm) x axis to energy domain (eV)
E = h*c./(eV*x.*10^(-9));

%% Fitting Equations
transition_energy = @(par)(E_GaAs + hbar^2*pi^2/(2*me*m0*(par(1)*1E-9)^2*(eV))...
    + hbar^2*pi^2/(2*mh*m0*(par(1)*1E-9)^2*(eV)));

%LSW model
%2D dos
urb_LSW1 =@(par,E) ((E<=transition_energy(par))  .*exp((E-transition_energy(par))/par(3))); %Urbach tail par(3) is characteristic width of peak function
B_QWT1 = @(par,E) ((E>transition_energy(par))  .*1) ;
%disorder term
Gauss=@(par,E) (exp(-((E-E(round(end/2))).^2)/(2*(abs(par(4))^2))));
%define boltzmann term
Boltz1=@(par,E)(exp((E-(par(5)+transition_energy(par)))/(k_B*par(2))));
%LSW
I_num1 = @(par,E)(E>transition_energy(par)).*E.^2.*(1-exp(B_QWT1(par,E))).*(1-2./(Boltz1(par,E).^0.5+1));
I_den1 = @(par,E)1./(Boltz1(par,E)-1);
I_frac1 = @(par,E)I_num1(par,E).* I_den1(par,E);
comb1 = @(par,E) (par(6).*urb_LSW1(par,E) - I_frac1(par,E));

%Planck-Einstein relationship for LSW (FINAL EQUATION)
I = @(par,E) par(7) .* (conv(Gauss(par,E),comb1(par,E),'same'));

%Array to store parameters
num_params = 7;
f_prms = zeros(regions_found,num_params);
rms_r = zeros(regions_found,1);

%Array for errors from Jacobian of fitting parameters
jac_l = zeros(regions_found,num_params);
jac_u = zeros(regions_found,num_params);

%Initialise arrays for Peak Intensity, average PL spectrum per region
%and LSW fitting curve
Intensity = zeros(regions_found,1);
avspec = zeros(regions_found,1024);
curve = zeros(regions_found,1024);
    %% Fitting Loops (lsqcurvefit)
    for i = regions
        % Find mask (region of an individual structure)
        mask = (w(:)==i);
            
        % Select the spectra for that region
        spec = flatdata(mask,:);
        % Average over the region
        avspec(i,:) = sum(spec,1);

        %outlier data points removal and normalise
        avspec(i,:) = filloutliers(avspec(i,:),'nearest','movmedian',20);
        avspec(i,:) = avspec(i,:)./trapz(settings.wl,avspec(i,:)); %normalise spectra
        
        %prm name:  L                            T                                gamma     sigma       Fermi_Energy    Urb_Scale          scale            ];
        %prms    = [prm(1)                       prm(2)                           prm(3)    prm(4)      prm(5)          parm(6)            prm(7)           ];
        LB       = [6                            300                              0.00      0.004       0.0             1.0                0                ];
        par0     = [7                            325                              0.02      0.010       0.1             1.5                max(avspec(i,:)) ];
        UB       = [8                            Inf                              0.10      0.300       0.2             2.5                2*max(avspec(i,:)) ];

        %[params,resnorm,residual,output,lambda,jacobian] = lsqcurvefit(I,par0,E,avspec(i,:),LB,UB,opts);
        [params,resnorm,~,~,~,jacobian] = lsqcurvefit(I,par0,E,avspec(i,:),LB,UB,opts);
        f_prms(i,:) = params;
        Intensity(i) = max(I(params,E));
        
        jac_l(i,:) = jacobian.lower;
        jac_u(i,:) = jacobian.upper;
        rms_r(i) = sqrt(resnorm); %RMS of Residual

        curve(i,:) = I(params,E);
        clear params
    end
    %% Sort Parameters
    QW_Width = struct('val',f_prms(:,1),'yneg',jac_l(:,1),'ypos',jac_u(:,1));
    Temp = struct('val',f_prms(:,2),'yneg',jac_l(:,2),'ypos',jac_u(:,2)); 
    Gamma = struct('val',f_prms(:,3),'yneg',jac_l(:,3),'ypos',jac_u(:,3));
    Sigma = struct('val',f_prms(:,4),'yneg',jac_l(:,4),'ypos',jac_u(:,4));
    E_f = struct('val',f_prms(:,5),'yneg',jac_l(:,5),'ypos',jac_u(:,5));

    E_g = zeros(regions_found,3);
    E_g(:,1) = arrayfun(transition_energy,f_prms(:,1));
    transition_energy_err = @(par) (hbar^2*pi^2*eV/m0).*(     (1./(me*(par(1)*1E-9).^3))  +   (1./(mh*(par(1)*1E-9).^3))   );
    E_g(:,2) = jac_l(:,1).*arrayfun(transition_energy_err,f_prms(:,1));
    E_g(:,3) = jac_u(:,1).*arrayfun(transition_energy_err,f_prms(:,1));
    
    cow = zeros(regions_found,3);
    cow(:,1) = (10^(9)*h*c) ./ (eV.*E_g(:,1)); %might need to add transition_energy here
    cow(:,2) = E_g(:,2).*(10^(9)*h*c) ./ (eV.*E_g(:,1));
    cow(:,3) = E_g(:,3).*(10^(9)*h*c) ./ (eV.*E_g(:,1));
    
    cow = struct('val',cow(:,1),'yneg',cow(:,2),'ypos',cow(:,3));
    E_g = struct('val',E_g(:,1),'yneg',E_g(:,2),'ypos',E_g(:,3));

    parameters = struct('QW_Width',QW_Width,'Transition_Energy',E_g,'Temperature',Temp,'Gamma',Gamma,'Sigma',Sigma,'Fermi_Energy',E_f,'Intensity',Intensity,'cow',cow,'rms_r',rms_r);
end
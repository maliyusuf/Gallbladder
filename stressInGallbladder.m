function [hGB sigma_max_p] = stressInGallbladder(PatientNumber, PercentChangeInhGB, PercentageOfLesion, plt)
% function stressInGallbladder(PatientNumber,change, plt)
%  plt = 0 for no plot, 1 for 3D plot
%  'change' is in percentage
%  The complete working function to plot graphs shown in the paper
%  Li - A Mechanical Model for CCK-Induced Acalculous Gallbladder Pain - 2011
%  Written by M. Ali Yousuf
%  Johns Hopkins
%  Ver 1.0, Feb 17, 2012
%  Ver 2.0 with partial ellipsoids, Feb 27, 2012

cd('E:\Dropbox\Research - Private Folders\utilities\Gallbladder');

[num,txt] = xlsread('data2010.xlsx');  % Read Table 1 of paper
% txt = PatientNo	k1	k2	D3 (mm)	EF(%)	Peak Pressure(mmHg)	Peak Total
% Stress(mmHg)	Peak passive stress (mmHg)	Peak active stress (mmHg)	
% te (min) texp (min)	Vexp (mL) 

I = PatientNumber;

hGB = 2.5;  % Gallbladder wall thickness in mm 
            % Reference: Engel - Gallbladder Wall Thickness - Sonographic
            % Accuracy and Relation to Disease - 1980 

hGB = (1+PercentChangeInhGB/100)*hGB;

ti = 5;     % Time in which the CKK reaches the GB, 5 min in the paper
            % Li - A Mechanical Model for CCK-Induced Acalculous Gallbladder Pain - 2011
te = num(I,10);    % GB emptying time in min, a subject dependent quantity, 
                   % Some exact values shown in Fig 4a, of 
            % Li - A Mechanical Model for CCK-Induced Acalculous Gallbladder Pain - 2011
            % It requres a special function to calculate in the model shown in paper
texp = num(I,11);  % specific time during emtptying at which Vexp is measured
pe = 11;    % gallbladder internal pressure while emptying, 11 mmHg in the paper
            % Ref: Majeed - Continuous in vivo manometry of the human gallbladder - 1994
pd = 6;     % pressure in the duodenum, 6 mm Hg in the paper: Li - A Mechanical 
            % Model for CCK-Induced Acalculous Gallbladder Pain - 2011 
C = 2.731;  % Wall compliance (the inverse of stiffness), in mL/mmHg, taken from the paper
            % Middelfart, H.V., Jensen, P., Hojgaard, L. and Funch-Jensen,
            % P., 1998, Pain patterns after distension of the gallbladder
            % in patients with acute cholecystis, Scandinavica Journal of
            % Gastroenterology, 33, 982–987.   
            %%%%% The above value will change for a Lesion. It should be
            %%%%% smaller for stiffer wall
%C = (1-PercentChangeInhGB/100)*C;

% Refer to Li - Correlation of mechanical factors and gallbladder pain - 2008
% Table 1, page 9 for a detialed list of the following variables for 51
% patients.  
% D1 = 20;      % Diameter along the first minor axis, in mm
% D2 = 25;      % Diameter along the second minor axis, in mm 
D3 = num(I,4)/5;  % Diameter along the major axis at time t, in mm
k1 = num(I,2);  % k1 = D1/D3 ;    
k2 = num(I,3);  % k2 = D2/D3;
                % The above values taken from page 791 of Li - A Mechanical Model
                % for CCK-Induced Acalculous Gallbladder Pain - 2011 
D1 = k1*D3/2;
D2 = k2*D3/2;

%Volume_mm3 = (4/3) * pi * D1/2 * D2/2 * D3/2; % Volume of the ellipsoid in mm^3
Volume_mm3 = (1/6) * pi * k1 * k2 * power(D3,3); % Volume of the ellipsoid in mm^3

V0 = Volume_mm3 / 1000; % Volume of the gallbladder AFTER refilling
                    % Converted to ml, using 1000 mm^3 = 1 mL
VB = V0;            % Volume of the gallbladder after CKK intravenous infusion, 
                    % at point B in the graph
                    % Li - A Mechanical Model for CCK-Induced Acalculous Gallbladder Pain - 2011
Vexp = num(I,12);  % Volume measured at a specific time 'texp' during emtying 

D1B = D1;   % Diameter along the first minor axis, in mm at time ti
D2B = D2;   % Diameter along the second minor axis, in mm at time ti
D3B = D3;   % Diameter along the major axis at time t, in mm at time ti

dCBD = 10; % Mean diameter of the bile duct, in mm
           % Ferris, D., and J. Vibert. The common bile duct: significance
           % of its diameter. Ann. Surg. 149(2):249–251, 1959. 
hCBD = 1;  % Mean thickness of the bile duct, in mm
           % Mahour, G., K. Wakim, et al. The common bile duct in man: its
           % diameter and circumference. Ann. Surg. 165(3): 415, 1967. 
pCBD = 35; % mmHg average pressure threshold for pain
           % Ref: Gaensler, E. Quantitative determination of the visceral 
           % pain threshold in man. J. Clin. Invest. 30:406, 1951. 

%*********************Some definitions****************************

Ve = 0.3*V0;    % Volume when the emptying is finsihed (patient dependent)

k1B = D1B/D3B ; k2B = D2B/D3B;

gamm = Ve - C * (pe - pd);
R = (texp - ti) / (C * log((V0-gamm)/(Vexp - gamm)));% Flow resistance
pB=pd + (pe - pd)*exp(te / (R*C));
PainThreshold = pCBD * dCBD / (2*hCBD);

%******************Start Calculating Stress***********************
%*****************************************************************

n = 50;     % The mesh size - higher values generate finer mesh
cut = PercentageOfLesion*2;   % Percentage of the GB conveted to lesion 
cut = (100-cut) / 100;   
cut = asin(cut);
%[xx,yy,zz]= ellipsoido(D1,D1,D1,n); % generates three
[xx,yy,zz]= cutellipsoido(D1,D1,D1,cut,n); 

sphi_a = []; stheta_a = []; tauthetaphi_a = [];
sphi_p = []; stheta_p = []; tauthetaphi_p = [];

% The correct setup
% theta = (-pi:2*pi/n:pi);
% phi   = (-pi/2 : 2*pi/2/n : cut)';

for theta = -pi : 2*pi/n : pi;
    s1_a = []; s2_a = []; s3_a = [];
    s1_p = []; s2_p = []; s3_p = [];
        for phi   = -pi/2 : 2*pi/2/n : cut;
            pp = pe;
            [sigma_theta_p, sigma_phi_p, tau_theta_phi_p] = fStress(pp, hGB, D3,  k1,  k2,  theta, phi);
            pp = pB - pe;
            [sigma_theta_a, sigma_phi_a, tau_theta_phi_a] = fStress(pp, hGB, D3B, k1B, k2B, theta, phi);
            s1_a = [s1_a; sigma_theta_a];
            s2_a = [s2_a; sigma_phi_a];
            s3_a = [s3_a; tau_theta_phi_a];
            
            s1_p = [s1_p; sigma_theta_p];
            s2_p = [s2_p; sigma_phi_p];
            s3_p = [s3_p; tau_theta_phi_p];
        end
       stheta_a = [stheta_a, s1_a];
       sphi_a = [sphi_a, s2_a];
       tauthetaphi_a = [tauthetaphi_a, s3_a];

       stheta_p = [stheta_p, s1_p];
       sphi_p = [sphi_p, s2_p];
       tauthetaphi_p = [tauthetaphi_p, s3_p];

end

%****Calculate Peak Stress and Compare with Threshold of Pain*****
%*****************************************************************
sigma_theta   = stheta_p       + stheta_a; 
sigma_phi     = sphi_p         + sphi_a; 
tau_theta_phi = tauthetaphi_p  + tauthetaphi_a;

sigma_max_t = max( max(max(sigma_theta))   , max(max(sigma_phi)  ) );
sigma_max_p = max( max(max(stheta_p   ))   , max(max(sphi_p)  ) );
sigma_max_a = max( max(max(stheta_a   ))   , max(max(sphi_a)  ) );

% home
% disp('              ');
% disp(['For Patient ', num2str(I),':']);
% disp('===============');
% disp(['With C = ', num2str(C),' and hGB = ',num2str(hGB),' we get (a change of ', num2str(change), ' percent)'])
% disp(['The peak pressure = ', num2str(pB),' mmHg, Compare paper value = ',num2str(num(I,6)),' mmHg'])
% disp(['The peak total stress = ',   num2str(round(sigma_max_t)),' mmHg, Compare paper value = ', num2str(num(I,7)),' mmHg'])
% disp(['The peak passive stress = ', num2str(round(sigma_max_p)),' mmHg, Compare paper value = ', num2str(num(I,8)),' mmHg'])
% disp(['The peak active stress = ',  num2str(round(sigma_max_a)),' mmHg, Compare paper value = ', num2str(num(I,9)),' mmHg'])
% disp(['Threshold of Pain = ',num2str(PainThreshold),' mmHg'])

schange = abs( 137.7884-sigma_max_p) * 100 /  137.7884;

% disp('              ');
% disp(['With C = ', num2str(C),' and hGB = ',num2str(hGB),' we get (a change of ', num2str(change), ...
%       ' percent), Peak stress = ', num2str(round(sigma_max_p)),' mmHg (a change of ', num2str(schange),'%)'])
%[PercentageOfLesion sigma_max_p schange]

%******************Start Plotting******************************
%**************************************************************

    if plt ==1
    plotCombined(xx, yy, zz, D1, D2, D3, sphi_p, stheta_p, tauthetaphi_p, sphi_a, stheta_a, tauthetaphi_a)
   % plotSeparate(xx, yy, zz, D1, D2, D3, sphi_p, stheta_p, tauthetaphi_p, sphi_a, stheta_a, tauthetaphi_a)
    end

end 

function [xxx,yyy,zzz] = gall(xr,yr,zr,n)
% This function is needed to plot the ellipsoid in spherical coordinates.
% That will allow us to chop off its head when needed.

xxx = []; yyy = []; zzz = [];

    for theta = -pi/2 : pi/n : pi/2;
        xx = []; yy = []; zz = [];
        
         for phi  = -pi : 2*pi/n : pi;
            x = xr*cos(phi)*cos(theta);
            y = yr*cos(phi)*sin(theta);
            z = zr*sin(phi)*ones(1,n+1);
            xx = [xx; x];
            yy = [yy; y];
            zz = [zz; z];
         end 
            xxx = [xxx; xx];
            yyy = [yyy; yy];
            zzz = [zzz; zz];   
     
    end 
 
end

function [xx,yy,zz]= ellipsoido(xr,yr,zr,n)
%ELLIPSOIDO Generate ellipsoid.
%   [X,Y,Z]=ELLIPSOID(XC,YC,ZC,XR,YR,ZR,N) generates three
%   (N+1)-by-(N+1) matrices so that SURF(X,Y,Z) produces an
%   ellipsoid with center (XC,YC,ZC) and radii XR, YR, ZR.

[x,y,z] = sphere(n);

xx = xr*x;
yy = yr*y;
zz = zr*z;

end 

function [xx,yy,zz]= cutellipsoido(xr,yr,zr,cut,n)
%ELLIPSOIDO Generate ellipsoid.
%   [X,Y,Z]=ELLIPSOID(XC,YC,ZC,XR,YR,ZR,N) generates three
%   (N+1)-by-(N+1) matrices so that SURF(X,Y,Z) produces an
%   ellipsoid with center (XC,YC,ZC) and radii XR, YR, ZR.
%PercentageOfLesion = 10;

theta = (-pi:2*pi/n:pi);
phi   = (-pi/2 : 2*pi/2/n : cut)';

x = cos(phi)*cos(theta);
y = cos(phi)*sin(theta);
z = sin(phi)*ones(1,n+1);
xx = xr*x;
yy = yr*y;
zz = zr*z;

end 
                                                            
function [sigma_theta_p, sigma_phi_p, tau_theta_phi_p] = fStress(pp, hGB, D3, k1, k2, theta, phi)    
% function [sigma_theta_p, sigma_phi_p, tau_theta_phi_p] = fStress(pp, hGB, D3, k1, k2, theta, phi)    
%  Three components of the passive stress 
sigma_theta_p   = pp * f_theta(hGB, D3, k1, k2, phi) * f_n(k1, k2, theta, phi);
sigma_phi_p     = pp * f_phi(hGB, D3, k1, k2, theta, phi) / f_n(k1, k2, theta, phi);
tau_theta_phi_p = pp * f_tau(hGB, D3, k1, k2, theta, phi); 
end

function Fn = f_n(k1, k2, theta, phi)
% function Fn = F_n(k1, k2, theta, phi)
%  One of the four functions that describe the instantaneous shape change of a gallbladder   
%  Originally from Novozhilov, V.V., 1964, Thin Shell Theory (Groningen: P.
%  Noordhoff Ltd), pp. 125–130. 
Fn = sqrt(k1^2 * cos(theta)^2 * cos(phi).^2 + k2^2 * cos(theta)^2 * sin(phi).^2 + sin(theta)^2) / ...
     sqrt(k1^2 * sin(phi).^2 + k2^2 * cos(phi).^2);
end % manually checked and corrected

function Fphi = f_phi(hGB, D3, k1, k2, theta, phi)
% function Fphi = F_phi(hGB, D3, k1, k2, theta, phi)
%  One of the four functions that describe the instantaneous shape change of a gallbladder   
%  Originally from Novozhilov, V.V., 1964, Thin Shell Theory (Groningen: P.
%  Noordhoff Ltd), pp. 125–130.p       
Fphi = D3 / (4*k1*k2*hGB)  *( k1^2 * k2^2 + (k1^2 + k2^2 - 2 * k1^2 * k2^2 ) * sin(theta)^2 ...
                                + (k1^2 - k2^2) * cos(theta)^2 * cos(2*phi) );
end % manually checked and corrected

function Ftau = f_tau(hGB, D3, k1, k2, theta, phi)
% function Ftau = F_tau(hGB, D3, k1, k2, theta, phi)
%  One of the four functions that describe the instantaneous shape change of a gallbladder   
%  Originally from Novozhilov, V.V., 1964, Thin Shell Theory (Groningen: P.
%  Noordhoff Ltd), pp. 125–130.p
Ftau = D3 / (4*k1*k2*hGB)  * (k1^2 - k2^2) * cos(theta) * sin(2*phi);
end  % manually checked and corrected

function Ftheta = f_theta(hGB, D3, k1, k2, phi)   % manually checked and corrected
% function Ftheta = F_theta(hGB, D3, k1, k2, phi)
%  One of the four functions that describe the instantaneous shape change of a gallbladder   
%  Originally from Novozhilov, V.V., 1964, Thin Shell Theory (Groningen: P.
%  Noordhoff Ltd), pp. 125–130.
Ftheta = (D3 * k1 * k2 / (4*hGB) ) * (1 - (k1^2 - k2^2) * cos(2*phi) / (k1^2 * k2^2)  );
end

function plotCombined(xx, yy, zz, D1, D2, D3, sphi_p, stheta_p, tauthetaphi_p, sphi_a, stheta_a, tauthetaphi_a)
        scrsz = get(0,'ScreenSize');
        fig = figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);
        
        subplot(3,2,1)
            surf(xx, yy, zz, sphi_p,'EdgeColor','none')
            title('\sigma_\phi^p (latitude)','FontSize',16)
            addPlotDetails(D1, D2, D3)
        subplot(3,2,3)
            surf(xx, yy, zz, stheta_p,'EdgeColor','none')
            title('\sigma_\theta^p (meridian)','FontSize',16)
            addPlotDetails(D1, D2, D3)
        subplot(3,2,5)
            surf(xx, yy, zz, tauthetaphi_p,'EdgeColor','none')
            title('\tau_\theta_\phi^p (in surface)','FontSize',16)
            addPlotDetails(D1, D2, D3)
            
        subplot(3,2,2)
            surf(xx, yy, zz, sphi_a,'EdgeColor','none')
            title('\sigma_\phi^a (latitude)','FontSize',16)
            addPlotDetails(D1, D2, D3)
        subplot(3,2,4)
            surf(xx, yy, zz, stheta_a,'EdgeColor','none')
            title('\sigma_\theta^a (meridian)','FontSize',16)
            addPlotDetails(D1, D2, D3)
        subplot(3,2,6)
            surf(xx, yy, zz, tauthetaphi_a,'EdgeColor','none')
            title('\tau_\theta_\phi^a (in surface)','FontSize',16)
            addPlotDetails(D1, D2, D3)
        pause;
        set(0,'ShowHiddenHandles','on')
        delete(get(0,'Children'))
    end

function plotSeparate(xx, yy, zz, D1, D2, D3, sphi_p, stheta_p, tauthetaphi_p, sphi_a, stheta_a, tauthetaphi_a)

            scrsz = get(0,'ScreenSize');
            fig = figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);  
            surf(xx, yy, zz, sphi_p,'EdgeColor','none')
            title('\sigma_\phi^p (latitude)','FontSize',16)
            addPlotDetails(D1, D2, D3)
            
            scrsz = get(0,'ScreenSize');
            fig = figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);           
        %    [size(xx) size(yy) size(zz) size(stheta_p)]
            surf(xx, yy, zz, stheta_p,'EdgeColor','none')
            title('\sigma_\theta^p (meridian)','FontSize',16)
            addPlotDetails(D1, D2, D3)

            scrsz = get(0,'ScreenSize');
            fig = figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);                       
            surf(xx, yy, zz, tauthetaphi_p,'EdgeColor','none')
            title('\tau_\theta_\phi^p (in surface)','FontSize',16)
            addPlotDetails(D1, D2, D3)

            scrsz = get(0,'ScreenSize');
            fig = figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);                       
            surf(xx, yy, zz, sphi_a,'EdgeColor','none')
            title('\sigma_\phi^a (latitude)','FontSize',16)
            addPlotDetails(D1, D2, D3)
   
            scrsz = get(0,'ScreenSize');
            fig = figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);                         
            surf(xx, yy, zz, stheta_a,'EdgeColor','none')
            title('\sigma_\theta^a (meridian)','FontSize',16)
            addPlotDetails(D1, D2, D3)
        
            scrsz = get(0,'ScreenSize');
            fig = figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);           
            surf(xx, yy, zz, tauthetaphi_a,'EdgeColor','none')
            title('\tau_\theta_\phi^a (in surface)','FontSize',16)
            addPlotDetails(D1, D2, D3)
   
        pause;
        set(0,'ShowHiddenHandles','on')
        delete(get(0,'Children'))
end

function addPlotDetails(D1, D2, D3)
h1 = xlabel('x');
h2 = ylabel('y');
h3 = zlabel('z');
% xlim([-(D1/2)*1.1 (D1/2)*1.1])
% ylim([-(D2/2)*1.1 (D2/2)*1.1])
% zlim([-(D3/2)*1.1 (D3/2)*1.1])
fsize = 10;
colorbar('location','EastOutside')
% set(gca,'fontsize',fsize) % increase the size
%set(gca,'Xtick',round([-D1/2 D1/2]),'Ytick',round([-D2/2 D2/2]),'Ztick',round([-D3/2 D3/2]) );
set(h1,'fontsize',fsize) % increase the size
set(h2,'fontsize',fsize) % increase the size	
set(h3,'fontsize',fsize) % increase the size	
% set(gca,'DataAspectRatio',[1 1 1])
% set(gca,'LooseInset',get(gca,'TightInset'));  % Remove excessive figure margins
%grid off
%set(gca,'color','none') 
end

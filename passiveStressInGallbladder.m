function Dat = passiveStressInGallbladder(PatientNumber, reduced_num, PercentChangeInhGB, PercentageOfLesion, PercentOfGallstones, rho, s1, s2, s3, plt)
% function Dat = passiveStressInGallbladder(PatientNumber, reduced_num, PercentChangeInhGB, PercentageOfLesion, PercentOfGallstones, rho, s1, s2, s3, plt)
%  plt = 0 for no plot, 1 for 3D plot
%  PercentChangeInhGB is the change in hGB
%  The complete working function to plot graphs for our paper
%  Dimensions scaling factors
%  s = [1/s1 1/s2 1/s3]; 
%  Written by M. Ali Yousuf
%  Johns Hopkins Medicine
%  Ver 1.0, Feb 29, 2012

cd('E:\Dropbox\Research - Private Folders\utilities\Gallbladder');

%num = xlsread('CombinedData.xlsx');  % Read Table 1 of paper
% txt = PatientNo	D1	D2	D3 (mm)	

I = PatientNumber;

hGB = 2.5;  % Gallbladder wall thickness in mm 
            % Reference: Engel - Gallbladder Wall Thickness - Sonographic
            % Accuracy and Relation to Disease - 1980 

hGB = (1+PercentChangeInhGB/100)*hGB;

pe = 11;    % gallbladder internal pressure while emptying, 11 mmHg in the paper
            % Ref: Majeed - Continuous in vivo manometry of the human gallbladder - 1994

D1 = reduced_num(I,2);  % Diameter along the major axis at time t, in mm
D2 = reduced_num(I,3);  % Diameter along the major axis at time t, in mm
D3 = reduced_num(I,4);  % Diameter along the major axis at time t, in mm

% D1 = reduced_num(I,2);  % Diameter along the major axis at time t, in mm
% D2 = reduced_num(I,4);  % Diameter along the major axis at time t, in mm
% D3 = reduced_num(I,3);  % Diameter along the major axis at time t, in mm

% D1 = reduced_num(I,3);  % Diameter along the major axis at time t, in mm
% D2 = reduced_num(I,2);  % Diameter along the major axis at time t, in mm
% D3 = reduced_num(I,4);  % Diameter along the major axis at time t, in mm

% D1 = reduced_num(I,3);  % Diameter along the major axis at time t, in mm
% D2 = reduced_num(I,4);  % Diameter along the major axis at time t, in mm
% D3 = reduced_num(I,2);  % Diameter along the major axis at time t, in mm

% D1 = reduced_num(I,4);  % Diameter along the major axis at time t, in mm
% D2 = reduced_num(I,2);  % Diameter along the major axis at time t, in mm
% D3 = reduced_num(I,3);  % Diameter along the major axis at time t, in mm
% 
% D1 = reduced_num(I,4);  % Diameter along the major axis at time t, in mm
% D2 = reduced_num(I,3);  % Diameter along the major axis at time t, in mm
% D3 = reduced_num(I,2);  % Diameter along the major axis at time t, in mm

D1 = D1 / s1;
D2 = D2 / s2;
D3 = D3 / s3;

k1 = D1/D3;
k2 = D2/D3;

%p = pe;

p = cal_p(D1, D2, D3, pe, PercentOfGallstones, rho);
%p = pe; 

%******************Start Calculating Stress***********************
n = 50;     % The mesh size - higher values generate finer mesh
cut = PercentageOfLesion*2;   % Percentage of the GB converted to lesion 
cut = (100-cut) / 100;   
cut = asin(cut);

plt2 = 0;
[xx,yy,zz]= cutellipsoid(D1,D2,D3,cut,n,plt2); 

sphi_p = []; stheta_p = [];

for theta = -pi : 2*pi/n : pi;
    s1_p = []; s2_p = []; 
        for phi   = -pi/2 : 2*pi/2/n : cut;
            [sigma_theta_p, sigma_phi_p, tau_theta_phi_p] = fStress(p, hGB, D3,  k1,  k2, theta, phi);         
            s1_p = [s1_p; sigma_theta_p];
            s2_p = [s2_p; sigma_phi_p];
        end

       stheta_p = [stheta_p, s1_p];
       sphi_p = [sphi_p, s2_p];
   
end

%*****************Calculate Peak Stress ***********************

sigma_max_p = max( max(max(stheta_p   ))   , max(max(sphi_p)  ) );

Dat = [PatientNumber, PercentageOfLesion, D1, D2, D3, hGB, sigma_max_p];

%******************Start Plotting******************************
%**************************************************************

    if plt ==1
    plotCombined(xx, yy, zz, D1, D2, D3, sphi_p, stheta_p,hGB);
%    save(['Plot_',num2str(s1),'_',num2str(s2),'_',num2str(s3),num2str(PercentChangeInhGB),'.x'],'xx','-ASCII')
    end

end 

function [xx,yy,zz]= cutellipsoid(xr,yr,zr,cut,n,plt2)
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

    if plt2 == 1
        surf(xx, yy, zz)
        pause
    end 

end 
                                                            
function [sigma_theta_p, sigma_phi_p, tau_theta_phi_p] = fStress(pp, hGB, D3, k1, k2, theta, phi)    
% function [sigma_theta_p, sigma_phi_p, tau_theta_phi_p] = fStress(pp, hGB, D3, k1, k2, theta, phi)    
%  Three components of the passive stress 
sigma_theta_p   = pp * f_theta(hGB, D3, k1, k2, phi) * f_n(k1, k2, theta, phi);
sigma_phi_p     = pp * f_phi(hGB, D3, k1, k2, theta, phi) / f_n(k1, k2, theta, phi);
tau_theta_phi_p = 0;  % this reduces the amount of work, as it is not needed for the paper
%tau_theta_phi_p = pp * f_tau(hGB, D3, k1, k2, theta, phi); 
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

function plotCombined(xx, yy, zz, D1, D2, D3, sphi_p, stheta_p, hGB)
        scrsz = get(0,'ScreenSize');
        fig = figure('OuterPosition',[1 5 scrsz(3)/3 scrsz(4)]);
        
        subplot(2,1,1)
            surf(xx, yy, zz, sphi_p,'EdgeColor','none')
            title(['\sigma_\phi^p, ',' h_G_B = ',num2str(hGB),' mm'], 'FontSize',16)
            addPlotDetails(D1, D2, D3, hGB)
        subplot(2,1,2)
            surf(xx, yy, zz, stheta_p,'EdgeColor','none')
            title(['\sigma_\theta^p, ',' h_G_B = ',num2str(hGB),' mm'],'FontSize',16)
            addPlotDetails(D1, D2, D3, hGB)
        
        pause;
        set(0,'ShowHiddenHandles','on')
        delete(get(0,'Children'))
    end

function addPlotDetails(D1, D2, D3, hGB)
h1 = xlabel('x');
h2 = ylabel('y');
h3 = zlabel('z');
fsize = 16;
%title(['Wall thickness, h_G_B = ',num2str(hGB),' mm'],'FontSize',fsize)
colorbar('location','EastOutside')
set(gca,'fontsize',fsize) % increase the size
set(h1,'fontsize',fsize) % increase the size
set(h2,'fontsize',fsize) % increase the size	
set(h3,'fontsize',fsize) % increase the size	
grid off
end

function simpleStress = CylinderSphereModel(D1, D2, D3, p, hGB)
% This model assumes a cylindrical body and ellipsoidal caps
%    
%              --------------------------
%           (                               )
%           (                               )
%              -------------------------- 
%
%       80% of D3 = central cylindrical part 
%       Volume of ellipsoidal caps = (4/3)*a*b * 0.1*c 
%
a = D1/2;
b = D2/2;
c = D3/2;
g = 9.8 * 10^3;       % mm/sec^2
rho = 2000 *(1/10^9); % density of gallstones in Kg/mm^3
                      % For comparison, rho = 1000 for water 
Volume_of_ellipsoidal_caps = (4/3)*a*b * 0.1*c; 
Volume_of_ellipsoide = (4/3)*a*b*c;
perct = 0.66;         % Percentage of gallstones in the gallbladder
VolumeOfStones = perct*Volume_of_ellipsoide;        % Volume of gallstones in mm^3
F = rho*g*VolumeOfStones;    % Lengths measured in mm

% ellipsoidal part
gi = (a-b)^2 / (a+b)^2;
Circum = pi*(a+b)*(1+  3*gi / (10+sqrt(4-3*gi))  );
Area = pi*a*b;   % in mm^2

Tot_p_ellipse = p + F/Area;    %Fluid pressure + pressure due to gallstones
ellipsicalStress = Tot_p_ellipse * Area / (hGB*Circum);

%cylindrical part
r = min(a,b);      % the side with smallest area will have largest pressure
Area = 2*r*D3;
Tot_p_cylinder = p + F/Area;
Sigma_h = Tot_p_cylinder*Area / (2*D3*hGB);                    % Hoop stress
%Sigma_l = Tot_p*Area / (hGB*Circum);		% Longitudinal stress, same as for
                                            %ellipse

simpleStress = max(ellipsicalStress, Sigma_h);

end 

function p = cal_p(D1, D2, D3, pe, PercentOfGallstones,rho)
% u = 1.6075;
% a = D1/2;
% b = D2/2;
% c = D3/2;

g = 9.8;              % m/sec^2
u = 1.6075;
a = D1/2/1000;        % Convert to m by dividing with 1000
b = D2/2/1000;
c = D3/2/1000;

SurfaceAreaOfEllipsoid = 4*pi*( (a^u*b^u + a^u*c^u + b^u*c^u ) /3     )^(1/u);
MinimumSurfaceArea = pi*a*b;           % This is what you see from the top, even if there is a contraction 
Volume_of_ellipsoid = (4/3)*a*b*c; 
frac = PercentOfGallstones/100;                % Fraction of gallstone volume in the gallbladder
VolumeOfStones = frac*Volume_of_ellipsoid;     % Volume of gallstones in mm^3
F = rho*g*VolumeOfStones;                      % Lengths measured in mm, Force in Newtons 
p_correction = F / MinimumSurfaceArea; % Divide by two as stones area on one side only
p_correction = p_correction / 133;             % Convert pascals into mmHg
p = pe + p_correction;                         % Fluid pressure + pressure due to gallstones                         
                                            
end
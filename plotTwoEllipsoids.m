function plotTwoEllipsoids

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);

n = 100; 
f1 = 2; f2 = 20.5;
D4 = 2; D5 = 2; D6 = 5;
D7 = 1; D8 = 1; D9 = 8;
% solving the two ellipsoids of revolution we get a quadratic with 
A = D4^2/D6^2 - D7^2/D9^2;
B = 2*f2*D4^2/D6^2 + 2*f1*D7^2/D9^2;
C = f2^2 * D4^2/D6^2 - f1^2 * D7^2/D9^2 - (D4/2)^2 + (D7/2)^2;
radi = B^2 - 4 * A * C
if radi < 0 
    disp('Ellipsoids do not intersect')
else
    z1 = (-B + sqrt(radi) ) / (2*A);
    z2 = (-B - sqrt(radi) ) / (2*A);
    cut = 0.1; %asin(z1/D6);
    [x, y, z] = cutTheUpperEllipsoid(0,0,f1,D4/2,D5/2,D6/2,cut,n);
    surfl(x, y, z)
    colormap(spring)
    hold on
    hold all
    [x, y, z] = cutTheLowerEllipsoid(0,0,f2,D7/2,D8/2,D9/2,cut,n);
    surfl(x, y, z)
    colormap(autumn)

    grid off
%     xlim([-(D4/2)*1.1 (D4/2)*1.1])
%     ylim([-(D5/2)*1.1 (D5/2)*1.1])
%     zlim([-(D6+D9)*1.1/2 (D6+D9)*1.1/2])
    alpha(0.5)   % Makes the ellipsoid transparent
    h1 = xlabel('x');
    h2 = ylabel('y');
    h3 = zlabel('z');

end

pause;
set(0,'ShowHiddenHandles','on')
delete(get(0,'Children'))
end 

function [xx,yy,zz]= cutTheUpperEllipsoid(xc, yc, zc, xr, yr, zr, cut, n)
%Based on ELLIPSOID which Generates ellipsoid.
%   [X,Y,Z]=ELLIPSOID(XC,YC,ZC,XR,YR,ZR,N) generates three
%   (N+1)-by-(N+1) matrices so that SURF(X,Y,Z) produces an
%   ellipsoid with center (XC,YC,ZC) and radii XR, YR, ZR.

theta = (-pi:2*pi/n:pi);
phi   = (-pi/2 : 2*pi/2/n : cut)';

x = cos(phi)*cos(theta);
y = cos(phi)*sin(theta);
z = sin(phi)*ones(1,n+1);

xx = xr*x + xc;
yy = yr*y + yc;
zz = zr*z + zc;

end 

function [xx,yy,zz]= cutTheLowerEllipsoid(xc, yc, zc, xr, yr, zr, cut, n)
%Based on ELLIPSOID which Generates ellipsoid.
%   [X,Y,Z]=ELLIPSOID(XC,YC,ZC,XR,YR,ZR,N) generates three
%   (N+1)-by-(N+1) matrices so that SURF(X,Y,Z) produces an
%   ellipsoid with center (XC,YC,ZC) and radii XR, YR, ZR.

theta = (-pi:2*pi/n:pi);
phi   = (cut : 2*pi/2/n : pi / 2)';

x = cos(phi)*cos(theta);
y = cos(phi)*sin(theta);
z = sin(phi)*ones(1,n+1);

xx = xr*x + xc;
yy = yr*y + yc;
zz = zr*z + zc;

end 
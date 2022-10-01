function plotCollapsedGallbladder(PercentageOfLesion)

cut = PercentageOfLesion*2;   % Percentage of the GB conveted to lesion 
cut = (100-cut) / 100;   
cut = asin(cut);

n = 50;
xr = 1;
yr = 1;
zr = 2;
[xx,yy,zz] = galloriginal(0,0,0,xr,yr,zr,n,cut);  % (xc,yc,zc,xr,yr,zr,n, alpha)
scrsz = get(0,'ScreenSize');
figure('Name','figure name','NumberTitle','off','OuterPosition',[1 scrsz(4)/7.5 scrsz(3) 0.9*scrsz(4)]);
surf(xx, yy, zz)
grid off
%addPlotDetails
end 

function [xx,yy,zz] = galloriginal(xc,yc,zc,xr,yr,zr,n, cut)

theta = (-pi:2*pi/n:pi);
phi   = (-pi/2 : 2*pi/2/n : cut)';

xx = xr*cos(phi)*cos(theta);
yy = yr*cos(phi)*sin(theta);
zz = zr*sin(phi)*ones(1,size(yy,2));

[size(theta) size(phi) size(xx) size(yy) size(zz)]

end
function addPlotDetails
h1 = xlabel(['xr']);
h2 = ylabel(['yr']);
h3 = zlabel(['zr']);
% xlim([-1 1]);
% ylim([-1 1]);
% zlim([-1 1]);
fsize = 16;
set(gca,'fontsize',fsize) % increase the size
set(h1,'fontsize',fsize) % increase the size
set(h2,'fontsize',fsize) % increase the size	
set(h3,'fontsize',fsize) % increase the size	
colorbar('location','EastOutside')
grid on
hold on
hold all
set(gca,'DataAspectRatio',[1 1 1])

pause;
set(0,'ShowHiddenHandles','on')
delete(get(0,'Children'))
end 

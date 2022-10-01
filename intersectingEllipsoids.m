function intersectingEllipsoids(f1)

D4 = 4; 
D5 = 4; 
D6 = 10;
f2 = 0;

D7 = 2; 
D8 = 2; 
D9 = 4;
%f1 = 3; %6; 

A = (D4/D6)^2 - (D7/D9)^2;
B = 2*f2*(D4/D6)^2 + 2*f1*(D7/D9)^2;
C = f2^2 * (D4/D6)^2 - f1^2 * (D7/D9)^2 - (D4/2)^2 + (D7/2)^2;

if A ~= 0
    radi = B^2 - 4 * A * C;
    if radi < 0; disp('No real roots'); end; 
    z1 = (-B + sqrt(radi) ) / (2*A);
    z2 = (-B - sqrt(radi) ) / (2*A);
    [subs(z1) subs(z2)]
else
    z1 = subs( (D9/2)^2 - (D6/2)^2 + f1^2 ) / (2 * f1);
%    z1 = C / B
end

[xx,yy,zz] = ellipsoid(0,0,f1,D4/2, D5/2, D6/2,100);  % (xc,yc,zc,xr,yr,zr,n)
surf(xx, yy, zz)
hold on

[xx,yy,zz] = ellipsoid(0,0,f2,D7/2, D8/2, D9/2,100);  % (xc,yc,zc,xr,yr,zr,n)
surf(xx, yy, zz)
alpha(0.2)   % Makes the ellipsoid transparent

%grid off
h1 = xlabel('x');
h2 = ylabel('y');
h3 = zlabel('z');
%set(gca,'DataAspectRatio',[1 1 1])
view([0,20,0])

pause;
set(0,'ShowHiddenHandles','on')
delete(get(0,'Children'))

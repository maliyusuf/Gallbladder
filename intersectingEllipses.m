function intersectingEllipses(f)
% function intersectingEllipses
%  Works perfectly fine
scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);

% x^2 / a1^2 +  y^2  /   b1^2 = 1
% x^2 / a2^2 + (y-f)^2 / b2^2 = 1
a1 = 8;
b1 = 6;
a2 = 3;
b2 = 3;
%f = 4.0;
syms x y
[x y] = solve(x^2 / a1^2 +  y^2 / b1^2 - 1, x^2 / a2^2 + (y-f)^2 / b2^2 - 1 );
xlimit = subs(x(1,1));

% xlimit = sqrt(r2^2 - (r1^2 - r2^2 - f^2)^2 / (4*f^2)  );

for x = -a1:0.1:-xlimit
    radi = 1 - x^2/a1^2;
    if radi < 0; disp('Ellipses do not intersect'); end 
    y = b1* sqrt(radi);
    plot(x,y,'b*')
    hold on
    plot(x,-y,'y*')
end
for x = -xlimit:0.1:+xlimit
    radi2 = 1 - x^2/a2^2;
    if radi < 0; disp('Ellipses do not intersect'); end 
    y = b2* sqrt(radi2) + f;
    plot(x,y,'r*')
    y2 = b1* sqrt(1 - x^2/a1^2);
    plot(x,-y2,'y*')
    hold on
end
for x = xlimit:0.1:a1
    radi = 1 - x^2/a1^2;
    if radi < 0; disp('Ellipses do not intersect'); end 
    y = b1* sqrt(radi);
    plot(x,y,'b*')
    plot(x,-y,'y*')
    hold on
end
grid on
%xlim([-a1*1.1, a1*1.1])
%ylim([-b1*1.1, b1+b2*1.1])
set(gca,'DataAspectRatio',[1 1 1])
pause;
set(0,'ShowHiddenHandles','on')
delete(get(0,'Children'))

end 
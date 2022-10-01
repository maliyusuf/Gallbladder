function intersectingCircles
% function intersectingCircles
% Works perfectly fine
scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);
r1 = 6;
r2 = 2;
f = 5;
% x^2 + y^2     = r1^2
% x^2 + (y-f)^2 = r2^2
syms x y
[x y] = solve(x^2 + y^2 - r1^2,x^2 + (y-f)^2 - r2^2);
xlimit = subs(x);
%xlimit = sqrt(r2^2 - (r1^2 - r2^2 - f^2)^2 / (4*f^2)  );

for x = -r1:0.1:-xlimit
    y = sqrt(r1^2 - x^2);
    plot(x,y,'b*')
    hold on
end
for x = -xlimit:0.1:+xlimit
    y = sqrt(r2^2 - x^2) + f;
    plot(x,y,'b*')
    hold on
end
for x = xlimit:0.1:r1
    y = sqrt(r1^2 - x^2);
    plot(x,y,'b*')
    hold on
end

pause;
set(0,'ShowHiddenHandles','on')
delete(get(0,'Children'))

end 
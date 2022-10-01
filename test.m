data =[   111    73    57    64   154   137   123   115   135   196   200   147   134   191    91    76    64    64    73    57 ...
    63   235   160    92   121   120    67    57    55   110   194   110    78   187    86    91    62   106    60   101 ...
   201    75   272   188    66    67   107   143   141    62   128    93   143   131   142   368   290   150   210   301 ...
   134   514   209   271   770   217    86    99    95   247   111    83   157   404   138   249   247   247   133   106 ...
   211   590    99   132   596   119   173   104];


%x = min(data)-25:50:max(data)+25;
%x = [75 125 175 225];
x = 25:25:800;
hist(data,x)
pause
a = histc(data,x);
bar(x,a)
xlim([0 max(data)])
h1 = xlabel('Stress (mmHg)');
h2 = ylabel('Number of Subjects');
fsize = 20;
set(gca,'fontsize',fsize) % increase the size
set(h1,'fontsize',fsize) % increase the size
set(h2,'fontsize',fsize) % increase the size
%plot(x,a)


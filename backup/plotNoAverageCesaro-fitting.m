clear all;
set(0, "defaultlinelinewidth", 4);

%md.dc = load("md-lag20000.dcNoAverageCesaro" );
%md3.dc = load("md3-lag20000.dcNoAverageCesaro" );

md.ec = load("lag1000000.ecNoAverageCesaro" );

%f1 = figure(1);
%clf;
%hold on;
%p1.na = plot(md.dc.timeLags, md.dc.dcNoAverageCesaro(:,1), "-1");
%p1.cl = plot(md.dc.timeLags, md.dc.dcNoAverageCesaro(:,2), "-2");
%p2.na = plot(md3.dc.timeLags, md3.dc.dcNoAverageCesaro(:,1), "-3");
%p2.cl = plot(md3.dc.timeLags, md3.dc.dcNoAverageCesaro(:,2), "-4");
%
%set([p1.na, p1.cl, p2.na, p2.cl], "linewidth", 4);
%
%title("Non-averaged Cesaro sum for diffusion coefficients of 1m NaCl solution");
%legend("{/Symbol D}t = 5 fs, Na+", "{/Symbol D}t = 5 fs, Cl-", "{/Symbol D}t = 1 fs, Na+", "{/Symbol D}t = 1 fs, Cl-");
%xlabel("partial sum upper limit (ps)");
%ylabel("Non-averaged partial sum (m^2/s)");
%axis([0,100, 0, 3e-7]);
%
%print('dcNoAverageCesaro.eps', '-deps', '-color');
%
%axis([0,0.1,0,2.5e-10]);
%print('dcNoAverageCesaro-zoom.eps', '-deps', '-color');

f1 = figure(1);
clf;
hold on;
p1 = plot(md.ec.timeLags, md.ec.ecTotalNoAverageCesaro, "-1");
p2 = plot(md.ec.timeLags, md.ec.ecAutocorrNoAverageCesaro(:,1), "-3");
p3 = plot(md.ec.timeLags, md.ec.ecAutocorrNoAverageCesaro(:,2), "-", "color", [0, 0.6, 0]);
p4 = plot(md.ec.timeLags, md.ec.ecCorrNoAverageCesaro(:,2)+md.ec.ecCorrNoAverageCesaro(:,3), "-4");

fitRange=[40001:80000];
slope1 = polyfit(md.ec.timeLags(fitRange), md.ec.ecTotalNoAverageCesaro(fitRange), 1)(1);
slope2 = polyfit(md.ec.timeLags(fitRange), md.ec.ecAutocorrNoAverageCesaro(:,1)(fitRange), 1)(1);
slope3 = polyfit(md.ec.timeLags(fitRange), md.ec.ecAutocorrNoAverageCesaro(:,2)(fitRange), 1)(1);
slope4 = polyfit(md.ec.timeLags(fitRange), (md.ec.ecCorrNoAverageCesaro(:,2)+md.ec.ecCorrNoAverageCesaro(:,3))(fitRange), 1)(1);

text(60,1200, num2str(slope1));
text(60,300, num2str(slope2));
text(60,700, num2str(slope3));
text(60,-50, num2str(slope4));

set(p1, "linewidth", 6);
set(p2, "linewidth", 6);
set(p3, "linewidth", 6);
set(p4, "linewidth", 6);

title("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution");
%legend("{/Symbol D}t = 1 fs");
legend("Total", "Auto Na+", "Auto Cl-", "Cross Na-Cl", "location", "north");
xlabel("partial sum upper limit (ps)");
ylabel("Non-averaged partial sum (ps*S/m)");
axis([0,80,-300, 2000]);

print('ecNoAverageCesaro.eps', '-deps', '-color');

%axis([0,0.1,0,1.5]);
%print('ecNoAverageCesaro-zoom.eps', '-deps', '-color');


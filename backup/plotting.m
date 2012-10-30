clear all;
set(0, "defaultlinelinewidth", 4);

md.dc = load("md-lag20000.dcCesaro" );
md3.dc = load("md3-lag20000.dcCesaro" );

md.ec = load("md-lag20000.ecCesaro" );
md3.ec = load("md3-lag20000.ecCesaro" );

f1 = figure(1);
clf;
hold on;
p1.na = plot(md.dc.timeLags, md.dc.dcCesaro(:,1), "-1");
p1.cl = plot(md.dc.timeLags, md.dc.dcCesaro(:,2), "-2");
p2.na = plot(md3.dc.timeLags, md3.dc.dcCesaro(:,1), "-3");
p2.cl = plot(md3.dc.timeLags, md3.dc.dcCesaro(:,2), "-4");

set([p1.na, p1.cl, p2.na, p2.cl], "linewidth", 4);

title("Cesaro sum for diffusion coefficients of 1m NaCl solution");
legend("{/Symbol D}t = 5 fs, Na+", "{/Symbol D}t = 5 fs, Cl-", "{/Symbol D}t = 1 fs, Na+", "{/Symbol D}t = 1 fs, Cl-");
xlabel("partial sum upper limit (ps)");
ylabel("average of partial sum (m^2/s)");
axis([0,100,1.5e-9, 3e-9]);

print('dcCesaro.eps', '-deps', '-color');


f2 = figure(2);
clf;
hold on;
p3 = plot(md.ec.timeLags, md.ec.ecTotalCesaro, "-1");
p4 = plot(md3.ec.timeLags, md3.ec.ecTotalCesaro, "-3");

set([p3, p4], "linewidth", 4);

title("Cesaro sum for electrical conductivity of 1m NaCl solution");
legend("{/Symbol D}t = 5 fs", "{/Symbol D}t = 1 fs");
xlabel("partial sum upper limit (ps)");
ylabel("average of partial sum (S/m)");
axis([0,100,8, 15]);

print('ecCesaro.eps', '-deps', '-color');


clear all
hold off;
md = load("md-lag20000.vCorr");
md3 = load("md3-lag20000.vCorr");

f3 = figure(3);
p1.na = plot(md.timeLags, md.vAutocorr{1}, "color", "red");
axis([0,10, -0.4,0.5]);
title("<v(0)v(t)> of Na+ in 1m NaCl solution");
legend("{/Symbol D}t = 5 fs");
xlabel("time lag (ps)");
ylabel("<v(0)v(t)> (nm^2/ps^2)");

print('corr.eps', '-deps', '-color');

f4 = figure(4);
p2.na = plot(md.timeLags, abs(md.vAutocorr{1}), "color", "red");
title("|<v(0)v(t)>| of Na+ in 1m NaCl solution");
legend("{/Symbol D}t = 5 fs");
xlabel("time lag (ps)");
ylabel("|<v(0)v(t)>| (nm^2/ps^2)");
axis([0,10, -0.4,0.5]);

print('abs-corr.eps', '-deps', '-color');

f5 = figure(5);
p3.na = loglog(md.timeLags, abs(md.vAutocorr{1}), "color", "blue");
title("loglog plot for |<v(0)v(t)>| of Na+ in 1m NaCl solution");
legend("{/Symbol D}t = 5 fs");
xlabel("time lag (ps)");
ylabel("|<v(0)v(t)>| (nm^2/ps^2)");
axis([1e-3, 10]);

print('log-abs-corr.eps', '-deps', '-color');


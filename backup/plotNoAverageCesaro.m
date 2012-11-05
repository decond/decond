clear all;
set(0, "defaultlinelinewidth", 4);

md.dc = load("md-lag20000.dcNoAverageCesaro" );
md3.dc = load("md3-lag20000.dcNoAverageCesaro" );

md.ec = load("md-lag20000.ecNoAverageCesaro" );
md3.ec = load("md3-lag20000.ecNoAverageCesaro" );

f1 = figure(1);
clf;
hold on;
p1.na = plot(md.dc.timeLags, md.dc.dcNoAverageCesaro(:,1), "-1");
p1.cl = plot(md.dc.timeLags, md.dc.dcNoAverageCesaro(:,2), "-2");
p2.na = plot(md3.dc.timeLags, md3.dc.dcNoAverageCesaro(:,1), "-3");
p2.cl = plot(md3.dc.timeLags, md3.dc.dcNoAverageCesaro(:,2), "-4");

set([p1.na, p1.cl, p2.na, p2.cl], "linewidth", 4);

title("Non-averaged Cesaro sum for diffusion coefficients of 1m NaCl solution");
legend("{/Symbol D}t = 5 fs, Na+", "{/Symbol D}t = 5 fs, Cl-", "{/Symbol D}t = 1 fs, Na+", "{/Symbol D}t = 1 fs, Cl-");
xlabel("partial sum upper limit (ps)");
ylabel("Non-averaged partial sum (m^2/s)");
axis([0,100, 0, 3e-7]);

print('dcNoAverageCesaro.eps', '-deps', '-color');

axis([0,0.1,0,2.5e-10]);
print('dcNoAverageCesaro-zoom.eps', '-deps', '-color');

f2 = figure(2);
clf;
hold on;
p3 = plot(md.ec.timeLags, md.ec.ecTotalNoAverageCesaro, "-1");
p4 = plot(md3.ec.timeLags, md3.ec.ecTotalNoAverageCesaro, "-3");

set([p3, p4], "linewidth", 4);

title("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution");
legend("{/Symbol D}t = 5 fs", "{/Symbol D}t = 1 fs");
xlabel("partial sum upper limit (ps)");
ylabel("Non-averaged partial sum (S/m)");
axis([0,100,0, 1200]);

print('ecNoAverageCesaro.eps', '-deps', '-color');

axis([0,0.1,0,1.5]);
print('ecNoAverageCesaro-zoom.eps', '-deps', '-color');


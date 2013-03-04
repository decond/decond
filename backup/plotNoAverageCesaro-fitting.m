#!/home/kmtu/bin/octave -qf

clear all;

if (nargin() < 1)
    error("Usage: $plotNoAverageCesaro.m <md#>")
endif

md_str = argv(){1}

set(0, "defaultlinelinewidth", 4);

%md.dc = load("md-lag20000.dcNoAverageCesaro" );
%md3.dc = load("md3-lag20000.dcNoAverageCesaro" );

md.ec = load(strcat("./md", md_str, "/lag1000000.ecNoAverageCesaro"));

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
%p3 = plot(md.ec.timeLags, md.ec.ecAutocorrNoAverageCesaro(:,2), "-", "color", [0, 0.6, 0]);
%p4 = plot(md.ec.timeLags, md.ec.ecCorrNoAverageCesaro(:,2)+md.ec.ecCorrNoAverageCesaro(:,3), "-4");
p4 = plot(md.ec.timeLags, md.ec.ecCorrNoAverageCesaro(:,1), "-2");
p5 = plot(md.ec.timeLags, md.ec.ecCorrNoAverageCesaro(:,2), "-4");
p6 = plot(md.ec.timeLags, md.ec.ecCorrNoAverageCesaro(:,3), "-5");
p7 = plot(md.ec.timeLags, md.ec.ecCorrNoAverageCesaro(:,4), "-", "color", [0, 0, 0]);

fitRange = [20, 40; 40, 60; 60, 80; 80, 100];
fitRange .*= 1000;
for r = [1:size(fitRange, 1)]
    slope(1,r) = polyfit(md.ec.timeLags(fitRange(r, 1):fitRange(r, 2)), md.ec.ecTotalNoAverageCesaro(fitRange(r, 1):fitRange(r, 2)), 1)(1);
    slope(2,r) = polyfit(md.ec.timeLags(fitRange(r, 1):fitRange(r, 2)), md.ec.ecAutocorrNoAverageCesaro(:,1)(fitRange(r, 1):fitRange(r, 2)), 1)(1);
    slope(3,r) = polyfit(md.ec.timeLags(fitRange(r, 1):fitRange(r, 2)), md.ec.ecAutocorrNoAverageCesaro(:,2)(fitRange(r, 1):fitRange(r, 2)), 1)(1);
    slope(4,r) = polyfit(md.ec.timeLags(fitRange(r, 1):fitRange(r, 2)), (md.ec.ecCorrNoAverageCesaro(:,1))(fitRange(r, 1):fitRange(r, 2)), 1)(1);
    slope(5,r) = polyfit(md.ec.timeLags(fitRange(r, 1):fitRange(r, 2)), (md.ec.ecCorrNoAverageCesaro(:,2))(fitRange(r, 1):fitRange(r, 2)), 1)(1);
    slope(6,r) = polyfit(md.ec.timeLags(fitRange(r, 1):fitRange(r, 2)), (md.ec.ecCorrNoAverageCesaro(:,3))(fitRange(r, 1):fitRange(r, 2)), 1)(1);
    slope(7,r) = polyfit(md.ec.timeLags(fitRange(r, 1):fitRange(r, 2)), (md.ec.ecCorrNoAverageCesaro(:,4))(fitRange(r, 1):fitRange(r, 2)), 1)(1);
endfor

left = 25;
h_space = 20;
top = 972;
v_space = 44;
for r = [1:size(fitRange, 1)]
    for i = [1:7]
        text(left + (r-1)*h_space, top - (i-1)*v_space, num2str(slope(i, r)));
    endfor
endfor

%set(p1, "linewidth", 3);
%set(p2, "linewidth", 3);
%set(p3, "linewidth", 3);
%set(p4, "linewidth", 3);

title(strcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution - md", md_str));
%legend("{/Symbol D}t = 1 fs");
legend("Total", "Auto Na+", "Auto Cl-", "Cross Na-Na", "Cross Na-Cl", "Cross Cl-Na", "Cross Cl-Cl", "location", "northwest");
xlabel("partial sum upper limit (ps)");
ylabel("Non-averaged partial sum (ps*S/m)");
axis([0,100,-600, 1000]);

print(strcat('ecNoAverageCesaro-', md_str, '.eps'), '-deps', '-color');

%axis([0,0.1,0,1.5]);
%print('ecNoAverageCesaro-zoom.eps', '-deps', '-color');


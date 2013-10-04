#!/home/kmtu/bin/octave -qf

clear all;
global numIonTypes;

if (nargin() < 3)
    error("Usage: $plotNoAverageCesaro.m <dataFilename> <skip> <dt>")
endif

dataFilename = argv(){1}
skip = str2num(argv(){2}) #skipped interval in vCorr data
deltaStep = str2num(argv(){3})

set(0, "defaultlinelinewidth", 4);

%md.dc = load("md-lag20000.dcNoAverageCesaro" );
%md3.dc = load("md3-lag20000.dcNoAverageCesaro" );

numMD = 1;

initialize = load(strcat("./md0/", dataFilename));
numIonTypes = size(initialize.ecAutocorrNoAverageCesaro, 2);
timeLags = initialize.timeLags;
clear initialize;

function index = zipIndexPair(idx1, idx2)
    global numIonTypes;
%    index = (idx1 - 1) * numIonTypes + idx2;
    # accept only the "upper-half" index pair, because cross-correlation should 
    # be the same for (i,j) and (j,i)
    if (idx1 > idx2)
      error("Error - zipIndexPair: idx1 > idx2");
    else
      index = (idx1 - 1) * numIonTypes + idx2 - (idx1-1);
    endif
endfunction

# packing and distributing data to a single md array
# md(fileIndex, corrIndex, corrData)
for n = [1:numMD]
    data = load(strcat("./md", num2str(n-1), "/", dataFilename));
    md(n, 1, :) = data.ecTotalNoAverageCesaro;
    for i = [1:numIonTypes]
        for j = [i:numIonTypes]
            if (i == j)
                md(n, 1+i, :) = data.ecAutocorrNoAverageCesaro(:, i);
            endif
            crossCorrIndex = zipIndexPair(i,j);
            md(n, 1+numIonTypes+crossCorrIndex, :) = data.ecCorrNoAverageCesaro(:, crossCorrIndex);
        endfor
    endfor
%    md_ave.ec.ecTotalNoAverageCesaro += md(n).ec.ecTotalNoAverageCesaro;
%    md_ave.ec.ecAutocorrNoAverageCesaro += md(n).ec.ecAutocorrNoAverageCesaro;
%    md_ave.ec.ecCorrNoAverageCesaro += md(n).ec.ecCorrNoAverageCesaro;
endfor

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

# md_ave(corrAve, corrIndex);
# md_std(corrStd, corrIndex);
md_ave = squeeze(mean(md, 1))';
md_std = squeeze(std(md, 0, 1))';
md_err = squeeze(std(md, 0, 1) ./ sqrt(numMD))';  # standard error 

fitRange = [20, 40; 40, 60; 60, 80; 80, 100]; #ps
fitRange *= floor(1000 / skip / deltaStep); #fs (frame)

# calculate slope for each segment of md_ave
for i = [1:size(md_ave, 2)]
    for r = [1:size(fitRange, 1)]
        slope(r,i) = polyfit(timeLags(fitRange(r, 1):fitRange(r, 2)), md_ave(fitRange(r, 1):fitRange(r, 2), i), 1)(1);
    endfor
endfor

if (numMD > 1)
  # evaluate the uncertainty in the slope of the fitting line
  # reference: Numerical Recipes Chapter 15.2 (p.656)
  for i = [1:size(md_ave, 2)]
      for r = [1:size(fitRange, 1)]
          rec_sig2 = 1 ./ (md_std(fitRange(r, 1):fitRange(r, 2), i) .^ 2);
          S(r, i) = sum(rec_sig2, 1);
          Sx(r, i) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)) .* rec_sig2, 1); 
          Sxx(r, i) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)).^2 .* rec_sig2, 1); 
          Sy(r, i) = sum(md_ave(fitRange(r, 1):fitRange(r, 2), i) .* rec_sig2, 1);
          Syy(r, i) = sum(md_ave(fitRange(r, 1):fitRange(r, 2), i).^2 .* rec_sig2, 1);
          Sxy(r, i) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)) .* md_ave(fitRange(r, 1):fitRange(r, 2), i) .* rec_sig2, 1); 
      endfor
  endfor
  Delta = S .* Sxx - Sx .* Sx;

  # output slope for double check
  slope_b = (S .* Sxy - Sx .* Sy) ./ Delta
  slopeSD = sqrt(S ./ Delta);
  # slopeSD(fitRange, corrIndex)
else
  slopeSD = 0;
endif

save(strcat(dataFilename, '-ave', num2str(numMD), '.fit'), "numIonTypes", "timeLags", "md_ave", "md_std", "md_err", "slope", "slopeSD");

%# calculate slope for each segment of each md
%fitRange .*= 1000;
%for n = [1:numMD]
%    for r = [1:size(fitRange, 1)]
%        md(n).totalSlope(r) = polyfit(md.ec.timeLags(fitRange(r, 1):fitRange(r, 2)), md(n).ec.ecTotalNoAverageCesaro(fitRange(r, 1):fitRange(r, 2)), 1)(1);
%        for i = [1:numIonTypes]
%            for j = [1:numIonTypes]
%                if (i == j)
%                    md(n).autoSlope(i,r) = polyfit(md.ec.timeLags(fitRange(r, 1):fitRange(r, 2)), md(n).ec.ecAutocorrNoAverageCesaro(fitRange(r, 1):fitRange(r, 2)), i)(1);
%                   endif
%                    md(n).crossSlope(zipIndexPair(i,j),r) = polyfit(md(n).ec.timeLags(fitRange(r, 1):fitRange(r, 2)), md(n).ec.ecCorrNoAverageCesaro(fitRange(r, 1):fitRange(r, 2), zipIndexPair(i,j))(1);
%            endfor
%        endfor
%    endfor
%endfor



#numPlots = 1 + numIonTypes + numIonTypes*numIonTypes;

# standard error for selected values
errFrameInterval = floor(20000 / skip / deltaStep);
errFrames = [1:errFrameInterval:length(timeLags)]';


f1 = figure(1);
clf;
hold on;
#errOffset = floor(200 / skip / deltaStep);
errOffset = 0

for i = [1:size(md_ave, 2)]
    drawFrames = errFrames .+ (i-1)*errOffset;
    drawFrames = drawFrames(drawFrames <= length(timeLags));
    if (i == 6)
        p(i) = plot(timeLags, md_ave(:,i), "-", "color", [0, 0.6, 0]);
        ebar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i));
        set(ebar, "color", [0, 0.6, 0], "linestyle", "none");
%    elseif (i == 7)
%#        p(i).plot = plot(timeLags, md_ave(:,i), "-", "color", [0, 0, 0]);
%        p(i).errorbar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i));
%        set(p(i).errorbar(1), "color", [0, 0, 0]);
    else
        plotFormat = strcat("-", num2str(i));
        errorbarFormat = strcat("~", num2str(i));
        ebar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i), errorbarFormat);
        set(ebar, "linestyle", "none");
        p(i) = plot(timeLags, md_ave(:,i), plotFormat);
    endif
endfor

legend(p, {"Total", "Auto Na+", "Auto Cl-", "Cross Na-Na", "Cross Na-Cl", "Cross Cl-Cl"}, "location", "northwest");

left = 25;
h_space = 20;
top = 1272;
v_space = 44;
for i = [1:size(md_ave, 2)]
    for r = [1:size(fitRange, 1)]
        text(left + (r-1)*h_space, top - (i-1)*v_space, num2str(slope(r, i)));
    endfor
endfor

%set(p1, "linewidth", 3);
%set(p2, "linewidth", 3);
%set(p3, "linewidth", 3);
%set(p4, "linewidth", 3);

%title(strcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution - md", "-ave"));
title(cstrcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution\n",\
             "data interval = ", num2str(skip), ", integration interval = ", num2str(deltaStep)));
%legend("{/Symbol D}t = 1 fs");
xlabel("partial sum upper limit (ps)");
ylabel("Non-averaged partial sum (ps*S/m)");
axis([0,100,-400, 1300]);

print(strcat(dataFilename, '.fit.ave-', num2str(numMD), '.eps'), '-deps', '-color');
hold off
%axis([0,0.1,0,1.5]);
%print('ecNoAverageCesaro-zoom.eps', '-deps', '-color');

##############################

# standard deviation for selected values
for i = [1:size(fitRange, 1)]
    stdFrames(i) = (fitRange(i, 1) + fitRange(i, 2)) / 2;
endfor

f2 = figure(2);
%clf;
hold on;
%stdOffset = floor(100 / skip / deltaStep);
stdOffset = 0;

for i = [1:size(md_ave, 2)]
    drawFrames = stdFrames .+ (i-1)*stdOffset;
    if (i == 6)
        ebar = errorbar(timeLags(drawFrames), slope(:,i), slopeSD(:,i));
        set(p(i).errorbar(1), "color", [0, 0.6, 0]);
%    elseif (i == 7)
%        p(i).errorbar = errorbar(timeLags(drawFrames), slope(:,i), slopeSD(:,i));
%        set(p(i).errorbar(1), "color", [0, 0, 0]);
    else
        plotFormat = strcat("-", num2str(i));
        errorbarFormat = strcat("~", num2str(i));
        ebar = errorbar(timeLags(drawFrames), slope(:,i), slopeSD(:,i), errorbarFormat);
    endif
endfor

legend("Total", "Auto Na+", "Auto Cl-", "Cross Na-Na", "Cross Na-Cl", "Cross Cl-Cl", "location", "northwest");

%left = 25;
%h_space = 20;
%top = 9.75;
%v_space = 0.41;
%for i = [1:size(md_ave, 2)]
%    for r = [1:size(fitRange, 1)]
%        text(left + (r-1)*h_space, top - (i-1)*v_space, num2str(slope(r, i)));
%    endfor
%endfor

%set(p1, "linewidth", 3);
%set(p2, "linewidth", 3);
%set(p3, "linewidth", 3);
%set(p4, "linewidth", 3);

%title(strcat("Electrical conductivity of 1m NaCl solution - md", "-ave"));
title(cstrcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution\n",\
             "data interval = ", num2str(skip), ", integration interval = ", num2str(deltaStep)));
%legend("{/Symbol D}t = 1 fs");
xlabel("partial sum upper limit (ps)");
ylabel("Electrical conductivity (S/m)");
axis([0,100,-5, 13]);

print(strcat(dataFilename, '.fit.slope-ave-', num2str(numMD), '.eps'), '-deps', '-color');

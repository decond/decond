#!/home/kmtu/bin/octave -qf

clear all;
global numIonTypes constant;
constant.kB = 1.3806488E-23; #(J K-1)
constant.beta = 1.0/(constant.kB*300); #(J-1)
constant.basicCharge = 1.60217646E-19; #(Coulomb)
constant.ps = 1.0E-12; #(s)
constant.nm = 1.0E-9; #(m)

if (nargin() < 3)
    error("Usage: $fitAveNoAverageCesaro-ND.m <dataFilename> <numMD> <dt (fs)> <intDelta>")
endif

dataFilename = argv(){1}
numMD = str2num(argv(){2})
dt = str2num(argv(){3}) #data interval in vCorr data (unit: fs)
intStep = str2num(argv(){4})

set(0, "defaultlinelinewidth", 4);

%md.dc = load("md-lag20000.dcNoAverageCesaro" );
%md3.dc = load("md3-lag20000.dcNoAverageCesaro" );

load(strcat("./md0/", dataFilename));
%numIonTypes
%timeLags

function index = zipIndexPair2(idx1, idx2)
    global numIonTypes;
%    index = (idx1 - 1) * numIonTypes + idx2;
    # accept only the "upper-half" index pair, because cross-correlation should 
    # be the same for (i,j) and (j,i)
    if (idx1 > idx2)
      error("Error - zipIndexPair2: idx1 > idx2");
    else
      index = (idx1 - 1) * numIonTypes + idx2 - (idx1-1);
    endif
endfunction

zz(1:numIonTypes) = 1;
for i = [1:numIonTypes]
  for j = [i:numIonTypes]
    zz(numIonTypes + zipIndexPair2(i,j)) = charge(i)*charge(j);
  endfor
endfor

numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;
md = zeros(length(timeLags), numIonTypes + numIonTypePairs, numMD);
# packing and distributing data to a single md array
for n = [1:numMD]
    data_tmp = load(strcat("./md", num2str(n-1), "/", dataFilename));
    md(:, 1:numIonTypes, n) = data_tmp.autoNDNoAverageCesaro;
    md(:, numIonTypes+1:end, n) = data_tmp.crossNDNoAverageCesaro;
    volume_data(n) = prod(data_tmp.cell);
endfor
clear("data_tmp");

md_total = sum(md.*zz, 2);
for i = [1:numIonTypes]
  for j = [i+1:numIonTypes]
    idx = numIonTypes + zipIndexPair2(i,j);
    md_total += zz(idx)*md(:, idx, :);
  endfor
endfor

# NDNoAveCesaro.ave(length(timeLags), corrIndex);
# NDNoAveCesaro.std(length(timeLags), corrIndex);
if (numMD > 1)
  NDNoAveCesaro.ave = mean(md, 3);
  NDNoAveCesaro.std = std(md, 0, 3);
  NDNoAveCesaro.err = NDNoAveCesaro.std ./ sqrt(numMD);  # standard error 

  NDNoAveCesaro_total.ave = mean(md_total, 3);
  NDNoAveCesaro_total.std = std(md_total, 0, 3);
  NDNoAveCesaro_total.err = NDNoAveCesaro_total.std ./ sqrt(numMD);  # standard error 

  volume.ave = mean(volume_data);
  volume.std = std(volume_data);
  volume.err = volume.std ./ sqrt(numMD);
else
  NDNoAveCesaro.ave = md;
  NDNoAveCesaro.std = NDNoAveCesaro.ave*0;
  NDNoAveCesaro.err = NDNoAveCesaro.ave*0;

  NDNoAveCesaro_total.ave = md_total;
  NDNoAveCesaro_total.std = NDNoAveCesaro_total.ave*0;
  NDNoAveCesaro_total.err = NDNoAveCesaro_total.ave*0;

  volume.ave = volume_data;
  volume.std = volume.ave*0;
  volume.err = volume.ave*0;
endif

fitRange = [20, 40; 40, 60; 60, 80; 80, 100]*10; #ps
fitRange *= floor(1000 / dt / intStep); #fs (frame)

# calculate slope for each segment of NDNoAveCesaro.ave to get ND
for i = [1:size(NDNoAveCesaro.ave, 2)]
    for r = [1:size(fitRange, 1)]
        ND(r,i) = polyfit(timeLags(fitRange(r, 1):fitRange(r, 2)), NDNoAveCesaro.ave(fitRange(r, 1):fitRange(r, 2), i), 1)(1);
    endfor
endfor
for r = [1:size(fitRange, 1)]
    ND_total(r) = polyfit(timeLags(fitRange(r, 1):fitRange(r, 2)), NDNoAveCesaro_total.ave(fitRange(r, 1):fitRange(r, 2), 1), 1)(1);
endfor

if (numMD > 1)
  # evaluate the uncertainty in the slope of the fitting line
  # reference: Numerical Recipes Chapter 15.2 (p.656)
  for i = [1:size(NDNoAveCesaro.ave, 2)]
      for r = [1:size(fitRange, 1)]
          rec_sig2 = 1 ./ (NDNoAveCesaro.std(fitRange(r, 1):fitRange(r, 2), i) .^ 2);
          S(r, i) = sum(rec_sig2, 1);
          Sx(r, i) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)) .* rec_sig2, 1); 
          Sxx(r, i) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)).^2 .* rec_sig2, 1); 
          Sy(r, i) = sum(NDNoAveCesaro.ave(fitRange(r, 1):fitRange(r, 2), i) .* rec_sig2, 1);
          Syy(r, i) = sum(NDNoAveCesaro.ave(fitRange(r, 1):fitRange(r, 2), i).^2 .* rec_sig2, 1);
          Sxy(r, i) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)) .* NDNoAveCesaro.ave(fitRange(r, 1):fitRange(r, 2), i) .* rec_sig2, 1); 
      endfor
  endfor
  Delta = S .* Sxx - Sx .* Sx;

  # output slope for double check
  slope_b = (S .* Sxy - Sx .* Sy) ./ Delta
  ND_err = sqrt(S ./ Delta);
  # ND_err(fitRange, corrIndex)

  for r = [1:size(fitRange, 1)]
      rec_sig2_total = 1 ./ (NDNoAveCesaro_total.std(fitRange(r, 1):fitRange(r, 2), 1) .^ 2);
      S_total(r) = sum(rec_sig2_total, 1);
      Sx_total(r) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)) .* rec_sig2_total, 1); 
      Sxx_total(r) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)).^2 .* rec_sig2_total, 1); 
      Sy_total(r) = sum(NDNoAveCesaro_total.ave(fitRange(r, 1):fitRange(r, 2), 1) .* rec_sig2_total, 1);
      Syy_total(r) = sum(NDNoAveCesaro_total.ave(fitRange(r, 1):fitRange(r, 2), 1).^2 .* rec_sig2_total, 1);
      Sxy_total(r) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)) .* NDNoAveCesaro_total.ave(fitRange(r, 1):fitRange(r, 2), 1) .* rec_sig2_total, 1); 
  endfor
  Delta_total = S_total .* Sxx_total - Sx_total .* Sx_total;
  ND_total_err = sqrt(S_total ./ Delta_total);

else
  ND_err = ND*0;
  ND_total_err = ND_total*0;
endif

constant.ND2EC = constant.beta * constant.basicCharge^2 / (volume.ave*(constant.nm^3));
save(strcat(dataFilename, '-ave', num2str(numMD), '.fit'), "constant", "charge",
     "numIonTypes", "numMol", "timestep", "timeLags", "volume", "NDNoAveCesaro", "ND", "ND_err",
     "ND_total", "ND_total_err", "zz");

#For drawing:
NDNoAveCesaro.ave .*= constant.ND2EC;
NDNoAveCesaro.err .*= constant.ND2EC;
ND .*= zz*constant.ND2EC;
ND_err .*= zz*constant.ND2EC;

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
errFrameInterval = floor(20000 / dt / intStep);
errFrames = [1:errFrameInterval:length(timeLags)]';


f1 = figure(1);
clf;
hold on;
#errOffset = floor(200 / dt / intStep);
errOffset = 0;

for i = [1:size(NDNoAveCesaro.ave, 2)]
    drawFrames = errFrames .+ (i-1)*errOffset;
    drawFrames = drawFrames(drawFrames <= length(timeLags));
    if (i == 6)
        p(i) = plot(timeLags, NDNoAveCesaro.ave(:,i), "-", "color", [0, 0.6, 0]);
        ebar = errorbar(timeLags(drawFrames), NDNoAveCesaro.ave(drawFrames,i), NDNoAveCesaro.err(drawFrames,i));
        set(ebar, "color", [0, 0.6, 0], "linestyle", "none");
%    elseif (i == 7)
%#        p(i).plot = plot(timeLags, md_ave(:,i), "-", "color", [0, 0, 0]);
%        p(i).errorbar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i));
%        set(p(i).errorbar(1), "color", [0, 0, 0]);
    else
        plotFormat = strcat("-", num2str(i));
        errorbarFormat = strcat("~", num2str(i));
        ebar = errorbar(timeLags(drawFrames), NDNoAveCesaro.ave(drawFrames,i), NDNoAveCesaro.err(drawFrames,i), errorbarFormat);
        set(ebar, "linestyle", "none");
        p(i) = plot(timeLags, NDNoAveCesaro.ave(:,i), plotFormat);
    endif
endfor

legend(p, {"Auto Na+", "Auto Cl-", "Cross Na-Na", "Cross Na-Cl", "Cross Cl-Cl"}, "location", "northwest");

left = 25;
h_space = 20;
top = 1272;
v_space = 44;
for i = [1:size(NDNoAveCesaro.ave, 2)]
    for r = [1:size(fitRange, 1)]
        text(left + (r-1)*h_space, top - (i-1)*v_space, num2str(ND(r, i)));
    endfor
endfor

%set(p1, "linewidth", 3);
%set(p2, "linewidth", 3);
%set(p3, "linewidth", 3);
%set(p4, "linewidth", 3);

%title(strcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution - md", "-ave"));
title(cstrcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution\n",
             "data interval = ", num2str(dt), ", integration interval = ", num2str(intStep)));
%legend("{/Symbol D}t = 1 fs");
xlabel("partial sum upper limit (ps)");
ylabel("Non-averaged partial sum (ps*S/m)");
axis([0,100,-400, 1300]);

print(strcat(dataFilename, '-ave', num2str(numMD), '.eps'), '-deps', '-color');
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
%stdOffset = floor(100 / dt / intStep);
stdOffset = 0;

for i = [1:size(NDNoAveCesaro.ave, 2)]
    drawFrames = stdFrames .+ (i-1)*stdOffset;
    if (i == 6)
        ebar = errorbar(timeLags(drawFrames), ND(:,i), ND_err(:,i));
        set(p(i).errorbar(1), "color", [0, 0.6, 0]);
%    elseif (i == 7)
%        p(i).errorbar = errorbar(timeLags(drawFrames), ND(:,i), ND_err(:,i));
%        set(p(i).errorbar(1), "color", [0, 0, 0]);
    else
        plotFormat = strcat("-", num2str(i));
        errorbarFormat = strcat("~", num2str(i));
        ebar = errorbar(timeLags(drawFrames), ND(:,i), ND_err(:,i), errorbarFormat);
    endif
endfor

legend("Auto Na+", "Auto Cl-", "Cross Na-Na", "Cross Na-Cl", "Cross Cl-Cl", "location", "northwest");

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
title(cstrcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution\n",
             "data interval = ", num2str(dt), ", integration interval = ", num2str(intStep)));
%legend("{/Symbol D}t = 1 fs");
xlabel("partial sum upper limit (ps)");
ylabel("Electrical conductivity (S/m)");
axis([0,100,-5, 13]);

print(strcat(dataFilename, '-slope-ave', num2str(numMD), '.eps'), '-deps', '-color');

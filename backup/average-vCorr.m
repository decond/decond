#!/home/kmtu/bin/octave -qf

clear all
format long

vCorrDataName = "lag1000000.vCorr";
numData = 10;
for i = [0:numData-1]
    md(i+1) = load(strcat("md", num2str(i), "/", vCorrDataName));
endfor

timestep = md(1).timestep;
charge = md(1).charge;
numAtoms = md(1).numAtoms;
timeLags = md(1).timeLags;

vAutocorr = cell(size(md(1).vAutocorr));
vCorr = cell(size(md(2).vCorr));

vAutocorr(:) = 0;
vCorr(:) = 0;

% sum
for n = [1:numData]
    for i = [1:length(vAutocorr)]
        vAutocorr{i} = vAutocorr{i} + md(n).vAutocorr{i};
        for j = [1:length(vAutocorr)]
            vCorr{i,j} = vCorr{i,j} + md(n).vCorr{i,j};
        endfor
    endfor
endfor

% average    
for i = [1:length(vAutocorr)]
    vAutocorr{i} = vAutocorr{i} / numData;
    for j = [1:length(vAutocorr)]
        vCorr{i,j} = vCorr{i,j} / numData;
    endfor
endfor

save(strcat(vCorrDataName, "Average"), "timestep", "charge", "numAtoms", "timeLags", "vAutocorr", "vCorr");


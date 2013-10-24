#!/home/kmtu/bin/octave -qf

clear all
format long

ecCesaroDataName = "lag1000000.ecNoAverageCesaro";
numData = 10;
for i = [0:numData-1]
    md(i+1) = load(strcat("md", num2str(i), "/", ecCesaroDataName));
endfor

timestep = md(1).timestep;
timeLags = md(1).timeLags;
numIonTypes = size(md(1).ecAutocorrNoAverageCesaro, 2);

fitRange=[80001:100000];

% read fitting values
for n = [1:numData]
    for i = [1:numIonTypes]
        ecTotal(n,1) = polyfit(md(n).timeLags(fitRange), md(n).ecTotalNoAverageCesaro(fitRange), 1)(1);
        ecAutocorr(n,i) = polyfit(md(n).timeLags(fitRange), md(n).ecAutocorrNoAverageCesaro(:,i)(fitRange), 1)(1);
        for j = [1:numIonTypes]
            idx = (i-1)*numIonTypes+ j;
            ecCorr(n,idx) = polyfit(md(n).timeLags(fitRange), md(n).ecCorrNoAverageCesaro(:,idx)(fitRange), 1)(1);
        endfor
    endfor
endfor

% average    
ecTotal_average = sum(ecTotal) / numData;
ecAutocorr_average = sum(ecAutocorr) / numData;
ecCorr_average = sum(ecCorr) / numData;

% rmsd
ecTotal_rmsd = sqrt(sum((ecTotal - ecTotal_average).^2) / numData);
ecAutocorr_rmsd = sqrt(sum((ecAutocorr - ecAutocorr_average).^2) / numData);
ecCorr_rmsd = sqrt(sum((ecCorr - ecCorr_average).^2) / numData);

%save(strcat(vCorrDataName, "Average"), "timestep", "charge", "numAtoms", "timeLags", "vAutocorr", "vCorr");
save(strcat(ecCesaroDataName, ".results"), "ecTotal", "ecTotal_average", "ecTotal_rmsd",\
                                        "ecAutocorr", "ecAutocorr_average", "ecAutocorr_rmsd",\
                                        "ecCorr", "ecCorr_average", "ecCorr_rmsd");


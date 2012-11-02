#!/home/kmtu/bin/octave -qf

clear all
format long

global kB beta basicCharge ps nm volume timestep t numIonTypes;
kB = 1.3806488E-23; #(J K-1)
beta = 1.0/(kB*300); #(J-1)
basicCharge = 1.60217646E-19; #(Coulomb)
ps = 1.0E-12; #(s)
nm = 1.0E-9; #(m)

if (nargin() < 2)
    error("Usage: $calCesaroEC.m <filename.vCorr> <maxLag -1=max> <systemVolume(nm^3)> <maxLagLowerBound>\n\
<filename> is used for both input and output\n\
<maxLag> and <maxLagLowerBound> are in units of number of frames (positive integers)");
else
    filename = argv(){1};
    maxLag = str2num(argv(){2});
    volume = str2num(argv(){3}) * (1.0E-9)**3; #(m3)
    maxLagLowerBound = str2num(argv(){4})

    extnamePos = rindex(filename, "."); #locate the position of the extension name
    baseFilename = filename(1:extnamePos-1);
endif

#.vCorr file contains timestep, charge(), numAtoms(), timeLags(), vAutocorr{}, and vCorr{}
load(filename);

numIonTypes = length(vAutocorr);
if (numIonTypes != length(vCorr))
    error(strcat("Numbers of ion types are inconsistent!\n\
vAutocorr: ", num2str(length(vAutocorr)), ", vCorr: ", num2str(length(vCorr))));
endif

if (maxLag < 0)
    maxLag = length(vAutocorr{1}) - 1;
endif
maxLag #for showing

t = [0:maxLag];

function ec = integrateEC(corrData, maxLag)
    global kB beta basicCharge ps nm volume timestep t;
    if (length(corrData) == 1 && maxLag > 0)
        #there is only one ion so no mutual-corr{i}{i}
        ec = 0;
    else
        ec = beta * basicCharge**2 / volume * trapz(t(1:maxLag+1)', corrData(1:maxLag+1)) * timestep * nm**2 / ps;
    endif
endfunction

function index = zipIndexPair(idx1, idx2)
    global numIonTypes;
    index = (idx1 - 1) * numIonTypes + idx2;
endfunction

for T = [maxLagLowerBound:maxLag]
    ecTotal(T) = 0;
    for i = [1:numIonTypes]
        ecAutocorr(T, i) = charge(i) * charge(i) * integrateEC(vAutocorr{i}, T);
        ecAutocorrCesaro(T, i) = sum(ecAutocorr(maxLagLowerBound:T, i)) / (T - maxLagLowerBound + 1);
        ecTotal(T) = ecTotal(T) + ecAutocorr(T, i);
        for j = [1:numIonTypes]
            ecCorr(T,zipIndexPair(i,j)) = charge(i) * charge(j) * integrateEC(vCorr{i,j}, T);
            ecCorrCesaro(T, zipIndexPair(i,j)) = sum(ecCorr(maxLagLowerBound:T, zipIndexPair(i,j))) / (T - maxLagLowerBound + 1);
            ecTotal(T) = ecTotal(T) + ecCorr(T, zipIndexPair(i,j));
        endfor
    endfor
    ecTotalCesaro(T) = sum(ecTotal(maxLagLowerBound:T)) / (T - maxLagLowerBound + 1);
endfor
ecAutocorrCesaro = ecAutocorrCesaro(maxLagLowerBound:end, :);
ecCorrCesaro = ecCorrCesaro(maxLagLowerBound:end, :);
ecTotalCesaro = ecTotalCesaro'(maxLagLowerBound:end, :);

timeLags = [maxLagLowerBound:maxLag]' * timestep;
save(strcat(baseFilename, ".ecPartialCesaro-", num2str(maxLagLowerBound)), "timestep", "timeLags", "ecTotalCesaro", "ecAutocorrCesaro", "ecCorrCesaro");


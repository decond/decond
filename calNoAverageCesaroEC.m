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
    error("Usage: $calCesaroEC.m <filename.vCorr> <maxLag -1=max> <systemVolume(nm^3)>\n\
where <filename> is used for both input and output: filename.vCorr and filename.xvg");
else
    filename = argv(){1};
    maxLag = str2num(argv(){2});
    volume = str2num(argv(){3}) * (1.0E-9)**3; #(m3)
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

%function ec = integrateEC(corrData, maxLag)
%    global kB beta basicCharge ps nm volume timestep t;
%    if (length(corrData) == 1 && maxLag > 0)
%        #there is only one ion so no mutual-corr{i}{i}
%        ec = 0;
%    else
%        ec = beta * basicCharge**2 / volume * trapz(t(1:maxLag+1)', corrData(1:maxLag+1)) * timestep * nm**2 / ps;
%    endif
%endfunction

function index = zipIndexPair(idx1, idx2)
    global numIonTypes;
    index = (idx1 - 1) * numIonTypes + idx2;
endfunction

%for T = [0:maxLag]
%    ecTotal(T+1) = 0;
%    for i = [1:numIonTypes]
%        ecAutocorr(T+1, i) = charge(i) * charge(i) * integrateEC(vAutocorr{i}, T);
%#        ecAutocorrNoAverageCesaro(T, i) = sum(ecAutocorr(1:T, i));
%        if (T > 0)
%            ecAutocorrNoAverageCesaro(T, i) = trapz([0:T]', ecAutocorr(1:T+1, i)) * timestep;
%        endif
%        ecTotal(T+1) = ecTotal(T+1) + ecAutocorr(T+1, i);
%        for j = [1:numIonTypes]
%            ecCorr(T+1,zipIndexPair(i,j)) = charge(i) * charge(j) * integrateEC(vCorr{i,j}, T);
%#            ecCorrNoAverageCesaro(T, zipIndexPair(i,j)) = sum(ecCorr(1:T, zipIndexPair(i,j)));
%            if (T > 0)
%                ecCorrNoAverageCesaro(T, zipIndexPair(i,j)) = trapz([0:T]', ecCorr(1:T+1, zipIndexPair(i,j))) * timestep;
%            endif
%            ecTotal(T+1) = ecTotal(T+1) + ecCorr(T+1, zipIndexPair(i,j));
%        endfor
%    endfor
%#    ecTotalNoAverageCesaro(T) = sum(ecTotal(1:T));
%    if (T > 0)
%        ecTotalNoAverageCesaro(T) = trapz([0:T], ecTotal(1:T+1)) * timestep;
%    endif
%endfor
%ecTotalNoAverageCesaro = ecTotalNoAverageCesaro';

function ec = cumIntegrateEC(corrData, maxLag)
    global kB beta basicCharge ps nm volume timestep t;
    if (length(corrData) == 1 && maxLag > 0)
        #there is only one ion so no mutual-corr{i}{i}
        ec = 0;
    else
        ec = beta * basicCharge**2 / volume * cumtrapz(t(1:maxLag+1)', corrData(1:maxLag+1)) * timestep * nm**2 / ps;
    endif
endfunction

ecTotal = zeros(maxLag+1, 1); 
for i = [1:numIonTypes]
    ecAutocorr(:, i) = charge(i) * charge(i) * cumIntegrateEC(vAutocorr{i}, maxLag);
    ecAutocorrNoAverageCesaro(:, i) = cumtrapz(ecAutocorr(:, i)) * timestep;
    ecTotal .+= ecAutocorr(:, i);
    for j = [1:numIonTypes]
        ecCorr(:, zipIndexPair(i,j)) = charge(i) * charge(j) * cumIntegrateEC(vCorr{i,j}, maxLag);
        ecCorrNoAverageCesaro(:, zipIndexPair(i,j)) = cumtrapz(ecCorr(:, zipIndexPair(i,j))) * timestep;
        ecTotal .+= ecCorr(:, zipIndexPair(i,j));
    endfor
endfor
ecTotalNoAverageCesaro = cumtrapz(ecTotal) * timestep;

timeLags = [0:maxLag]' * timestep;
save(strcat(baseFilename, ".ecNoAverageCesaro-cum"), "timestep", "timeLags", "ecTotalNoAverageCesaro", "ecAutocorrNoAverageCesaro", "ecCorrNoAverageCesaro");

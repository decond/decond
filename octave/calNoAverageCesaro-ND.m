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
    error("Usage: $calCesaroEC.m <filename.vCorr> <maxLag -1=max (step)> [deltaStep] \n\
where <filename> is used for both input and output: filename.vCorr and filename.xvg,\n\
optional deltaStep is for integrating the vCorr every deltaStep, default is 1.");
else
    filename = argv(){1};
    maxLag = str2num(argv(){2});
    extnamePos = rindex(filename, "."); #locate the position of the extension name
    baseFilename = filename(1:extnamePos-1);
    if (nargin() > 2)
        deltaStep = round(str2num(argv(){3}))
        if (deltaStep <= 0)
            deltaStep = 1;
        endif
    else
        deltaStep = 1
    endif
endif

#.vCorr file contains timestep, charge(), numAtom(), timeLags(), autoCorr(), and crossCorr(), cell()
load(filename);
volume = prod(cell)*(1.0E-9)**3; #(m3)

numIonTypes = size(autoCorr, 2);
numIonTypes_check = size(crossCorr, 2);
if (numIonTypes**2 != numIonTypes_check)
    error(strcat("Numbers of ion types are inconsistent!\n\
autoCorr: ", num2str(numIonTypes), ", crossCorr: ", num2str(numIonTypes_check)));
endif

# modifying data according to deltaStep if it is greater than 1
if (deltaStep > 1)
    timestep *= deltaStep;
    timeLags = timeLags(1:deltaStep:end);
    for i = [1:numIonTypes]
        autoCorr_tmp(:,i) = autoCorr(1:deltaStep:end, i);
    endfor
    for i = [1:numIonTypes**2]
        crossCorr_tmp(:,i) = crossCorr(1:deltaStep:end, i);
    endfor
    autoCorr = autoCorr_tmp;
    crossCorr = crossCorr_tmp;
    clear("autoCorr_tmp");
    clear("crossCorr_tmp");
    if (maxLag > 0)
        maxLag = floor(maxLag / deltaStep);
    endif
endif

if (maxLag < 0)
    maxLag = size(autoCorr, 1) - 1;
endif

#for showing
timestep
maxLag 

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

function index = zipIndexPair2(idx1, idx2)
    global numIonTypes;
%    index = (idx1 - 1) * numIonTypes + idx2;
    # accept only the "upper-half" index pair, because cross-correlation should 
    # be the same for (i,j) and (j,i)
    if (idx1 > idx2)
      error("Error in zipIndexPair2: idx1 > idx2");
    else
      index = (idx1 - 1) * numIonTypes + idx2 - (idx1-1);
    endif
endfunction

%for T = [0:maxLag]
%    ecTotal(T+1) = 0;
%    for i = [1:numIonTypes]
%        ecAutocorr(T+1, i) = charge(i) * charge(i) * integrateEC(autoCorr{i}, T);
%#        ecAutocorrNoAverageCesaro(T, i) = sum(ecAutocorr(1:T, i));
%        if (T > 0)
%            ecAutocorrNoAverageCesaro(T, i) = trapz([0:T]', ecAutocorr(1:T+1, i)) * timestep;
%        endif
%        ecTotal(T+1) = ecTotal(T+1) + ecAutocorr(T+1, i);
%        for j = [1:numIonTypes]
%            ecCorr(T+1,zipIndexPair(i,j)) = charge(i) * charge(j) * integrateEC(crossCorr{i,j}, T);
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

t = [0:maxLag]';

function ec = cumIntegrateEC(corrData, maxLag)
    global kB beta basicCharge ps nm volume timestep t;
    if (length(corrData) == 1 && maxLag > 0)
        #there is only one ion so no mutual-corr{i}{i}
        ec = 0;
    else
        ec = beta * basicCharge**2 / volume * cumtrapz(t(1:maxLag+1), corrData(1:maxLag+1)) * timestep * nm**2 / ps;
    endif
endfunction

function ND = cumIntegrateND(corrData, maxLag)
    global ps nm timestep t;
    if (length(corrData) == 1 && maxLag > 0)
        #there is only one ion so no mutual-sdcorr
        ND = 0;
    else
        ND = cumtrapz(t(1:maxLag+1), corrData(1:maxLag+1)) * timestep * nm**2 / ps;
    endif
endfunction

%ecTotal = zeros(maxLag+1, 1); 
%for i = [1:numIonTypes]
%    ecAutocorr(:, i) = charge(i) * charge(i) * cumIntegrateEC(autoCorr{i}, maxLag);
%    ecAutocorrNoAverageCesaro(:, i) = cumtrapz(ecAutocorr(:, i)) * timestep;
%    ecTotal .+= ecAutocorr(:, i);
%    for j = [1:numIonTypes]
%        ecCorr(:, zipIndexPair(i,j)) = charge(i) * charge(j) * cumIntegrateEC(crossCorr{i,j}, maxLag);
%        ecCorrNoAverageCesaro(:, zipIndexPair(i,j)) = cumtrapz(ecCorr(:, zipIndexPair(i,j))) * timestep;
%        ecTotal .+= ecCorr(:, zipIndexPair(i,j));
%    endfor
%endfor
for i = [1:numIonTypes]
    autoND(:, i) = cumIntegrateND(autoCorr(:,i), maxLag);
    autoNDNoAverageCesaro(:, i) = cumtrapz(autoND(:, i)) * timestep;
%    ecTotal .+= autoD(:, i);
    for j = [i:numIonTypes]
        idx(1) = zipIndexPair(i, j);
        idx(2) = zipIndexPair(j, i);
        idx2 = zipIndexPair2(i, j);
        if (idx(1) != idx(2))
          ND = (crossCorr(:, idx(1)) + crossCorr(:, idx(2))) / 2;
        else
          ND = crossCorr(:, idx(1));
        endif
        crossND(:, idx2) = cumIntegrateND(ND, maxLag);
        crossNDNoAverageCesaro(:, idx2) = cumtrapz(crossND(:, idx2)) * timestep;
%        if (i == j)
%          ecTotal .+= ecCorr(:, zipIndexPair2(i,j));
%        else
%          ecTotal .+= ecCorr(:, zipIndexPair2(i,j)).*2;
%        endif
    endfor
endfor
%ecTotalNoAverageCesaro = cumtrapz(ecTotal) * timestep;

timeLags = [0:maxLag]' * timestep;
save(strcat(baseFilename, ".NDNoAverageCesaro-dt-", num2str(deltaStep)), "charge", "numIonTypes", "numAtom", "cell", "timestep", "timeLags", "autoNDNoAverageCesaro", "crossNDNoAverageCesaro");

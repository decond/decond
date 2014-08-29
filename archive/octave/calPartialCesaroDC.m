#!/home/kmtu/bin/octave -qf

clear all
format long

global kB beta basicCharge ps nm volume timestep t;
kB = 1.3806488E-23; #(J K-1)
beta = 1.0/(kB*300); #(J-1)
basicCharge = 1.60217646E-19; #(Coulomb)
ps = 1.0E-12; #(s)
nm = 1.0E-9; #(m)

if (nargin() < 3)
    error("Usage: $calPartialCesaroDC.m <filename.vCorr> <maxLag -1=max> <maxLagLowerBound>\n\
<filename> is used for both input and output\n\
<maxLag> and <maxLagLowerBound> are in units of number of frames (positive integers)");
else
    filename = argv(){1};
    maxLag = str2num(argv(){2});
    maxLagLowerBound = str2num(argv(){3})

    extnamePos = rindex(filename, "."); #locate the position of the extension name
    baseFilename = filename(1:extnamePos-1);
endif

#.vCorr file contains timestep, charge(), numAtoms(), timeLags(), vAutocorr{}, and vCorr{}
load(filename, "timestep", "numAtoms", "vAutocorr");

numIonTypes = length(vAutocorr);

if (maxLag < 0)
    maxLag = length(vAutocorr{1}) - 1;
endif
maxLag #for showing
t = [0:maxLag]';

function dc = integrateDC(corrData, maxLag)
    global kB beta basicCharge ps nm volume timestep t;
    if (length(corrData) == 1 && maxLag > 0)
        #there is only one ion so no mutual-corr{i}{i}
        dc = 0;
    else
        dc = trapz(t(1:maxLag+1), corrData(1:maxLag+1)) * timestep * nm**2 / ps;
    endif
endfunction


for T = [maxLagLowerBound:maxLag]
    for i = [1:numIonTypes]
        dc(T, i) = integrateDC(vAutocorr{i} / numAtoms(i), T);
        dcCesaro(T,i) = sum(dc(maxLagLowerBound:T, i)) / (T - maxLagLowerBound + 1);
    endfor
endfor

dcCesaro = dcCesaro(maxLagLowerBound:end, :)

timeLags = [maxLagLowerBound:maxLag]' * timestep;
save(strcat(baseFilename, ".dcPartialCesaro-", num2str(maxLagLowerBound)), "timestep", "timeLags", "dcCesaro");


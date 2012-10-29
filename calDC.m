#!/home/kmtu/bin/octave -qf

clear all
format long

global ps nm timestep t maxLag numAtoms;

ps = 1.0E-12; #(s)
nm = 1.0E-9; #(m)

if (nargin() < 2)
    error("Usage: $calDC.m <filename.vCorr> <maxLag -1=max>\n\
where <filename> is used for both input and output: filename.vCorr and filename.ec");
else
    filename = argv(){1};
    extnamePos = rindex(filename, "."); #locate the position of the extension name
    baseFilename = filename(1:extnamePos-1);
    maxLag = str2num(argv(){2});
endif

#.vCorr file contains timestep, charge{}, vAutocorr{}, and vCorr{}
load(filename, "timestep", "numAtoms", "vAutocorr");

numIonTypes = length(vAutocorr);

if (maxLag < 0)
    maxLag = length(vAutocorr{1}) - 1;
endif
maxLag #for showing

t = [0:maxLag];

function dc = integrateDC(corrData)
    global ps nm timestep t maxLag numAtoms;
    if (length(corrData) == 1 && maxLag > 0)
        #there is only one ion so no mutual-corr{i}{i}
        dc = 0;
    else
        dc = trapz(t', corrData(1:maxLag+1)) * timestep * nm**2 / ps;
    endif
endfunction

for i = [1:numIonTypes]
    dc(i) = integrateDC(vAutocorr{i} / numAtoms(i));
endfor

dc = dc'
save(strcat(baseFilename, ".dc"), "dc");


#!/home/kmtu/local/octave-3.6.2/bin/octave -qf

clear all
format long
pkg load signal;

ps = 1.0E-12; #(s)
nm = 1.0E-9; #(m)


if (nargin() < 2)
    error("Usage: $calculateEC.m <velData.binary> <outFileName> [<maxLag>]")
else
    vFileName = argv(){1};
    outFileName = argv(){2};
    if (nargin() > 2)
        maxLag = str2num(argv(){3}) # in the unit of frame number
    else
        maxLag = -1;
    endif
endif

data = readGmx2Matlab(vFileName);

numFrames = data.num_frames
timeStep = data.time_step
if (maxLag < 0)
    maxLag = idivide(numFrames, 2)
endif

## data.trajectory(atoms, dimension, frames) 
## squeeze(data.trajectory(:,1,:))' = (frames, atoms) in x-dimension
corrData = 0
for dim = [1:3]
    vData = squeeze(data.trajectory(:,dim,:))'; #(nm / ps)
    for i = [1:data.num_atoms]
        if dim == 1 && i == 1
            [tmp, t] = xcorr(vData(:,i), maxLag, "unbiased");
            corrData = corrData + tmp;
        else
            corrData = corrData + xcorr(vData(:,i), maxLag, "unbiased");
        endif
    endfor
endfor

t = t(maxLag + 1:end);
corrData = corrData(maxLag + 1:end) / data.num_atoms;
diffConst = trapz(t', corrData) * timeStep * nm**2 / ps / 3

size(t)
#numTimeLag
numFrames

save(strcat(outFileName, ".corr"), "corrData");
save(strcat(outFileName, ".diff"), "diffConst");


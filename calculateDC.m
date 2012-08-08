#!/home1/kmtu/local/octave-3.6.2/bin/octave -qf

clear all
format long
pkg load signal;

ps = 1.0E-12 #(s)
nm = 1.0E-9 #(m)


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

## x-dimension
## data.trajectory(atoms, dimension, frames) 
## squeeze(data.trajectory(:,1,:))' = (frames, atoms) in x-dimension
vData = squeeze(data.trajectory(:,1,:))'; #(nm / ps)
for i = [1:data.num_atoms]
    [corrData(:,i), t] = xcorr(vData(:,i), maxLag, "unbiased");
endfor

t = t(maxLag + 1:end);
corrData = corrData(maxLag + 1:end,:);
corrData = sum(corrData, 2) / data.num_atoms;
diffConst = trapz(t', corrData) * timeStep * nm**2 / ps

size(t)
#numTimeLag
numFrames

save(strcat(outFileName, ".corr"), "corrData");
save(strcat(outFileName, ".diff"), "diffConst");

################
#{
for i = indexes
    fileName = strcat(baseFileName, num2str(i, "%02d"))
    data = load(fileName)
    if (!exist("xdata"))
        xdata= data(:,1)
        ydata = data(:,2)
    else
#        if (size(xdata) != size(data) || xdata != data(:,1))
#            error("The xdata indexes of file %s is different", fileName)
#        endif
    # don't check the consistency of xdata
    # use the first xdata
        minSize = min(size(ydata, 1), size(data, 1))
        ydata(1:minSize, end+1) = data(1:minSize, 2)
    endif
endfor

if (nargin() > 3)
    switch (mode)
    case "even"
        ydata_std(1:2:j) = 0
    case "odd"
        ydata_std(2:2:j) = 0
    case "custom"
        #find the closest xdata, both its value and index
        for i = [1:length(errbarXPos)]
            [m, mi] = min(abs(errbarXPos(i) - xdata))
            errbarIndexes(i) = mi
        endfor
        #assign values of OTHER indexes as 0
        for i = [1:length(xdata)]
            if (all(errbarIndexes != i))
                ydata_std(i) = 0
            endif
        endfor
    otherwise
        error("Unknown mode: %s", mode)
    endswitch
endif

out = [xdata, ydata_mean, ydata_std]
outFileName = strcat("rdf_std_", num2str(ibegin), "-", num2str(iend))
save(outFileName,"out")
#}

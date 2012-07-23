#!/home1/kmtu/local/octave-3.6.2/bin/octave -qf

clear all
format long
pkg load signal;

kB = 1.3806488E-23 #(J K-1)
beta = 1.0/(kB*300) #(J-1)
basicCharge = 1.60217646E-19 #(Coulomb)
ps = 1.0E-12 #(s)
nm = 1.0E-9 #(m)

volume = (3.10913E-9)**3 #(m3)

if (nargin() < 5)
    error("Usage: $calculateEC.m <velData1.binary> <charge1> <velData2.binary> \
           <charge2> <outFileName> [<maxLag>]")
else
    v1FileName = argv(){1};
    charge1 = str2num(argv(){2});
    v2FileName = argv(){3};
    charge2 = str2num(argv(){4});
    outFileName = argv(){5};
    if (nargin() > 5)
        maxLag = str2num(argv(){6}) # in the unit of frame number
    else
        maxLag = -1;
    endif
endif

data1 = readGmx2Matlab(v1FileName);
data2 = readGmx2Matlab(v2FileName);

if (data1.num_frames != data2.num_frames)
    error(strcat("Numbers of frames are different between ", v1FileName, " and ", v2FileName, "."))
endif
numFrames = data1.num_frames
timeStep = data1.time_step
if (maxLag < 0)
    maxLag = numFrames - 1
endif

## x-dimension
## data.trajectory(atoms, dimension, frames) 
## squeeze(data1.trajectory(:,1,:))' = (frames, atoms) in x-dimension
jData1 = charge1*squeeze(data1.trajectory(:,1,:))' * nm / ps;
jData2 = charge2*squeeze(data2.trajectory(:,1,:))' * nm / ps;
[corrData, t] = xcorr([jData1, jData2], maxLag, "unbiased");
numIons = data1.num_atoms + data2.num_atoms

#numTimeLag = 2 * numFrames - 1;
t = t(maxLag + 1:end);
corrData = corrData(maxLag + 1:end,:);
corrData = sum(corrData, 2);
elecCond = beta * basicCharge**2 / volume * trapz(t', corrData) * timeStep * ps

size(t)
#numTimeLag
numFrames

save(strcat(outFileName, ".corr"), "corrData");
save(strcat(outFileName, ".elec"), "elecCond");

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

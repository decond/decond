#!/home/kmtu/bin/octave -qf

clear all
format long
pkg load signal;

kB = 1.3806488E-23; #(J K-1)
beta = 1.0/(kB*300); #(J-1)
basicCharge = 1.60217646E-19; #(Coulomb)
ps = 1.0E-12; #(s)
nm = 1.0E-9; #(m)


if (nargin() < 4)
    error("Usage: $calculateEC.m <velData.binary> <charge> <outFileName> <systemVolume(nm^3)> [<maxLag>]")
else
    vFileName = argv(){1};
    charge = str2num(argv(){2});
    outFileName = argv(){3};
    volume = str2num(argv(){4}) * (1.0E-9)**3; #(m3)
    if (nargin() > 4)
        maxLag = str2num(argv(){5}) # in the unit of frame number
    else
        maxLag = -1;
    endif
endif

data = readGmx2Matlab(vFileName);

numFrames = data.num_frames
timeStep = data.time_step
if (maxLag < 0)
    maxLag = numFrames - 1
endif

## data.trajectory(atoms, dimension, frames) 
## squeeze(data1.trajectory(:,1,:))' = (frames, atoms) in x-dimension
autocorrData = 0;
corrData = 0;
for dim = [1:3]
    puts(cstrcat("dim = ", num2str(dim), "\n"));
    if data.num_atoms == 1
        #don't transpose because after squeeze it becomes (frames, 1) directly
        jData = charge*squeeze(data.trajectory(:,dim,:)); #(nm / ps)
    else
        jData = charge*squeeze(data.trajectory(:,dim,:))'; #(nm / ps)
    endif
    
    #calculate each diffusion term separately
    for i = [1:data.num_atoms]
        #self-diffusion (autocorrelation) for data
        puts(cstrcat("calculating autocorrData of dimension ", num2str(dim), " for atom ", num2str(i), "\n"));
        if dim == 1 && i == 1
            #get t
            [autocorrData, t] = xcorr(jData(:,i), maxLag, "unbiased");
        else
            autocorrData = autocorrData + xcorr(jData(:,i), maxLag, "unbiased");
        endif
        for j = [1:data.num_atoms]
            if i != j
                #mutual-diffusion for data1-data1
                puts(cstrcat("calculating corrData of dimension ", num2str(dim), " for atom-pair(", num2str(i), ", ", num2str(j), ")\n"));
                corrData = corrData + xcorr(jData(:,i), jData(:,j), maxLag, "unbiased");
            endif
        endfor
    endfor
endfor

corrDataTotal = autocorrData + corrData;

t = t(maxLag + 1:end);

corrDataTotal = corrDataTotal(maxLag + 1:end);
total = beta * basicCharge**2 / volume * trapz(t', corrDataTotal) * timeStep * nm**2 / ps / 3;
save(strcat(outFileName, ".corrDataTotal"), "corrDataTotal");

if data.num_atoms > 1
    autocorrData = autocorrData(maxLag + 1:end);
    corrData = corrData(maxLag + 1:end);

    autocorr = beta * basicCharge**2 / volume * trapz(t', autocorrData) * timeStep * nm**2 / ps / 3;
    corr = beta * basicCharge**2 / volume * trapz(t', corrData) * timeStep * nm**2 / ps / 3;

    save(strcat(outFileName, ".econd"), "total", "autocorr", "corr");
    save(strcat(outFileName, ".autocorrData"), "autocorrData");
    save(strcat(outFileName, ".corrData"), "corrData");
else
    save(strcat(outFileName, ".econd"), "total");
endif

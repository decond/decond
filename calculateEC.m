#!/home1/kmtu/local/octave-3.6.2/bin/octave -qf

clear all
format long
pkg load signal;

kB = 1.3806488E-23; #(J K-1)
beta = 1.0/(kB*300); #(J-1)
basicCharge = 1.60217646E-19; #(Coulomb)
ps = 1.0E-12; #(s)
nm = 1.0E-9; #(m)


if (nargin() < 6)
    error("Usage: $calculateEC.m <velData1.binary> <charge1> <velData2.binary> \
<charge2> <outFileName> <systemVolume(nm^3)> [<maxLag>]")
else
    v1FileName = argv(){1};
    charge1 = str2num(argv(){2});
    v2FileName = argv(){3};
    charge2 = str2num(argv(){4});
    outFileName = argv(){5};
    volume = str2num(argv(){6}) * (1.0E-9)**3; #(m3)
    if (nargin() > 6)
        maxLag = str2num(argv(){7}) # in the unit of frame number
    else
        maxLag = -1;
    endif
endif

data1 = readGmx2Matlab(v1FileName);
data2 = readGmx2Matlab(v2FileName);

if (data1.num_frames != data2.num_frames)
    error(cstrcat("Numbers of frames are different between ", v1FileName, " and ", v2FileName, "."))
endif
numFrames = data1.num_frames
timeStep = data1.time_step
if (maxLag < 0)
    maxLag = numFrames - 1
endif

## data.trajectory(atoms, dimension, frames) 
## squeeze(data1.trajectory(:,1,:))' = (frames, atoms) in x-dimension
autocorrData1 = 0;
autocorrData2 = 0;
corrData11 = 0;
corrData22 = 0;
corrData12 = 0;
corrData21 = 0;
for dim = [1:3]
    puts(cstrcat("dim = ", num2str(dim), "\n"));
    jData1 = charge1*squeeze(data1.trajectory(:,dim,:))'; #(nm / ps)
    jData2 = charge2*squeeze(data2.trajectory(:,dim,:))'; #(nm / ps)
    
#    #calculate total diffusion directly
#    puts("calculating corrData\n");
#    if dim == 1
#        [corrData, t] = xcorr([jData1, jData2], maxLag, "unbiased");
#    else
#        corrData = corrData + xcorr([jData1, jData2], maxLag, "unbiased");
#    endif

    #calculate each diffusion term separately
    for i = [1:data1.num_atoms]
        #self-diffusion (autocorrelation) for data1
        puts(cstrcat("calculating autocorrData1 of dimension ", num2str(dim), " for atom ", num2str(i), "\n"));
        if dim == 1 && i == 1
            #get t
            [autocorrData1, t] = xcorr(jData1(:,i), maxLag, "unbiased");
        else
            autocorrData1 = autocorrData1 + xcorr(jData1(:,i), maxLag, "unbiased");
        endif
        for j = [1:data1.num_atoms]
            if i != j
                #mutual-diffusion for data1-data1
                puts(cstrcat("calculating corrData11 of dimension ", num2str(dim), " for atom-pair(", num2str(i), ", ", num2str(j), ")\n"));
                corrData11 = corrData11 + xcorr(jData1(:,i), jData1(:,j), maxLag, "unbiased");
            endif
        endfor
        for j = [1:data2.num_atoms]
            #mutual-diffusion for data1-data2
            puts(cstrcat("calculating corrData12 of dimension ", num2str(dim), " for atom-pair(", num2str(i), ", ", num2str(j), ")\n"));
            corrData12 = corrData12 + xcorr(jData1(:,i), jData2(:,j), maxLag, "unbiased");
        endfor
    endfor
    for i = [1:data2.num_atoms]
        #self-diffusion (autocorrelation) for data2
        puts(cstrcat("calculating autocorrData2 of dimension ", num2str(dim), " for atom ", num2str(i), "\n"));
        autocorrData2 = autocorrData2 + xcorr(jData2(:,i), maxLag, "unbiased");
        for j = [1:data2.num_atoms]
            if i != j
                #mutual-diffusion for data2-data2
                puts(cstrcat("calculating corrData22 of dimension ", num2str(dim), " for atom-pair(", num2str(i), ", ", num2str(j), ")\n"));
                corrData22 = corrData22 + xcorr(jData2(:,i), jData2(:,j), maxLag, "unbiased");
            endif
        endfor
        for j = [1:data1.num_atoms]
            #mutual-diffusion for data2-data1
            puts(cstrcat("calculating corrData21 of dimension ", num2str(dim), " for atom-pair(", num2str(i), ", ", num2str(j), ")\n"));
            corrData21 = corrData21 + xcorr(jData2(:,i), jData1(:,j), maxLag, "unbiased");
        endfor
    endfor
endfor

corrDataTotal = autocorrData1 + autocorrData2 + corrData11 + corrData22 + corrData12 + corrData21;


#numTimeLag = 2 * numFrames - 1;

t = t(maxLag + 1:end);

#corrData = corrData(maxLag + 1:end,:);
#corrData = sum(corrData, 2);
#elecCond = beta * basicCharge**2 / volume * trapz(t', corrData) * timeStep * nm**2 / ps / 3;

corrDataTotal = corrDataTotal(maxLag + 1:end);
autocorrData1 = autocorrData1(maxLag + 1:end);
autocorrData2 = autocorrData2(maxLag + 1:end);
corrData11 = corrData11(maxLag + 1:end);
corrData22 = corrData22(maxLag + 1:end);
corrData12 = corrData12(maxLag + 1:end);
corrData21 = corrData21(maxLag + 1:end);
corrData1221 = corrData12 + corrData21;


total = beta * basicCharge**2 / volume * trapz(t', corrDataTotal) * timeStep * nm**2 / ps / 3;
autocorr1 = beta * basicCharge**2 / volume * trapz(t', autocorrData1) * timeStep * nm**2 / ps / 3;
autocorr2 = beta * basicCharge**2 / volume * trapz(t', autocorrData2) * timeStep * nm**2 / ps / 3;
corr11 = beta * basicCharge**2 / volume * trapz(t', corrData11) * timeStep * nm**2 / ps / 3;
corr22 = beta * basicCharge**2 / volume * trapz(t', corrData22) * timeStep * nm**2 / ps / 3;
corr12 = beta * basicCharge**2 / volume * trapz(t', corrData12) * timeStep * nm**2 / ps / 3;
corr21 = beta * basicCharge**2 / volume * trapz(t', corrData21) * timeStep * nm**2 / ps / 3;
corr1221 = corr12 + corr21;

#save(strcat(outFileName, ".econd"), "elecCond");
#save(strcat(outFileName, ".corr"), "corrData");

save(strcat(outFileName, ".econd"), "total", "autocorr1", "autocorr2", "corr11", "corr22", "corr1221", "corr12", "corr21");
save(strcat(outFileName, ".corrDataTotal"), "corrDataTotal");
save(strcat(outFileName, ".autocorrData1"), "autocorrData1");
save(strcat(outFileName, ".autocorrData2"), "autocorrData2");
save(strcat(outFileName, ".corrData11"), "corrData11");
save(strcat(outFileName, ".corrData22"), "corrData22");
save(strcat(outFileName, ".corrData12"), "corrData12");
save(strcat(outFileName, ".corrData21"), "corrData21");
save(strcat(outFileName, ".corrData1221"), "corrData1221");


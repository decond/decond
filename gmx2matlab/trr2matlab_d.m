%% trr2matlab.m  
% by Evan Arthur, University of Michigan, October 2011

% Matlab outputs trajectories in a relatively consistent format that is
% fundamentally challanging and inefficient for Matlab to read directly. 
% This program translates most trr files from recent versions of Gromacs  
% into binary files that can be read quickly and efficiently into Matlab
% via readGmx2Matlab.m. 
% readGmx2Matlab.m is a sibling program that reads the output of 
% this program. Currently only coordinates, velocities, and forces are 
% output. If I get requests for other outputs (box dimensions, lambda 
% parameters, etc) I'll fit it in. 

% Requirements: 
%    - Gromacs trr trajectory file (GMX trn format)
%         tested on version 4.5 and later
%    - readGmx2Matlab.m  (reads the output from this script)
%    - Free RAM: not much. Less than 500 kb for most simulations.
%    - Free Hard Disk: between 1 to 2 times the .trr you input.
%         By default the entire trajectory is copied and reformatted. It
%         takes the output, converts it into a usable format, and then it
%         rewrites the output with useful information. Temp files are 
%         removed after all calculations are done, so at most the
%         trajectory is just duplicated in a cleaner format.

% Limitations:
%    - Broken trr files. If there is a broken frame, it probably should be
%         removed before inputting the trajectory.

% Inputs:
%     - path to trr file (required, must be first input)
%     - path to output file (optional)
%         if none given, a default name is chosen (such as 'xdata.binary')
%     - 'x' (optional)
%         outputs xyz atomic coordinates
%     - 'v' (optional)
%         outputs xyz of instantaneous velocities
%     - 'f' (optional)
%         outputs xyz of atomic forces

% Outputs:
%     - xyz data 
%         output either by default or if 'x' option is given
%         default name is 'xdata.binary'
%     - velocity data
%         output either by default or if 'v' option is given
%         default name is 'vdata.binary'
%     - force data
%         output either by default or if 'f' option is given
%         default name is 'fdata.binary'

% Example inputs and outputs:
%     trr2matlab ('traj.trr')
%           outputs all atomic coordinates, velocities, and forces as files
%           'xdata.binary', 'vdata.binary', and 'fdata.binary'
%     trr2matlab ('traj.trr', 'x', 'f')
%           outputs all atomic coordinates and forces as files
%           'xdata.binary' and 'fdata.binary' (velocity data is not output)
%     trr2matlab ('traj.trr', 'x')
%           outputs only atomic coordinates as file 'xdata.binary' 
%           (velocity and force data are not output)
%     trr2matlab ('traj.trr', 'f', 'proteinA')
%           outputs only atomic forces as file 'proteinA_xdata.binary' 
%           (velocity and coordinates data are not output)

% notes on Single/Double Precision:
%     This program detects the precision of the trr file automatically. If
%     the detection fails, write in the input 'single' or 'double'. The
%     program outputs garbage or fails spectacularly when this is not done
%     properly. 
%     
%     Since single precision exceeds the margin of error for most analyses,
%     this program only outputs data as single-precision numbers. If I get
%     requests for the excrutiatingly accurate double-precision output, I
%     will put in an option for it.
%     

%% set file inputs and outputs

function trr2matlab( trajFile, varargin )

tic;

% check if input exists
if exist(trajFile, 'file')
    fprintf('\n    Reading trajectory %s... \n', trajFile);
else
    fprintf('\n    Trajectory file not found, exiting program \n');
    return;
end

% decide if output name was in input, and output frequency
outputX = 'xdata.binary';
outputV = 'vdata.binary';
outputF = 'fdata.binary';

noX = varargin(~strcmp(varargin, 'x'));
noV = noX(~strcmp(noX, 'v'));
noF = noV(~strcmp(noV, 'f'));

outputEveryIntFrames = 50;
for n = 1:numel(noF)
    if isfloat(cell2mat(noF(n))) == true
        outputEveryIntFrames = floor(cell2mat(noF(n)));
    end
    if ischar(cell2mat(noF(n))) == true
        outputName = cell2mat(noF(n));
        outputX = [outputName '_xdata.binary'];
        outputV = [outputName '_vdata.binary'];
        outputF = [outputName '_fdata.binary'];
    end
end
if outputEveryIntFrames < 1
    outputEveryIntFrames = 1;
end
fprintf('outputting data every %g frames\n', outputEveryIntFrames);

% decide what files to output with booleans
writeX = false;
writeV = false;
writeF = false;

if sum(strcmp(varargin, 'x')) > 0
    writeX = true;
end
if sum(strcmp(varargin, 'v')) > 0
    writeV = true;
end
if sum(strcmp(varargin, 'f')) > 0
    writeF = true;
end
% by default, write out everything
if writeX + writeV + writeF == 0
    writeX = true;
    writeV = true;
    writeF = true;
end

%% determine single/double precision

% user override
if sum(strcmp(varargin, 'single')) > 0
    fileType = 'single';
end
if sum(strcmp(varargin, 'double')) > 0
    fileType = 'double';
end

% detect single/double
if sum(strcmp(varargin, 'single')) + sum(strcmp(varargin, 'double')) == 0
    
    OPENfile = fopen(trajFile, 'r', 'b');
    fseek(OPENfile, 0 ,'bof');
    precisionTest = fread(OPENfile, [1 9], 'int32');
    
    if precisionTest(9) == 36
        fprintf('single precision detected\n');
        fileType = 'single';
    elseif precisionTest(9) == 72
        fprintf('double precision detected\n');
        fileType = 'double';
    else
        fprintf('no precision dectected, defaulting to single precision\n');
        fileType = 'single';
    end
    fclose(OPENfile);
    
end

%% parse binary from trr file

% open input trr
OPENfile = fopen(trajFile, 'r', 'b');
fseek(OPENfile, 0 ,'bof');

% remove semi-processed trajectories in same folder
system('rm -f tempx tempv tempf');

%prepare output files
WRITEfile_X = fopen('tempx', 'w', 'b');
    fwrite(WRITEfile_X, 0, 'int32');
    fwrite(WRITEfile_X, 0, 'int32');
    fwrite(WRITEfile_X, 0, 'double');

WRITEfile_V = fopen('tempv', 'w', 'b');
    fwrite(WRITEfile_V, 0, 'int32');
    fwrite(WRITEfile_V, 0, 'int32');
    fwrite(WRITEfile_V, 0, 'double');

WRITEfile_F = fopen('tempf', 'w', 'b');
    fwrite(WRITEfile_F, 0, 'int32');
    fwrite(WRITEfile_F, 0, 'int32');
    fwrite(WRITEfile_F, 0, 'double');

% initialize local variables - throws error for traj > 500,000 frames
maxNumFrames = 500000; % increase this if more than this # frames in sim!
frame_num = 1;
coord_frame = 1;
coord_time = zeros([1, maxNumFrames]);
veloc_frame = 1;
veloc_time = zeros([1, maxNumFrames]);
force_frame = 1;
force_time = zeros([1, maxNumFrames]);

data_present_words = {'coordinate  ', 'velocity  ', 'force  '};

% fix offset for end of file
empty_spacing = fread(OPENfile, [1], '*char');

while ~feof(OPENfile)
    
    intro_words = fread(OPENfile, [1 51], '*char');
    data_present = fread(OPENfile, [1 3], 'int32');
    num_atoms = fread(OPENfile, 1, 'int32');
    frame_step = fread(OPENfile, [1 2], 'int32');
    frame_time = fread(OPENfile, 1, fileType);
    box_params = fread(OPENfile, [1 10], fileType);
    
    % display statistics in periodic intervals
    if ~mod(frame_num, outputEveryIntFrames)
    fprintf(['\nopening file "%s",\n        ...%s\n' ...
         '    frame number %g contains %g atoms,\n'...
         '    and is located at time step %u\n' ...
         '    frame is set at %g ps from beginning of simulation\n' ...
         '    box dimensions are %g by %g by %g nm\n' ...
         '  %sdata present in this frame\n' ], ...
         trajFile, intro_words(12:25), frame_num, num_atoms, frame_step(1), ...
         frame_time, box_params(2), box_params(6), box_params(10), ... 
         cell2mat(data_present_words(data_present ~= 0)));
    end
     
    % coordinate data
    if (data_present(1) ~= 0)
        
        coord_frame = coord_frame + 1;
        coord_time(frame_num) = frame_time;
        
        coord_XYZ = fread(OPENfile, [1 (num_atoms * 3)], fileType);
        if writeX == true
            fwrite(WRITEfile_X, coord_XYZ, 'double');
        end
    end
    
    % velocity data
    if (data_present(2) ~= 0)
        
        veloc_frame = veloc_frame + 1;
        veloc_time(frame_num) = frame_time;
        
        velocity_XYZ = fread(OPENfile, [1 (num_atoms * 3)], fileType);
        if writeV == true
            fwrite(WRITEfile_V, velocity_XYZ, 'double');
        end
    end
    
    % force data
    if (data_present(3) ~= 0)
        
        force_frame = force_frame + 1;
        force_time(frame_num) = frame_time;
        
        force_XYZ = fread(OPENfile, [1 (num_atoms * 3)], fileType);
        if writeF == true
            fwrite(WRITEfile_F, force_XYZ, 'double');
        end
    end
    
    frame_num = frame_num + 1;
    % fix offset for end of file
    empty_spacing = fread(OPENfile, [1], '*char');
    
end

fclose(OPENfile);


%% finish processing binary output, delete temporary files

% xyz coord output
if writeX == true
    fprintf('\ncorrecting intro binaries of xyz data %s\n', outputX);
    
    coord_frame = coord_frame - 1;
    coord_incriment = coord_time(2) - coord_time(1);
    
    frewind(WRITEfile_X);
    fwrite(WRITEfile_X, num_atoms, 'int32');
    fwrite(WRITEfile_X, coord_frame, 'int32');
    fwrite(WRITEfile_X, coord_incriment, 'double');
    
    fclose(WRITEfile_X);
    system(['mv -f tempx ' outputX]);
else
    fclose(WRITEfile_X);
    system('rm -f tempx');
end


% velocity output
if writeV == true
    fprintf('\ncorrecting intro binaries of velocity data %s\n', outputV);
    
    veloc_frame = veloc_frame - 1;
    veloc_incriment = veloc_time(2) - veloc_time(1);
    
    frewind(WRITEfile_V);
    fwrite(WRITEfile_V, num_atoms, 'int32');
    fwrite(WRITEfile_V, veloc_frame, 'int32');
    fwrite(WRITEfile_V, veloc_incriment, 'double');
    
    fclose(WRITEfile_V);
    system(['mv -f tempv ' outputV]);
else
    fclose(WRITEfile_V);
    system('rm -f tempv');
end

% force output
if writeF == true
    fprintf('\ncorrecting intro binaries of force data %s\n', outputF);
    
    force_frame = force_frame - 1;
    force_incriment = force_time(2) - force_time(1);
    
    frewind(WRITEfile_F);
    fwrite(WRITEfile_F, num_atoms, 'int32');
    fwrite(WRITEfile_F, force_frame, 'int32');
    fwrite(WRITEfile_F, force_incriment, 'double');
    
    fclose(WRITEfile_F);
    system(['mv -f tempf ' outputF]);
else
    fclose(WRITEfile_F);
    system('rm -f tempf');
end

timeSpent = toc;
fprintf('\n%g seconds spent processing trajectory %s\n\n', timeSpent, trajFile);

end

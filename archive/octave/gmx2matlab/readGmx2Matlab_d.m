%% readGmx2Matlab.m  
% by Evan Arthur, University of Michigan, October 2011

% Matlab outputs trajectories in a relatively consistent format that is
% fundamentally challanging and inefficient for Matlab to read directly. 
% This program turns the output from trr2matlab.m into matricies for other
% programs to read. These are by default in a ".binary format". The matrix
% has introductory code, and the trajectory. There are options to read only
% a small portion of the trajectory with a starting frame and ending frame
% option. Skipping frames during the reading process (say, to read in every 
% other frame), is not implimented. If I get requests, I will add it.

% Requirements: 
%    - binary file from trr2matlab.m program
%    - Free RAM: a little more than the size of the binaries being read. 
%         10,000 atoms * 3 axes * 1000 frames * 4 bytes = 160 mb (single precision)
%    - Free Hard Disk: none

% Inputs:
%     - path to binary file (required, must be first input)
%     - start frame (optional)
%         integer, starts reading at this point
%     - end frame (optional)
%         integer, stops reading at this point
%
%     please note! if only one numeric input is given, this program
%     interprets that as the "end frame" and begins reading from the first
%     frame of the simulation

% Outputs:
%     - trajectory matrix
%         this is output as "coodData.trajectory" in the file
%         this is a 3D matrix is made of the trajectory with the format
%           (atom number ; xyz value ; frame number)
%     - information of trajectory 
%         coodData.num_atoms has number of atoms
%         coodData.num_frames has number of frames
%         coodData.time_step has the time incriment between frames

% Example inputs and outputs:
%     [coodData] = readGmx2Matlab('xdata.binary')
%         - makes a 3D matrix (coodData.trajectory) of entire coordinate trajectory
%     [coodData] = readGmx2Matlab('vdata.binary', 1000)
%         - makes a 3D matrix (coodData.trajectory) of velocity trajectory from frames 1 to 1000
%     [coodData] = readGmx2Matlab('fdata.binary', 1000, 2000)
%         - makes a 3D matrix (coodData.trajectory) of force trajectory from frames 1000 to 2000
%
%     [coodData] = readGmx2Matlab('xdata.binary');
%     trajectory = coodData.trajectory(:,:,1:2:end);
%     for n = 1:size(trajectory,3)
%       plot3(trajectory(:,1,n), trajectory(:,2,n), trajectory(:,3,n),'.');
%       drawnow;
%       pause(0.2);
%     end
%        - plot out every other frame of trajectory as a 3D figure

function [coodData] = readGmx2Matlab_d(coordFile, start_frame, end_frame)

if exist(coordFile, 'file')
    fprintf('\n    Reading %s... \n', coordFile);
else
    fprintf('\n    Trajectory file %s not found, exiting program\n', coordFile);
    
    % return null data if trajectory not found, then exit program
    coodData.trajectory = 0;
    coodData.num_atoms = 0;
    coodData.num_frames = 0;
    coodData.time_step = 0;
    
    return;
end

% open file: read-only, big-endian
file_ID = fopen(coordFile, 'r', 'b');

% extract statistics about simulation
coodData.num_atoms = fread(file_ID, 1, 'int32');
coodData.num_frames = fread(file_ID, 1, 'int32');
coodData.time_step = fread(file_ID, 1, 'double');

% scrub input options
if nargin == 3
    if (end_frame > coodData.num_frames)
        end_frame = coodData.num_frames;
    end
    
    if (start_frame > end_frame)
        fprintf(['    starting frame %g is not in trajectory\n' ...
                 ' ...exiting\n'], start_frame);
        return;
    end
    
    fprintf('    beginning from %g frame \n', start_frame);
    fprintf('    processing simulation until frame %g \n', end_frame);
elseif nargin == 2
    end_frame = start_frame;
    start_frame = 1;
    if (end_frame > coodData.num_frames)
        end_frame = coodData.num_frames;
    end
    fprintf('    processing simulation until frame %g \n', end_frame);
else 
    fprintf('    processing all frames from simulation \n');
    start_frame = 1;
    end_frame = coodData.num_frames;
end

% find first frame to read (3 axes * 8 hex sectors per float)
fseek(file_ID, (start_frame - 1)*3*8*coodData.num_atoms, 'cof');
% stream in data
rawdata = fread(file_ID, [3 coodData.num_atoms*end_frame], 'double');
fclose(file_ID);

% structure data as XYZ * atoms * frames 3D matrix
coodData.trajectory = permute(reshape(rawdata,3,coodData.num_atoms,[]),[2 1 3]);

fprintf('Done reading file \n');

end

function buffer_recon_SPEN(connection) % 'Buffer' is a colloquial name for Recon Data.
disp("Matlab SPEN buffer reconstruction running.")

%% Add to path location of reconstruct_SPEN function
%addpath(genpath([pwd filesep 'SPEN_Gadgetron_Reco']));
addpath(genpath('/home/mygadg/Code/Gadg_SPEN_reco_Matlab/gadgetron4.1_Matlab_functions/+SPEN/SPEN_Gadgetron_Reco'));

Parameters = struct;

%   connection.next(); % Discard the first acquisition - it's noise data.
try
    while true
        next_acquisition = @connection.next;
        [Parameters] = gadgetron.SPEN.reconstruct_SPEN(next_acquisition,connection,Parameters);
    end
catch ME
    if ~strcmp(ME.identifier, 'Connection:noNextItem'), rethrow(ME); end
end
end

function buffer_recon_SPEN(connection) % 'Buffer' is a colloquial name for Recon Data.
disp("Matlab SPEN buffer reconstruction running.")

%% Add to path location of reconstruct_SPEN function
addpath('/home/mygadg/Code/Gadg_SPEN_reco_Matlab/gadgetron4.1_functions');
Parameters = struct;

%   connection.next(); % Discard the first acquisition - it's noise data.
try
    while true
        next_acquisition = @connection.next;
        [Parameters] = reconstruct_SPEN(next_acquisition,connection,Parameters);
    end
catch ME
    if ~strcmp(ME.identifier, 'Connection:noNextItem'), rethrow(ME); end
end
end

function buffer_recon_SPEN_diff(connection) % 'Buffer' is a colloquial name for Recon Data.
disp("Matlab SPEN buffer reconstruction running.")

Parameters = struct;
%   connection.next(); % Discard the first acquisition - it's noise data.
try
    while true
        next_acquisition = @connection.next;
        [Parameters] = gadgetron.SPEN.reconstruct_SPEN_diff(next_acquisition,connection,Parameters);
    end
catch ME
    if ~strcmp(ME.identifier, 'Connection:noNextItem'), rethrow(ME); end
end
end

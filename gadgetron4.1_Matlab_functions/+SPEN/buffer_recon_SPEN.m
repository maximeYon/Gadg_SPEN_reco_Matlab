function buffer_recon_SPEN(connection) % 'Buffer' is a colloquial name for Recon Data.
disp("Matlab SPEN buffer reconstruction running.")

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

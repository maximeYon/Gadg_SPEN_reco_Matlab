%% This function perform SPEN reconstruction in case of :
% The configuration is sent to the external process. It's left pretty empty here.

% Listen for incoming connections on port 18000. Call 'handle_connection'
% when a connection is established.
gadgetron.external.listen(18000, @handle_connection);

function handle_connection(connection)
disp("handle_connection was called.")
Parameters = struct;
%   connection.next(); % Discard the first acquisition - it's noise data.
try
    while true
        next_acquisition = @connection.next;
        [Parameters] = gadgetron.SPEN.reconstruct_SPEN(next_acquisition,connection,Parameters);
%         [Parameters] = gadgetron.SPEN.reconstruct_SPEN_diff(next_acquisition,connection,Parameters);
    end
catch ME
    if ~strcmp(ME.identifier, 'Connection:noNextItem'), rethrow(ME); end
end
end


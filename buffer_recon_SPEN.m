function buffer_recon_SPEN(connection) % 'Buffer' is a colloquial name for Recon Data.
disp("Matlab SPEN buffer reconstruction running.")

%% Add to path location of reconstruct_SPEN function
addpath('/home/mygadg/Documents/MATLAB/SPEN_Gadgetron_4');
Parameters = struct;

%% determine the number of repetitions, not ideal solution....
Nsegments = round(connection.header.userParameters.userParameterLong(1,2).value/connection.header.encoding.encodingLimits.kspace_encoding_step_1.maximum); 
Nset = (connection.header.encoding.encodingLimits.set.maximum+1)/Nsegments;

%% loop over the repetitions while keeping the constant parameter (in the Parameters structure)
for ind = 1:Nset
    next_acquisition = @connection.next;
    [Parameters] = reconstruct_SPEN(next_acquisition,connection,Parameters);
end

end

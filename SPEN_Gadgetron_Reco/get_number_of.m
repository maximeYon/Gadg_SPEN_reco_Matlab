function [ number_of_slices, number_of_coils  , number_of_repetitions, number_of_contrasts, number_of_phase, number_of_average , number_of_segments, number_of_sets] = get_number_of( hdr )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here


try
    number_of_slices = hdr.encoding.encodingLimits.slice.maximum + 1;
catch
    number_of_slices = 1;
end

try
    number_of_coils = hdr.acquisitionSystemInformation.receiverChannels;
catch
    number_of_coils = 1;
end

try
    number_of_repetitions = hdr.encoding.encodingLimits.repetition.maximum + 1;
catch
    number_of_repetitions = 1;
end

try
    number_of_contrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 ;
catch
    number_of_contrasts = 1;
end

try
    number_of_phase = hdr.encoding.encodingLimits.phase.maximum + 1 ;
catch
    number_of_phase = 1;
end

try
    number_of_average = hdr.encoding.encodingLimits.average.maximum + 1 ;
catch
    number_of_average = 1;
end
 
try
    number_of_segments = hdr.encoding.encodingLimits.segment.maximum + 1 ;
catch
    number_of_segments = 1;
end

try
    number_of_sets = hdr.encoding.encodingLimits.set.maximum + 1 ;
catch
    number_of_sets = 1;
end

if (number_of_contrasts>1)
    number_of_contrasts
     disp('number_of_contrasts>1');
     return;
end

if (number_of_phase>1)
    number_of_phase
    disp('number_of_phase>1');
    return; 
end

if (number_of_average>1)
    number_of_average
    disp('number_of_average>1');
    return; 
end


end


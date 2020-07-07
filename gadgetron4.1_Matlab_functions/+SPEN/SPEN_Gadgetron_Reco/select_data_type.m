function out = select_data_type(data)
    switch class(data)
        case 'int16'
            out = gadgetron.types.Image.SHORT;
        case 'uint16'
            out = gadgetron.types.Image.USHORT;
        case 'int32'
            out = gadgetron.types.Image.INT;
        case 'uint32'
            out = gadgetron.types.Image.UINT;
        case 'single'
            if isreal(data)
                out = gadgetron.types.Image.FLOAT;
            else
                out = gadgetron.types.Image.CXFLOAT;
            end
        case 'double'
            if isreal(data)
                out = gadgetron.types.Image.DOUBLE;
            else
                out = gadgetron.types.Image.CXDOUBLE;
            end
        otherwise
            error("Unsupported image data type: %s", class(data))
    end
end
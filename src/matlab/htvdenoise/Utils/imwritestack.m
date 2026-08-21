function imwritestack(stack, filename,bit)

im = Tiff(filename, 'w');

infostruct.ImageLength = size(stack, 1);
infostruct.ImageWidth = size(stack, 2);
infostruct.Photometric = Tiff.Photometric.MinIsBlack;
if bit==8
    infostruct.BitsPerSample =8;
elseif bit ==16
    infostruct.BitsPerSample =16;
end
infostruct.SampleFormat = Tiff.SampleFormat.UInt;
infostruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

if bit==8
    for k = 1:size(stack, 3)
    im.setTag(infostruct)
    im.write(uint8(stack(:, :, k)));
    im.writeDirectory();
    end
elseif bit ==16
    for k = 1:size(stack, 3)
    im.setTag(infostruct)
    im.write(uint16(stack(:, :, k)));
    im.writeDirectory();
    end
end


im.close();
end
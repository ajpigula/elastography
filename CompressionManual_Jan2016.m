function track = CompressionManual_Jan2016(patient,msl,visit)

warning off

[swp, frc] = getSweep(patient,msl,visit);

peak = find(frc(:,2) == max(frc(:,2)));
swp = swp(:,:,1:peak);

height = size(swp,1);
width = size(swp,2);
snaps = size(swp,3);

startFrame = 1;
track = [];

for frame = startFrame:3:snaps
    
    figure(1);
    [x,y] = getpts(1);
    
    track(frame,:) = [x y];
    
end

end

function sweep = plotTrack(sweep, track)

for j = 1:size(sweep,3)
    sweep(round(track(j,1)),:,j) = 0;
    sweep([round(track(j,1))-4:round(track(j,1))+4],round(track(j,2)),j) = 0;
end

end
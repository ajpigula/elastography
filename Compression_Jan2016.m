
function Compression_Jan2016(patient,msl,visit)

warning off

[swp, frc] = getSweep(patient,msl,visit);

peak = find(frc(:,2) == max(frc(:,2)));
swp = swp(:,:,1:peak);

height = size(swp,1);
width = size(swp,2);
snaps = size(swp,3);

startFrame = 1;

xx = [0 100];
while range(xx) >= 5
    [xx,yy] = getStarterPoints(swp,startFrame);
end

s = [yy(1) xx(1)];
b = [yy(2) xx(2)];

subQ = trackMe(swp,s,startFrame);
boneTop = trackMe(swp,b,startFrame);

if length(xx) == 3
    i = [yy(3) xx(3)];
    interMus = trackMe(swp,i,startFrame);
    swp = plotTrack(swp,interMus);
    save(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016/',patient,'_',msl,'_',visit,'.mat'),'subQ','boneTop','interMus')

else
    save(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016/',patient,'_',msl,'_',visit,'.mat'),'subQ','boneTop')   
end
    
swp = plotTrack(swp,boneTop);
swp = plotTrack(swp,subQ);

implay(swp)

end


function [xx,yy] = getStarterPoints(swp,startFrame)

figure(1); imshow(swp(:,:,startFrame));
[xx,yy] = getpts(1);
close(1);
xx = round(xx);
yy = round(yy);

end


function [sweep,sweepForces] = getSweep(patient,msl,visit)

sweep = [];
sweepForces = [];

load('/home/sisir/code/data/DMD_Complete_20Aug2015/directory_18Oct2015.mat');
direct = s;
clear s;

bas = '/home/sisir/code/data/DMD_Complete_20Aug2015/Sweeps_18Oct2015/';

x = findstr([direct.patient], patient);
x = (x+4)/5;

for i = 1:length(x)
    if abs(x(i)-round(x(i))) < .01
        if isequal(direct(x(i)).visit,visit)
            
            load(strcat(bas,patient,'_',visit,'/',direct(x(i)).filer,'_',msl,'_sweep.mat'));
            load(strcat(bas,patient,'_',visit,'/',direct(x(i)).filer,'_',msl,'_forcedata.mat'));
            
        end 
    end   
end

if isempty(sweep)
    error('Did not load files properly');
end

end


function track = trackMe(sweep,seed,startFrame)

warning off

height = size(sweep,1);
width = size(sweep,2);
snaps = size(sweep,3);

track = seed;


% for frame = startFrame+1:5
for frame = startFrame+1:snaps

    lowest = zeros(10,5,2);
    
    for dx = 1:5
        for dy = 1:10
            
            template = sweep(seed(1)-dy:seed(1)+dy,seed(2)-dx:seed(2)+dx,frame-1);
            
            if frame > startFrame+1
                lastSeed = track(frame-2,:);
                lastTemplate = sweep(lastSeed(1)-dy:lastSeed(1)+dy,lastSeed(2)-dx:lastSeed(2)+dx,frame-2);
                template = (template+lastTemplate)/2;
            end
            
            scores = ones(21,5); % dy dx i j
            
            for i = -10:10
                for j = -2:2
                    
%                     fprintf('\n%d %d %d %d %d',i,j,dy,dx,frame)
                    seedOffset = seed + [i j];
                    roi = sweep(seedOffset(1)-dy:seedOffset(1)+dy,seedOffset(2)-dx:seedOffset(2)+dx,frame);
                    scores(i+11,j+3) = mean2( (template-roi).^2 );
                                        
                end
            end
            
            idx = find(scores == min(min(scores)));
            [a,b] = ind2sub(size(scores),idx);
            a = max(a);
            b = max(b);
            
            lowest(dy,dx,:) = [a-11 b-3];
            
        end
    end
    
    pickMe = median(reshape(lowest,1,10*5,2),2);
    seed = track(frame-1,:)+transpose(squeeze(round(pickMe)));
    
    track(frame,:) = seed;
    
end

end


function sweep = plotTrack(sweep, track)

for j = 1:size(sweep,3)
    sweep(round(track(j,1)),:,j) = 0;
    sweep([round(track(j,1))-4:round(track(j,1))+4],round(track(j,2)),j) = 0;
end

end
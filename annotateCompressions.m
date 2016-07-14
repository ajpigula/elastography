
function [boneTop, interMus, subQ, sweep, sweepForces] = annotateCompressions(patient,visit,muscle)

clear paths
clear j

basePath = '/home/sisir/code/data/DMD_Complete_20Aug2015/Sweeps_22Mar2016/';
path = strcat(basePath,patient,'_',visit,'/');
d = dir(path);

paths = {};

if findstr(muscle, 'XDoNotUse')
    for j = 3:length(d)
        if findstr(d(j).name,strcat(muscle,'_forcedata.mat'))
            paths{end+1} = d(j).name;
            paths{end+1} = d(j+1).name;
            
        end
    end
else 
    for j = 3:length(d)
        if isempty(findstr(d(j).name, 'XDoNotUse'))
            if findstr(d(j).name,muscle)
                paths{end+1} = d(j).name;
            end
        end
    end
end

load(fullfile(path,paths{1}))
load(fullfile(path,paths{2}))
    

peak = find(sweepForces(:,2) == max(sweepForces(:,2)),1);

[xx,yy] = getStarterPoints(sweep(:,:,1));
s = [yy(1) xx(1)];
b = [yy(2) xx(2)];
i = [yy(3) xx(3)];

for z = 1:size(sweep,3)
    sm_sweep(:,:,z) = imfilter(sweep(:,:,z),ones(3,3)/9);
end

subQ = trackMe(sweep(:,:,1:peak),s);
boneTop = trackMe(sweep(:,:,1:peak),b);
interMus = trackMe(sweep(:,:,1:peak),i);

% save(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016_quads/',patient,'_',muscle,'_',visit,'.mat'),'subQ','boneTop','interMus','sweepForces');

swp = sweep;
for k = 1:length(subQ)
    swp([subQ(k,1) boneTop(k,1) interMus(k,1)],:,k) = 0;
end

implay(swp(:,:,1:peak))

end


function [xx,yy] = getStarterPoints(swp,startFrame)

figure(1); imshow(swp);
[xx,yy] = getpts(1);
close(1);
xx = round(xx);
yy = round(yy);

end


function track = trackMe(sweep,seed)

warning off

height = size(sweep,1);
width = size(sweep,2);
snaps = size(sweep,3);

track = seed;


for frame = 2:snaps

    lowest = zeros(10,5,2);
    
    for dx = 1:5
        for dy = 1:10
            
            template = sweep(seed(1)-dy:seed(1)+dy,seed(2)-dx:seed(2)+dx,frame-1);
            
            if frame > 2
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
                    scores(i+11,j+3) = mse(template,roi);
                                        
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


function mutualInformation = mutualInfo(im1,im2)

im1 = round(255*im1);
im2 = round(255*im2);

indrow = double(im1(:)) + 1;
indcol = double(im2(:)) + 1;

jointHistogram = accumarray([indrow indcol], 1);
jointProb = jointHistogram / numel(indrow);
indNoZero = jointHistogram ~= 0;
jointProb1DNoZero = jointProb(indNoZero);
jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));

entropy1 = entropy(im1/255);
entropy2 = entropy(im2/255);

mutualInformation = entropy1 + entropy2 - jointEntropy;

end


function MSE = mse(im1,im2)

MSE = mean2((im1-im2).^2);

end














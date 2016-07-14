% function CompressionAnalysis
function Q = CompressionAnalysis_expanded

clc
fprintf('\n\n');

global quadsData
quadsData = struct;
quadsData(1).patient = '00000';


functionals = getFunc;

d = dir('/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016_quads');
d = d(3:end);

%Calculate and store data
for i = 1:length(d)
    
    [topMus,botMus,wholMus,frc] = getMuscleThick(d(i).name);
        
    E(:,:,1) = getYoungMod(topMus,frc);
    E(:,:,2) = getYoungMod(botMus,frc);
    E(:,:,3) = getYoungMod(wholMus,frc); 
    
    createEntry(d(i).name,functionals,E);
    
end

quadsData = quadsData(2:end);
Q = quadsData;

YoungsModulus;



end


function YoungsModulus(r1,r2)

global quadsData

% fprintf('\nDoes Youngs modulus differentiate DMD/normal? (t-test)\n')
% checkStats(.001);

% plotBox
% vsAge
% plotBoxAndAge
% vs6MW
% vsThingy


% longitudinal('01005','rx')
% longitudinal('01009','rx')
% longitudinal('01011','rx')
% longitudinal('01014','rx')
% longitudinal('01016','rx')
% longitudinal('01023','rx')
longitudinal('01026','d',[244/255 131/255 55/255])

% 
% longitudinal('02003','bo')
longitudinal('02012','o',[173/255 216/255 235/255])
% longitudinal('02017','bo')
% longitudinal('02019','bo')
% longitudinal('02021','bo')
% longitudinal('02022','bo')
% longitudinal('02026','bo')

% vsAmb

% legend('DMD','DMD fit','Healthy','Healthy fit')

end


function data = getFunc

filename = '/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016_info/DataSummary_Compression_Functionals.csv';
delimiter = ';';
formatSpec = '%d%s%s%s%s%f%s%f%s%s%s%s%s%d%s%f%s%s%s%s%s%s';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'headerLines', 1, 'ReturnOnError', false);
fclose(fileID);

data = struct;
for i = 1:length(dataArray{1})
    d = dataArray{1}(i);
    data(i).patient = d;
    data(i).visit = dataArray{4}{i};
    data(i).age = dataArray{6}(i);
    data(i).ambulation = dataArray{14}(i);
    data(i).MW6 = dataArray{16}(i);
end


end


function [topMus,botMus,wholMus,frc] = getMuscleThick(filename)

load(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016_quads/',filename))

sweepForces = sweepForces(1:length(subQ),2);

idx1_5 = find(abs(sweepForces-1.5) == min(abs(sweepForces-1.5)));
idx10 = find(abs(sweepForces-10) == min(abs(sweepForces-10)));

idx1_5 = idx1_5(1);
idx10 = idx10(1);

frc = sweepForces([idx1_5 idx10]);

topMus = interMus([idx1_5 idx10],1) - subQ([idx1_5 idx10],1);
botMus = boneTop([idx1_5 idx10],1) - interMus([idx1_5 idx10],1);
wholMus = boneTop([idx1_5 idx10],1) - subQ([idx1_5 idx10],1);

end


function E = getYoungMod(x,f)

A = .013*.047; % Could make bone size

E = (-x(1) * diff(f) ) / (A* diff(x));
E = E/1000;


if isinf(E)
    E = NaN;
end

if E < 0
    E = NaN;
end

end


function createEntry(filename,funcs,E)

global quadsData

ID1 = strsplit(filename,'.');
ID1 = ID1{1};
ID2 = strsplit(ID1,'_');

walk_6min = NaN;
for j = 1:length(funcs)
    if isequal(funcs(j).patient, str2num(ID2{1}))
        if isequal(funcs(j).visit,ID2{end})
            age = funcs(j).age;
            amb = funcs(j).ambulation;
            walk_6min = funcs(j).MW6;
        end
    end    
end


if length(ID2) == 4
    mus = strcat(ID2{2},'_',ID2{3});
else
    mus = ID2{2};
end


quadsData(end+1).patient = ID2{1};
quadsData(end).visit = ID2{end};
quadsData(end).muscle = mus;
quadsData(end).age = age;
quadsData(end).MW6 = walk_6min;
quadsData(end).ambulatory = amb;
quadsData(end).E = squeeze(E);


end


function checkStats(alpha)

global quadsData

dmd = [];
norm = [];
for i = 1:length(quadsData)
    if strcmp(quadsData(i).patient(2), '1')
        dmd(i,:) = quadsData(i).E;
        
    elseif strcmp(quadsData(i).patient(2), '2')
        norm(i,:) = quadsData(i).E;
    end
end

[t_lin_q_top, p_lin_q_top] = ttest2(dmd(:,1),norm(:,1),'Alpha',alpha,'Vartype','unequal');
fprintf('--------------------------------------------');
fprintf('\nQuads, upper muscle, alpha = %f\n',alpha);
t_lin_q_top
p_lin_q_top

[t_lin_q_bottom, p_lin_q_bottom] = ttest2(dmd(:,2),norm(:,2),'Alpha',alpha,'Vartype','unequal');
fprintf('--------------------------------------------');
fprintf('\nQuads, lower muscle, alpha = %f\n',alpha);
t_lin_q_bottom
p_lin_q_bottom

[t_lin_q_whole, p_lin_q_whole] = ttest2(dmd(:,3),norm(:,3),'Alpha',alpha,'Vartype','unequal');
fprintf('--------------------------------------------');
fprintf('\nQuads, whole muscle, alpha = %f\n',alpha);
t_lin_q_bottom
p_lin_q_bottom


% fprintf('--------------------------------------------');
% fprintf('\nQuads, whole muscle, DMD    : mean = %f, stdev = %f',m_qd,s_qd);
% fprintf('\nQuads, whole muscle, normal : mean = %f, stdev = %f',m_qn,s_qn);
% fprintf('\nFor alpha = %f\n',alpha);
% fprintf('t = %d, p = %f\n',t_lin_q_whole,p_lin_q_whole)


end


function plotBox

global quadsData

dmd = [];
norm = [];
for i = 1:length(quadsData)
    if strcmp(quadsData(i).patient(2), '1')
        dmd(end+1,:) = quadsData(i).E;
        
    elseif strcmp(quadsData(i).patient(2), '2')
        norm(end+1,:) = quadsData(i).E;
    end
end

types = {};
for i = 1:length(norm)    
    types{i} = 'Healthy';
end
for i = 1+length(norm):length(dmd)+length(norm)
    types{i} = 'DMD';
end

squeeze(norm(:,1))
squeeze(dmd(:,1))

figure
boxplot([squeeze(norm(:,1));squeeze(dmd(:,1))],types,'Whisker',2);
ylabel('Youngs modulus (kPa)')
title('Upper muscle')

% figure; hold on
% plot(squeeze(norm(:,1)),'rx')
% plot(squeeze(dmd(:,1)),'bo')
% title('Upper muscle')
% ylabel('Youngs modulus (kPa)')

figure
boxplot([squeeze(norm(:,2));squeeze(dmd(:,2))],types);
ylabel('Avg. Youngs modulus (kPa)')
title('Lower muscle')

% 
% figure; hold on
% plot(squeeze(norm(:,2)),'rx')
% plot(squeeze(dmd(:,2)),'bo')
% title('Lower muscle')

figure
boxplot([squeeze(norm(:,3));squeeze(dmd(:,3))],types);
ylabel('Youngs modulus (kPa)')
title('Whole muscle')

% figure; hold on
% plot(squeeze(norm(:,3)),'rx')
% plot(squeeze(dmd(:,3)),'bo')
% title('Whole muscle')


end


function plotBoxAndAge

global quadsData

dmd = [];
norm = [];
age_dmd = [];
age_norm = [];
for i = 1:length(quadsData)
    if strcmp(quadsData(i).patient(2), '1')
        dmd(end+1,:) = quadsData(i).E;
        age_dmd(end+1) = quadsData(i).age;
        
    elseif strcmp(quadsData(i).patient(2), '2')
        norm(end+1,:) = quadsData(i).E;
        age_norm(end+1) = quadsData(i).age;
    end
end

types = {};
for i = 1:length(norm)    
    types{i} = 'Healthy';
end
for i = 1+length(norm):length(dmd)+length(norm)
    types{i} = 'DMD';
end



figure;
subplot(1,3,1);
boxplot([squeeze(norm(:,1));squeeze(dmd(:,1))],types,'Whisker',5);
ylim([0 500])
ylabel('Youngs modulus (kPa)')
set(gca,'fontsize', 16,'XTickLabelRotation',50);
set(findobj('Tag','Box'),'Color','k')
set(findobj('Tag','Median'),'Color',[244/255 131/255 55/255],'LineWidth',1.5)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
box = get(gca,'position');

subplot(1,3,2:3); hold on;
plot(age_norm,squeeze(norm(:,1)),'o','MarkerFaceColor',[173/255 216/255 235/255],'MarkerEdgeColor','k')
plot(age_dmd,squeeze(dmd(:,1)),'d','MarkerFaceColor',[244/255 131/255 55/255],'MarkerEdgeColor','k')
ylim([0 500])
pl = get(gca,'position');
pl([2 4]) = box([2 4]);
set(gca,'fontsize',16, 'YTickLabel','', 'position', pl);
xlabel('Age (years)')
legend('Healthy','DMD')
set(suptitle('Rectus femoris'),'FontSize',20)




figure;
subplot(1,3,1);
boxplot([squeeze(norm(:,2));squeeze(dmd(:,2))],types,'Whisker',5);
ylim([0 200])
ylabel('Youngs modulus (kPa)')
set(gca,'fontsize', 16,'XTickLabelRotation',50);
set(findobj('Tag','Box'),'Color','k')
set(findobj('Tag','Median'),'Color',[244/255 131/255 55/255],'LineWidth',1.5)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
box = get(gca,'position');

subplot(1,3,2:3); hold on;
plot(age_norm,squeeze(norm(:,2)),'o','MarkerFaceColor',[173/255 216/255 235/255],'MarkerEdgeColor','k')
plot(age_dmd,squeeze(dmd(:,2)),'d','MarkerFaceColor',[244/255 131/255 55/255],'MarkerEdgeColor','k')
ylim([0 200])
pl = get(gca,'position');
pl([2 4]) = box([2 4]);
set(gca,'fontsize',16, 'YTickLabel','', 'position', pl);
xlabel('Age (years)')
legend('DMD','Healthy')
set(suptitle('Vastus intermedius'),'FontSize',20)




figure;
subplot(1,3,1);
boxplot([squeeze(norm(:,3));squeeze(dmd(:,3))],types,'Whisker',5);
ylim([0 300])
ylabel('Youngs modulus (kPa)')
set(gca,'fontsize', 16,'XTickLabelRotation',50);
set(findobj('Tag','Box'),'Color','k')
set(findobj('Tag','Median'),'Color',[244/255 131/255 55/255],'LineWidth',1.5)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
box = get(gca,'position');

subplot(1,3,2:3); hold on;
plot(age_norm,squeeze(norm(:,3)),'o','MarkerFaceColor',[173/255 216/255 235/255],'MarkerEdgeColor','k')
plot(age_dmd,squeeze(dmd(:,3)),'d','MarkerFaceColor',[244/255 131/255 55/255],'MarkerEdgeColor','k')
ylim([0 300])
pl = get(gca,'position');
pl([2 4]) = box([2 4]);
set(gca,'fontsize',16, 'YTickLabel','', 'position', pl);
xlabel('Age (years)')
legend('DMD','Healthy')
set(suptitle('Rectus femoris + vastus intermedius'),'FontSize',20)




end


function vsAge

global quadsData

for i = 1:length(quadsData)
    if strcmp(quadsData(i).patient(2), '1')
        ages_d(i) = quadsData(i).age;
        es_d(i,:) = quadsData(i).E;
    elseif strcmp(quadsData(i).patient(2), '2')
        ages_n(i) = quadsData(i).age;
        es_n(i,:) = quadsData(i).E;
    end
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Youngs modulus correlations with age in DMD\n'); 
          
C = corrcoef(es_d(:,1),ages_d);
fprintf('Upper muscle: %.2f\n',C(1,2));          
C = corrcoef(es_d(:,2),ages_d);
fprintf('Lower muscle: %.2f\n',C(1,2));           
C = corrcoef(es_d(:,3),ages_d);
fprintf('Whole muscle: %.2f\n',C(1,2));   


figure; hold on
plot(ages_d,es_d(:,1),'rx','LineWidth',1.5)
plot(ages_n,es_n(:,1),'bo','LineWidth',1.5)
legend('DMD','Normal')
title('Rectus femoris')
ylabel('Youngs modulus (kPa)')
xlabel('Age (years)')

figure; hold on
plot(ages_d,es_d(:,2),'rx','LineWidth',1.5)
plot(ages_n,es_n(:,2),'bo','LineWidth',1.5)
legend('DMD','Normal')
title('Vastus intermedius')
ylabel('Youngs modulus (kPa)')
xlabel('Age (years)')

figure; hold on
plot(ages_d,es_d(:,3),'rx','LineWidth',1.5)
plot(ages_n,es_n(:,3),'bo','LineWidth',1.5)
legend('DMD','Normal')
title('Rectus femoris + vastus intermedius')
ylabel('Youngs modulus (kPa)')
xlabel('Age (years)')


end


function vsThingy

global quadsData

for i = 1:length(quadsData)
    if strcmp(quadsData(i).patient(2), '1')
        thingy_d(i) = quadsData(i).age/quadsData(i).MW6;
        es_d(i,:) = quadsData(i).E;
    elseif strcmp(quadsData(i).patient(2), '2')
        thingy_n(i) = quadsData(i).age/quadsData(i).MW6;
        es_n(i,:) = quadsData(i).E;
    end
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Youngs modulus correlations with thingy in DMD\n'); 
          
C = corrcoef(es_d(:,1),thingy_d);
fprintf('Upper muscle: %.2f\n',C(1,2));          
C = corrcoef(es_d(:,2),thingy_d);
fprintf('Lower muscle: %.2f\n',C(1,2));           
C = corrcoef(es_d(:,3),thingy_d);
fprintf('Whole muscle: %.2f\n',C(1,2));   


figure; hold on
plot(thingy_d,es_d(:,1),'rx','LineWidth',1.5)
plot(thingy_n,es_n(:,1),'bo','LineWidth',1.5)
legend('DMD','Normal')
title('Upper muscle')
ylabel('Youngs modulus (kPa)')
xlabel('Thingy (m/yr)')

figure; hold on
plot(thingy_d,es_d(:,2),'rx','LineWidth',1.5)
plot(thingy_n,es_n(:,2),'bo','LineWidth',1.5)
legend('DMD','Normal')
title('Lower muscle')
ylabel('Youngs modulus (kPa)')
xlabel('Thingy (m/yr)')

figure; hold on
plot(thingy_d,es_d(:,3),'rx','LineWidth',1.5)
plot(thingy_n,es_n(:,3),'bo','LineWidth',1.5)
legend('DMD','Normal')
title('Whole muscle')
ylabel('Youngs modulus (kPa)')
xlabel('Thingy (m/yr)')


end


function vs6MW

global quadsData

for i = 1:length(quadsData)
    if strcmp(quadsData(i).patient(2), '1')
        mw6_d(i) = quadsData(i).MW6;
        es_d(i,:) = quadsData(i).E;
    elseif strcmp(quadsData(i).patient(2), '2')
        mw6_n(i) = quadsData(i).MW6;
        es_n(i,:) = quadsData(i).E;
    end
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Youngs modulus correlations with 6MW in DMD\n'); 
          
C = corrcoef(es_d(:,1),mw6_d);
fprintf('Upper muscle: %.2f\n',C(1,2));          
C = corrcoef(es_d(:,2),mw6_d);
fprintf('Lower muscle: %.2f\n',C(1,2));           
C = corrcoef(es_d(:,3),mw6_d);
fprintf('Whole muscle: %.2f\n',C(1,2));   


figure; hold on
plot(mw6_d,es_d(:,1),'rx','LineWidth',1.5)
title('Rectus femoris')
ylabel('Youngs modulus (kPa)')
xlabel('6MWT (m)')

figure; hold on
plot(mw6_d,es_d(:,2),'rx','LineWidth',1.5)
title('Vastus intermedius')
ylabel('Youngs modulus (kPa)')
xlabel('6MWT (m)')

figure; hold on
plot(mw6_d,es_d(:,3),'rx','LineWidth',1.5)
title('Rectus femoris + vastus intermedius')
ylabel('Youngs modulus (kPa)')
xlabel('6MWT (m)')


end


function vsAmb

global quadsData

zz = [];
oo = [];
tt = [];
for i = 1:length(quadsData)
    if strcmp(quadsData(i).patient(2), '1')
        
        if quadsData(i).ambulatory == 0
            zz(end+1,:) = quadsData(i).E;

        elseif quadsData(i).ambulatory == 1
            oo(end+1,:) = quadsData(i).E;

        elseif quadsData(i).ambulatory == 2
            tt(end+1,:) = quadsData(i).E;

        end   
    end
end

types = {};
for j = 1:size(tt,1)
    types{j} = 'Fully ambulatory';
end
for j = 1:size(oo,1)
    types{end+1} = 'Partially Ambulatory';
end
for j = 1:size(zz,1)
    types{end+1} = 'Nonambulatory';
end

figure
boxplot([tt(:,1)' oo(:,1)' zz(:,1)'],types)
title('Upper muscle')
ylabel('Youngs modulus (kPa)')

figure
boxplot([tt(:,2)' oo(:,2)' zz(:,2)'],types)
title('Lower muscle')
ylabel('Youngs modulus (kPa)')

figure
boxplot([tt(:,3)' oo(:,3)' zz(:,3)'],types)
title('Whole muscle')
ylabel('Youngs modulus (kPa)')

end


function longitudinal(patient,plotStyle,col)


global quadsData

age = [];
mw6 = [];
es = [];

for i = 1:length(quadsData)
    if strcmp(quadsData(i).patient, patient)
        age(end+1) = quadsData(i).age;
        mw6(end+1) = quadsData(i).MW6;
        es(end+1,:) = quadsData(i).E;
    end
end

agerange = [min(age)-.5 max(age)+.5];

figure(1); hold on
plot(age,es(:,1),plotStyle,'MarkerFaceColor',col,'MarkerEdgeColor','k')
xlabel('Age (years)')
ylabel('Youngs modulus (kPa)')

p = polyfit(age',es(:,1),1);
plot(agerange,p(1)*agerange+p(2),'--','Color',col)
set(gca,'fontsize', 16)
title('Rectus femoris','FontSize',20)


figure(2); hold on
plot(age,es(:,2),plotStyle,'MarkerFaceColor',col,'MarkerEdgeColor','k')
xlabel('Age (years)')
ylabel('Youngs modulus (kPa)')

p = polyfit(age',es(:,2),1);
plot(agerange,p(1)*agerange+p(2),'--','Color',col)
set(gca,'fontsize', 16)
title('Vastus intermedius','FontSize',20)


figure(3); hold on
plot(age,es(:,3),plotStyle,'MarkerFaceColor',col,'MarkerEdgeColor','k')
xlabel('Age (years)')
ylabel('Youngs modulus (kPa)')

p = polyfit(age',es(:,3),1);
plot(agerange,p(1)*agerange+p(2),'--','Color',col)
set(gca,'fontsize', 16)
title('Rectus femoris + vastus intermedius','FontSize',20)

% figure; hold on
% plot(mw6,es(:,1),plotStyle,'LineWidth',2)
% title('Upper muscle vs 6MW')
% xlabel('6MWT (m)')
% ylabel('Youngs modulus (kPa)')
% 
% figure; hold on
% plot(mw6,es(:,2),plotStyle,'LineWidth',2)
% title('Lower muscle vs 6MW')
% xlabel('6MWT (m)')
% ylabel('Youngs modulus (kPa)')
% 
% figure; hold on
% plot(mw6,es(:,3),plotStyle,'LineWidth',2)
% title('Whole muscle vs 6MW')
% xlabel('6MWT (m)')
% ylabel('Youngs modulus (kPa)')


end


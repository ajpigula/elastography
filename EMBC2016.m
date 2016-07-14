

function EMBC2016

clear
warning off

load('/home/sisir/code/dmd/directories/EMBC2016.mat')

basePath = '/home/sisir/code/data/DMD_Complete_20Aug2015/';
runID1 = 'Sweeps_22Mar2016/';
runID_q = 'Compression_Jan2016_quads/';
runID_b = 'Compression_Jan2016_biceps/';

labels = {};
spring_dmd_q = zeros(23,3);
spring_norm_q = zeros(20,3);

for i = 1:length(d)
    stable_i = i;
    
    %% Quads
    
    folder = strcat(d(i).name,'_',d(i).q_visit);
    tree = dir(fullfile(basePath,runID1,folder));
    
    idx = [];
    for j = 3:length(tree)
        p = findstr(tree(j).name,'_');
        thisMuscle = tree(j).name(p(1)+1:p(end)-1);
        if length(d(i).q_muscle) > length(thisMuscle)
            continue
        end
        if strcmp(thisMuscle(1:length(d(i).q_muscle)),d(i).q_muscle)
            load(fullfile(basePath,runID1,folder,tree(j).name))
        end
    end
    
    pth_q = fullfile( basePath,runID_q, strcat(d(i).name,'_',d(i).q_muscle,'_',d(i).q_visit) );
    load(pth_q);
    
    rect = interMus(:,1)-subQ(:,1);
    vast = boneTop(:,1)-interMus(:,1);
    both = boneTop(:,1)-subQ(:,1);
        
    [k_rect,k_vast,k_both] = linearSpringConstant(rect,vast,both,sweepForces(:,2));
    [E_rect,E_vast,E_both] = YoungsModulus(rect,vast,both,sweepForces(:,2));
    [Cn_rect,Cn_vast,Cn_both] = Nonlinear(rect,vast,both,sweepForces(:,2));
    
    if stable_i < 24
        
        spring_dmd_q(stable_i,:) = [k_rect k_vast k_both];
        youngs_dmd_q(stable_i,:) = [E_rect E_vast E_both];
        expo_dmd_q(stable_i,:) = [Cn_rect(1) Cn_vast(1) Cn_both(1)];
        strength_dmd_q(stable_i,:) = [exp(Cn_rect(2)) exp(Cn_vast(2)) exp(Cn_both(2))];
        
        labels{20+stable_i}= 'DMD';
        
    elseif stable_i > 23
                
        spring_norm_q(stable_i-23,:) = [k_rect k_vast k_both];
        youngs_norm_q(stable_i-23,:) = [E_rect E_vast E_both];
        expo_norm_q(stable_i-23,:) = [Cn_rect(1) Cn_vast(1) Cn_both(1)];
        strength_norm_q(stable_i-23,:) = [exp(Cn_rect(2)) exp(Cn_vast(2)) exp(Cn_both(2))];
        
        labels{stable_i-23}= 'Normal';
        
    end
    
    
    %% Biceps
  
    folder = strcat(d(i).name,'_',d(i).b_visit);
    tree = dir(fullfile(basePath,runID1,folder));
    
    idx = [];
    for j = 3:length(tree)
        p = findstr(tree(j).name,'_');
        thisMuscle = tree(j).name(p(1)+1:p(end)-1);
        if length(d(i).b_muscle) > length(thisMuscle)
            continue
        end
        if strcmp(thisMuscle(1:length(d(i).b_muscle)),d(i).b_muscle)
            load(fullfile(basePath,runID1,folder,tree(j).name))
        end
    end
    
    pth_b = fullfile( basePath,runID_b, strcat(d(i).name,'_',d(i).b_muscle,'_',d(i).b_visit) );
    load(pth_b);
    
    bicep = interMus(:,1)-subQ(:,1);
    brach = boneTop(:,1)-interMus(:,1);
    both = boneTop(:,1)-subQ(:,1);
        
    [k_bicep,k_brach,k_both] = linearSpringConstant(bicep,brach,both,sweepForces(:,2));
    [E_bicep,E_brach,E_both] = YoungsModulus(bicep,brach,both,sweepForces(:,2));
    [Cn_bicep,Cn_brach,Cn_both] = Nonlinear(bicep,brach,both,sweepForces(:,2));
    
    if stable_i < 24
        
        spring_dmd_b(stable_i,:) = [k_rect k_brach k_both];
        youngs_dmd_b(stable_i,:) = [E_rect E_brach E_both];
        expo_dmd_b(stable_i,:) = [Cn_rect(1) Cn_brach(1) Cn_both(1)];
        strength_dmd_b(stable_i,:) = [exp(Cn_rect(2)) exp(Cn_brach(2)) exp(Cn_both(2))];
        
    elseif stable_i > 23
                
        spring_norm_b(stable_i-23,:) = [k_bicep k_brach k_both];
        youngs_norm_b(stable_i-23,:) = [E_bicep E_brach E_both];
        expo_norm_b(stable_i-23,:) = [Cn_bicep(1) Cn_brach(1) Cn_both(1)];
        strength_norm_b(stable_i-23,:) = [exp(Cn_bicep(2)) exp(Cn_brach(2)) exp(Cn_both(2))];
        
    end
    
end

% Print statistics

alpha = .01;

fprintf('\n')
[h,p] = ttest2(spring_norm_b(:,3),spring_dmd_b(:,3),'Alpha',alpha,'Vartype','unequal');
fprintf('\nSpring constant')
fprintf('\nBiceps, healthy: Mean = %.4f, Stdev = %.4f',mean(spring_norm_b(:,3)),std(spring_norm_b(:,3)))
fprintf('\nBiceps, DMD: Mean = %.4f, Stdev = %.4f',mean(spring_dmd_b(:,3)),std(spring_dmd_b(:,3)))
fprintf('\np = %.4f',p)

fprintf('\n')
[h,p] = ttest2(spring_norm_q(:,3),spring_dmd_q(:,3),'Alpha',alpha,'Vartype','unequal');
fprintf('\nSpring constant')
fprintf('\nQuads, healthy: Mean = %.4f, Stdev = %.4f',mean(spring_norm_q(:,3)),std(spring_norm_q(:,3)))
fprintf('\nQuads, DMD: Mean = %.4f, Stdev = %.4f',mean(spring_dmd_q(:,3)),std(spring_dmd_q(:,3)))
fprintf('\np = %.4f',p)


fprintf('\n')
[h,p] = ttest2(strength_norm_b(:,3),strength_dmd_b(:,3),'Alpha',alpha,'Vartype','unequal');
fprintf('\nStrength coefficient')
fprintf('\nBiceps, healthy: Mean = %.4f, Stdev = %.4f',mean(strength_norm_b(:,3)),std(strength_norm_b(:,3)))
fprintf('\nBiceps, DMD: Mean = %.4f, Stdev = %.4f',mean(strength_dmd_b(:,3)),std(strength_dmd_b(:,3)))
fprintf('\np = %.4f',p)

fprintf('\n')
[h,p] = ttest2(strength_norm_q(:,3),strength_dmd_q(:,3),'Alpha',alpha,'Vartype','unequal');
fprintf('\nStrength coefficient')
fprintf('\nQuads, healthy: Mean = %.4f, Stdev = %.4f',mean(strength_norm_q(:,3)),std(strength_norm_q(:,3)))
fprintf('\nQuads, DMD: Mean = %.4f, Stdev = %.4f',mean(strength_dmd_q(:,3)),std(strength_dmd_q(:,3)))
fprintf('\np = %.4f',p)


fprintf('\n')
[h,p] = ttest2(expo_norm_b(:,3),expo_dmd_b(:,3),'Alpha',alpha,'Vartype','unequal');
fprintf('\nNonlinear strain exponent')
fprintf('\nBiceps, healthy: Mean = %.4f, Stdev = %.4f',mean(expo_norm_b(:,3)),std(expo_norm_b(:,3)))
fprintf('\nBiceps, DMD: Mean = %.4f, Stdev = %.4f',mean(expo_dmd_b(:,3)),std(expo_dmd_b(:,3)))
fprintf('\np = %.4f',p)

fprintf('\n')
[h,p] = ttest2(expo_norm_q(:,3),expo_dmd_q(:,3),'Alpha',alpha,'Vartype','unequal');
fprintf('\nNonlinear strain exponent')
fprintf('\nQuads, healthy: Mean = %.4f, Stdev = %.4f',mean(expo_norm_q(:,3)),std(expo_norm_q(:,3)))
fprintf('\nQuads, DMD: Mean = %.4f, Stdev = %.4f',mean(expo_dmd_q(:,3)),std(expo_dmd_q(:,3)))
fprintf('\np = %.4f',p)


fprintf('\n')
[h,p] = ttest2(youngs_norm_b(:,3),youngs_dmd_b(:,3),'Alpha',alpha,'Vartype','unequal');
fprintf('\nYoungs modulus')
fprintf('\nBiceps, healthy: Mean = %.4f, Stdev = %.4f',mean(youngs_norm_b(:,3)),std(youngs_norm_b(:,3)))
fprintf('\nBiceps, DMD: Mean = %.4f, Stdev = %.4f',mean(youngs_dmd_b(:,3)),std(youngs_dmd_b(:,3)))
fprintf('\np = %.4f',p)

fprintf('\n')
[h,p] = ttest2(youngs_norm_q(:,3),youngs_dmd_q(:,3),'Alpha',alpha,'Vartype','unequal');
fprintf('\nYoungs modulus')
fprintf('\nQuads, healthy: Mean = %.4f, Stdev = %.4f',mean(youngs_norm_q(:,3)),std(youngs_norm_q(:,3)))
fprintf('\nQuads, DMD: Mean = %.4f, Stdev = %.4f',mean(youngs_dmd_q(:,3)),std(youngs_dmd_q(:,3)))
fprintf('\np = %.8f',p)


% Figures

figure;
boxplot([spring_norm_b(:,3);spring_dmd_b(:,3)],labels)
title('Spring constant')
set(gca,'fontsize',12,'LineWidth',1.5);
set(findobj('Tag','Box'),'Color','k','LineWidth',1.5)
set(findobj('Tag','Median'),'Color','b','LineWidth',1)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
ylabel('N/m')
pause(1);


figure; boxplot([strength_norm_b(:,3);strength_dmd_b(:,3)],labels)
title('Strength coefficient')
set(gca,'fontsize',12,'LineWidth',1.5);
set(findobj('Tag','Box'),'Color','k','LineWidth',1.5)
set(findobj('Tag','Median'),'Color','b','LineWidth',1)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
ylabel('kPa')
pause(1);


figure; boxplot([expo_norm_b(:,3);expo_dmd_b(:,3)],labels)
title('Nonlinear strain exponent')
set(gca,'fontsize',12,'LineWidth',1.5);
set(findobj('Tag','Box'),'Color','k','LineWidth',1.5)
set(findobj('Tag','Median'),'Color','b','LineWidth',1)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
ylabel('')
pause(1);


figure; boxplot([youngs_norm_b(:,3);youngs_dmd_b(:,3)],labels)
title('Youngs modulus')
set(gca,'fontsize',12,'LineWidth',1.5);
set(findobj('Tag','Box'),'Color','k','LineWidth',1.5)
set(findobj('Tag','Median'),'Color','b','LineWidth',1)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
ylabel('kPa')
pause(1);


figure; boxplot([spring_norm_q(:,3);spring_dmd_q(:,3)],labels)
title('Spring constant')
set(gca,'fontsize',12,'LineWidth',1.5);
set(findobj('Tag','Box'),'Color','k','LineWidth',1.5)
set(findobj('Tag','Median'),'Color','b','LineWidth',1)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
ylabel('N/m')
pause(1);


figure; boxplot([strength_norm_q(:,3);strength_dmd_q(:,3)],labels)
title('Strength coefficient')
set(gca,'fontsize',12,'LineWidth',1.5);
set(findobj('Tag','Box'),'Color','k','LineWidth',1.5)
set(findobj('Tag','Median'),'Color','b','LineWidth',1)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
ylabel('kPa')
pause(1);


figure; boxplot([expo_norm_q(:,3);expo_dmd_q(:,3)],labels)
title('Nonlinear strain exponent')
set(gca,'fontsize',12,'LineWidth',1.5);
set(findobj('Tag','Box'),'Color','k','LineWidth',1.5)
set(findobj('Tag','Median'),'Color','b','LineWidth',1)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
ylabel('')
pause(1);


figure; boxplot([youngs_norm_q(:,3);youngs_dmd_q(:,3)],labels)
title('Youngs modulus')
set(gca,'fontsize',12,'LineWidth',1.5);
set(findobj('Tag','Box'),'Color','k','LineWidth',1.5)
set(findobj('Tag','Median'),'Color','b','LineWidth',1)
set(findobj('Tag','Outliers'),'MarkerEdgeColor','k')
ylabel('kPa')
pause(1);



fprintf('\n')

end

function [k_upper k_lower k_both] = linearSpringConstant(upper,lower,both,f)

lower_force = 1.5;

upper = upper*.05/411; %50mm height = 411 pixels
lower = lower*.05/411; %50mm height = 411 pixels
both = both*.05/411; %50mm height = 411 pixels

start = find( abs(f-lower_force) == min(abs(f(1:30)-lower_force)) );
start = start(1);
peak = find(f == max(f));
    
k_upper = -(f(peak)-f(start)) / (upper(peak)-upper(start));
k_lower = -(f(peak)-f(start)) / (lower(peak)-lower(start));
k_both = -(f(peak)-f(start)) / (both(peak)-both(start));


end
function [E_upper E_lower E_both] = YoungsModulus(upper,lower,both,f)

A = .013*.047; % Could make bone size


lower_force = 1.5;

start = find( abs(f-lower_force) == min(abs(f(1:30)-lower_force)) );
start = start(1);
peak = find(f == max(f));
    
E_upper = upper(start) * (f(peak)-f(start)) / (1000 * A * (upper(start)-upper(peak))); %kPa
E_lower = lower(start) * (f(peak)-f(start)) / (1000 * A * (lower(start)-lower(peak)));
E_both = both(start) * (f(peak)-f(start)) / (1000 * A * (both(start)-both(peak)));

end
function [Cn_upper Cn_lower Cn_both] = Nonlinear(upper,lower,both,f)

A = .013*.047; % Could make bone size


stress = (f-f(1)) / (1000*A); %kPa
strain_upper = (upper(1)-upper) /upper(1);
strain_lower = (lower(1)-lower) /lower(1);
strain_both = (both(1)-both) /both(1);

stress_upper = stress(find(strain_upper));
strain_upper = strain_upper(find(strain_upper));
stress_lower = stress(find(strain_lower));
strain_lower = strain_lower(find(strain_lower));
stress_both = stress(find(strain_both));
strain_both = strain_both(find(strain_both));

Cn_upper = polyfit(log(strain_upper(1:end)),log(stress_upper(1:end)),1);
Cn_upper = real(Cn_upper);

Cn_lower = polyfit(log(strain_lower(1:end)),log(stress_lower(1:end)),1);
Cn_lower = real(Cn_lower);

Cn_both = polyfit(log(strain_both(1:end)),log(stress_both(1:end)),1);
Cn_both = real(Cn_both);

end



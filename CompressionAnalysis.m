% function CompressionAnalysis
function [Qd,Qn,Bd,Bn] = CompressionAnalysis
    
% dmd_B = 23;
% norm_B = 21;
% dmd_Q = 33;
% norm_Q = 28;

clc
fprintf('\n\n\n');

global quadsData_dmd
quadsData_dmd = struct;
quadsData_dmd(1).patient = '00000';
global quadsData_norm
quadsData_norm = struct;
quadsData_norm(1).patient = '00000';
global bicepsData_dmd
bicepsData_dmd = struct;
bicepsData_dmd(1).patient = '00000';
global bicepsData_norm
bicepsData_norm = struct;
bicepsData_norm(1).patient = '00000';


functionals = getFunc;

d = dir('/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016');
d = d(3:end-1);

%Calculate and store data
for i = 1:length(d)

    [topMus,botMus,wholMus,frc] = getMuscleThick(d(i).name);
    [topMus_c,botMus_c,wholMus_c,frc_c] = getMuscleThick_complete(d(i).name);
    
    k_top = getSpringConstant(topMus,frc);
    k_bottom = getSpringConstant(botMus,frc);
    k_whole = getSpringConstant(wholMus,frc);
    
    E_top = getYoungMod(topMus,frc);
    E_bottom = getYoungMod(botMus,frc);
    E_whole = getYoungMod(wholMus,frc); 
    
    [a_top,b_top] = getNonLinear(topMus,frc);
    [a_bottom,b_bottom] = getNonLinear(botMus,frc);
    [a_whole,b_whole] = getNonLinear(wholMus,frc);
    
    mb_top = doHollomons(topMus_c,frc_c);
    mb_bottom = doHollomons(botMus_c,frc_c);
    mb_whole = doHollomons(wholMus_c,frc_c);
    
%     k = getNewSpringConstant(wholMus,frc);
    
    createEntry(d(i).name,functionals,k_top,k_bottom,k_whole,E_top,E_bottom,E_whole,a_top,b_top,a_bottom,b_bottom,a_whole,b_whole,mb_top,mb_bottom,mb_whole);
    
end

quadsData_dmd = quadsData_dmd(2:end);
quadsData_norm = quadsData_norm(2:end);
bicepsData_dmd = bicepsData_dmd(2:end);
bicepsData_norm = bicepsData_norm(2:end);

Qd = quadsData_dmd;
Qn = quadsData_norm;
Bd = bicepsData_dmd;
Bn = bicepsData_norm;

% nonLinear;
% linearSpringConstant;
% Hollomons;
YoungsModulus;



end


function nonLinear

global bicepsData_dmd
global bicepsData_norm
global quadsData_dmd
global quadsData_norm

for i = 1:length(bicepsData_dmd)
    bd_top_b(:,:,i) = bicepsData_dmd(i).b_top;
    bd_bottom_b(:,:,i) = bicepsData_dmd(i).b_bottom;
    bd_whole_b(:,:,i) = bicepsData_dmd(i).b_whole;
    bd_top_a(:,:,i) = bicepsData_dmd(i).a_top;
    bd_bottom_a(:,:,i) = bicepsData_dmd(i).a_bottom;
    bd_whole_a(:,:,i) = bicepsData_dmd(i).a_whole;
end
for i = 1:length(bicepsData_norm)
    bn_top_b(:,:,i) = bicepsData_norm(i).b_top;
    bn_bottom_b(:,:,i) = bicepsData_norm(i).b_bottom;
    bn_whole_b(:,:,i) = bicepsData_norm(i).b_whole;
    bn_top_a(:,:,i) = bicepsData_norm(i).a_top;
    bn_bottom_a(:,:,i) = bicepsData_norm(i).a_bottom;
    bn_whole_a(:,:,i) = bicepsData_norm(i).a_whole;
end
for i = 1:length(quadsData_dmd)
    qd_top_b(:,:,i) = quadsData_dmd(i).b_top;
    qd_bottom_b(:,:,i) = quadsData_dmd(i).b_bottom;
    qd_whole_b(:,:,i) = quadsData_dmd(i).b_whole;
    qd_top_a(:,:,i) = quadsData_dmd(i).a_top;
    qd_bottom_a(:,:,i) = quadsData_dmd(i).a_bottom;
    qd_whole_a(:,:,i) = quadsData_dmd(i).a_whole;
end
for i = 1:length(quadsData_norm)
    qn_top_b(:,:,i) = quadsData_norm(i).b_top;
    qn_bottom_b(:,:,i) = quadsData_norm(i).b_bottom;
    qn_whole_b(:,:,i) = quadsData_norm(i).b_whole;
    qn_top_a(:,:,i) = quadsData_norm(i).a_top;
    qn_bottom_a(:,:,i) = quadsData_norm(i).a_bottom;
    qn_whole_a(:,:,i) = quadsData_norm(i).a_whole;
end

fprintf('\nDoes nonlinear _a_ differentiate DMD/normal? (t-test)\n')
checkStats(.001,bd_top_a,bd_bottom_a,bd_whole_a, bn_top_a,bn_bottom_a,bn_whole_a, qd_top_a,qd_bottom_a,qd_whole_a, qn_top_a,qn_bottom_a,qn_whole_a);
plotBox(bd_top_a,bn_top_a)
title('Biceps, upper muscle, a')
plotBox(bd_whole_a,bn_whole_a)
title('Biceps, whole muscle, a')
plotBox(bd_bottom_a,bn_bottom_a)
title('Biceps, lower muscle, a')
plotBox(qd_top_a,qn_top_a)
title('Quads, upper muscle, a')
plotBox(qd_whole_a,qn_whole_a)
title('Quads, whole muscle, a')
plotBox(qd_bottom_a,qn_bottom_a)
title('Quads, lower muscle, a')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _a_ correlate with age?\n');
fprintf('Biceps, upper muscle, DMD');
corrVsAge_nonlinear(bd_top_a);
fprintf('Biceps, lower muscle, DMD');
corrVsAge_nonlinear(bd_bottom_a);
fprintf('Biceps, whole muscle, DMD');
corrVsAge_nonlinear(bd_whole_a);
fprintf('Quads, upper muscle, DMD');
corrVsAge_nonlinear(qd_top_a);
fprintf('Quads, lower muscle, DMD');
corrVsAge_nonlinear(qd_bottom_a);
fprintf('Quads, whole muscle, DMD');
corrVsAge_nonlinear(qd_whole_a);

plotVsAge(bd_top_a,bn_top_a);
title('Biceps, upper muscle');
ylabel('a'); xlabel('Age');
plotVsAge(bd_bottom_a,bn_bottom_a);
title('Biceps, lower muscle');
ylabel('a'); xlabel('Age');
plotVsAge(bd_whole_a,bn_whole_a);
title('Biceps, whole muscle');
ylabel('a'); xlabel('Age');
plotVsAge(qd_top_a,qn_top_a);
title('Quads, upper muscle');
ylabel('a'); xlabel('Age');
plotVsAge(qd_bottom_a,qn_bottom_a);
title('Quads, lower muscle');
ylabel('a'); xlabel('Age');
plotVsAge(qd_whole_a,qn_whole_a);
title('Quads, whole muscle');
ylabel('a'); xlabel('Age');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _a_ correlate with 6MW?\n');
fprintf('Quads, upper muscle, DMD');
corrVs6MW_nonlinear(qd_top_a);
fprintf('Quads, lower muscle, DMD');
corrVs6MW_nonlinear(qd_bottom_a);
fprintf('Quads, whole muscle, DMD');
corrVs6MW_nonlinear(qd_whole_a);

plotVs6MW(qd_top_a);
title('Quads, upper muscle');
ylabel('a'); xlabel('6MW');
plotVs6MW(qd_bottom_a);
title('Quads, lower muscle');
ylabel('a'); xlabel('6MW');
plotVs6MW(qd_whole_a);
title('Quads, whole muscle');
ylabel('a'); xlabel('6MW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _a_ correlate with 10m run?\n');
fprintf('Quads, upper muscle, DMD');
corrVs10m_nonlinear(qd_top_a);
fprintf('Quads, lower muscle, DMD');
corrVs10m_nonlinear(qd_bottom_a);
fprintf('Quads, whole muscle, DMD');
corrVs10m_nonlinear(qd_whole_a);

plotVs10m(qd_top_a);
title('Quads, upper muscle');
ylabel('a'); xlabel('10m time');
plotVs10m(qd_bottom_a);
title('Quads, lower muscle');
ylabel('a'); xlabel('10m time');
plotVs10m(qd_whole_a);
title('Quads, whole muscle');
ylabel('a'); xlabel('10m time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _a_ correlate with supine-to-stand?\n');
fprintf('Quads, upper muscle, DMD');
corrVsSupStand_nonlinear(qd_top_a);
fprintf('Quads, lower muscle, DMD');
corrVsSupStand_nonlinear(qd_bottom_a);
fprintf('Quads, whole muscle, DMD');
corrVsSupStand_nonlinear(qd_whole_a);

plotVsSupStand(qd_top_a);
title('Quads, upper muscle');
ylabel('a'); xlabel('Sup/stand time');
plotVsSupStand(qd_bottom_a);
title('Quads, lower muscle');
ylabel('a'); xlabel('Sup/stand time');
plotVsSupStand(qd_whole_a);
title('Quads, whole muscle');
ylabel('a'); xlabel('Sup/stand time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _a_ correlate with NSAA?\n');
fprintf('Quads, upper muscle, DMD');
corrVsNSAA_nonlinear(qd_top_a);
fprintf('Quads, lower muscle, DMD');
corrVsNSAA_nonlinear(qd_bottom_a);
fprintf('Quads, whole muscle, DMD');
corrVsNSAA_nonlinear(qd_whole_a);

plotVsNSAA(qd_top_a);
title('Quads, upper muscle');
ylabel('a'); xlabel('NSAA');
plotVsNSAA(qd_bottom_a);
title('Quads, lower muscle');
ylabel('a'); xlabel('NSAA');
plotVsNSAA(qd_whole_a);
title('Quads, whole muscle');
ylabel('a'); xlabel('NSAA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _a_ correlate with HHD (quads)?\n');
fprintf('Quads, upper muscle, DMD');
corrVsHHD_q_nonlinear(qd_top_a);
fprintf('Quads, lower muscle, DMD');
corrVsHHD_q_nonlinear(qd_bottom_a);
fprintf('Quads, whole muscle, DMD');
corrVsHHD_q_nonlinear(qd_whole_a);

plotVsHHD_q(qd_top_a);
title('Quads, upper muscle');
ylabel('a'); xlabel('HHD, quads');
plotVsHHD_q(qd_bottom_a);
title('Quads, lower muscle');
ylabel('a'); xlabel('HHD, quads');
plotVsHHD_q(qd_whole_a);
title('Quads, whole muscle');
ylabel('a'); xlabel('HHD, quads');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _a_ correlate with HHD (biceps)?\n');
fprintf('Biceps, upper muscle, DMD');
corrVsHHD_b_nonlinear(bd_top_a);
fprintf('Biceps, lower muscle, DMD');
corrVsHHD_b_nonlinear(bd_bottom_a);
fprintf('Biceps, whole muscle, DMD');
corrVsHHD_b_nonlinear(bd_whole_a);

plotVsHHD_b(bd_top_a);
title('Biceps, upper muscle');
ylabel('a'); xlabel('HHD, biceps');
plotVsHHD_b(bd_bottom_a);
title('Biceps, lower muscle');
ylabel('a'); xlabel('HHD, biceps');
plotVsHHD_b(bd_whole_a);
title('Biceps, whole muscle');
ylabel('a'); xlabel('HHD, biceps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _a_ correlate with isometric pushup?\n');
fprintf('Biceps, upper muscle, DMD');
corrVsIsoPush_nonlinear(bd_top_a);
fprintf('Biceps, lower muscle, DMD');
corrVsIsoPush_nonlinear(bd_bottom_a);
fprintf('Biceps, whole muscle, DMD');
corrVsIsoPush_nonlinear(bd_whole_a);

plotVsIsoPush(bd_top_a);
title('Biceps, upper muscle');
ylabel('a'); xlabel('Isometric pushup');
plotVsIsoPush(bd_bottom_a);
title('Biceps, lower muscle');
ylabel('a'); xlabel('Isometric pushup');
plotVsIsoPush(bd_whole_a);
title('Biceps, whole muscle');
ylabel('a'); xlabel('Isometric pushup');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










fprintf('\nDoes nonlinear _b_ differentiate DMD/normal? (t-test)\n')
checkStats(.001,bd_top_b,bd_bottom_b,bd_whole_b, bn_top_b,bn_bottom_b,bn_whole_b, qd_top_b,qd_bottom_b,qd_whole_b, qn_top_b,qn_bottom_b,qn_whole_b);
plotBox(bd_top_b,bn_top_b)
title('Biceps, upper muscle, b')
plotBox(bd_whole_b,bn_whole_b)
title('Biceps, whole muscle, b')
plotBox(bd_bottom_b,bn_bottom_b)
title('Biceps, lower muscle, b')
plotBox(qd_top_b,qn_top_b)
title('Quads, upper muscle, b')
plotBox(qd_whole_b,qn_whole_b)
title('Quads, whole muscle, b')
plotBox(qd_bottom_b,qn_bottom_b)
title('Quads, lower muscle, b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _b_ correlate with age?\n');
fprintf('Biceps, upper muscle, DMD');
corrVsAge_nonlinear(bd_top_b);
fprintf('Biceps, lower muscle, DMD');
corrVsAge_nonlinear(bd_bottom_b);
fprintf('Biceps, whole muscle, DMD');
corrVsAge_nonlinear(bd_whole_b);
fprintf('Quads, upper muscle, DMD');
corrVsAge_nonlinear(qd_top_b);
fprintf('Quads, lower muscle, DMD');
corrVsAge_nonlinear(qd_bottom_b);
fprintf('Quads, whole muscle, DMD');
corrVsAge_nonlinear(qd_whole_b);

plotVsAge(bd_top_b,bn_top_b);
title('Biceps, upper muscle');
ylabel('b'); xlabel('Age');
plotVsAge(bd_bottom_b,bn_bottom_b);
title('Biceps, lower muscle');
ylabel('b'); xlabel('Age');
plotVsAge(bd_whole_b,bn_whole_b);
title('Biceps, whole muscle');
ylabel('b'); xlabel('Age');
plotVsAge(qd_top_b,qn_top_b);
title('Quads, upper muscle');
ylabel('b'); xlabel('Age');
plotVsAge(qd_bottom_b,qn_bottom_b);
title('Quads, lower muscle');
ylabel('b'); xlabel('Age');
plotVsAge(qd_whole_b,qn_whole_b);
title('Quads, whole muscle');
ylabel('b'); xlabel('Age');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _b_ correlate with 6MW?\n');
fprintf('Quads, upper muscle, DMD');
corrVs6MW_nonlinear(qd_top_b);
fprintf('Quads, lower muscle, DMD');
corrVs6MW_nonlinear(qd_bottom_b);
fprintf('Quads, whole muscle, DMD');
corrVs6MW_nonlinear(qd_whole_b);

plotVs6MW(qd_top_b);
title('Quads, upper muscle');
ylabel('b'); xlabel('6MW');
plotVs6MW(qd_bottom_b);
title('Quads, lower muscle');
ylabel('b'); xlabel('6MW');
plotVs6MW(qd_whole_b);
title('Quads, whole muscle');
ylabel('b'); xlabel('6MW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _b_ correlate with 10m run?\n');
fprintf('Quads, upper muscle, DMD');
corrVs10m_nonlinear(qd_top_b);
fprintf('Quads, lower muscle, DMD');
corrVs10m_nonlinear(qd_bottom_b);
fprintf('Quads, whole muscle, DMD');
corrVs10m_nonlinear(qd_whole_b);

plotVs10m(qd_top_b);
title('Quads, upper muscle');
ylabel('b'); xlabel('10m time');
plotVs10m(qd_bottom_b);
title('Quads, lower muscle');
ylabel('b'); xlabel('10m time');
plotVs10m(qd_whole_b);
title('Quads, whole muscle');
ylabel('b'); xlabel('10m time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _b_ correlate with supine-to-stand?\n');
fprintf('Quads, upper muscle, DMD');
corrVsSupStand_nonlinear(qd_top_b);
fprintf('Quads, lower muscle, DMD');
corrVsSupStand_nonlinear(qd_bottom_b);
fprintf('Quads, whole muscle, DMD');
corrVsSupStand_nonlinear(qd_whole_b);

plotVsSupStand(qd_top_b);
title('Quads, upper muscle');
ylabel('b'); xlabel('Sup/stand time');
plotVsSupStand(qd_bottom_b);
title('Quads, lower muscle');
ylabel('b'); xlabel('Sup/stand time');
plotVsSupStand(qd_whole_b);
title('Quads, whole muscle');
ylabel('b'); xlabel('Sup/stand time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _b_ correlate with NSAA?\n');
fprintf('Quads, upper muscle, DMD');
corrVsNSAA_nonlinear(qd_top_b);
fprintf('Quads, lower muscle, DMD');
corrVsNSAA_nonlinear(qd_bottom_b);
fprintf('Quads, whole muscle, DMD');
corrVsNSAA_nonlinear(qd_whole_b);

plotVsNSAA(qd_top_b);
title('Quads, upper muscle');
ylabel('b'); xlabel('NSAA');
plotVsNSAA(qd_bottom_b);
title('Quads, lower muscle');
ylabel('b'); xlabel('NSAA');
plotVsNSAA(qd_whole_b);
title('Quads, whole muscle');
ylabel('b'); xlabel('NSAA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _b_ correlate with HHD (quads)?\n');
fprintf('Quads, upper muscle, DMD');
corrVsHHD_q_nonlinear(qd_top_b);
fprintf('Quads, lower muscle, DMD');
corrVsHHD_q_nonlinear(qd_bottom_b);
fprintf('Quads, whole muscle, DMD');
corrVsHHD_q_nonlinear(qd_whole_b);

plotVsHHD_q(qd_top_b);
title('Quads, upper muscle');
ylabel('b'); xlabel('HHD, quads');
plotVsHHD_q(qd_bottom_b);
title('Quads, lower muscle');
ylabel('b'); xlabel('HHD, quads');
plotVsHHD_q(qd_whole_b);
title('Quads, whole muscle');
ylabel('b'); xlabel('HHD, quads');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _b_ correlate with HHD (biceps)?\n');
fprintf('Biceps, upper muscle, DMD');
corrVsHHD_b_nonlinear(bd_top_b);
fprintf('Biceps, lower muscle, DMD');
corrVsHHD_b_nonlinear(bd_bottom_b);
fprintf('Biceps, whole muscle, DMD');
corrVsHHD_b_nonlinear(bd_whole_b);

plotVsHHD_b(bd_top_b);
title('Biceps, upper muscle');
ylabel('b'); xlabel('HHD, biceps');
plotVsHHD_b(bd_bottom_b);
title('Biceps, lower muscle');
ylabel('b'); xlabel('HHD, biceps');
plotVsHHD_b(bd_whole_b);
title('Biceps, whole muscle');
ylabel('b'); xlabel('HHD, biceps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes nonlinear _b_ correlate with isometric pushup?\n');
fprintf('Biceps, upper muscle, DMD');
corrVsIsoPush_nonlinear(bd_top_b);
fprintf('Biceps, lower muscle, DMD');
corrVsIsoPush_nonlinear(bd_bottom_b);
fprintf('Biceps, whole muscle, DMD');
corrVsIsoPush_nonlinear(bd_whole_b);

plotVsIsoPush(bd_top_b);
title('Biceps, upper muscle');
ylabel('b'); xlabel('Isometric pushup');
plotVsIsoPush(bd_bottom_b);
title('Biceps, lower muscle');
ylabel('b'); xlabel('Isometric pushup');
plotVsIsoPush(bd_whole_b);
title('Biceps, whole muscle');
ylabel('b'); xlabel('Isometric pushup');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
function linearSpringConstant

global bicepsData_dmd
global bicepsData_norm
global quadsData_dmd
global quadsData_norm

for i = 1:length(bicepsData_dmd)
    bd_top(:,:,i) = bicepsData_dmd(i).k_top;
    bd_bottom(:,:,i) = bicepsData_dmd(i).k_bottom;
    bd_whole(:,:,i) = bicepsData_dmd(i).k_whole;
end
for i = 1:length(bicepsData_norm)
    bn_top(:,:,i) = bicepsData_norm(i).k_top;
    bn_bottom(:,:,i) = bicepsData_norm(i).k_bottom;
    bn_whole(:,:,i) = bicepsData_norm(i).k_whole;
end
for i = 1:length(quadsData_dmd)
    qd_top(:,:,i) = quadsData_dmd(i).k_top;
    qd_bottom(:,:,i) = quadsData_dmd(i).k_bottom;
    qd_whole(:,:,i) = quadsData_dmd(i).k_whole;
end
for i = 1:length(quadsData_norm)
    qn_top(:,:,i) = quadsData_norm(i).k_top;
    qn_bottom(:,:,i) = quadsData_norm(i).k_bottom;
    qn_whole(:,:,i) = quadsData_norm(i).k_whole;
end



fprintf('\nDoes spring constant differentiate DMD/normal? (t-test)\n')
checkStats(.01,bd_top,bd_bottom,bd_whole(r1,r2,:), bn_top,bn_bottom,bn_whole(r1,r2,:), qd_top,qd_bottom,qd_whole(r1,r2,:), qn_top,qn_bottom,qn_whole(r1,r2,:));
% checkStats(.01,bd_top,bd_bottom,bd_whole, bn_top,bn_bottom,bn_whole, qd_top,qd_bottom,qd_whole, qn_top,qn_bottom,qn_whole);
% plotBox(bd_top(r1,r2,:),bn_top(r1,r2,:))
% title('Biceps, upper muscle, n')
plotBox(bd_whole(r1,r2,:),bn_whole(r1,r2,:))
% title('Spring constant, biceps, whole muscle, 3-10N')
% plotBox(bd_bottom(r1,r2,:),bn_bottom(r1,r2,:))
% title('Biceps, lower muscle, n')
% plotBox(qd_top(r1,r2,:),qn_top(r1,r2,:))
% title('Quads, upper muscle, n')
plotBox(qd_whole(r1,r2,:),qn_whole(r1,r2,:))
% title('Spring constant, quads, whole muscle, 3-10N')
% plotBox(qd_bottom(r1,r2,:),qn_bottom(r1,r2,:))
% title('Quads, lower muscle, n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes spring constant correlate with age?\n');
% fprintf('Biceps, upper muscle, DMD');
% corrVsAge(bd_top);
% fprintf('Biceps, lower muscle, DMD');
% corrVsAge(bd_bottom);
% fprintf('Biceps, whole muscle, DMD');
% corrVsAge(bd_whole);
% fprintf('Quads, upper muscle, DMD');
% corrVsAge(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVsAge(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVsAge(qd_whole);

% plotVsAge(bd_top(r1,r2,:),bn_top(r1,r2,:));
% title('Biceps, upper muscle, 1-10N');
% ylabel('Spring constant'); xlabel('Age');
% plotVsAge(bd_bottom(r1,r2,:),bn_bottom(r1,r2,:));
% title('Biceps, lower muscle');
% ylabel('Spring constant'); xlabel('Age');
% plotVsAge(bd_whole(r1,r2,:),bn_whole(r1,r2,:));
% title('Biceps, whole muscle');
% ylabel('Spring constant'); xlabel('Age');
% plotVsAge(qd_top(r1,r2,:),qn_top(r1,r2,:));
% title('Quads, upper muscle');
% ylabel('Spring constant'); xlabel('Age');
% plotVsAge(qd_bottom(r1,r2,:),qn_bottom(r1,r2,:));
% title('Quads, lower muscle, 2-10N');
% ylabel('Spring constant'); xlabel('Age');
% plotVsAge(qd_whole(r1,r2,:),qn_whole(r1,r2,:));
% title('Quads, whole muscle');
% ylabel('Spring constant'); xlabel('Age');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes spring constant correlate with 6MW?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVs6MW(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVs6MW(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVs6MW(qd_whole);
% 
% plotVs6MW(qd_top(r1,r2,:));
% title('Quads, upper muscle');
% ylabel('Spring constant'); xlabel('6MW');
% plotVs6MW(qd_bottom(r1,r2,:));
% title('Quads, lower muscle');
% ylabel('Spring constant'); xlabel('6MW');
% plotVs6MW(qd_whole(r1,r2,:));
% title('Quads, whole muscle, 2-9N');
% ylabel('Spring constant'); xlabel('6MW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes spring constant correlate with 10m run?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVs10m(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVs10m(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVs10m(qd_whole);
% 
% plotVs10m(qd_top(r1,r2,:));
% title('Quads, upper muscle');
% ylabel('Spring constant'); xlabel('10m time');
% plotVs10m(qd_bottom(r1,r2,:));
% title('Quads, lower muscle, 1-10N');
% ylabel('Spring constant'); xlabel('10m time');
% plotVs10m(qd_whole(r1,r2,:));
% title('Quads, whole muscle');
% ylabel('Spring constant'); xlabel('10m time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes spring constant correlate with supine-to-stand?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsSupStand(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVsSupStand(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVsSupStand(qd_whole);
% % 
% plotVsSupStand(qd_top(r1,r2,:));
% title('Quads, upper muscle, 2-10N');
% ylabel('Spring constant'); xlabel('Sup/stand time');
% plotVsSupStand(qd_bottom(r1,r2,:));
% title('Quads, lower muscle');
% ylabel('Spring constant'); xlabel('Sup/stand time');
% plotVsSupStand(qd_whole(r1,r2,:));
% title('Quads, whole muscle');
% ylabel('Spring constant'); xlabel('Sup/stand time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes spring constant correlate with NSAA?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsNSAA(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVsNSAA(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVsNSAA(qd_whole);
% 
% plotVsNSAA(qd_top(r1,r2,:));
% title('Quads, upper muscle');
% ylabel('Spring constant'); xlabel('NSAA');
% plotVsNSAA(qd_bottom(r1,r2,:));
% title('Quads, lower muscle');
% ylabel('Spring constant'); xlabel('NSAA');
% plotVsNSAA(qd_whole(r1,r2,:));
% title('Quads, whole muscle');
% ylabel('Spring constant'); xlabel('NSAA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes spring constant correlate with HHD (quads)?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsHHD_q(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVsHHD_q(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVsHHD_q(qd_whole);
% 
% plotVsHHD_q(qd_top(r1,r2,:));
% title('Quads, upper muscle');
% ylabel('Spring constant'); xlabel('HHD, quads');
% plotVsHHD_q(qd_bottom(r1,r2,:));
% title('Quads, lower muscle');
% ylabel('Spring constant'); xlabel('HHD, quads');
% plotVsHHD_q(qd_whole(r1,r2,:));
% title('Quads, whole muscle');
% ylabel('Spring constant'); xlabel('HHD, quads');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
function Hollomons

global bicepsData_dmd
global bicepsData_norm
global quadsData_dmd
global quadsData_norm

for i = 1:length(bicepsData_dmd)
    bd_top_n(:,:,i) = bicepsData_dmd(i).n_top;
    bd_top_K(:,:,i) = bicepsData_dmd(i).K_top;
    bd_bottom_n(:,:,i) = bicepsData_dmd(i).n_bottom;
    bd_bottom_K(:,:,i) = bicepsData_dmd(i).K_bottom;
    bd_whole_n(:,:,i) = bicepsData_dmd(i).n_whole;
    bd_whole_K(:,:,i) = bicepsData_dmd(i).K_whole;
end
for i = 1:length(bicepsData_norm)
    bn_top_n(:,:,i) = bicepsData_norm(i).n_top;
    bn_top_K(:,:,i) = bicepsData_norm(i).K_top;
    bn_bottom_n(:,:,i) = bicepsData_norm(i).n_bottom;
    bn_bottom_K(:,:,i) = bicepsData_norm(i).K_bottom;
    bn_whole_n(:,:,i) = bicepsData_norm(i).n_whole;
    bn_whole_K(:,:,i) = bicepsData_norm(i).K_whole;
end
for i = 1:length(quadsData_dmd)
    qd_top_n(:,:,i) = quadsData_dmd(i).n_top;
    qd_top_K(:,:,i) = quadsData_dmd(i).K_top;
    qd_bottom_n(:,:,i) = quadsData_dmd(i).n_bottom;
    qd_bottom_K(:,:,i) = quadsData_dmd(i).K_bottom;
    qd_whole_n(:,:,i) = quadsData_dmd(i).n_whole;
    qd_whole_K(:,:,i) = quadsData_dmd(i).K_whole;
end
for i = 1:length(quadsData_norm)
    qn_top_n(:,:,i) = quadsData_norm(i).n_top;
    qn_top_K(:,:,i) = quadsData_norm(i).K_top;
    qn_bottom_n(:,:,i) = quadsData_norm(i).n_bottom;
    qn_bottom_K(:,:,i) = quadsData_norm(i).K_bottom;
    qn_whole_n(:,:,i) = quadsData_norm(i).n_whole;
    qn_whole_K(:,:,i) = quadsData_norm(i).K_whole;
end




fprintf('\nDoes nonlinear _n_ differentiate DMD/normal? (t-test)\n')
checkStats(.01,bd_top_n,bd_bottom_n,bd_whole_n, bn_top_n,bn_bottom_n,bn_whole_n, qd_top_n,qd_bottom_n,qd_whole_n, qn_top_n,qn_bottom_n,qn_whole_n);

% plotBox(bd_top_n,bn_top_n)
% title('Biceps, upper muscle, n')
plotBox(bd_whole_n,bn_whole_n)
% title('Strain hardening exponent, biceps, whole muscle')
% plotBox(bd_bottom_n,bn_bottom_n)
% title('Biceps, lower muscle, n')
% plotBox(qd_top_n,qn_top_n)
% title('Quads, upper muscle, n')
plotBox(qd_whole_n,qn_whole_n)
% title('Strain hardening exponent, quads, whole muscle')
% plotBox(qd_bottom_n,qn_bottom_n)
% title('Quads, lower muscle, n')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _n_ correlate with age?\n');
% fprintf('Biceps, upper muscle, DMD');
% corrVsAge_nonlinear(bd_top_n);
% fprintf('Biceps, lower muscle, DMD');
% corrVsAge_nonlinear(bd_bottom_n);
% fprintf('Biceps, whole muscle, DMD');
% corrVsAge_nonlinear(bd_whole_n);
% fprintf('Quads, upper muscle, DMD');
% corrVsAge_nonlinear(qd_top_n);
% fprintf('Quads, lower muscle, DMD');
% corrVsAge_nonlinear(qd_bottom_n);
% fprintf('Quads, whole muscle, DMD');
% corrVsAge_nonlinear(qd_whole_n);

% plotVsAge(bd_top_n,bn_top_n);
% title('Biceps, upper muscle');
% ylabel('n'); xlabel('Age');
% plotVsAge(bd_bottom_n,bn_bottom_n);
% title('Biceps, lower muscle');
% ylabel('n'); xlabel('Age');
% plotVsAge(bd_whole_n,bn_whole_n);
% title('Biceps, whole muscle');
% ylabel('n'); xlabel('Age');
% plotVsAge(qd_top_n,qn_top_n);
% title('Quads, upper muscle');
% ylabel('n'); xlabel('Age');
% plotVsAge(qd_bottom_n,qn_bottom_n);
% title('Quads, lower muscle');
% ylabel('n'); xlabel('Age');
% plotVsAge(qd_whole_n,qn_whole_n);
% title('Quads, whole muscle');
% ylabel('n'); xlabel('Age');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _n_ correlate with 6MW?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVs6MW_nonlinear(qd_top_n);
% fprintf('Quads, lower muscle, DMD');
% corrVs6MW_nonlinear(qd_bottom_n);
% fprintf('Quads, whole muscle, DMD');
% corrVs6MW_nonlinear(qd_whole_n);
% 
% plotVs6MW(qd_top_n);
% title('Quads, upper muscle');
% ylabel('n'); xlabel('6MW');
% plotVs6MW(qd_bottom_n);
% title('Quads, lower muscle');
% ylabel('n'); xlabel('6MW');
% plotVs6MW(qd_whole_n);
% title('Quads, whole muscle');
% ylabel('n'); xlabel('6MW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _n_ correlate with 10m run?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVs10m_nonlinear(qd_top_n);
% fprintf('Quads, lower muscle, DMD');
% corrVs10m_nonlinear(qd_bottom_n);
% fprintf('Quads, whole muscle, DMD');
% corrVs10m_nonlinear(qd_whole_n);
% 
% plotVs10m(qd_top_n);
% title('Quads, upper muscle');
% ylabel('n'); xlabel('10m time');
% plotVs10m(qd_bottom_n);
% title('Quads, lower muscle');
% ylabel('n'); xlabel('10m time');
% plotVs10m(qd_whole_n);
% title('Quads, whole muscle');
% ylabel('n'); xlabel('10m time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _n_ correlate with supine-to-stand?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsSupStand_nonlinear(qd_top_n);
% fprintf('Quads, lower muscle, DMD');
% corrVsSupStand_nonlinear(qd_bottom_n);
% fprintf('Quads, whole muscle, DMD');
% corrVsSupStand_nonlinear(qd_whole_n);
% 
% plotVsSupStand(qd_top_n);
% title('Quads, upper muscle');
% ylabel('n'); xlabel('Sup/stand time');
% plotVsSupStand(qd_bottom_n);
% title('Quads, lower muscle');
% ylabel('n'); xlabel('Sup/stand time');
% plotVsSupStand(qd_whole_n);
% title('Quads, whole muscle');
% ylabel('n'); xlabel('Sup/stand time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _n_ correlate with NSAA?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsNSAA_nonlinear(qd_top_n);
% fprintf('Quads, lower muscle, DMD');
% corrVsNSAA_nonlinear(qd_bottom_n);
% fprintf('Quads, whole muscle, DMD');
% corrVsNSAA_nonlinear(qd_whole_n);
% 
% plotVsNSAA(qd_top_n);
% title('Quads, upper muscle');
% ylabel('n'); xlabel('NSAA');
% plotVsNSAA(qd_bottom_n);
% title('Quads, lower muscle');
% ylabel('n'); xlabel('NSAA');
% plotVsNSAA(qd_whole_n);
% title('Quads, whole muscle');
% ylabel('n'); xlabel('NSAA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _n_ correlate with HHD (quads)?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsHHD_q_nonlinear(qd_top_n);
% fprintf('Quads, lower muscle, DMD');
% corrVsHHD_q_nonlinear(qd_bottom_n);
% fprintf('Quads, whole muscle, DMD');
% corrVsHHD_q_nonlinear(qd_whole_n);
% 
% plotVsHHD_q(qd_top_n);
% title('Quads, upper muscle');
% ylabel('n'); xlabel('HHD, quads');
% plotVsHHD_q(qd_bottom_n);
% title('Quads, lower muscle');
% ylabel('n'); xlabel('HHD, quads');
% plotVsHHD_q(qd_whole_n);
% title('Quads, whole muscle');
% ylabel('n'); xlabel('HHD, quads');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





fprintf('\nDoes nonlinear _K_ differentiate DMD/normal? (t-test)\n')
checkStats(.01,bd_top_K,bd_bottom_K,bd_whole_K, bn_top_K,bn_bottom_K,bn_whole_K, qd_top_K,qd_bottom_K,qd_whole_K, qn_top_K,qn_bottom_K,qn_whole_K);
% plotBox(bd_top_K,bn_top_K)
% title('Biceps, upper muscle, K')
plotBox(bd_whole_K,bn_whole_K)
% title('Strength index, biceps, whole muscle')
% plotBox(bd_bottom_K,bn_bottom_K)
% title('Biceps, lower muscle, K')
% plotBox(qd_top_K,qn_top_K)
% title('Quads, upper muscle, K')
plotBox(qd_whole_K,qn_whole_K)
% title('Strength index, quads, whole muscle')
% plotBox(qd_bottom_K,qn_bottom_K)
% title('Quads, lower muscle, K')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _K_ correlate with age?\n');
% fprintf('Biceps, upper muscle, DMD');
% corrVsAge_nonlinear(bd_top_K);
% fprintf('Biceps, lower muscle, DMD');
% corrVsAge_nonlinear(bd_bottom_K);
% fprintf('Biceps, whole muscle, DMD');
% corrVsAge_nonlinear(bd_whole_K);
% fprintf('Quads, upper muscle, DMD');
% corrVsAge_nonlinear(qd_top_K);
% fprintf('Quads, lower muscle, DMD');
% corrVsAge_nonlinear(qd_bottom_K);
% fprintf('Quads, whole muscle, DMD');
% corrVsAge_nonlinear(qd_whole_K);

% plotVsAge(bd_top_K,bn_top_K);
% title('Biceps, upper muscle');
% ylabel('K'); xlabel('Age');
% plotVsAge(bd_bottom_K,bn_bottom_K);
% title('Biceps, lower muscle');
% ylabel('K'); xlabel('Age');
% plotVsAge(bd_whole_K,bn_whole_K);
% title('Biceps, whole muscle');
% ylabel('K'); xlabel('Age');
% plotVsAge(qd_top_K,qn_top_K);
% title('Quads, upper muscle');
% ylabel('K'); xlabel('Age');
% plotVsAge(qd_bottom_K,qn_bottom_K);
% title('Quads, lower muscle');
% ylabel('K'); xlabel('Age');
% plotVsAge(qd_whole_K,qn_whole_K);
% title('Quads, whole muscle');
% ylabel('K'); xlabel('Age');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _K_ correlate with 6MW?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVs6MW_nonlinear(qd_top_K);
% fprintf('Quads, lower muscle, DMD');
% corrVs6MW_nonlinear(qd_bottom_K);
% fprintf('Quads, whole muscle, DMD');
% corrVs6MW_nonlinear(qd_whole_K);
% 
% plotVs6MW(qd_top_K);
% title('Quads, upper muscle');
% ylabel('K'); xlabel('6MW');
% plotVs6MW(qd_bottom_K);
% title('Quads, lower muscle');
% ylabel('K'); xlabel('6MW');
% plotVs6MW(qd_whole_K);
% title('Quads, whole muscle');
% ylabel('K'); xlabel('6MW');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _K_ correlate with 10m run?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVs10m_nonlinear(qd_top_K);
% fprintf('Quads, lower muscle, DMD');
% corrVs10m_nonlinear(qd_bottom_K);
% fprintf('Quads, whole muscle, DMD');
% corrVs10m_nonlinear(qd_whole_K);
% 
% plotVs10m(qd_top_K);
% title('Quads, upper muscle');
% ylabel('K'); xlabel('10m time');
% plotVs10m(qd_bottom_K);
% title('Quads, lower muscle');
% ylabel('K'); xlabel('10m time');
% plotVs10m(qd_whole_K);
% title('Quads, whole muscle');
% ylabel('K'); xlabel('10m time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _K_ correlate with supine-to-stand?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsSupStand_nonlinear(qd_top_K);
% fprintf('Quads, lower muscle, DMD');
% corrVsSupStand_nonlinear(qd_bottom_K);
% fprintf('Quads, whole muscle, DMD');
% corrVsSupStand_nonlinear(qd_whole_K);
% 
% plotVsSupStand(qd_top_K);
% title('Quads, upper muscle');
% ylabel('K'); xlabel('Sup/stand time');
% plotVsSupStand(qd_bottom_K);
% title('Quads, lower muscle');
% ylabel('K'); xlabel('Sup/stand time');
% plotVsSupStand(qd_whole_K);
% title('Quads, whole muscle');
% ylabel('K'); xlabel('Sup/stand time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _K_ correlate with NSAA?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsNSAA_nonlinear(qd_top_K);
% fprintf('Quads, lower muscle, DMD');
% corrVsNSAA_nonlinear(qd_bottom_K);
% fprintf('Quads, whole muscle, DMD');
% corrVsNSAA_nonlinear(qd_whole_K);
% 
% plotVsNSAA(qd_top_K);
% title('Quads, upper muscle');
% ylabel('K'); xlabel('NSAA');
% plotVsNSAA(qd_bottom_K);
% title('Quads, lower muscle');
% ylabel('K'); xlabel('NSAA');
% plotVsNSAA(qd_whole_K);
% title('Quads, whole muscle');
% ylabel('K'); xlabel('NSAA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes nonlinear _K_ correlate with HHD (quads)?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsHHD_q_nonlinear(qd_top_K);
% fprintf('Quads, lower muscle, DMD');
% corrVsHHD_q_nonlinear(qd_bottom_K);
% fprintf('Quads, whole muscle, DMD');
% corrVsHHD_q_nonlinear(qd_whole_K);
% 
% plotVsHHD_q(qd_top_K);
% title('Quads, upper muscle');
% ylabel('K'); xlabel('HHD, quads');
% plotVsHHD_q(qd_bottom_K);
% title('Quads, lower muscle');
% ylabel('K'); xlabel('HHD, quads');
% plotVsHHD_q(qd_whole_K);
% title('Quads, whole muscle');
% ylabel('K'); xlabel('HHD, quads');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
function YoungsModulus

r1 = 7;
r2 = 9;


global bicepsData_dmd
global bicepsData_norm
global quadsData_dmd
global quadsData_norm

for i = 1:length(bicepsData_dmd)
    bd_top(:,:,i) = bicepsData_dmd(i).E_top;
    bd_bottom(:,:,i) = bicepsData_dmd(i).E_bottom;
    bd_whole(:,:,i) = bicepsData_dmd(i).E_whole;
end
for i = 1:length(bicepsData_norm)
    bn_top(:,:,i) = bicepsData_norm(i).E_top;
    bn_bottom(:,:,i) = bicepsData_norm(i).E_bottom;
    bn_whole(:,:,i) = bicepsData_norm(i).E_whole;
end
for i = 1:length(quadsData_dmd)
    qd_top(:,:,i) = quadsData_dmd(i).E_top;
    qd_bottom(:,:,i) = quadsData_dmd(i).E_bottom;
    qd_whole(:,:,i) = quadsData_dmd(i).E_whole;
    
    E_dmd(:,:,i) = quadsData_dmd(i).MW6;
end
for i = 1:length(quadsData_norm)
    qn_top(:,:,i) = quadsData_norm(i).E_top;
    qn_bottom(:,:,i) = quadsData_norm(i).E_bottom;
    qn_whole(:,:,i) = quadsData_norm(i).E_whole;
    
    E_norm(:,:,i) = quadsData_norm(i).MW6;
end


fprintf('\nDoes Youngs modulus differentiate DMD/normal? (t-test)\n')
% checkStats(.01,bd_top,bd_bottom,bd_whole, bn_top,bn_bottom,bn_whole, qd_top,qd_bottom,qd_whole, qn_top,qn_bottom,qn_whole);
% checkStats(.01,bd_top,bd_bottom,bd_whole(r1,r2,:), bn_top,bn_bottom,bn_whole(r1,r2,:), qd_top,qd_bottom,qd_whole(r1,r2,:), qn_top,qn_bottom,qn_whole(r1,r2,:));
% plotBox(bd_top(6,9,:),bn_top(6,9,:))
% title('Biceps, upper muscle, 1-10N')
% plotBox(bd_whole(r1,r2,:),bn_whole(r1,r2,:))
% title('Youngs modulus, biceps, whole muscle, 3-10N')
% plotBox(bd_bottom(r1,r2,:),bn_bottom(r1,r2,:))
% title('Biceps, lower muscle, E')
% plotBox(qd_top(r1,r2,:),qn_top(r1,r2,:))
% figure(1)
% title('Quads, upper muscle, E')
plotBox(qd_bottom(r1,r2,:),qn_bottom(r1,r2,:))
% figure(2)
title({'Youngs modulus of quadriceps';'Lower muscle'})
plotBox(qd_whole(r1,r2,:),qn_whole(r1,r2,:))
title({'Youngs modulus of quadriceps';'Whole muscle'})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes Young's modulus correlate with age?\n');
% fprintf('Biceps, upper muscle, DMD');
% corrVsAge(bd_top);
% fprintf('Biceps, lower muscle, DMD');
% corrVsAge(bd_bottom);
% fprintf('Biceps, whole muscle, DMD');
% corrVsAge(bd_whole);
% fprintf('Quads, upper muscle, DMD');
% corrVsAge(qd_top);
% fprintf('Quads, lower muscle, DMD');
corrVsAge(qd_bottom);
fprintf('Quads, whole muscle, DMD');
% corrVsAge(qd_whole);

% plotVsAge(bd_top(9,9,:),bn_top(9,9,:));
% title('Biceps, upper muscle, 1-10N');
% ylabel('Young's modulus'); xlabel('Age');
% plotVsAge(bd_bottom(3,3,:),bn_bottom(3,3,:));
% title('Biceps, lower muscle, 1-10N');
% ylabel('Young's modulus'); xlabel('Age');
% plotVsAge(bd_whole(r1,r2,:),bn_whole(r1,r2,:));
% title('Biceps, whole muscle');
% ylabel('Young's modulus'); xlabel('Age');
% plotVsAge(qd_top(r1,r2,:),qn_top(r1,r2,:));
% title('Quads, upper muscle, E(3-10) vs age');
% ylabel('Youngs modulus'); xlabel('Age');
plotVsAge(qd_bottom(r1,r2,:),qn_bottom(r1,r2,:));
title({'Youngs modulus of quadriceps over time';'Lower muscle'});
ylabel('Youngs modulus (kPa)'); xlabel('Age (days)');
% plotVsAge(qd_whole(r1,r2,:),qn_whole(r1,r2,:));
% title('Quads, whole muscle, E(3-10) vs age');
% ylabel('Youngs modulus'); xlabel('Age');

% plot 6MW vs age too

plotVsAge(E_dmd,E_norm);
title('6MWT vs age');
ylabel('6MW (m)'); xlabel('Age');
legend('Healthy','DMD')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n\n\n\n\n\n\n\n');
fprintf('\nDoes Youngs modulus correlate with 6MW?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVs6MW(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVs6MW(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVs6MW(qd_whole);

% plmailotVs6MW(qd_top(r1,r2,:));
% title('Quads, upper muscle, E(3-10) vs 6MW');
% ylabel('Youngs modulus'); xlabel('6MW');
plotVs6MW(qd_bottom(r1,r2,:));
title({'Youngs modulus of quadriceps';'Lower muscle'})
ylabel('Youngs modulus (kPa)'); xlabel('6MW (m)');
% plotVs6MW(qd_whole(r1,r2,:));
% title('Quads, whole muscle, E(3-10) vs 6MW');
% ylabel('Youngs modulus'); xlabel('6MW');
% axis([200 500 3.5*10^4 10^5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes Young's modulus correlate with 10m run?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVs10m(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVs10m(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVs10m(qd_whole);

% plotVs10m(qd_top(r1,r2,:));
% title('Quads, upper muscle, E(3-10) vs 10m run');
% ylabel('Youngs modulus'); xlabel('10m time');
% plotVs10m(qd_bottom(r1,r2,:));
% title('Quads, lower muscle, E(3-10) vs 10m run');
% ylabel('Youngs modulus'); xlabel('10m time');
% plotVs10m(qd_whole(r1,r2,:));
% title('Quads, whole muscle, E(3-10) vs 10m run');
% ylabel('Youngs modulus'); xlabel('10m time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes Young's modulus correlate with supine-to-stand?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsSupStand(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVsSupStand(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVsSupStand(qd_whole);

% plotVsSupStand(qd_top(r1,r2,:));
% title('Quads, upper muscle, E(3-10) vs supine-to-stand');
% ylabel('Youngs modulus'); xlabel('Sup/stand time');
% plotVsSupStand(qd_bottom(r1,r2,:));
% title('Quads, lower muscle, E(3-10) vs supine-to-stand');
% ylabel('Youngs modulus'); xlabel('Sup/stand time');
% plotVsSupStand(qd_whole(r1,r2,:));
% title('Quads, whole muscle, E(3-10) vs supine-to-stand');
% ylabel('Youngs modulus'); xlabel('Sup/stand time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes Young's modulus correlate with NSAA?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsNSAA(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVsNSAA(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVsNSAA(qd_whole);

% plotVsNSAA(qd_top(r1,r2,:));
% title('Quads, upper muscle, E(3-10) vs NSAA');
% ylabel('Youngs modulus'); xlabel('NSAA');
% plotVsNSAA(qd_bottom(r1,r2,:));
% title('Quads, lower muscle, E(3-10) vs NSAA');
% ylabel('Youngs modulus'); xlabel('NSAA');
% plotVsNSAA(qd_whole(r1,r2,:));
% title('Quads, whole muscle, E(3-10) vs NSAA');
% ylabel('Youngs modulus'); xlabel('NSAA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fprintf('\n\n\n\n\n\n\n\n\n');
% fprintf('\nDoes Young's modulus correlate with HHD (quads)?\n');
% fprintf('Quads, upper muscle, DMD');
% corrVsHHD_q(qd_top);
% fprintf('Quads, lower muscle, DMD');
% corrVsHHD_q(qd_bottom);
% fprintf('Quads, whole muscle, DMD');
% corrVsHHD_q(qd_whole);

% plotVsHHD_q(qd_top(r1,r2,:));
% title('Quads, upper muscle, E(3-10) vs HHD');
% ylabel('Youngs modulus'); xlabel('HHD, quads');
% plotVsHHD_q(qd_bottom(r1,r2,:));
% title('Quads, upper muscle, E(3-10) vs HHD');
% ylabel('Youngs modulus'); xlabel('HHD, quads');
% plotVsHHD_q(qd_whole(r1,r2,:));
% title('Quads, upper muscle, E(3-10) vs HHD');
% ylabel('Youngs modulus'); xlabel('HHD, quads');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end


function corrVs6MW(dmd)

global quadsData_dmd

scores = zeros(9,9);


for i = 1:length(quadsData_dmd)
    times_6MW(i) = quadsData_dmd(i).MW6;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with 6MW\n'); 


for i = 1:9
    for j = i:9
        imageBiomarkers = squeeze(dmd(i,j,:));            
        C = corrcoef(imageBiomarkers(~isnan(times_6MW)),times_6MW(~isnan(times_6MW)));
        scores(i,j) = C(1,2);   
    end
end

scores

fprintf('\nMax is %f', scores(find(abs(scores) == max(max(abs(scores))))) )
fprintf('\nMean is %f\n\n',mean(scores(scores~=0)));

end
function corrVs6MW_nonlinear(dmd)

global quadsData_dmd

score = [];

for i = 1:length(quadsData_dmd)
    times_6MW(i) = quadsData_dmd(i).MW6;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with 6MW\n'); 


imageBiomarkers = squeeze(dmd);
C = corrcoef(imageBiomarkers(~isnan(times_6MW)),times_6MW(~isnan(times_6MW)));
score = C(1,2)


end
function corrVs10m(dmd)

global quadsData_dmd

scores = zeros(9,9);


for i = 1:length(quadsData_dmd)
    times_10m(i) = quadsData_dmd(i).m10;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with 10m\n'); 


for i = 1:9
    for j = i:9
        imageBiomarkers = squeeze(dmd(i,j,:));            
        C = corrcoef(imageBiomarkers(~isnan(times_10m)),times_10m(~isnan(times_10m)));
        scores(i,j) = C(1,2);   
    end
end

scores

fprintf('\nMax is %f', scores(find(abs(scores) == max(max(abs(scores))))) )
fprintf('\nMean is %f\n\n',mean(scores(scores~=0)));

end
function corrVs10m_nonlinear(dmd)

global quadsData_dmd

score = [];

for i = 1:length(quadsData_dmd)
    times_10m(i) = quadsData_dmd(i).m10;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with 10m\n'); 


imageBiomarkers = squeeze(dmd);
C = corrcoef(imageBiomarkers(~isnan(times_10m)),times_10m(~isnan(times_10m)));
score = C(1,2)


end
function corrVsSupStand(dmd)

global quadsData_dmd

scores = zeros(9,9);


for i = 1:length(quadsData_dmd)
    times_sup(i) = quadsData_dmd(i).supStand;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with supine-to-stand\n');


for i = 1:9
    for j = i:9
        imageBiomarkers = squeeze(dmd(i,j,:));            
        C = corrcoef(imageBiomarkers(~isnan(times_sup)),times_sup(~isnan(times_sup)));
        scores(i,j) = C(1,2);   
    end
end

scores

fprintf('\nMax is %f', scores(find(abs(scores) == max(max(abs(scores))))) )
fprintf('\nMean is %f\n\n',mean(scores(scores~=0)));

end
function corrVsSupStand_nonlinear(dmd)

global quadsData_dmd

score = [];

for i = 1:length(quadsData_dmd)
    times_ss(i) = quadsData_dmd(i).supStand;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with supine-to-stand\n'); 


imageBiomarkers = squeeze(dmd);
C = corrcoef(imageBiomarkers(~isnan(times_ss)),times_ss(~isnan(times_ss)));
score = C(1,2)


end
function corrVsNSAA(dmd)

global quadsData_dmd

scores = zeros(9,9);


for i = 1:length(quadsData_dmd)
    nsas(i) = quadsData_dmd(i).NSAA;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with NSAA\n');


for i = 1:9
    for j = i:9
        imageBiomarkers = squeeze(dmd(i,j,:));            
        C = corrcoef(imageBiomarkers(~isnan(nsas)),nsas(~isnan(nsas)));
        scores(i,j) = C(1,2);   
    end
end

scores

fprintf('\nMax is %f', scores(find(abs(scores) == max(max(abs(scores))))) )
fprintf('\nMean is %f\n\n',mean(scores(scores~=0)));

end
function corrVsNSAA_nonlinear(dmd)

global quadsData_dmd

score = [];

for i = 1:length(quadsData_dmd)
    nsas(i) = quadsData_dmd(i).NSAA;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with NSAA\n'); 


imageBiomarkers = squeeze(dmd);
C = corrcoef(imageBiomarkers(~isnan(nsas)),nsas(~isnan(nsas)));
score = C(1,2)


end
function corrVsHHD_q(dmd)

global quadsData_dmd

scores = zeros(9,9);


for i = 1:length(quadsData_dmd)
    hhd_qs(i) = quadsData_dmd(i).HHD_quads;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with HHD\n');


for i = 1:9
    for j = i:9
        imageBiomarkers = squeeze(dmd(i,j,:));            
        C = corrcoef(imageBiomarkers(~isnan(hhd_qs)),hhd_qs(~isnan(hhd_qs)));
        scores(i,j) = C(1,2);   
    end
end

scores

fprintf('\nMax is %f', scores(find(abs(scores) == max(max(abs(scores))))) )
fprintf('\nMean is %f\n\n',mean(scores(scores~=0)));

end
function corrVsHHD_q_nonlinear(dmd)

global quadsData_dmd

score = [];

for i = 1:length(quadsData_dmd)
    hhd_qs(i) = quadsData_dmd(i).HHD_quads;
end

fprintf('---------------------------------------------');
fprintf('\nQuads: Biomarker correlations with HHD\n'); 


imageBiomarkers = squeeze(dmd);
C = corrcoef(imageBiomarkers(~isnan(hhd_qs)),hhd_qs(~isnan(hhd_qs)));
score = C(1,2)


end
function corrVsHHD_b(dmd)

global bicepsData_dmd

scores = zeros(9,9);


for i = 1:length(bicepsData_dmd)
    hhd_bs(i) = bicepsData_dmd(i).HHD_biceps;
end

fprintf('---------------------------------------------');
fprintf('\nBiceps: Biomarker correlations with HHD\n');


for i = 1:9
    for j = i:9
        imageBiomarkers = squeeze(dmd(i,j,:));            
        C = corrcoef(imageBiomarkers(~isnan(hhd_bs)),hhd_bs(~isnan(hhd_bs)));
        scores(i,j) = C(1,2);   
    end
end

scores

fprintf('\nMax is %f', scores(find(abs(scores) == max(max(abs(scores))))) )
fprintf('\nMean is %f\n\n',mean(scores(scores~=0)));

end
function corrVsHHD_b_nonlinear(dmd)

global bicepsData_dmd

score = [];

for i = 1:length(bicepsData_dmd)
    hhd_bs(i) = bicepsData_dmd(i).HHD_biceps;
end

fprintf('---------------------------------------------');
fprintf('\nBiceps: Biomarker correlations with HHD\n'); 


imageBiomarkers = squeeze(dmd);
C = corrcoef(imageBiomarkers(~isnan(hhd_bs)),hhd_bs(~isnan(hhd_bs)));
score = C(1,2)


end
function corrVsIsoPush(dmd)

global bicepsData_dmd

scores = zeros(9,9);


for i = 1:length(bicepsData_dmd)
    push(i) = bicepsData_dmd(i).isoPush;
end

fprintf('---------------------------------------------');
fprintf('\nBiceps: Biomarker correlations with isometric pushup\n');


for i = 1:9
    for j = i:9
        imageBiomarkers = squeeze(dmd(i,j,:));            
        C = corrcoef(imageBiomarkers(~isnan(push)),push(~isnan(push)));
        scores(i,j) = C(1,2);   
    end
end

scores

fprintf('\nMax is %f', scores(find(abs(scores) == max(max(abs(scores))))) )
fprintf('\nMean is %f\n\n',mean(scores(scores~=0)));

end
function corrVsIsoPush_nonlinear(dmd)

global bicepsData_dmd

score = [];

for i = 1:length(bicepsData_dmd)
    push(i) = bicepsData_dmd(i).isoPush;
end

fprintf('---------------------------------------------');
fprintf('\nBiceps: Biomarker correlations with isometric pushup\n'); 


imageBiomarkers = squeeze(dmd);
C = corrcoef(imageBiomarkers(~isnan(push)),push(~isnan(push)));
score = C(1,2)


end
function corrVsAge(dmd)

global bicepsData_dmd
global quadsData_dmd

scores = zeros(9,9);

if size(dmd,3) == 21 %biceps
    for i = 1:length(bicepsData_dmd)
        ages(i) = bicepsData_dmd(i).age;
    end
    
    fprintf('---------------------------------------------');
    fprintf('\nBiceps: Biomarker correlations with age\n');    

elseif size(dmd,3) == 33 %quads
    for i = 1:length(quadsData_dmd)
        ages(i) = quadsData_dmd(i).age;
    end
    
    fprintf('---------------------------------------------');
    fprintf('\nQuads: Biomarker correlations with age\n'); 
end


for i = 1:9
    for j = i:9
        imageBiomarkers = squeeze(dmd(i,j,:));            
        C = corrcoef(imageBiomarkers,ages);
        scores(i,j) = C(1,2);   
    end
end

scores

fprintf('\nMax is %f', scores(find(abs(scores) == max(max(abs(scores))))) )
fprintf('\nMean is %f\n\n',mean(scores(scores~=0)));


end
function corrVsAge_nonlinear(dmd)

global bicepsData_dmd
global quadsData_dmd

score = [];

if length(dmd) == 21 %biceps
    for i = 1:length(bicepsData_dmd)
        ages(i) = bicepsData_dmd(i).age;
    end
    
    fprintf('---------------------------------------------');
    fprintf('\nBiceps: Biomarker correlations with age\n');    

elseif length(dmd) == 31 %quads
    for i = 1:length(quadsData_dmd)
        ages(i) = quadsData_dmd(i).age;
    end
    
    fprintf('---------------------------------------------');
    fprintf('\nQuads: Biomarker correlations with age\n'); 
end

imageBiomarkers = squeeze(dmd);
C = corrcoef(imageBiomarkers,ages);
score = C(1,2);   

score


end


function plotVs6MW(dmd)

dmd = squeeze(dmd);
global quadsData_dmd

figure('Position',[600,600,400,300]); 
hold on;

toLeg = {};

for i = 1:length(quadsData_dmd)
    if ~isnan(quadsData_dmd(i).MW6)
%         if isequal(quadsData_dmd(i).patient,'01003')
%             plot(quadsData_dmd(i).MW6,dmd(i),'bo')
%             toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - ineligible');
        if isequal(quadsData_dmd(i).patient,'01012')
            plot(quadsData_dmd(i).MW6,dmd(i),'ko')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - scoliosis');
        elseif isequal(quadsData_dmd(i).patient,'01002')
            plot(quadsData_dmd(i).MW6,dmd(i),'bd')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - nonambulatory');
        elseif isequal(quadsData_dmd(i).patient,'01004')
            plot(quadsData_dmd(i).MW6,dmd(i),'bd')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - nonambulatory');
        elseif isequal(quadsData_dmd(i).patient,'01022')
            plot(quadsData_dmd(i).MW6,dmd(i),'bd')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - nonambulatory');
        elseif isequal(quadsData_dmd(i).patient,'01039')
            plot(quadsData_dmd(i).MW6,dmd(i),'bd')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - nonambulatory');
        else
            plot(quadsData_dmd(i).MW6,dmd(i),'ro')
            toLeg{end+1} = quadsData_dmd(i).patient;
        end
    end
end

% legend(toLeg)
set(gca,'fontsize', 12);

end
function plotVs10m(dmd)

dmd = squeeze(dmd);
global quadsData_dmd

figure; hold on;

toLeg = {};

for i = 1:length(quadsData_dmd)
    if ~isnan(quadsData_dmd(i).m10)
%         if isequal(quadsData_dmd(i).patient,'01003')
%             plot(quadsData_dmd(i).m10,dmd(i),'bo')
%             toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - ineligible');
        if isequal(quadsData_dmd(i).patient,'01012')
            plot(quadsData_dmd(i).m10,dmd(i),'ko')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - scoliosis');
        elseif isequal(quadsData_dmd(i).patient,'01002')
            plot(quadsData_dmd(i).m10,dmd(i),'bd')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - nonambulatory');
        elseif isequal(quadsData_dmd(i).patient,'01004')
            plot(quadsData_dmd(i).m10,dmd(i),'bd')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - nonambulatory');
        elseif isequal(quadsData_dmd(i).patient,'01022')
            plot(quadsData_dmd(i).m10,dmd(i),'bd')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - nonambulatory');
        elseif isequal(quadsData_dmd(i).patient,'01039')
            plot(quadsData_dmd(i).m10,dmd(i),'bd')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - nonambulatory');
        else
            plot(quadsData_dmd(i).m10,dmd(i),'ro')
            toLeg{end+1} = quadsData_dmd(i).patient;
        end
    end
end

legend(toLeg)
set(gca,'fontsize', 12);

end
function plotVsIsoPush(dmd)

dmd = squeeze(dmd);
global bicepsData_dmd

figure; hold on;

toLeg = {};

for i = 1:length(bicepsData_dmd)
    if ~isnan(bicepsData_dmd(i).isoPush)
%         if isequal(bicepsData_dmd(i).patient,'01003')
%             plot(bicepsData_dmd(i).isoPush,dmd(i),'bo')
%             toLeg{end+1} = strcat(bicepsData_dmd(i).patient,' - ineligible');
        if isequal(bicepsData_dmd(i).patient,'01012')
            plot(bicepsData_dmd(i).isoPush,dmd(i),'ko')
            toLeg{end+1} = strcat(bicepsData_dmd(i).patient,' - scoliosis');
        else
            plot(bicepsData_dmd(i).isoPush,dmd(i),'ro')
            toLeg{end+1} = bicepsData_dmd(i).patient;
        end
    end
end

% legend(toLeg)
set(gca,'fontsize', 12);

end
function plotVsSupStand(dmd)

dmd = squeeze(dmd);
global quadsData_dmd

figure; hold on;

toLeg = {};

for i = 1:length(quadsData_dmd)
    if ~isnan(quadsData_dmd(i).supStand)
%         if isequal(quadsData_dmd(i).patient,'01003')
%             plot(quadsData_dmd(i).supStand,dmd(i),'bo')
%             toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - ineligible');
        if isequal(quadsData_dmd(i).patient,'01012')
            plot(quadsData_dmd(i).supStand,dmd(i),'ko')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - scoliosis');
        else
            plot(quadsData_dmd(i).supStand,dmd(i),'ro')
            toLeg{end+1} = quadsData_dmd(i).patient;
        end
    end
end

% legend(toLeg)
set(gca,'fontsize', 12);

end
function plotVsNSAA(dmd)

dmd = squeeze(dmd);
global quadsData_dmd

figure; hold on;

toLeg = {};

for i = 1:length(quadsData_dmd)
    if ~isnan(quadsData_dmd(i).NSAA)
%         if isequal(quadsData_dmd(i).patient,'01003')
%             plot(quadsData_dmd(i).NSAA,dmd(i),'bo')
%             toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - ineligible');
        if isequal(quadsData_dmd(i).patient,'01012')
            plot(quadsData_dmd(i).NSAA,dmd(i),'ko')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - scoliosis');
        else
            plot(quadsData_dmd(i).NSAA,dmd(i),'ro')
            toLeg{end+1} = quadsData_dmd(i).patient;
        end
    end
end

% legend(toLeg)
set(gca,'fontsize', 12);

end
function plotVsHHD_b(dmd)

dmd = squeeze(dmd);
global bicepsData_dmd

figure; hold on;

toLeg = {};

for i = 1:length(bicepsData_dmd)
    if ~isnan(bicepsData_dmd(i).HHD_biceps)
%         if isequal(bicepsData_dmd(i).patient,'01003')
%             plot(bicepsData_dmd(i).HHD_biceps,dmd(i),'bo')
%             toLeg{end+1} = strcat(bicepsData_dmd(i).patient,' - ineligible');
        if isequal(bicepsData_dmd(i).patient,'01012')
            plot(bicepsData_dmd(i).HHD_biceps,dmd(i),'ko')
            toLeg{end+1} = strcat(bicepsData_dmd(i).patient,' - scoliosis');
        else
            plot(bicepsData_dmd(i).HHD_biceps,dmd(i),'ro')
            toLeg{end+1} = bicepsData_dmd(i).patient;
        end
    end
end

% legend(toLeg)
set(gca,'fontsize', 12);

end
function plotVsHHD_q(dmd)

dmd = squeeze(dmd);
global quadsData_dmd

figure; hold on;

toLeg = {};

for i = 1:length(quadsData_dmd)
    if ~isnan(quadsData_dmd(i).HHD_quads)
%         if isequal(quadsData_dmd(i).patient,'01003')
%             plot(quadsData_dmd(i).HHD_quads,dmd(i),'bo')
%             toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - ineligible');
        if isequal(quadsData_dmd(i).patient,'01012')
            plot(quadsData_dmd(i).HHD_quads,dmd(i),'ko')
            toLeg{end+1} = strcat(quadsData_dmd(i).patient,' - scoliosis');
        else
            plot(quadsData_dmd(i).HHD_quads,dmd(i),'ro')
            toLeg{end+1} = quadsData_dmd(i).patient;
        end
    end
end

% legend(toLeg)
set(gca,'fontsize', 12);

end
function plotVsAge(dmd,norm)

% Will need to update this to reflect accurate pop size

global bicepsData_dmd
global bicepsData_norm
global quadsData_dmd
global quadsData_norm


if length(norm) == 17

    for i = 1:length(bicepsData_dmd)
        ages_d(i) = bicepsData_dmd(i).age;
    end

    for i = 1:length(bicepsData_norm)
        ages_n(i) = bicepsData_norm(i).age;
    end
      
elseif length(norm) == 28

    for i = 1:length(quadsData_dmd)
        ages_d(i) = quadsData_dmd(i).age;
    end

    for i = 1:length(quadsData_norm)
        ages_n(i) = quadsData_norm(i).age;
    end

end


figure('Position',[600,600,400,300]); hold on;
plot(ages_n,squeeze(norm),'kx')
plot(ages_d,squeeze(dmd),'ro')

legend('Healthy','DMD')


set(gca,'fontsize', 12);

end
function plotBox(dmd,norm)


types = {};
for i = 1:28
    types{i} = 'Healthy';
end
for i = 29:61
    types{i} = 'DMD';
end
figure('Position',[600,600,400,300])
boxplot([transpose(squeeze(norm)) transpose(squeeze(dmd))],types);
set(gca,'fontsize', 12);
ylabel('Avg. Youngs modulus (kPa)')


end


function checkStats(alpha,bd_top,bd_bottom,bd_whole, bn_top,bn_bottom,bn_whole, qd_top,qd_bottom,qd_whole, qn_top,qn_bottom,qn_whole)

% [t_lin_b_top, p_lin_b_top] = ttest2(bd_top,bn_top,'Alpha',alpha,'dim',3,'Vartype','unequal');
% fprintf('---------------------------------------------');
% fprintf('\nBiceps, upper muscle, alpha = %f\n',alpha);
% t_lin_b_top
% p_lin_b_top
% 
% [t_lin_b_bottom, p_lin_b_bottom] = ttest2(bd_bottom,bn_bottom,'Alpha',alpha,'dim',3,'Vartype','unequal');
% fprintf('---------------------------------------------');
% fprintf('\nBiceps, lower muscle, alpha = %f\n',alpha);
% t_lin_b_bottom
% p_lin_b_bottom

m_bd = mean(bd_whole(bd_whole~=0));
s_bd = std(bd_whole(bd_whole~=0));
m_bn = mean(bn_whole(bn_whole~=0));
s_bn = std(bn_whole(bn_whole~=0));
[t_lin_b_whole, p_lin_b_whole] = ttest2(bd_whole,bn_whole,'Alpha',alpha,'dim',3,'Vartype','unequal');
fprintf('---------------------------------------------');
fprintf('\nBiceps, whole muscle, DMD    : mean = %f, stdev = %f',m_bd,s_bd);
fprintf('\nBiceps, whole muscle, normal : mean = %f, stdev = %f',m_bn,s_bn);
fprintf('\nFor alpha = %f\n',alpha);
fprintf('t = %d, p = %f\n',t_lin_b_whole,p_lin_b_whole)


% [t_lin_q_top, p_lin_q_top] = ttest2(qd_top,qn_top,'Alpha',alpha,'dim',3,'Vartype','unequal');
% fprintf('--------------------------------------------');
% fprintf('\nQuads, upper muscle, alpha = %f\n',alpha);
% t_lin_q_top
% p_lin_q_top
% 
% [t_lin_q_bottom, p_lin_q_bottom] = ttest2(qd_bottom,qn_bottom,'Alpha',alpha,'dim',3,'Vartype','unequal');
% fprintf('--------------------------------------------');
% fprintf('\nQuads, lower muscle, alpha = %f\n',alpha);
% t_lin_q_bottom
% p_lin_q_bottom

m_qd = mean(qd_whole(qd_whole~=0));
s_qd = std(qd_whole(qd_whole~=0));
m_qn = mean(qn_whole(qn_whole~=0));
s_qn = std(qn_whole(qn_whole~=0));
[t_lin_q_whole, p_lin_q_whole] = ttest2(qd_whole,qn_whole,'Alpha',alpha,'dim',3,'Vartype','unequal');
fprintf('--------------------------------------------');
fprintf('\nQuads, whole muscle, DMD    : mean = %f, stdev = %f',m_qd,s_qd);
fprintf('\nQuads, whole muscle, normal : mean = %f, stdev = %f',m_qn,s_qn);
fprintf('\nFor alpha = %f\n',alpha);
fprintf('t = %d, p = %f\n',t_lin_q_whole,p_lin_q_whole)


end


function [a,b] = getNonLinear(x,f)

x = x/x(1);
k = x./f;
p = polyfit(log(f),log(k),1);
b = p(1);
a = exp(p(2));

end
function mb = doHollomons(x,f)

A = .013*.047; % Could make bone size
stress = (f-f(1)) /A;
strain = (x(1)-x) /x(1);

stress = stress(find(strain));
strain = strain(find(strain));

mb = polyfit(log(strain(1:end)),log(stress(1:end)),1);
mb = real(mb);



end
function k = getSpringConstant(x,f)

x = x*.05/411; %50mm height = 411 pixels

for i = 1:9 % difference in N between points
    for j = 2:10 % larger force of the two
        if j-i >= 1
            
            if x(j) == x(j-i)
                k(i,j-1) = 3;
            else
                k(i,j-1) = -(f(j)-f(j-i)) / (x(j)-x(j-i));
            end
        
        end
    end
end


end
function E = getYoungMod(x,f)

A = .013*.047; % Could make bone size

for i = 1:9 % difference in N between points
    for j = 2:10 % larger force of the two
        if j-i >= 1

%             str(i,j-1) = (x(j-i)-x(j)) / (x(j-i) * (f(j)-f(j-i)) );
            E(i,j-1) = (x(j-i) * (f(j)-f(j-i)) ) / (A* (x(j-i)-x(j)) );
        
        end
    end
    
end

end
function data = getFunc

filename = '/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016_output/QED_Functional.csv';
delimiter = ';';
formatSpec = '%f%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'headerLines', 1, 'ReturnOnError', false);
fclose(fileID);

data = struct;
for i = 1:length(dataArray{1})
    data(i).patient = dataArray{1}(i);
    data(i).visit = dataArray{2}{i};
    data(i).HHD_biceps = dataArray{4}(i);
    data(i).HHD_quads = dataArray{8}(i);
    data(i).isoPush = dataArray{10}(i);
    data(i).NSAA = dataArray{11}(i);
    data(i).supStand = dataArray{14}(i);
    data(i).m10 = dataArray{15}(i);
    data(i).MW6 = dataArray{17}(i);    
end

clear filename

end


function createEntry(filename,funcs,k_top,k_bottom,k_whole,E_top,E_bottom,E_whole,a_top,b_top,a_bottom,b_bottom,a_whole,b_whole,mb_top,mb_bottom,mb_whole)

global bicepsData_dmd
global bicepsData_norm
global quadsData_dmd
global quadsData_norm

load('/home/sisir/code/data/DMD_Complete_20Aug2015/directory_18Oct2015.mat');

ID1 = strsplit(filename,'.');
ID1 = ID1{1};
ID2 = strsplit(ID1,'_');

for i = 1:length(s)
    if isequal(s(i).patient, ID2{1})
        if isequal(s(i).visit, ID2{3})
            age = datenum(s(i).DoV)-datenum(s(i).DoB);
        end
    end
end

walk_6min = NaN;
walk_10m = NaN;
pushup = NaN;
HHD_bcps = NaN;
HHD_qds = NaN;
standup = NaN;
nsa = NaN;
for j = 1:length(funcs)
    if isequal(funcs(j).patient, str2num(ID2{1}))
        if isequal(funcs(j).visit, ID2{3})
            walk_10m = funcs(j).m10;
            walk_6min = funcs(j).MW6;
            pushup = funcs(j).isoPush;
            HHD_bcps = funcs(j).HHD_biceps;
            HHD_qds = funcs(j).HHD_quads;
            standup = funcs(j).supStand;
            nsa = funcs(j).NSAA;
        end
    end
    
end



switch ID1([2 7])
    case '1b'
        bicepsData_dmd(end+1).patient = ID2{1};
        bicepsData_dmd(end).visit = ID2{3};
        bicepsData_dmd(end).muscle = ID2{2};
        bicepsData_dmd(end).age = age;
        bicepsData_dmd(end).isoPush = pushup;
        bicepsData_dmd(end).HHD_biceps = HHD_bcps;
        bicepsData_dmd(end).k_top = k_top;
        bicepsData_dmd(end).k_bottom = k_bottom;
        bicepsData_dmd(end).k_whole = k_whole;
        bicepsData_dmd(end).E_top = E_top;
        bicepsData_dmd(end).E_bottom = E_bottom;
        bicepsData_dmd(end).E_whole = E_whole;
        bicepsData_dmd(end).a_top = a_top;
        bicepsData_dmd(end).b_top = b_top;
        bicepsData_dmd(end).a_bottom = a_bottom;
        bicepsData_dmd(end).b_bottom = b_bottom;
        bicepsData_dmd(end).a_whole = a_whole;
        bicepsData_dmd(end).b_whole = b_whole;
        bicepsData_dmd(end).n_top = mb_top(1);
        bicepsData_dmd(end).K_top = exp(mb_top(2));
        bicepsData_dmd(end).n_bottom = mb_bottom(1);
        bicepsData_dmd(end).K_bottom = exp(mb_bottom(2));
        bicepsData_dmd(end).n_whole = mb_whole(1);
        bicepsData_dmd(end).K_whole = exp(mb_whole(2));
        
    case '2b'
        bicepsData_norm(end+1).patient = ID2{1};
        bicepsData_norm(end).visit = ID2{3};
        bicepsData_norm(end).muscle = ID2{2};
        bicepsData_norm(end).age = age;
        bicepsData_norm(end).isoPush = pushup;
        bicepsData_norm(end).HHD_biceps = HHD_bcps;
        bicepsData_norm(end).k_top = k_top;
        bicepsData_norm(end).k_bottom = k_bottom;
        bicepsData_norm(end).k_whole = k_whole;
        bicepsData_norm(end).E_top = E_top;
        bicepsData_norm(end).E_bottom = E_bottom;
        bicepsData_norm(end).E_whole = E_whole;
        bicepsData_norm(end).a_top = a_top;
        bicepsData_norm(end).b_top = b_top;
        bicepsData_norm(end).a_bottom = a_bottom;
        bicepsData_norm(end).b_bottom = b_bottom;
        bicepsData_norm(end).a_whole = a_whole;
        bicepsData_norm(end).b_whole = b_whole;
        bicepsData_norm(end).n_top = mb_top(1);
        bicepsData_norm(end).K_top = exp(mb_top(2));
        bicepsData_norm(end).n_bottom = mb_bottom(1);
        bicepsData_norm(end).K_bottom = exp(mb_bottom(2));
        bicepsData_norm(end).n_whole = mb_whole(1);
        bicepsData_norm(end).K_whole = exp(mb_whole(2));
        
    case '1q'
        quadsData_dmd(end+1).patient = ID2{1};
        quadsData_dmd(end).visit = ID2{3};
        quadsData_dmd(end).muscle = ID2{2};
        quadsData_dmd(end).age = age;
        quadsData_dmd(end).MW6 = walk_6min;
        quadsData_dmd(end).m10 = walk_10m;
        quadsData_dmd(end).supStand = standup;
        quadsData_dmd(end).NSAA = nsa;
        quadsData_dmd(end).HHD_quads = HHD_qds;
        quadsData_dmd(end).k_top = k_top;
        quadsData_dmd(end).k_bottom = k_bottom;
        quadsData_dmd(end).k_whole = k_whole;
        quadsData_dmd(end).E_top = E_top;
        quadsData_dmd(end).E_bottom = E_bottom;
        quadsData_dmd(end).E_whole = E_whole;
        quadsData_dmd(end).a_top = a_top;
        quadsData_dmd(end).b_top = b_top;
        quadsData_dmd(end).a_bottom = a_bottom;
        quadsData_dmd(end).b_bottom = b_bottom;
        quadsData_dmd(end).a_whole = a_whole;
        quadsData_dmd(end).b_whole = b_whole;
        quadsData_dmd(end).n_top = mb_top(1);
        quadsData_dmd(end).K_top = exp(mb_top(2));
        quadsData_dmd(end).n_bottom = mb_bottom(1);
        quadsData_dmd(end).K_bottom = exp(mb_bottom(2));
        quadsData_dmd(end).n_whole = mb_whole(1);
        quadsData_dmd(end).K_whole = exp(mb_whole(2));
        
    case '2q'
        quadsData_norm(end+1).patient = ID2{1};
        quadsData_norm(end).visit = ID2{3};
        quadsData_norm(end).muscle = ID2{2};
        quadsData_norm(end).age = age;
        quadsData_norm(end).MW6 = walk_6min;
        quadsData_norm(end).m10 = walk_10m;
        quadsData_norm(end).supStand = standup;
        quadsData_norm(end).NSAA = nsa;
        quadsData_norm(end).HHD_quads = HHD_qds;
        quadsData_norm(end).k_top = k_top;
        quadsData_norm(end).k_bottom = k_bottom;
        quadsData_norm(end).k_whole = k_whole;
        quadsData_norm(end).E_top = E_top;
        quadsData_norm(end).E_bottom = E_bottom;
        quadsData_norm(end).E_whole = E_whole;
        quadsData_norm(end).a_top = a_top;
        quadsData_norm(end).b_top = b_top;
        quadsData_norm(end).a_bottom = a_bottom;
        quadsData_norm(end).b_bottom = b_bottom;
        quadsData_norm(end).a_whole = a_whole;
        quadsData_norm(end).b_whole = b_whole;
        quadsData_norm(end).n_top = mb_top(1);
        quadsData_norm(end).K_top = exp(mb_top(2));
        quadsData_norm(end).n_bottom = mb_bottom(1);
        quadsData_norm(end).K_bottom = exp(mb_bottom(2));
        quadsData_norm(end).n_whole = mb_whole(1);
        quadsData_norm(end).K_whole = exp(mb_whole(2));
        
end



end


function [topMus,botMus,wholMus,frc] = getMuscleThick(filename)

load(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016/',filename))

sweepForces = sweepForces(1:length(subQ),2);

for i = 1:10
    holdme = find(abs(sweepForces-i) == min(abs(sweepForces-i)));
    idx(i) = holdme(1);
end

frc = sweepForces(idx);

topMus = interMus(idx,1) - subQ(idx,1);
botMus = boneTop(idx,1) - interMus(idx,1);
wholMus = boneTop(idx,1) - subQ(idx,1);

end
function [topMus,botMus,wholMus,frc] = getMuscleThick_complete(filename)

load(strcat('/home/sisir/code/data/DMD_Complete_20Aug2015/Compression_Jan2016/',filename))

frc = sweepForces(1:length(subQ),2);

topMus = interMus(:,1) - subQ(:,1);
botMus = boneTop(:,1) - interMus(:,1);
wholMus = boneTop(:,1) - subQ(:,1);

end
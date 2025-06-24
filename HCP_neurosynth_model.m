function HCP_neurosynth_kjhnotes_aparc(task,yourpath,yourhcppath)
yourpath='/fs/ess/PAS1305/ajuna';
task='reading';
yourhcppath='/fs/ess/PDE0021';
load([yourpath '/HCP_CIFTI_INDS.mat'])
  % these are the indicies for the cortex & subortical regions


load([ yourpath '/APARCconn_meanconn.mat' ]) %%
% this is every voxel to the mean signal of another region, already averaged across all individuals.

% load volume -% this is the 2mm volume, this is for fitting subcortical
% stuff
mri = MRIread([yourpath '/inputs_query_new/' task '_2.nii.gz']);

% initialize variables to iterate on
Mdl = {};
predcorr = [];
R2 = [];
b = {};
beta2mean = [];
beta2max0 = [];
beta2min0 = [];
beta2extrem0 = [];
beta2max = [];
beta2min = [];
beta2extrem = [];
ssnames = {};

%% left hemi
%  these are vertices on the left surface.
bBrain = gifti([ yourhcppath '/100307/MNINonLinear/fsaverage_LR32k/100307.L.atlasroi.32k_fs_LR.shape.gii']);
bBrain = logical(bBrain.cdata);
% use this to subset our functional data (cog domain)

% load the functional task data (this is the neuroquery map in the same space as the HCP data)
func = gifti([yourpath '/surf_query/' task '.L.func.gii']); % edit the path
func = func.cdata(bBrain); % include cortex not the median wall
    % func and conn should be the same size!
conn = c(indsStart(1):indsEnd(1),startsWith(regions,'L'));

% initialize the variables
meanfunc = [];
maxfunc = [];
minfunc = [];
extremfunc = [];

for i=1:35 %left hemi
      if i==4 || i==39; continue; end % skip region 4 because we dont have corpus callosum data
    tmp = func(atlas(atlas>=1 & atlas<=35)==i); %left hemi
    meanfunc(end+1) = mean(tmp(tmp~=0));
    maxfunc(end+1) = max(tmp);
    minfunc(end+1) = min(tmp);
    [extremfunc(end+1), ind] = max(abs(tmp));
    extremfunc(end) = extremfunc(end)*sign(tmp(ind));
end


for i=1:length(structure)
    if startsWith(structure{i},'L') && ~endsWith(structure{i},'CORTEX')
        tmp = mri.vol(volinds{i});
        meanfunc(end+1) = mean(tmp(tmp~=0));
        maxfunc(end+1) = max(tmp);
        minfunc(end+1) = min(tmp);
        [extremfunc(end+1), ind] = max(abs(tmp));
        extremfunc(end) = extremfunc(end)*sign(tmp(ind));
    end
end

atlas_lh=atlas(atlas>=1 & atlas<=35);

% do model for LH cortial regions
for iSS = 1:35
    % do the prediction one region at a time
    if iSS==0 || iSS==4 || iSS==39; continue;end

    % get connectivity info for all voxels in region iSS
    connSS = conn(atlas_lh==iSS,:);
    connSS(:,4)=[]; % remove 4 because it is the corpus callosum which we dont have data for
    funcSS = func(atlas_lh==iSS,:);
    % do the hyperparameter optimization and fit the model
    Mdl{end+1} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('ShowPlots',false,'MaxObjectiveEvaluations',100));

    % uses the model to predict the activation in each voxel
    pred = predict(Mdl{end},zscore(connSS));
    % how well correlated is voxel prediction based on NeuroQuery Map activation
    predcorr(end+1) = corr(pred,funcSS);

    % calculate model fit
    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(end+1) = 1-(SS_err / SS_tot);

    ssnames{end+1} = ['L' num2str(iSS)];
end

%combine striatum into one unit
tmp = mri.vol([volinds{3};volinds{8};volinds{18}]);
bSS = tmp~=0; % this should never be the case
if sum(bSS)>=0
    funcSS = tmp(bSS);
    % combine the accumbens, caudate and putamen
    conn = c([indsStart(3):indsEnd(3) indsStart(8):indsEnd(8) indsStart(18):indsEnd(18)],startsWith(regions,'L'));
    connSS = conn(bSS,:);
    connSS(:,4)=[]; % remove the corpus callosum

    % do the hyperparameter optimization and fit the model
    Mdl{end+1} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('ShowPlots',false,'MaxObjectiveEvaluations',100));

    % predict activation based on connectivity
    pred = predict(Mdl{end},zscore(connSS));
    predcorr(end+1) = corr(pred,funcSS);
    % this isnt predicting with new data -- this is just getting residuals
    % by ftting the same data with same model
    % residuals tells us how good our fit is for a region - error

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(end+1) = 1-(SS_err / SS_tot); % compare variance of residuals to the variance in original data

    ssnames{end+1} = 'L_STRIATUM'; % dont really need that, its in the atlas
end

% do the other subcortical regions
for i = [16 20 5 14 10 12] % already did 3 8 18 for striatum; the L Cortex and brainstem are excluded
    tmp = mri.vol(volinds{i});
    bSS = tmp~=0;
    if sum(bSS)<0; continue; end

    funcSS = tmp(bSS);
    conn = c(indsStart(i):indsEnd(i),startsWith(regions,'L'));
    connSS = conn(bSS,:);
    connSS(:,4)=[];

    Mdl{end+1} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('ShowPlots',false,'MaxObjectiveEvaluations',100));

    pred = predict(Mdl{end},zscore(connSS));
    predcorr(end+1) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(end+1) = 1-(SS_err / SS_tot);


    b{end+1} = Mdl{end}.Beta;
    beta2mean(end+1) = corr(b{end},meanfunc','r','p');
    beta2max0(end+1) = corr(b{end},maxfunc','r','p');
    beta2min0(end+1) = corr(b{end},minfunc','r','p');
    beta2extrem0(end+1) = corr(b{end},extremfunc','r','p');
    tmp = maxfunc; tmp(tmp==0) = nan; beta2max(end+1) = corr(b{end},tmp','r','p');
    tmp = minfunc; tmp(tmp==0) = nan; beta2min(end+1) = corr(b{end},tmp','r','p');
    tmp = extremfunc; tmp(tmp==0) = nan; beta2extrem(end+1) = corr(b{end},tmp','r','p');

    ssnames{end+1} = structure{i};
end

%% right hemi
bBrain = gifti([yourhcppath '/100307/MNINonLinear/fsaverage_LR32k/100307.R.atlasroi.32k_fs_LR.shape.gii']);
bBrain = logical(bBrain.cdata);

func = gifti([yourpath '/surf_query/' task '.R.func.gii']);
func = func.cdata(bBrain);

conn = c(indsStart(2):indsEnd(2),startsWith(regions,'R'));

%
meanfunc = [];
maxfunc = [];
minfunc = [];
extremfunc = [];

for i=36:70 %right hemi
    if i==4 || i==39; continue; end
    tmp = func(atlas(atlas>35 & atlas<=70)==i); %right hemi
    meanfunc(end+1) = mean(tmp(tmp~=0));
    maxfunc(end+1) = max(tmp);
    minfunc(end+1) = min(tmp);
    [extremfunc(end+1), ind] = max(abs(tmp));
    extremfunc(end) = extremfunc(end)*sign(tmp(ind));
end

%subcortical
for i=1:length(structure)
    if startsWith(structure{i},'R') && ~endsWith(structure{i},'CORTEX')
        tmp = mri.vol(volinds{i});
        meanfunc(end+1) = mean(tmp(tmp~=0));
        maxfunc(end+1) = max(tmp);
        minfunc(end+1) = min(tmp);
        [extremfunc(end+1), ind] = max(abs(tmp));
        extremfunc(end) = extremfunc(end)*sign(tmp(ind));
    end
end

%
atlas_rh=atlas(atlas>35 & atlas<=70);
for iSS = 36:70
    if iSS==0 || iSS==4 || iSS==39; continue;end

    connSS = conn(atlas_rh==iSS,:);
    funcSS = func(atlas_rh==iSS,:);
    connSS(:,4)=[];

    Mdl{end+1} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('ShowPlots',false,'MaxObjectiveEvaluations',100));

    pred = predict(Mdl{end},zscore(connSS));
    predcorr(end+1) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(end+1) = 1-(SS_err / SS_tot);


    b{end+1} = Mdl{end}.Beta;
    beta2mean(end+1) = corr(b{end},meanfunc','r','p');
    beta2max0(end+1) = corr(b{end},maxfunc','r','p');
    beta2min0(end+1) = corr(b{end},minfunc','r','p');
    beta2extrem0(end+1) = corr(b{end},extremfunc','r','p');
    tmp = maxfunc; tmp(tmp==0) = nan; beta2max(end+1) = corr(b{end},tmp','r','p');
    tmp = minfunc; tmp(tmp==0) = nan; beta2min(end+1) = corr(b{end},tmp','r','p');
    tmp = extremfunc; tmp(tmp==0) = nan; beta2extrem(end+1) = corr(b{end},tmp','r','p');

    ssnames{end+1} = ['R' num2str(iSS)];
end

%combine striatum into one unit
tmp = mri.vol([volinds{4};volinds{9};volinds{19}]);
bSS = tmp~=0;
if sum(bSS)>=0
    funcSS = tmp(bSS);
    conn = c([indsStart(4):indsEnd(4) indsStart(9):indsEnd(9) indsStart(19):indsEnd(19)],startsWith(regions,'R'));

    connSS = conn(bSS,:);
    connSS(:,4)=[];

    Mdl{end+1} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('ShowPlots',false,'MaxObjectiveEvaluations',100));

    pred = predict(Mdl{end},zscore(connSS));
    predcorr(end+1) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(end+1) = 1-(SS_err / SS_tot);

    ssnames{end+1} = 'R_STRIATUM';
end

for i = [17 21 6 15 11 13]
    tmp = mri.vol(volinds{i});
    bSS = tmp~=0;
    if sum(bSS)<0; continue; end

    funcSS = tmp(bSS);
    conn = c(indsStart(i):indsEnd(i),startsWith(regions,'R'));
    connSS = conn(bSS,:);
    connSS(:,4)=[];

    Mdl{end+1} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',struct('ShowPlots',false,'MaxObjectiveEvaluations',100));

    pred = predict(Mdl{end},zscore(connSS));
    predcorr(end+1) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(end+1) = 1-(SS_err / SS_tot);

    ssnames{end+1} = structure{i};
    ssnames_true = [regions(1:3);regions(5:35);ssnames(35:41)';regions(36:38);regions(40:70);ssnames(76:82)'];
    % fix the names so that they are in the same order as the data
end

clear SS* bBrain bSS c* func* i iSS ind m* ss tmp z
if ~exist([yourpath '/models_query_aparc/'])
    unix(['mkdir ' yourpath '/models_query_aparc/'])
end
save([yourpath '/models_query_aparc/' task '_wpred.mat'])

function HCP_neurosynth_model_perm(task,numperm,yourhcppath,yourpath)
%% do permuataions on CF models

load([yourhcppath '/HCP_CIFTI_INDS.mat'])
    % load indices because we saved data in the same order it came in
    % vertices and subcortical stored together, this just tells us what is what

load([yourpath '/APARCconn_meanconn.mat']) %%
    % this is every voxel to the mean signal of another region.
    % already averaged across all individuals.

Model=load([yourpath '/models_query_aparc/' task '_wpred.mat']);
    % load the previous model information, so that we can use the parameters
    % for the permutations


mri = MRIread([yourpath '/inputs_query_new/' task '_2.nii.gz']);
    % load volume -% this is the 2mm volume, this is for fitting subcortical

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

%% Set up your data
    % load LH vertices, load functional data, subset connectivity by left cortex
    %then get mean, max and min function in each of the LH regions.
%  these are vertices on the left surface.
bBrain = gifti([yourhcppath '100307/MNINonLinear/fsaverage_LR32k/100307.L.atlasroi.32k_fs_LR.shape.gii']);
bBrain = logical(bBrain.cdata);
    % use this to subset our functional data (cog domain)


func = gifti([yourpath '/surf_query/' task '.L.func.gii']); % edit the path
func = func.cdata(bBrain); % include cortex not the median wall
% func and conn are the same size!

conn = c(indsStart(1):indsEnd(1),startsWith(regions,'L'));
    % subset the connectivity data by all the cortex stuff

% initialize information about func activation in each region LH
meanfunc = [];
maxfunc = [];
minfunc = [];
extremfunc = [];

% calc mean, max, min for each region
for i=1:35 %left hemi cortex
      if i==4 || i==39; continue; end % skip region 4 because we dont have corpus callosum data
    tmp = func(atlas(atlas>=1 & atlas<=35)==i); %left hemi
    meanfunc(end+1) = mean(tmp(tmp~=0));
    maxfunc(end+1) = max(tmp);
    minfunc(end+1) = min(tmp);
    [extremfunc(end+1), ind] = max(abs(tmp));
    extremfunc(end) = extremfunc(end)*sign(tmp(ind));
end


for i=1:length(structure) %subcortical
    if startsWith(structure{i},'L') && ~endsWith(structure{i},'CORTEX')
        tmp = mri.vol(volinds{i});
        meanfunc(end+1) = mean(tmp(tmp~=0));
        maxfunc(end+1) = max(tmp);
        minfunc(end+1) = min(tmp);
        [extremfunc(end+1), ind] = max(abs(tmp));
        extremfunc(end) = extremfunc(end)*sign(tmp(ind));
    end
end

%% do model for LH cortex
% subset the data for LH only
atlas_lh=atlas(atlas>=1 & atlas<=35);
for iSS = 1:35
    % do the prediction one region at a time
    if iSS==0 || iSS==4 || iSS==39; continue;end
    % b/c we skip 4, we have to look in the right part of the previously
    % run model by subtracting 1 for all numbers >4
    if iSS>4
        mSS=iSS-1;
    elseif iSS==4
        mSS=NaN;
    else
        mSS=iSS;
    end
    outtmp(iSS)=mSS; % just to check that I did it correctly


    % just give me the connectivity info for all voxels in region 1
    connSS = conn(atlas_lh==iSS,:);
    % here is all the voxels in region iSS, and their connectivity to all
    % 44 LH regions
    connSS(:,4)=[];

    % get the functional activation for all the voxels in region iSS
    funcSS_tmp = func(atlas_lh==iSS,:);
    % shuffle funcSS then each voxel connectivity is paired with new funcSS

    % randomly permute that funcSS variable
    for p=1:numperm
     permorder=randperm(size(funcSS_tmp,1));
     funcSS=funcSS_tmp(permorder');
     % run the model with the permuted data, using the hyperparameter
     % previously optimized
     Mdl{p,iSS} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','Epsilon',Model.Mdl{1,mSS}.Epsilon,'ResponseTransform',Model.Mdl{1,mSS}.ResponseTransform,'Learner',Model.Mdl{1,mSS}.Learner,'Bias',Model.Mdl{1,mSS}.Bias,'Lambda',Model.Mdl{1,mSS}.Lambda );

    % uses the model to predict the activation in each voxel
    pred = predict(Mdl{p,iSS},zscore(connSS));
    % how well correlated is my voxel prediction based on true activation
    predcorr(p,iSS) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(p,iSS) = 1-(SS_err / SS_tot);


    b{p,iSS} = Mdl{p,iSS}.Beta;
    beta2mean(p,iSS) = corr(b{p,iSS},meanfunc','r','p');
    beta2max0(p,iSS) = corr(b{p,iSS},maxfunc','r','p');
    beta2min0(p,iSS) = corr(b{p,iSS},minfunc','r','p');
    beta2extrem0(p,iSS) = corr(b{p,iSS},extremfunc','r','p');
    tmp = maxfunc; tmp(tmp==0) = nan; beta2max(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = minfunc; tmp(tmp==0) = nan; beta2min(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = extremfunc; tmp(tmp==0) = nan; beta2extrem(p,iSS) = corr(b{p,iSS},tmp','r','p');

    end
        ssnames{iSS} = ['L' num2str(iSS)];
end

%% remove 4 b/c we skipped it
Mdl(:,4)=[];
ssnames(4)=[];
predcorr(:,4)=[];
R2(:,4)=[];
beta2extrem(:,4)=[];
beta2max(:,4)=[];
beta2mean(:,4)=[];
beta2min(:,4)=[];
b(:,4)=[];

%% run model for LH Striatum
%combine striatum into one unit
tmp = mri.vol([volinds{3};volinds{8};volinds{18}]);
bSS = tmp~=0;
if sum(bSS)>=0
    funcSS_tmp = tmp(bSS);
    conn = c([indsStart(3):indsEnd(3) indsStart(8):indsEnd(8) indsStart(18):indsEnd(18)],startsWith(regions,'L'));
    connSS = conn(bSS,:);
    connSS(:,4)=[];

    iSS=size(Mdl,2)+1;

    for p=1:numperm
     permorder=randperm(size(funcSS_tmp,1));
     funcSS=funcSS_tmp(permorder');

    Mdl{p,iSS} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','Epsilon',Model.Mdl{1,iSS}.Epsilon,'ResponseTransform',Model.Mdl{1,iSS}.ResponseTransform,'Learner',Model.Mdl{1,iSS}.Learner,'Bias',Model.Mdl{1,iSS}.Bias,'Lambda',Model.Mdl{1,iSS}.Lambda );
    pred = predict(Mdl{p,iSS},zscore(connSS));
    predcorr(p,iSS) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(p,iSS) = 1-(SS_err / SS_tot);


    b{p,iSS} = Mdl{p,iSS}.Beta;
    beta2mean(p,iSS) = corr(b{p,iSS},meanfunc','r','p');
    beta2max0(p,iSS) = corr(b{p,iSS},maxfunc','r','p');
    beta2min0(p,iSS) = corr(b{p,iSS},minfunc','r','p');
    beta2extrem0(p,iSS) = corr(b{p,iSS},extremfunc','r','p');
    tmp = maxfunc; tmp(tmp==0) = nan; beta2max(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = minfunc; tmp(tmp==0) = nan; beta2min(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = extremfunc; tmp(tmp==0) = nan; beta2extrem(p,iSS) = corr(b{p,iSS},tmp','r','p');
    end

    ssnames{end+1} = 'L_STRIATUM'; % dont really need that, its in the atlas
end

%% Run the model for the LH subcortical

mdl_lh=size(Mdl,2);
lhsub=[16 20 5 14 10 12];
for i = 1:length(lhsub) % already did 3 8 18 for striatum; the L Cortex and brainstem are excluded
    tmp = mri.vol(volinds{lhsub(i)});
    bSS = tmp~=0;
    if sum(bSS)<0; continue; end
    iSS=mdl_lh+i;
    funcSS_tmp = tmp(bSS);
    conn = c(indsStart(lhsub(i)):indsEnd(lhsub(i)),startsWith(regions,'L'));
    connSS = conn(bSS,:);
    connSS(:,4)=[];

     for p=1:numperm
     permorder=randperm(size(funcSS_tmp,1));
     funcSS=funcSS_tmp(permorder');
    Mdl{p,iSS} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','Epsilon',Model.Mdl{1,iSS}.Epsilon,'ResponseTransform',Model.Mdl{1,iSS}.ResponseTransform,'Learner',Model.Mdl{1,iSS}.Learner,'Bias',Model.Mdl{1,iSS}.Bias,'Lambda',Model.Mdl{1,iSS}.Lambda );

    pred = predict(Mdl{p,iSS},zscore(connSS));
    predcorr(p,iSS) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(p,iSS) = 1-(SS_err / SS_tot);


    b{p,iSS} = Mdl{p,iSS}.Beta;
    beta2mean(p,iSS) = corr(b{p,iSS},meanfunc','r','p');
    beta2max0(p,iSS) = corr(b{p,iSS},maxfunc','r','p');
    beta2min0(p,iSS) = corr(b{p,iSS},minfunc','r','p');
    beta2extrem0(p,iSS) = corr(b{p,iSS},extremfunc','r','p');
    tmp = maxfunc; tmp(tmp==0) = nan; beta2max(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = minfunc; tmp(tmp==0) = nan; beta2min(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = extremfunc; tmp(tmp==0) = nan; beta2extrem(p,iSS) = corr(b{p,iSS},tmp','r','p');
    end
    size(Mdl)
    ssnames{iSS} = structure{lhsub(i)};
end


%% right hemi

%% prep the variables like we did for the LH
bBrain = gifti([yourhcppath '/100307/MNINonLinear/fsaverage_LR32k/100307.R.atlasroi.32k_fs_LR.shape.gii');
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

%% run the model for the RH cortex
atlas_rh=atlas(atlas>35 & atlas<=70);
rhval=[36:70];
lhmdl_all=size(Mdl,2);
for i = 1:length(rhval)
    irh=rhval(i);
    if irh==0 || irh==4 || irh==39; continue;end

    % have to set mSS
     if irh>39
        mSS=irh-1;
    elseif irh==39
        mSS=NaN;
    else
        mSS=irh;
     end


    connSS = conn(atlas_rh==irh,:);
    funcSS_tmp = func(atlas_rh==irh,:);
    connSS(:,4)=[];
    iSS=lhmdl_all+i;
     for p=1:numperm
     permorder=randperm(size(funcSS_tmp,1));
     funcSS=funcSS_tmp(permorder');
     Mdl{p,iSS} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','Epsilon',Model.Mdl{1,mSS}.Epsilon,'ResponseTransform',Model.Mdl{1,mSS}.ResponseTransform,'Learner',Model.Mdl{1,mSS}.Learner,'Bias',Model.Mdl{1,mSS}.Bias,'Lambda',Model.Mdl{1,mSS}.Lambda );
    pred = predict(Mdl{p,iSS},zscore(connSS));
    predcorr(p,iSS) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(p,iSS) = 1-(SS_err / SS_tot);


    b{p,iSS} = Mdl{p,iSS}.Beta;
    beta2mean(p,iSS) = corr(b{p,iSS},meanfunc','r','p');
    beta2max0(p,iSS) = corr(b{p,iSS},maxfunc','r','p');
    beta2min0(p,iSS) = corr(b{p,iSS},minfunc','r','p');
    beta2extrem0(p,iSS) = corr(b{p,iSS},extremfunc','r','p');
    tmp = maxfunc; tmp(tmp==0) = nan; beta2max(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = minfunc; tmp(tmp==0) = nan; beta2min(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = extremfunc; tmp(tmp==0) = nan; beta2extrem(p,iSS) = corr(b{p,iSS},tmp','r','p');
    end

    ssnames{end+1} = ['R' num2str(iSS)];
end

%% we dont have brainstem info, so remove it -- its all blanks anyway
Mdl(:,45)=[];
predcorr(:,45)=[];
R2(:,45)=[];
beta2extrem(:,45)=[];
beta2max(:,45)=[];
beta2mean(:,45)=[];
beta2min(:,45)=[];
b(:,45)=[];

%% RH striatum model
%combine striatum into one unit
tmp = mri.vol([volinds{4};volinds{9};volinds{19}]);
bSS = tmp~=0;
if sum(bSS)>=0
    funcSS_tmp = tmp(bSS);
     conn = c([indsStart(4):indsEnd(4) indsStart(9):indsEnd(9) indsStart(19):indsEnd(19)],startsWith(regions,'R'));

    connSS = conn(bSS,:);
    connSS(:,4)=[];

        iSS=size(Mdl,2)+1;
     for p=1:numperm
     permorder=randperm(size(funcSS_tmp,1));
     funcSS=funcSS_tmp(permorder');
    Mdl{p,iSS} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','Epsilon',Model.Mdl{1,iSS}.Epsilon,'ResponseTransform',Model.Mdl{1,iSS}.ResponseTransform,'Learner',Model.Mdl{1,iSS}.Learner,'Bias',Model.Mdl{1,iSS}.Bias,'Lambda',Model.Mdl{1,iSS}.Lambda );
    pred = predict(Mdl{p,iSS},zscore(connSS));
    predcorr(p,iSS) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(p,iSS) = 1-(SS_err / SS_tot);


    b{p,iSS} = Mdl{p,iSS}.Beta;
    beta2mean(p,iSS) = corr(b{p,iSS},meanfunc','r','p');
    beta2max0(p,iSS) = corr(b{p,iSS},maxfunc','r','p');
    beta2min0(p,iSS) = corr(b{p,iSS},minfunc','r','p');
    beta2extrem0(p,iSS) = corr(b{p,iSS},extremfunc','r','p');
    tmp = maxfunc; tmp(tmp==0) = nan; beta2max(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = minfunc; tmp(tmp==0) = nan; beta2min(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = extremfunc; tmp(tmp==0) = nan; beta2extrem(p,iSS) = corr(b{p,iSS},tmp','r','p');
    end


    ssnames{end+1} = 'R_STRIATUM';
end
%% RH subcortical
% rh subcortical indices
    rhsub=[17 21 6 15 11 13];
   mdlsize=size(Mdl,2);
for i = 1:length(rhsub)
    irh=rhsub(i);
    tmp = mri.vol(volinds{irh});
    bSS = tmp~=0;
    if sum(bSS)<0; continue; end

    funcSS_tmp = tmp(bSS);
    conn = c(indsStart(irh):indsEnd(irh),startsWith(regions,'R'));
    connSS = conn(bSS,:);
    connSS(:,4)=[];
    iSS=mdlsize+i;
   for p=1:numperm
     permorder=randperm(size(funcSS_tmp,1));
     funcSS=funcSS_tmp(permorder');
    Mdl{p,iSS} = fitrlinear(zscore(connSS),funcSS,'Regularization','ridge','Epsilon',Model.Mdl{1,iSS}.Epsilon,'ResponseTransform',Model.Mdl{1,iSS}.ResponseTransform,'Learner',Model.Mdl{1,iSS}.Learner,'Bias',Model.Mdl{1,iSS}.Bias,'Lambda',Model.Mdl{1,iSS}.Lambda );
    pred = predict(Mdl{p,iSS},zscore(connSS));
    predcorr(p,iSS) = corr(pred,funcSS);

    SS_tot = sum((funcSS - mean(funcSS)).^2);
    SS_err = sum((funcSS - pred).^2);
    R2(p,iSS) = 1-(SS_err / SS_tot);


    b{p,iSS} = Mdl{p,iSS}.Beta;
    beta2mean(p,iSS) = corr(b{p,iSS},meanfunc','r','p');
    beta2max0(p,iSS) = corr(b{p,iSS},maxfunc','r','p');
    beta2min0(p,iSS) = corr(b{p,iSS},minfunc','r','p');
    beta2extrem0(p,iSS) = corr(b{p,iSS},extremfunc','r','p');
    tmp = maxfunc; tmp(tmp==0) = nan; beta2max(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = minfunc; tmp(tmp==0) = nan; beta2min(p,iSS) = corr(b{p,iSS},tmp','r','p');
    tmp = extremfunc; tmp(tmp==0) = nan; beta2extrem(p,iSS) = corr(b{p,iSS},tmp','r','p');
    end
    ssnames{end+1} = structure{irh};
end

clear SS* bBrain bSS c* func* i iSS ind m* ss tmp z
if ~exist([yourpath 'models_query_aparc_perm/'])
    unix(['mkdir ' yourpath 'models_query_aparc_perm/'])
end
save([yourpath 'models_query_aparc_perm/' task '_' mat2str(numperm) 'permutations_wpred.mat'])

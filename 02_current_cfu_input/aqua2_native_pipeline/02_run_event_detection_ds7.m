%% AQuA2 native 2D event detection on DS7 patch slice Z06
% One MATLAB client with a controlled local parpool.
% Output is isolated from the previous DS2 and patch results.

close all;
clc;
clearvars;

z = str2double(getenv('AQUA_Z'));
if ~isfinite(z) || z ~= round(z) || z < 0 || z > 12
    error('AQUA_Z must be an integer in 0:12');
end
z = round(z);

runStart = tic;
fprintf('RUN_START %s\n', datestr(now, 31));
fprintf('SLICE_Z %02d\n', z);

% AQuA2 contains parfor sections in preprocessing and active-region
% detection. Use a bounded pool rather than launching independent MATLAB
% processes, so the input movie is not duplicated unnecessarily.
delete(gcp('nocreate'));
try
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
catch
end
poolWorkers = str2double(getenv('AQUA_POOL_WORKERS'));
if ~isfinite(poolWorkers) || poolWorkers < 1, poolWorkers = 8; end
pool = parpool('local', poolWorkers);
fprintf('PARPOOL_WORKERS %d\n', pool.NumWorkers);
maxNumCompThreads(4);

% Paths can be supplied by AQUA_BASE_DIR and AQUA2_DIR for batch runs.
% Current server input/output layout (directly runnable here):
baseDir = '/home/cyf/wbi/wbi_code/experiments/_trash/roi_ca_dependency/aqua_validation/';
% Portable placeholder (uncomment and edit for another machine):
% baseDir = '/path/to/your/aqua_validation/';
envBase = getenv('AQUA_BASE_DIR');
if ~isempty(envBase), baseDir = [envBase filesep]; end
sliceDir = [baseDir 'output/slices/dff/'];
outDir = [baseDir 'output/slices/results_ds7_all_parallel8/'];
mkdir(outDir);

startupDir = getenv('AQUA2_DIR');
if isempty(startupDir), startupDir = [baseDir '../AQuA2/']; end
cd(startupDir);
startup;

f1 = sprintf('slice_Z%02d.mat', z);
outFile = sprintf('%sslice_Z%02d_ds7_AQuA2_parallel8.mat', outDir, z);
fprintf('INPUT %s\n', fullfile(sliceDir, f1));
fprintf('OUTPUT %s\n', outFile);

% Same native AQuA2 parameter family as the existing DS7 2D pipeline.
opts.singleChannel = true;
opts.regMaskGap = 5;
opts.frameRate = 2.0;
opts.spatialRes = 1.0;

opts.registrateCorrect = 1;
opts.bleachCorrect = 1;
opts.medSmo = 0;
opts.smoXY = 0.5;

opts.thrARScl = 3;
opts.minDur = 5;
opts.minSize = 10;
opts.maxSize = inf;
opts.circularityThr = 0;
opts.spaMergeDist = 0;

opts.needTemp = 1;
opts.seedSzRatio = 0.01;
opts.sigThr = 3.5;
opts.maxDelay = 0.6;
opts.needRefine = 0;
opts.needGrow = 0;

opts.needSpa = 1;
opts.sourceSzRatio = 0.01;
opts.sourceSensitivity = 8;
opts.whetherExtend = 1;

opts.detectGlo = 0;
opts.ignoreTau = 1;
opts.propMetric = 0;
opts.networkFeatures = 0;

opts.gtwSmo = 0.2;
opts.ratio = 0.5;
opts.cut = 800;
opts.movAvgWin = 25;
opts.minShow1 = 0.2;
opts.correctTrend = 1;
opts.propthrmin = 0.5;
opts.propthrstep = 0.1;
opts.propthrmax = 0.5;
opts.compress = 0;
opts.gapExt = 5;
opts.TPatch = 20;
opts.maxSpaScale = 7;
opts.minSpaScale = 3;

opts.varEst = 0.02;
opts.fgFluo = 0;
opts.bgFluo = 0;
opts.northx = 0;
opts.northy = 1;

bd = containers.Map;
bd('None') = [];

sliceStart = tic;
try
    fprintf('PREP_START %s\n', datestr(now, 31));
    [datOrg1, datOrg2, opts] = burst.prep1(sliceDir, f1, sliceDir, [], [], opts);
    [H, W, ~, T] = size(datOrg1);
    fprintf('DATA H=%d W=%d T=%d\n', H, W, T);
    opts.singleChannel = isempty(datOrg2);
    opts.sz = size(datOrg1);
    evtSpatialMask = true(opts.sz(1:3));

    fprintf('BASELINE_START %s\n', datestr(now, 31));
    [dF1, opts] = pre.baselineRemoveAndNoiseEstimation(datOrg1, opts, evtSpatialMask, 1, []);
    opts.maxdF1 = min(100, max(dF1(:)));

    fprintf('AR_START %s\n', datestr(now, 31));
    arLst1 = act.acDetect(dF1, opts, evtSpatialMask, 1, []);
    fprintf('ACTIVE_REGIONS %d\n', numel(arLst1));

    fprintf('SE_START %s\n', datestr(now, 31));
    opts.step = 0.5;
    [seLst1, subEvtLst1, seLabel1, majorInfo1, opts, ~, ~, ~] = ...
        se.seDetection(dF1, datOrg1, arLst1, opts, []);
    fprintf('SUPER_EVENTS %d\n', numel(seLst1));

    fprintf('EVENT_START %s\n', datestr(now, 31));
    [riseLst1, datR1, evt1, ~] = evt.se2evtTop(dF1, seLst1, subEvtLst1, ...
        seLabel1, majorInfo1, opts, []);
    fprintf('EVENTS %d\n', numel(evt1));

    fprintf('FEATURE_START %s\n', datestr(now, 31));
    opts.stdMapOrg = opts.stdMapOrg1;
    opts.maxValueDat = opts.maxValueDat1;
    opts.minValueDat = opts.minValueDat1;
    opts.tempVarOrg = opts.tempVarOrg1;
    opts.correctPars = opts.correctPars1;
    [fts1, dffMat1, dMat1, dffAlignedMat1] = fea.getFeaturesTop(datOrg1, evt1, opts, []);

    res.maxVal = opts.maxValueDat1;
    res.opts = opts;
    res.datOrg1 = datOrg1;
    res.evt1 = evt1;
    res.fts1 = fts1;
    res.dffMat1 = dffMat1;
    res.dMat1 = dMat1;
    res.dffAlignedMat1 = dffAlignedMat1;
    res.riseLst1 = riseLst1;
    res.dF1 = dF1;
    res.seLst1 = seLst1;
    res.stg.post = 1;
    res.stg.detect = 1;
    res.bd = bd;

    ov = containers.Map('UniformValues', 0);
    ov('None') = [];
    ov1 = ui.over.getOv([], evt1, opts.sz, datR1, 1);
    ov1.name = 'Events';
    ov1.colorCodeType = {'Random'};
    ov('Events_Red') = ov1;
    res.ov = ov;

    fprintf('SAVE_START %s\n', datestr(now, 31));
    save(outFile, 'res', '-v7.3');
    fprintf('SAVED %s\n', outFile);
    fprintf('RUN_END %s\n', datestr(now, 31));
    fprintf('ELAPSED_SEC %.3f\n', toc(runStart));
catch ME
    fprintf(2, 'RUN_ERROR %s\n', ME.message);
    for i = 1:numel(ME.stack)
        fprintf(2, '  at %s line %d\n', ME.stack(i).file, ME.stack(i).line);
    end
    rethrow(ME);
end

delete(pool);

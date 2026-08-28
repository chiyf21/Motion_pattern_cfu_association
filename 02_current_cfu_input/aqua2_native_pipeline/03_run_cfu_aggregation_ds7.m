%% Native AQuA2 CFU aggregation for DS7 patch Z06, min 5 events
% Reuses the completed DS7 event result; does not rerun event detection.

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

% Paths can be supplied by AQUA_BASE_DIR and AQUA2_DIR for batch runs.
% Current server input/output layout (directly runnable here):
baseDir = '/home/cyf/wbi/wbi_code/experiments/_trash/roi_ca_dependency/aqua_validation/';
% Portable placeholder (uncomment and edit for another machine):
% baseDir = '/path/to/your/aqua_validation/';
envBase = getenv('AQUA_BASE_DIR');
if ~isempty(envBase), baseDir = [envBase filesep]; end
aquaDir = getenv('AQUA2_DIR');
if isempty(aquaDir), aquaDir = [baseDir '../AQuA2/']; end
resultFile = sprintf('%soutput/slices/results_ds7_all_parallel8/slice_Z%02d_ds7_AQuA2_parallel8.mat', ...
    baseDir, z);
outDir = [baseDir 'output/slices/cfu_ds7_all_ot030_min5_group5/'];
outFile = sprintf('%sslice_Z%02d_ds7_native_CFU_ot030_min5_group5.mat', outDir, z);
if ~isfolder(outDir)
    mkdir(outDir);
end

cd(aquaDir);
startup;

fprintf('INPUT %s\n', resultFile);
fprintf('OUTPUT %s\n', outFile);

try
    S = load(resultFile, 'res');
    res = S.res;
    nEvents = numel(res.evt1);
    fprintf('EVENTS %d\n', nEvents);

    res.fts1.curve.dffMaxFrame = nan(1, nEvents);

    cfuOpts.cfuDetect.overlapThr1 = 0.30;
    cfuOpts.cfuDetect.overlapThr2 = 0.50;
    cfuOpts.cfuDetect.minNumEvt1 = 5;
    cfuOpts.cfuDetect.minNumEvt2 = 3;
    cfuOpts.cfuAnalysis.maxDist = 10;
    cfuOpts.cfuAnalysis.shift = 0;
    cfuOpts.cfuGroup.pValueThr = 1e-5;
    cfuOpts.cfuGroup.cfuNumThr = 5;

    fprintf('CFU_START %s\n', datestr(now, 31));
    [cfuInfo1, cfuInfo2] = cfu.CFUdetectScript(res, cfuOpts);
    fprintf('CFUS %d\n', size(cfuInfo1,1));

    fprintf('RELATION_START %s\n', datestr(now, 31));
    cfuRelation = cfu.calAllDependencyScript(cfuInfo1, cfuInfo2, cfuOpts);
    fprintf('RELATIONS %d\n', size(cfuRelation,1));

    fprintf('GROUP_START %s\n', datestr(now, 31));
    cfuGroupInfo = cfu.groupCFUscript(cfuInfo1, cfuInfo2, cfuRelation, cfuOpts);
    fprintf('GROUPS %d\n', size(cfuGroupInfo,1));

    datPro = util.normalize01(mean(single(res.datOrg1), 4));
    elapsedSeconds = toc(runStart);
    save(outFile, 'cfuInfo1', 'cfuInfo2', 'cfuRelation', 'cfuGroupInfo', ...
        'cfuOpts', 'datPro', 'nEvents', 'elapsedSeconds', '-v7.3');
    fprintf('RUN_END %s\n', datestr(now, 31));
    fprintf('ELAPSED_SEC %.3f\n', elapsedSeconds);
    fprintf('SAVED %s\n', outFile);
catch ME
    fprintf(2, 'RUN_ERROR %s\n', ME.message);
    for i = 1:numel(ME.stack)
        fprintf(2, '  at %s line %d\n', ME.stack(i).file, ME.stack(i).line);
    end
    rethrow(ME);
end

delete(pool);

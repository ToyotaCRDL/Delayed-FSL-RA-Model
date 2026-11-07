function check_with_python(what)
% CHECK_WITH_PYTHON  Compare the revised MATLAB scripts with the results computed in Python.
%
% The results in result1/, result2/ and python/results/ were computed with the Python implementation
% (python/run_all.py). This function checks that the MATLAB scripts give the same values.
%
%   check_with_python('probe')
%       Runs the probe (all sets in stimulus1_probe/param.mat, six speeds) with the revised settings and
%       compares the spike trains with stimulus1_probe/result2/probe_spiketrain_param_*_ps_*.csv.
%
%   For the scanning stimuli, run the main script first (in its folder), then call this function from the
%   same MATLAB session; it reads the variables the script left in the base workspace. A full run over the
%   150 sets takes long at dt = 0.0125 ms, so the set loop can be restricted, e.g.
%   "for mlp2 = [1 87 92 118]" in place of "for mlp2 = 1:size(param,1)". Sets that were not run are skipped.
%
%   check_with_python('braille')   after stimulus2_braille/mainBrailleRA.m   (spike trains of the last set run)
%   check_with_python('dot')       after stimulus3_dot_texture/mainDotTexRA.m (FR and TV RMSE)
%   check_with_python('delay_dot') after stimulus3_dot_texture_delay/mainDotTexRA.m
%   check_with_python('single')    after stimulus4_single_dot/mainSingleDotRA.m (RW)
%   check_with_python('delay_single') after stimulus4_single_dot_delay/mainSingleDotRA.m
%   check_with_python('modi_RW')   after stimulus4_single_dot_modi_RW/mainSingleDotRA.m (spike width, Fig. S6)
%
% The noise analyses (*_noise) use random numbers and are not compared.

root = fileparts(mfilename('fullpath'));
res = fullfile(root, '..', 'python', 'results', 'unitA');
switch what
    case 'probe'
        check_probe(root);
    case 'braille'
        check_braille(root);
    case 'dot'
        S = readtable(fullfile(res, 'summary.csv'), 'VariableNamingRule', 'preserve');
        report('FR RMSE', base_var('setFiringRateRMSE'), S.fr_rmse);
        report('TV RMSE', base_var('setTempoVariRMSE'), S.tv_rmse);
    case 'single'
        S = readtable(fullfile(res, 'summary.csv'), 'VariableNamingRule', 'preserve');
        report('RW (4.1, 2.6, 1.6 mm)', base_var('setRespWidth'), [S.('rw_4-1'), S.('rw_2-6'), S.('rw_1-6')]);
    case 'modi_RW'
        S = readtable(fullfile(res, 'summary.csv'), 'VariableNamingRule', 'preserve');
        report('spike width (4.1, 2.6, 1.6 mm)', base_var('spikeWidth'), [S.('sw_4-1'), S.('sw_2-6'), S.('sw_1-6')]);
    case 'delay_dot'
        D = readtable(fullfile(res, 'delay_control.csv'), 'VariableNamingRule', 'preserve');
        report('FR RMSE (0-30 ms)', base_var('setFiringRateRMSE'), D.fr_rmse);
        report('TV RMSE (0-30 ms)', base_var('setTempoVariRMSE'), D.tv_rmse);
    case 'delay_single'
        D = readtable(fullfile(res, 'delay_control.csv'), 'VariableNamingRule', 'preserve');
        rw = [D.('dRW_4-1'), D.('dRW_2-6'), D.('dRW_1-6')] + [5.13, 4.84, 3.88];
        report('RW (0-30 ms)', base_var('setRespWidth'), rw);
    otherwise
        error('unknown check: %s', what);
end
end


function check_probe(root)
global R C tau M beta lamda v_thr ref u_thr %#ok<NUSED>
R = 66.11e6; C = 1.01e-6; tau = 66.77; M = 1.0;
theta = 0.13; dt = 0.0125; u_thr = theta * dt;   % as in mainProbe.m
here = pwd;
cleanup = onCleanup(@() cd(here));
cd(fullfile(root, 'stimulus1_probe'));
load('param.mat', 'param');
data = xlsread(fullfile('finger61-1', 'f61-1.xlsx'), 'sheet1');
speeds = {'9-5', '14-8', '20-6', '34-0', '46-5', '60-6'};
start = 0.00; finish = 200.0; point = round((finish-start)/dt) + 1;
xq = linspace(start, finish, point);
nOK = 0; nBad = 0; worst = 0;
for mlp = 1:size(param, 1)
    beta = param(mlp,1); lamda = param(mlp,2); ref = param(mlp,3); v_thr = param(mlp,4);
    for d = 1:6
        F = griddedInterpolant(data(:,7)*1000, data(:,d), 'spline');
        dataMain = zeros(point, 4);
        dataMain(:,1) = F(xq); dataMain(:,2) = xq; dataMain(:,3) = dt; dataMain(:,4) = linspace(start, finish/1000, point);
        [~, spikeTime] = spike(dataMain);
        m = spikeTime(spikeTime > 0);
        p = readmatrix(fullfile('result2', sprintf('probe_spiketrain_param_%d_ps_%s.csv', mlp, speeds{d})));
        p = p(p > 0);
        if numel(m) == numel(p) && all(abs(m(:) - p(:)) < 1e-6)
            nOK = nOK + 1;
        else
            nBad = nBad + 1;
            fprintf('  set %d, %s mm/s: MATLAB %d spikes, Python %d spikes\n', mlp, speeds{d}, numel(m), numel(p));
        end
        if numel(m) == numel(p) && ~isempty(m); worst = max(worst, max(abs(m(:) - p(:)))); end
    end
end
fprintf('probe spike trains: %d match, %d differ (largest time difference among equal counts %.3g ms)\n', nOK, nBad, worst);
end


function check_braille(root)
allSpikeTime = base_var('allSpikeTime');
mlp2 = base_var('mlp2');
p = readmatrix(fullfile(root, 'stimulus2_braille', 'result1', sprintf('braille_sep_param_%d.csv', mlp2)));
nOK = 0; nBad = 0;
for i = 1:51
    m = allSpikeTime{1275 + i}; m = m(m > 0);
    q = p(i, :); q = q(q > 0);
    if numel(m) == numel(q) && all(abs(m(:) - q(:)) < 1e-6)
        nOK = nOK + 1;
    else
        nBad = nBad + 1;
    end
end
fprintf('Braille, set %d: %d afferents match, %d differ\n', mlp2, nOK, nBad);
end


function v = base_var(name)
if ~evalin('base', sprintf('exist(''%s'', ''var'')', name))
    error('variable %s not found: run the main script of this folder first', name);
end
v = evalin('base', name);
end


function report(label, m, p)
% m: MATLAB values (rows of sets or periods that were not run are 0); p: Python values
n = min(size(m, 1), size(p, 1));
m = m(1:n, :); p = p(1:n, :);
run = any(m ~= 0 | isnan(m), 2);
d = abs(m(run, :) - p(run, :));
nanDiff = xor(isnan(m(run, :)), isnan(p(run, :)));
fprintf('%s: %d rows compared, largest |MATLAB - Python| = %.3g, NaN in only one of them: %d\n', ...
    label, sum(run), max([d(~isnan(d)); 0]), sum(nanDiff(:)));
end

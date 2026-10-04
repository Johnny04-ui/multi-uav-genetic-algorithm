function run_two_stage_ga(random_seed)
% =========================================================================
%  名称：run_two_stage_ga
%  功能：两阶段遗传算法（GA）优化 5 架 UAV × 每架 3 枚烟雾弹（共 40 基因）
%        以“33 条视线同时被任一存活云遮蔽”的总时长为硬指标进行优化
%       （软指标为各导弹被遮蔽比例的时间积分）。
%
%  关键设定（与题设一致）：
%    R=10，tdown=3，T_eff=20，v_missile=300，g=9.80665
%    “云落地仍有效”：z<0 时钳到 z=0（仅竖直下沉，不再水平移动）
%    **重要：评估时间与起爆时刻均被限制在 ≤ 67 s**
%
%  两阶段策略：
%    阶段1：步长较大、搜索更“大胆”（dt=1.0）
%    阶段2：在阶段1 可行解附近“精修”（dt=0.01）
%
%  输出：
%    - 终评估（硬/软指标）
%    - Excel（参数与结果）
%    - 遮蔽时间线图
% =========================================================================

% 如果调用者没有指定随机种子，就使用 42。
% 固定种子可以让同一份代码更容易复现和排查问题。
if nargin < 1
    random_seed = 42;
end

clc; close all;

% 结果统一写入仓库的 results 文件夹，避免散落在当前工作目录。
repo_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(repo_dir, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% ========================================================================
%% 一、场景与物理常量（不改题设）
%% ========================================================================
% 三枚导弹初始坐标（直线匀速飞向原点）
M1_ini = [20000, 0, 2000];
M2_ini = [19000, 600, 2100];
M3_ini = [18000, -600, 1900];

% 五架 UAV 初始坐标
FY_ini = [17800,    0, 1800;
          12000, 1400, 1400;
           6000,-3000,  700;
          11000, 2000, 1800;
          13000,-2000, 1300];

% 真实目标 11 个采样点（固定）
P = [0,193,10; 0,207,10; 7,200,10; -7,200,10;
     0,207,5; 7,200,5; -7,200,5; 7,200,0;
     -7,200,0; 0,207,0; 0,193,0];

% 物理常量
g           = 9.80665;   % 重力加速度
R           = 10;        % 云等效半径
v_missile   = 300;       % 导弹速度
tdown       = 3;         % 云竖直下沉速度
T_eff       = 20;        % 单枚云有效时长

% 导弹直线飞向原点：预定义 3 个“位置函数句柄”
M_funcs = { ...
    @(t) M1_ini - v_missile .* t .* (M1_ini ./ norm(M1_ini)); ...
    @(t) M2_ini - v_missile .* t .* (M2_ini ./ norm(M2_ini)); ...
    @(t) M3_ini - v_missile .* t .* (M3_ini ./ norm(M3_ini)) ...
};

%% ========================================================================
%% 二、基因（染色体）结构与边界
%% ========================================================================
% 每架 UAV 8 基因：[theta, v, t1, f1, t2, f2, t3, f3]
U = 5;                 % UAV 数
genes_per_uav = 8;     % 每架 8 个基因
nVars = U * genes_per_uav;

% 边界（按题设）
theta_min = 0;      theta_max = 2*pi;    % 航向角
v_min     = 70;     v_max     = 140;     % 平飞速度 m/s
t_min     = 0;      t_max     = 67;      % 投放相对时刻上限（**重要：总时轴上限**）
f_min     = 0.1;    f_max     = 30;      % 引信时间
min_sep   = 1.0;                        % 投放最小间隔（工程约束）

[lb, ub] = build_bounds(U, genes_per_uav, theta_min, theta_max, v_min, v_max, t_min, t_max, f_min, f_max);

%% ========================================================================
%% 三、评估时间步（两阶段）
%% ========================================================================
dt_stage1 = 1.0;   % 阶段1：粗搜索，每 1 秒检查一次
dt_stage2 = 0.01;  % 阶段2：细

%% ========================================================================
%% 四、GA 参数（两阶段）
%% ========================================================================
rng(random_seed, 'twister');  % 固定随机数生成方式，便于复现

% 阶段1：探索性强
stage1.pop_size    = 80;
stage1.max_gen     = 120;
stage1.elite_num   = 2;
stage1.tour_k      = 3;
stage1.p_crossover = 0.9;
stage1.p_mutation  = 0.35;
stage1.p_gene_mut  = 0.12;
stage1.noise_frac  = 0.10;   % 高斯噪声相对边界长度的比例
stage1.print_int   = 10;
stage1.lambda_soft = 0.6;    % 软目标权重

% 阶段2：精修
stage2.pop_size    = 60;
stage2.max_gen     = 100;
stage2.elite_num   = 3;
stage2.tour_k      = 3;
stage2.p_crossover = 0.9;
stage2.p_mutation  = 0.25;
stage2.p_gene_mut  = 0.08;
stage2.noise_frac  = 0.06;
stage2.print_int   = 10;
stage2.lambda_soft = 0.3;

% 阶段2 初始种群从阶段1 挑选前 K（优先硬指标>0，否则按软指标）
topK_for_stage2 = 12;

%% ========================================================================
%% 五、构造启发式“种子”个体（可行性增强）
%% ========================================================================
% 思路：令起爆点靠近采样点簇 (0,200,z≈10/5/0)，提升“33 条同时遮蔽”的可能性
seed_list = build_seed_population(U, genes_per_uav, FY_ini, ...
    v_min, v_max, t_min, t_max, f_min, f_max, min_sep);
ind40 = [0.0983, 135.3515, 0.0000, 0.0000, 1.0000, 0.0000, 32.7378, 23.0453, 3.1686, 77.4467, 46.5928, 20.4534, 67.8979, 5.8249, 69.0293, 23.6537, 2.6554, 125.7669, 6.2269, 48.1678, 55.1915, 25.1582, 66.7618, 22.8242, 4.9656, 111.1352, 57.6155, 23.0455, 66.7042, 0.9008, 70.0000, 18.6917, 2.6693, 135.8007, 3.9803, 36.8251, 61.2821, 46.9894, 70.0000, 13.0898];
seed_list = [seed_list; ind40];

%% ========================================================================
%% 六、阶段1：大胆搜索（dt=0.1）
%% ========================================================================
fprintf('\n=== 阶段 1：大胆搜索（dt=%.2f）===\n', dt_stage1);
[best1, best1_stats, pop1, fits1, aux1] = run_ga_stage( ...
    seed_list, stage1, ...
    U, genes_per_uav, lb, ub, ...
    M_funcs, FY_ini, P, R, g, tdown, T_eff, dt_stage1, ...
    theta_min, theta_max, v_min, v_max, t_min, t_max, f_min, f_max, min_sep);

fprintf('阶段1结束 | 最佳硬目标(33条) = %.4f s | 软目标 = %.4f s | 组合适应度 = %.4f\n', ...
    best1_stats.glob_cov, best1_stats.soft_cov, best1_stats.fitness);

%% ========================================================================
%% 七、阶段2：精修（dt=0.01）
%% ========================================================================
fprintf('\n=== 阶段 2：精修（dt=%.2f）===\n', dt_stage2);

% 从阶段1 结果中挑前 K 个体（优先硬指标>0），并在其附近加小扰动
seeds2 = select_top_for_stage2(pop1, aux1, topK_for_stage2, lb, ub, stage2.noise_frac);

[best2, best2_stats] = run_ga_stage( ...
    seeds2, stage2, ...
    U, genes_per_uav, lb, ub, ...
    M_funcs, FY_ini, P, R, g, tdown, T_eff, dt_stage2, ...
    theta_min, theta_max, v_min, v_max, t_min, t_max, f_min, f_max, min_sep);

fprintf('阶段2结束 | 最佳硬目标(33条) = %.4f s | 软目标 = %.4f s | 组合适应度 = %.4f\n', ...
    best2_stats.glob_cov, best2_stats.soft_cov, best2_stats.fitness);

best_ind = best2;
print_best_genes(best2, U, genes_per_uav);

%% ========================================================================
%% 八、终评估 + 导出 + 绘图（dt=阶段2）
%% ========================================================================
[glob_cov, cov_per_missile, cov_missile_vecs, cov_global_vec, tvec, soft_cov] = ...
    fitness_eval_33LOS_grounded(best_ind, M_funcs, FY_ini, P, R, g, tdown, T_eff, dt_stage2);

fprintf('\n【终评估】硬目标(33条同时) = %.4f s | 软目标 = %.4f s \n', glob_cov, soft_cov);

[release_pts, explode_pts] = decode_and_compute_points(best_ind, U, genes_per_uav, FY_ini, g);

result_file = fullfile(results_dir, 'latest_run.xlsx');
export_results_to_excel(result_file, best_ind, U, genes_per_uav, release_pts, explode_pts, ...
    cov_per_missile, glob_cov);

plot_results_GA(cov_missile_vecs, cov_global_vec, tvec, glob_cov, cov_per_missile);

% 将图保存下来，方便不安装 MATLAB 的读者直接查看结果。
figure_file = fullfile(results_dir, 'latest_run.png');
exportgraphics(gcf, figure_file, 'Resolution', 160);
fprintf('Figure saved to %s\n', figure_file);

end
% ============================ 主函数结束 ================================



%% ========================================================================
%% 子函数：单阶段 GA 主循环
%% ========================================================================
function [best_ind, best_stats, final_pop, final_fits, aux] = run_ga_stage( ...
    seed_list, ga, ...                 % 种子个体 + GA 参数结构体
    U, genes_per_uav, lb, ub, ...      % UAV 数量、基因数、基因上下界
    M_funcs, FY_ini, P, R, g, tdown, T_eff, dt, ...  % 物理参数 & 评估用函数
    theta_min, theta_max, v_min, v_max, t_min, t_max, f_min, f_max, min_sep)

% 说明：
% - 初始种群由给定 seed_list 填充，剩下个体随机生成
% - 适应度 = 硬指标（33条同时遮挡时间） + lambda_soft * 软指标
% - 遗传循环：精英保留 → 锦标赛选择 → 交叉 → 变异 → 修复 → 评估 → 更新最优

nVars = numel(lb);   % 每个个体的基因数（40 维）

%% ---------------- 1) 初始化种群 ----------------
pop = zeros(ga.pop_size, nVars);    % 种群矩阵（行=个体，列=基因）
S = size(seed_list, 1);             % 种子数量
useS = min(S, ga.pop_size);         % 实际能用的种子个数
if useS > 0
    pop(1:useS, :) = seed_list(1:useS, :);   % 前 useS 个直接用种子
end
for i = useS+1:ga.pop_size
    pop(i,:) = sample_individual(U, genes_per_uav, lb, ub); % 随机初始化
end
for i=1:ga.pop_size
    pop(i,:) = repair_individual(pop(i,:), ...              % 修复基因（含 t+f≤67）
        theta_min, theta_max, v_min, v_max, ...
        t_min, t_max, f_min, f_max, min_sep);
end

%% ---------------- 2) 初评估 ----------------
fits  = zeros(ga.pop_size,1);  % 适应度向量
gvals = zeros(ga.pop_size,1);  % 硬目标（33条完全遮挡时间）
svals = zeros(ga.pop_size,1);  % 软目标（部分遮挡程度）

for i=1:ga.pop_size
    [gvals(i), ~, ~, ~, ~, svals(i)] = ...
        fitness_eval_33LOS_grounded(pop(i,:), M_funcs, FY_ini, ...
                                    P, R, g, tdown, T_eff, dt);
    fits(i) = gvals(i) + ga.lambda_soft * svals(i);
end

[best_fit, idx] = max(fits);
best_ind = pop(idx,:);
best_stats.fitness  = best_fit;
best_stats.glob_cov = gvals(idx);
best_stats.soft_cov = svals(idx);

if ga.print_int>0
    fprintf('Gen %3d | Best=%.4f | Mean=%.4f | (Hard=%.4f, Soft=%.4f)\n', ...
        0, best_fit, mean(fits), best_stats.glob_cov, best_stats.soft_cov);
    print_best_genes(best_ind, U, genes_per_uav);
end

%% ---------------- 3) 遗传循环 ----------------
for gen = 1:ga.max_gen
    % ---- 3.1 精英保留 ----
    [~, order] = sort(fits, 'descend');
    new_pop = zeros(size(pop));
    new_pop(1:ga.elite_num,:) = pop(order(1:ga.elite_num),:);

    % ---- 3.2 生成后代 ----
    ptr = ga.elite_num + 1;
    while ptr <= ga.pop_size
        p1 = pop(tournament_select(fits, ga.tour_k),:);
        p2 = pop(tournament_select(fits, ga.tour_k),:);
        c1 = p1; c2 = p2;

        if rand < ga.p_crossover
            alpha = rand(1, nVars);
            c1 = alpha .* p1 + (1-alpha) .* p2;
            c2 = alpha .* p2 + (1-alpha) .* p1;
        end

        if rand < ga.p_mutation
            mask = rand(1,nVars) < ga.p_gene_mut;
            span = (ub - lb);
            noise = randn(1,nVars) .* (ga.noise_frac * span);
            c1(mask) = c1(mask) + noise(mask);
        end
        if rand < ga.p_mutation
            mask = rand(1,nVars) < ga.p_gene_mut;
            span = (ub - lb);
            noise = randn(1,nVars) .* (ga.noise_frac * span);
            c2(mask) = c2(mask) + noise(mask);
        end

        % 修复 & 边界裁剪（含 t+f ≤ 67）
        c1 = max(min(c1, ub), lb);
        c2 = max(min(c2, ub), lb);
        c1 = repair_individual(c1, theta_min, theta_max, v_min, v_max, ...
                               t_min, t_max, f_min, f_max, min_sep);
        c2 = repair_individual(c2, theta_min, theta_max, v_min, v_max, ...
                               t_min, t_max, f_min, f_max, min_sep);

        new_pop(ptr,:) = c1; ptr = ptr + 1;
        if ptr <= ga.pop_size
            new_pop(ptr,:) = c2; ptr = ptr + 1;
        end
    end

    % ---- 3.3 评估新种群 ----
    pop = new_pop;
    for i=1:ga.pop_size
        [gvals(i), ~, ~, ~, ~, svals(i)] = ...
            fitness_eval_33LOS_grounded(pop(i,:), M_funcs, FY_ini, ...
                                        P, R, g, tdown, T_eff, dt);
        fits(i) = gvals(i) + ga.lambda_soft * svals(i);
    end

    % ---- 3.4 更新最优个体 ----
    [cur_best, idx] = max(fits);
    if cur_best > best_fit
        best_fit = cur_best;
        best_ind = pop(idx,:);
        best_stats.fitness  = best_fit;
        best_stats.glob_cov = gvals(idx);
        best_stats.soft_cov = svals(idx);
    end

    if ga.print_int>0 && (mod(gen, ga.print_int)==0 || gen==ga.max_gen)
        fprintf('Gen %3d | Best=%.4f | Mean=%.4f | (Hard=%.4f, Soft=%.4f)\n', ...
            gen, best_fit, mean(fits), best_stats.glob_cov, best_stats.soft_cov);
        print_best_genes(best_ind, U, genes_per_uav);
    end
end

final_pop = pop;
final_fits = fits;
aux.gvals = gvals;
aux.svals = svals;
end



%% ========================================================================
%% 子函数：锦标赛选择
%% ========================================================================
function idx = tournament_select(fits, k)
n = numel(fits);
cands = randi(n, [k,1]);
[~, j] = max(fits(cands));
idx = cands(j);
end



%% ========================================================================
%% 评估函数：落实“落地仍有效”的 33 条同时遮蔽判定（**评估≤67 s**）
%% ========================================================================
function [glob_cov, cov_per_missile, cover_missile_vecs, cover_global_vec, tvec, soft_cov] = ...
    fitness_eval_33LOS_grounded(ind, M_funcs, FY_ini, P, R, g, tdown, T_eff, dt)

% --- 时间硬上限（场景总时轴） ---
SIM_TMAX = 67;   % **关键：评估窗口不会超过 67 s**

U = size(FY_ini,1);
genes_per_uav = 8;
numMiss = numel(M_funcs);

% ---------- 0) 解码 15 朵云：C0 与 t_exp ----------
expl_times   = zeros(U,3);   % 起爆时刻
exp_centers0 = zeros(U,3,3); % 起爆点 C0
for u=1:U
    base = (u-1)*genes_per_uav;
    th   = ind(base+1);
    v    = ind(base+2);
    t_rels  = [ind(base+3), ind(base+5), ind(base+7)];
    t_fuses = [ind(base+4), ind(base+6), ind(base+8)];
    dir = [cos(th), sin(th), 0];
    for k=1:3
        rel = FY_ini(u,:) + v * t_rels(k) * dir;                      % 投放点
        C0  = rel + v * dir * t_fuses(k) - [0,0,0.5*g*(t_fuses(k)^2)]; % 起爆点
        exp_centers0(u,k,:) = C0;
        expl_times(u,k)     = t_rels(k) + t_fuses(k);                  % 起爆时刻
    end
end

% ---------- 1) 评估时间轴（**截断到 ≤67 s**） ----------
t_start = max(0, min(expl_times(:)));
t_end_nominal = max(expl_times(:)) + T_eff;
t_end = min(t_end_nominal, SIM_TMAX);   % 关键：不超过 67

% 对齐到时间栅格，避免边界误差
t_start = ceil(t_start/dt)*dt;
t_end   = floor(t_end/dt)*dt;

if isempty(t_start) || t_end < t_start
    tvec = t_start; Tn = 1;
    cover_missile_vecs = false(numMiss, Tn);
    cover_global_vec   = false(1, Tn);
    cov_per_missile    = zeros(1, numMiss);
    glob_cov = 0; soft_cov = 0;
    return;
end

tvec = t_start:dt:t_end;
Tn   = numel(tvec);

% ---------- 2) 逐时刻遮蔽判定 ----------
cover_missile_vecs = false(numMiss, Tn);  % 每枚导弹：11 条是否同时被挡
cover_global_vec   = false(1, Tn);        % 是否 33 条同时
cov_per_missile    = zeros(1, numMiss);   % 每导弹 11 条同时遮蔽时长
soft_cov           = 0;                   % 软指标积分

for ti = 1:Tn
    t_abs = tvec(ti);

    % 2.1 收集存活云中心（落地仍参与遮蔽）
    live_centers = [];
    for u=1:U
        for k=1:3
            s = t_abs - expl_times(u,k);
            if s < 0 || s > T_eff, continue; end
            C0 = squeeze(exp_centers0(u,k,:))';
            cz = C0(3) - tdown * s;  if cz < 0, cz = 0; end
            live_centers(end+1,:) = [C0(1), C0(2), cz]; %#ok<AGROW>
        end
    end

    if isempty(live_centers)
        cover_global_vec(ti) = false;
        soft_cov = soft_cov + 0 * dt;
        continue;
    end

    % 2.2 对每枚导弹，统计 11 条视线被挡数量
    fracs = zeros(1,numMiss);
    all_ok = true;
    for m=1:numMiss
        Mpos = M_funcs{m}(t_abs);
        blocked_cnt = 0;
        for pi = 1:size(P,1)
            A = Mpos; B = P(pi,:);
            blocked = false;
            for ci=1:size(live_centers,1)
                if line_segment_sphere_intersect(A, B, live_centers(ci,:), R)
                    blocked = true; break;
                end
            end
            if blocked, blocked_cnt = blocked_cnt + 1; end
        end
        fracs(m) = blocked_cnt / size(P,1);
        cover_missile_vecs(m,ti) = (blocked_cnt == size(P,1));
        if ~cover_missile_vecs(m,ti), all_ok = false; end
    end

    % 2.3 全局 33 条同时
    cover_global_vec(ti) = all_ok;

    % 2.4 软指标积分
    soft_cov = soft_cov + min(fracs) * dt;
end

% ---------- 3) 统计时长 ----------
cov_per_missile = sum(cover_missile_vecs,2)' * dt;
glob_cov        = sum(cover_global_vec) * dt;
end



%% ========================================================================
%% 几何判定：线段-球体相交
%% ========================================================================
function tf = line_segment_sphere_intersect(A,B,C,R)
d = B - A;
a = dot(d,d);
if a == 0
    tf = (norm(A - C) <= R);
    return;
end
f = A - C;
b = 2*dot(f,d);
c = dot(f,f) - R^2;
disc = b^2 - 4*a*c;
if disc < 0
    tf = false;
    return;
end
t1 = (-b - sqrt(disc)) / (2*a);
t2 = (-b + sqrt(disc)) / (2*a);
tf = ((t1 >= 0 && t1 <= 1) || (t2 >= 0 && t2 <= 1));
end



%% ========================================================================
%% 工具：构建基因边界
%% ========================================================================
function [lb, ub] = build_bounds(U, genes_per_uav, theta_min, theta_max, v_min, v_max, t_min, t_max, f_min, f_max)
nVars = U*genes_per_uav;
lb = zeros(1,nVars);
ub = zeros(1,nVars);
for u = 1:U
    base = (u-1)*genes_per_uav;
    lb(base+1) = theta_min; ub(base+1) = theta_max; % theta
    lb(base+2) = v_min;     ub(base+2) = v_max;     % v
    lb(base+3) = t_min;     ub(base+3) = t_max;     % t1（投放）
    lb(base+4) = f_min;     ub(base+4) = f_max;     % f1（引信）
    lb(base+5) = t_min;     ub(base+5) = t_max;     % t2
    lb(base+6) = f_min;     ub(base+6) = f_max;     % f2
    lb(base+7) = t_min;     ub(base+7) = t_max;     % t3
    lb(base+8) = f_min;     ub(base+8) = f_max;     % f3
end
end



%% ========================================================================
%% 工具：随机采样一个个体（初始确保 t1<t2<t3）
%% ========================================================================
function x = sample_individual(U, genes_per_uav, lb, ub)
nVars = U*genes_per_uav;
x = lb + rand(1,nVars) .* (ub - lb);
for u=1:U
    base = (u-1)*genes_per_uav;
    t = [x(base+3), x(base+5), x(base+7)];
    t = sort(t);
    x(base+3)=t(1); x(base+5)=t(2); x(base+7)=t(3);
end
end



%% ========================================================================
%% 工具：个体修复（**含起爆时刻≤67 s 的强制约束**）
%% ========================================================================
function ind = repair_individual(ind, theta_min, theta_max, v_min, v_max, t_min, t_max, f_min, f_max, min_sep)
% 约束意图：
% 1) t_k（投放相对时刻）∈ [t_min, t_max] 但要给最小引信留出空间：t_k ≤ t_max - f_min
% 2) 引信 f_k ∈ [f_min, f_max] 且满足 t_k + f_k ≤ t_max（=67）
% 3) t1 < t2 < t3 且相邻间隔 ≥ min_sep
genes_per_uav = 8;
U = numel(ind)/genes_per_uav;

for u=1:U
    base = (u-1)*genes_per_uav;

    % 角度与标量边界
    ind(base+1) = mod(ind(base+1), 2*pi);                           % theta
    ind(base+2) = min(max(ind(base+2), v_min), v_max);              % v
    ind(base+4) = min(max(ind(base+4), f_min), f_max);              % f1
    ind(base+6) = min(max(ind(base+6), f_min), f_max);              % f2
    ind(base+8) = min(max(ind(base+8), f_min), f_max);              % f3

    % ---------- 投放时间：先考虑“必须留出最小引信时间”的上限 ----------
    t_release_max = t_max - f_min;              % 保证最少还能引爆且不超过 67
    t = [ind(base+3), ind(base+5), ind(base+7)];
    t = max(t, t_min);
    t = min(t, t_release_max);                  % 关键：t_k ≤ 67 - f_min
    t = sort(t);

    % 最小间隔
    for j=2:3
        if t(j)-t(j-1) < min_sep
            t(j) = t(j-1) + min_sep;
        end
    end
    % 如果因间隔挤到上界，整体回退
    if t(3) > t_release_max
        shift = t(3)-t_release_max;
        t = t - shift;
        if t(1) < t_min
            t(1) = t_min;
            t(2) = max(t(2), t(1)+min_sep);
            t(3) = max(t(3), t(2)+min_sep);
            if t(3) > t_release_max
                % 在线性铺开时也不能超过 t_release_max
                t = linspace(t_min, t_release_max, 3);
            end
        end
    end
    ind(base+3)=t(1); ind(base+5)=t(2); ind(base+7)=t(3);

    % ---------- 引信：确保 t_k + f_k ≤ 67 ----------
    for kk = 1:3
        t_idx = base + 2 + (kk-1)*2;  % 3,5,7
        f_idx = base + 3 + (kk-1)*2;  % 4,6,8
        % 先保证投放时刻自身合法（冗余保护）
        ind(t_idx) = min(max(ind(t_idx), t_min), t_release_max);
        % 可允许的最大引信：不让起爆超过 67
        max_f_allowed = max(f_min, t_max - ind(t_idx) - 1e-9);
        ind(f_idx) = min(max(ind(f_idx), f_min), min(f_max, max_f_allowed));
    end
end
end



%% ========================================================================
%% 启发式种子：让起爆点靠近采样簇 (0,200,z≈10/5/0)（**种子也遵守起爆≤67**）
%% ========================================================================
function seeds = build_seed_population(U, genes_per_uav, FY_ini, ...
    v_min, v_max, t_min, t_max, f_min, f_max, min_sep)

target_xy = [0,200];
target_z  = [10, 5, 0];
base_v = 0.5*(v_min+v_max);

seed_main = zeros(1, U*genes_per_uav);
for u=1:U
    base = (u-1)*genes_per_uav;
    du   = target_xy - FY_ini(u,1:2);
    th   = atan2(du(2), du(1));
    dist = norm(du);
    Ttot = dist / base_v;          % 水平到位总时长（粗估）

    z0 = FY_ini(u,3);
    fuses = zeros(1,3);
    for k=1:3
        zt = target_z(min(k, numel(target_z)));
        tf = sqrt(max(0, 2*max(0,z0-zt)/9.80665));
        fuses(k) = min(max(tf, f_min), f_max);
    end

    % 关键：t_rel ≤ 67 - f（保证起爆≤67）
    trels = max(t_min, Ttot - fuses);
    trels = min(trels, t_max - fuses);
    trels = max(trels, t_min);
    trels = sort(trels);
    for j=2:3
        if trels(j)-trels(j-1) < min_sep
            trels(j) = trels(j-1) + min_sep;
        end
    end
    % 防溢出
    trels = min(trels, t_max - f_min);

    seed_main(base+(1:8)) = [th, base_v, trels(1), fuses(1), trels(2), fuses(2), trels(3), fuses(3)];
end

% 生成若干扰动版种子（增强多样性）
J = 8;                  % 共 1+J 个种子
seeds = zeros(J+1, U*genes_per_uav);
seeds(1,:) = seed_main;

for j=1:J
    x = seed_main;
    noise = randn(size(x)) .* 0.05;
    for u=1:U
        base = (u-1)*genes_per_uav;
        x(base+1) = x(base+1) + noise(base+1)*0.3;      % theta
        x(base+2) = x(base+2) + noise(base+2)*5.0;      % v
        x(base+3) = x(base+3) + abs(noise(base+3))*3.0; % t1
        x(base+4) = x(base+4) + abs(noise(base+4))*1.5; % f1
        x(base+5) = x(base+5) + abs(noise(base+5))*3.0; % t2
        x(base+6) = x(base+6) + abs(noise(base+6))*1.5; % f2
        x(base+7) = x(base+7) + abs(noise(base+7))*3.0; % t3
        x(base+8) = x(base+8) + abs(noise(base+8))*1.5; % f3
    end
    % 统一修复（含 t+f ≤ 67）
    x = repair_individual(x, 0, 2*pi, v_min, v_max, t_min, t_max, f_min, f_max, min_sep);
    seeds(j+1,:) = x;
end
end



%% ========================================================================
%% 阶段2 初始种群选择：优先硬指标>0，并在附近加小扰动
%% ========================================================================
function seeds2 = select_top_for_stage2(pop1, aux1, K, lb, ub, noise_frac)
n = size(pop1,1);
hard = aux1.gvals;
soft = aux1.svals;

idx_hard = find(hard>0);
idx_easy = setdiff(1:n, idx_hard);
order = [idx_hard(:).' idx_easy(:).'];

[~,ord2] = sortrows([-double(hard(order)>0).', -soft(order).']);
order = order(ord2);

take = min(K, numel(order));
core = pop1(order(1:take), :);

seeds2 = zeros(take, size(pop1,2));
for i=1:take
    x = core(i,:);
    span = (ub - lb);
    noise = randn(size(x)) .* (noise_frac * span);
    x = max(min(x + noise, ub), lb);
    % 修复（含 t+f ≤ 67）
    x = repair_individual(x, 0, 2*pi, lb(2), ub(2), lb(3), ub(3), lb(4), ub(4), 1.0);
    seeds2(i,:) = x;
end
end



%% ========================================================================
%% 解码：从最优基因生成三次投放点与起爆点（用于导出展示）
%% ========================================================================
function [release_pts, explode_pts] = decode_and_compute_points(ind, U, genes_per_uav, FY_ini, g)
release_pts = zeros(U,3,3);
explode_pts = zeros(U,3,3);
for u=1:U
    base = (u-1)*genes_per_uav;
    th = ind(base+1);
    v  = ind(base+2);
    t_rels  = [ind(base+3), ind(base+5), ind(base+7)];
    t_fuses = [ind(base+4), ind(base+6), ind(base+8)];
    dir = [cos(th), sin(th), 0];
    for k=1:3
        rel  = FY_ini(u,:) + v * t_rels(k) * dir;
        expc = rel + v * dir * t_fuses(k) - [0,0,0.5*g*(t_fuses(k)^2)];
        release_pts(u,k,:) = rel;
        explode_pts(u,k,:) = expc;
    end
end
end



%% ========================================================================
%% 导出 Excel：参数表 + 每导弹 11 条时长 + 全局 33 条时长
%% ========================================================================
function export_results_to_excel(xlsxName, best_ind, U, genes_per_uav, release_pts, explode_pts, cov_per_missile, glob_cov)
theta_vec = zeros(U,1);
theta_deg = zeros(U,1);
v_vec = zeros(U,1);
for u=1:U
    theta_vec(u) = best_ind((u-1)*genes_per_uav+1);
    theta_deg(u) = rad2deg(theta_vec(u));
    v_vec(u)     = best_ind((u-1)*genes_per_uav+2);
end

Tparams = table((1:U)', theta_vec, theta_deg, v_vec, ...
    'VariableNames', {'UAV','theta_rad','theta_deg','v_mps'});

for k=1:3
    tcol=zeros(U,1); fcol=zeros(U,1);
    relx=zeros(U,1); rely=zeros(U,1); relz=zeros(U,1);
    expx=zeros(U,1); expy=zeros(U,1); expz=zeros(U,1);
    for u=1:U
        base=(u-1)*genes_per_uav;
        tcol(u) = best_ind(base + 2 + (k-1)*2);
        fcol(u) = best_ind(base + 3 + (k-1)*2);
        relx(u) = release_pts(u,k,1); rely(u) = release_pts(u,k,2); relz(u) = release_pts(u,k,3);
        expx(u) = explode_pts(u,k,1); expy(u) = explode_pts(u,k,2); expz(u) = explode_pts(u,k,3);
    end
    Tparams = [Tparams, table(tcol, fcol, relx, rely, relz, expx, expy, expz, ...
        'VariableNames', {sprintf('t%d_rel',k), sprintf('f%d',k), ...
                          sprintf('rel%d_x',k), sprintf('rel%d_y',k), sprintf('rel%d_z',k), ...
                          sprintf('exp%d_x',k), sprintf('exp%d_y',k), sprintf('exp%d_z',k)})];
end

Tcov_m = table((1:3)', cov_per_missile', 'VariableNames', {'Missile','coverage_s_11lines'});
Tsum   = table(glob_cov, 'VariableNames', {'global_coverage_s_33lines'});

writetable(Tparams, xlsxName, 'Sheet','params');
writetable(Tcov_m,  xlsxName, 'Sheet','per_missile_11lines');
writetable(Tsum,    xlsxName, 'Sheet','global_33lines');
fprintf('Results saved to %s\n', xlsxName);
end



%% ========================================================================
%% 绘图：全局 33 条 & 各导弹 11 条 时间线
%% ========================================================================
function plot_results_GA(coverage_vecs_m, coverage_vec_global, tvec, glob_cov, cov_per_missile)
figure('Position', [100, 100, 1200, 820]);

subplot(2,2,1);
stairs(tvec, double(coverage_vec_global), 'LineWidth', 2);
ylim([-0.1, 1.1]); xlabel('时间 (s)'); ylabel('全局遮蔽(33条)');
title(sprintf('全局遮蔽时间线 (总计: %.2f s)', glob_cov)); grid on;

for m = 1:3
    subplot(2,2,m+1);
    stairs(tvec, double(coverage_vecs_m(m,:)), 'LineWidth', 1.5);
    ylim([-0.1, 1.1]); xlabel('时间 (s)'); ylabel(sprintf('M%d (11条同时)', m));
    title(sprintf('导弹 M%d 遮蔽 (总计: %.2f s)', m, cov_per_missile(m))); grid on;
end

sgtitle('GA 两阶段结果 - 33 条同时遮蔽判定（落地云仍有效；评估≤67 s）');
end



%% ========================================================================
%% 打印：当前最优 1×40 基因向量与按 UAV 可读参数
%% ========================================================================
function print_best_genes(best_ind, U, genes_per_uav)
n = numel(best_ind);
fprintf('>> Best genes (1x%d):\n   [', n);
for i = 1:n
    if i < n
        fprintf('%.6f, ', best_ind(i));
    else
        fprintf('%.6f', best_ind(i));
    end
end
fprintf(']\n');

for u = 1:U
    base = (u-1)*genes_per_uav;
    th = best_ind(base+1); v = best_ind(base+2);
    t1 = best_ind(base+3); f1 = best_ind(base+4);
    t2 = best_ind(base+5); f2 = best_ind(base+6);
    t3 = best_ind(base+7); f3 = best_ind(base+8);
    fprintf(['   UAV%02d: th=%.6f rad (%.2f deg), v=%.3f m/s, ', ...
             't1=%.3f, f1=%.3f, t2=%.3f, f2=%.3f, t3=%.3f, f3=%.3f\n'], ...
             u, th, rad2deg(th), v, t1, f1, t2, f2, t3, f3);
end
end

# off_2270 系列在 dell 服务器上的测试汇总

更新日期：2026-04-17

## 1. 汇总范围

本文件只统计在 dell 服务器上已经确认存在目录的 `off_2270` 相关测试。远端检查范围目前覆盖：

- `guoyang@202.119.37.201`
  - 根目录：`/share/home/guoyang/liyihua/outer_corona_relaxation`
- `guoyang2@202.119.37.201`
  - 根目录：`/share/home/guoyang2/liyihua/outer_corona_relaxation`
- `guoyang2@202.119.37.251`
  - 本次扫描 `/share/home/guoyang2/liyihua` 下未返回 `off_2270` 相关目录

本次在远端确认到的测试目录共 7 个：

| 账号 | 测试名 | 远端目录 | 状态 |
|---|---|---|---|
| `guoyang` | `off_2270` | `/share/home/guoyang/liyihua/outer_corona_relaxation/off_2270` | 已多次运行，当前已有完整输出 |
| `guoyang` | `off_2270_cfl08` | `/share/home/guoyang/liyihua/outer_corona_relaxation/off_2270_cfl08` | 目录已建，但未发现正式提交结果 |
| `guoyang` | `off_2270_highstream/B2_ghost_antisym_mom1` | `/share/home/guoyang/liyihua/outer_corona_relaxation/off_2270_highstream/B2_ghost_antisym_mom1` | 正在运行 |
| `guoyang2` | `off_2270_initwind` | `/share/home/guoyang2/liyihua/outer_corona_relaxation/off_2270_initwind` | 正在运行 |
| `guoyang2` | `off_2270_cfl08` | `/share/home/guoyang2/liyihua/outer_corona_relaxation/off_2270_cfl08` | 目录已建，但未发现正式提交结果 |
| `guoyang2` | `off_2270_lts_32_linde` | `/share/home/guoyang2/liyihua/outer_corona_relaxation/off_2270_lts_32_linde` | 已运行并崩溃 |
| `guoyang2` | `off_2270_highstream/I1_ppbound` | `/share/home/guoyang2/liyihua/outer_corona_relaxation/off_2270_highstream/I1_ppbound` | 正在运行 |

说明：

- 本地目录里还存在 `off_2270_energy_split`、`off_2270_initwind_B0split`、`off_2270_highstream/B1_mom_soft`、`off_2270_highstream/I1_ppbound_B0split_ab` 等分支，但这次在上述 dell 远端扫描中没有确认到对应运行目录，因此不作为“dell 服务器现有 tests”纳入正文。

## 2. 测试对比总表

| 账号 | 测试 | 初始条件 | 边界条件 | 相对上一个测试的改动 | 运行结果 |
|---|---|---|---|---|---|
| `guoyang` | `off_2270` | 从 `../relaxation_230502_rho1/data/polytropic0070.dat` fresh start；磁场来自 `initial/OFF_combined_lmax10_q1.008_nr600.bin`；速度/密度/压强用 Parker 风背景初始化 | `r=min` 为 `special`；`r=max` 为 `special`；`theta` 两端 `pole`；`phi` 两端 `periodic` | 基线案例 | 已跑到 `t=2.0542`，已有 `off0000.dat` 到 `off0202.dat` |
| `guoyang` | `off_2270_cfl08` | 与 `off_2270` 相同，仍从 `polytropic0070.dat` fresh start | 与 `off_2270` 相同，`mod_usr.t` 未改 | `courantpar: 0.5 -> 0.8`；核心数提高到 `1280`；保持 fresh start | 未发现 `out` / `data/out.txt` / 数据输出，视为未正式运行 |
| `guoyang2` | `off_2270_initwind` | OFF 磁场 + polytropic Parker 背景；当前活跃配置输出到 `data/amr_probe*` | 内边界仍用背景重置且 `mom(:)=0`；外边界仍是零梯度外流型；磁场仍 line-tied + 外推 | 相对公共基线 `off_2270`，改为 polytropic 初始化/背景；活跃诊断版同时把 `ditregrid` 改到 `1`、`dtsave_dat` 改到 `0.01`、关闭 `saveprim/autoconvert` | 正在运行；截至 2026-04-17 已到 `t=1.2189`，已有 `amr_probe0000.dat` 到 `amr_probe0121.dat` |
| `guoyang2` | `off_2270_cfl08` | 与基线 `off_2270` 相同，仍从 `polytropic0070.dat` fresh start | 与基线 `off_2270` 相同 | `courantpar: 0.5 -> 0.8`；计划 `1280` 核 | 未发现正式运行输出 |
| `guoyang2` | `off_2270_lts_32_linde` | 从 `../off_2270/data/off0000.dat` restart；OFF 磁场仍从 `initial/OFF_combined...bin` 读取 | 边界总体仍沿用基线 Parker + `mom(:)=0` + line-tied/offflow 组合 | 相对基线 `off_2270`，切换 `typedivbfix: glm -> linde`，开启 `local_timestep=.true.`，并把自定义 AMR 改成全局 `J/B` 判据 | 已崩溃；在 `step=31`、`t≈1.6847e-4` 出现 `Time step too small` 和 `small value encountered` |
| `guoyang + guoyang2` | `off_2270_highstream（B2_ghost_antisym_mom1 / I1_ppbound）` | 两者都属于 highstream 定向修正分支，均保留 OFF 初始磁场并以 polytropic Parker 类型背景为基础；`B2` 用统一背景 fresh start，`I1` 进一步加载 `pp_inner_bc_rho_Tiso_p85.bin` 做 2D map 初始化 | 两者都保留 `r=max` 外流型 special、`theta=pole`、`phi=periodic`；差别集中在 `r=min`：`B2` 用反对称 `mom(1)` ghost，`I1` 直接写入 `mom(1)=rho*vr` 的 PKSW/PP 注入 | `B2` 相对 `off_2270` 一次同时改背景风、径向动量 ghost 和 AMR 极区抑制；`I1` 相对 `off_2270_initwind` 继续把底边界从 `mom(:)=0` 改成 map 驱动的 Parker 注入 | `B2` 当前跑到 `t=0.84143`、`dt≈3.2141e-6`；`I1` 当前跑到 `t=1.2388`、`dt≈9.7663e-6`，推进更快 |

## 3. 各测试详细记录

### 3.1 `off_2270`

远端目录：

- `/share/home/guoyang/liyihua/outer_corona_relaxation/off_2270`

#### 1. 初始条件

- 初始重启文件来源是 `../relaxation_230502_rho1/data/polytropic0070.dat`。
- 磁场由二进制文件 `./initial/OFF_combined_lmax10_q1.008_nr600.bin` 读入。
- `mod_usr.t` 中使用 Parker solar wind 背景来填充初始 `rho`、`mom(1)` 和 `p`。
- 当前目录中的活跃 `amrvac.par` 已经是续跑版本，设置为从 `./data/off0147.dat` 继续。

#### 2. 边界条件

- 内边界 `r=min`：
  - `rho` 和 `p` 用 Parker 背景重置。
  - `mom(:)=0`，也就是径向和切向动量在内边界都被强制为零。
  - 磁场第一层 ghost 固定为 OFF 场，后续 ghost 做外推。
- 外边界 `r=max`：
  - `rho`、`p`、`mom(1)` 采用零梯度外推。
  - `mom(2:3)=0`。
  - `mag(2:3)=0`。
  - `Br` 按 `r^2 Br` 守恒方式外推。
- 角向边界：
  - `theta` 两端为 `pole`。
  - `phi` 两端为 `periodic`。

#### 3. 比起上一个模拟修改的地方

- 这是当前 `off_2270` 系列在 dell 上能确认到的公共基线案例，可以视为后续 `cfl08`、`B2`、`initwind`、`I1` 的参照版本。
- 它的特征是：
  - 使用 isothermal/Parker 型背景风初始化；
  - 内边界径向动量直接钉为零；
  - 使用 `glm` 散度清理；
  - 球坐标拉伸网格，`domain_nx=(600,96,192)`，`refine_max_level=2`。

#### 4. 模拟结果（跑的时间）

- 已确认至少经历了以下几段运行：
  - 第 1 段：`2026-04-01 10:59:34` 到 `2026-04-01 17:43:02`，`640` 核，LSF 记录的 wall-clock 为 `24202 sec`。
  - 第 2 段前有两次 5 秒级短启动失败/退出：
    - `2026-04-10 01:48:30` 到 `2026-04-10 01:48:35`
    - `2026-04-10 02:10:07` 到 `2026-04-10 02:10:12`
  - 第 2 段有效续跑：`2026-04-10 02:11:06` 到 `2026-04-12 23:58:12`，`960` 核，wall-clock 为 `251226 sec`。
- 当前 `out` 文件尾部显示，案例已经推进到：
  - 迭代步：`it=1917204`
  - 物理时间：`t=2.0542`
  - 时间步长：`dt≈1.06e-6`
- 当前数据文件已到 `off0202.dat`，说明至少保存了 203 个快照。
- 这个案例的主要问题不是立即崩溃，而是 `dt` 很小，推进速度偏慢。

### 3.2 `off_2270_cfl08`（`guoyang`）

远端目录：

- `/share/home/guoyang/liyihua/outer_corona_relaxation/off_2270_cfl08`

#### 1. 初始条件

- 仍然从 `../relaxation_230502_rho1/data/polytropic0070.dat` fresh start。
- `base_filename='data/cfl08'`。
- `reset_grid=.true.`，说明希望按当前网格重新开始，而不是沿用旧网格直接续跑。
- `mod_usr.t` 与 `off_2270` 保持一致，因此磁场初始化和 Parker 背景构造逻辑不变。

#### 2. 边界条件

- 与 `off_2270` 完全一致：
  - `r=min` 为 Parker 背景 + 动量置零 + 磁场固定/外推；
  - `r=max` 为外流型 special；
  - `theta` 为 `pole`；
  - `phi` 为 `periodic`。

#### 3. 比起上一个模拟修改的地方

- 相比 `off_2270`，这个分支只做了参数层面的加速尝试：
  - `courantpar` 从 `0.5` 提到 `0.8`。
  - 计划核心数从基线常见的 `640/960` 提到 `1280`。
  - 使用 fresh start，而不是接着 `off_2270` 的中间结果续跑。
- 物理模型、网格尺寸和 `mod_usr.t` 都没有改。

#### 4. 模拟结果（跑的时间）

- 当前远端目录里没有发现：
  - `out`
  - `data/out.txt`
  - `data/*.dat`
- `RUNLOG.md` 里也明确写了“截至 2026-04-01 尚未提交”。
- 因此这个案例目前更适合归类为“已准备但未正式运行”的参数测试，暂时没有可总结的运行耗时和物理结果。

### 3.3 `off_2270_highstream`（`B2_ghost_antisym_mom1` + `I1_ppbound` 合并对照）

远端目录：

- `guoyang`：`/share/home/guoyang/liyihua/outer_corona_relaxation/off_2270_highstream/B2_ghost_antisym_mom1`
- `guoyang2`：`/share/home/guoyang2/liyihua/outer_corona_relaxation/off_2270_highstream/I1_ppbound`

#### 1. 初始条件

- 这两个案例都属于 `off_2270_highstream` 高速流修正链条，目标都是针对近底边界异常高速流做定向数值/边界修正。
- 共同点：
  - 都保留 OFF 初始磁场。
  - 都已经脱离原始 `off_2270` 的等温/旧 Parker 基线，转向 polytropic Parker 类型背景。
  - 都不是单纯调 CFL，而是直接修改了内边界/AMR/初始化逻辑。
- 区别：
  - `B2_ghost_antisym_mom1`：
    - 相对 `off_2270` 直接起分支。
    - 用统一的 polytropic Parker 背景 fresh start。
    - 当前磁场路径仍指向 OFF 初始场文件。
  - `I1_ppbound`：
    - 其上游是 `off_2270_initwind`。
    - 在 polytropic 背景基础上，再加载 `pp_inner_bc_rho_Tiso_p85.bin` 的 `rhob/Tiso` 二维 map。
    - 初始化逐点按 `(theta,phi)` 解本地 Parker 风，不再是统一标量背景。

#### 2. 边界条件

- 两者共同保留的部分：
  - `r=max` 都是外流型 `special` 外边界。
  - `theta` 两端都为 `pole`。
  - `phi` 两端都为 `periodic`。
  - 外边界上都保留：`rho/p/mom(1)` 外推、`mom(2:3)=0`、`mag(2:3)=0`、`Br` 做 `r^2 Br` 守恒外推。
- 两者真正分叉的地方在 `r=min`：
  - `B2_ghost_antisym_mom1`：
    - `rho/p` 用 polytropic Parker 背景重置。
    - `mom(1)` 改为反对称 ghost 填充，并考虑拉伸网格距离比。
    - `mom(2:3)=0`。
    - 磁场第一层 ghost 固定 OFF 场，其余外推。
  - `I1_ppbound`：
    - `rho/p` 按每个 `(theta,phi)` 的 `rhob_map` / `Tiso_map` 重置。
    - `mom(1)` 不再置零，也不做反对称 ghost，而是直接写成 `rho*vr` 的 PKSW/PP 注入。
    - `mom(2:3)=0`。
    - 磁场同样保持第一层 ghost 固定 OFF 场、其余外推。

#### 3. 比起上一个模拟修改的地方

- `B2_ghost_antisym_mom1` 可以看作 highstream 路线的第一层边界修正：
  - 相对基线 `off_2270`，它一次同时改了 3 件事情：
    - 背景风从原 Parker/isothermal 切到 polytropic Parker wind。
    - 内边界径向动量从“强制为 0”改成“反对称 ghost 延拓”。
    - AMR 判据改成更简单的单层 `J/B` 判据，并在极区附近直接禁止 refine。
  - 运行控制上还把：
    - `dtsave_dat: 0.01 -> 0.02`
    - `ditregrid: 8 -> 1`
- `I1_ppbound` 则是 highstream 路线在另一条服务器链上的第二层边界修正：
  - 它不是直接从 `off_2270` 分出去，而是在 `off_2270_initwind` 基础上继续改。
  - 相对 `off_2270_initwind`，核心变化是：
    - 从“统一 polytropic 背景”切到“按 `(theta,phi)` map 指定的局地底边界背景”。
    - 把内边界从 `mom(:)=0` 改成允许 Parker 径向动量注入，即 `mom(1)=rho*vr`。
    - 初始化也同步采用同一套 map，保证初值和内边界一致。
- 合并来看，这两个 test 实际上是在试两种不同的 highstream 处置方向：
  - `B2`：保守地改 ghost 构造和 AMR。
  - `I1`：更激进地改底边界输入本身，让内边界直接注入更物理的 PKSW/PP 径向流。

#### 4. 模拟结果（跑的时间）

- `B2_ghost_antisym_mom1` 当前状态：
  - 服务器：`guoyang@202.119.37.201`
  - 作业号：`202996`
  - 作业名：`off_highVr`
  - 提交时间：`2026-04-13 16:09`
  - 核心数：`960`
  - 当前推进到：
    - `it=258058`
    - `t=0.84143`
    - `dt≈3.2141e-6`
  - 已有快照：`off_anti_sym0000.dat` 到 `off_anti_sym0042.dat`
  - 备注：`smoke_1step.out` 曾以 `MPI_ABORT` 结束，但正式长程计算仍在推进。
- `I1_ppbound` 当前状态：
  - 服务器：`guoyang2@202.119.37.201`
  - 作业号：`203967`
  - 作业名：`PP_strategy`
  - 提交时间：`2026-04-15 19:40:29`
  - 核心数：`960`
  - 当前推进到：
    - `it=118247`
    - `t=1.2388`
    - `dt≈9.7663e-6`
  - 已有快照：`data_p85/off0000.dat` 到 `data_p85/off0082.dat`
- 合并判断：
  - `B2` 已经比基线 `off_2270` 更快，但仍属于“改 ghost + 改 AMR”的保守修正。
  - `I1` 的 `dt` 明显大于 `B2`，说明“直接给内边界注入 PKSW/PP 径向动量”的方案在推进效率上更占优。
  - 也因此，这两个案例应当被看作同一条 `off_2270_highstream` 实验链中的两种策略，而不是完全独立的 tests。

### 3.4 `off_2270_initwind`

远端目录：

- `/share/home/guoyang2/liyihua/outer_corona_relaxation/off_2270_initwind`

#### 1. 初始条件

- 当前活跃运行配置是诊断版 `amrvac.par`，输出前缀为 `data/amr_probe`。
- 磁场从 `../off_2270_lts_32_linde/initial/OFF_combined_lmax10_q1.008_nr600.bin` 读取。
- 流体背景使用 `parker_solar_wind_polytropic` 和 `set_background_polytropic_state`，也就是多方 Parker 风初始化，而不是旧的等温 Parker 写法。
- 当前运行参数同时具有明显的“诊断化”特征：
  - `saveprim=.false.`
  - `autoconvert=.false.`
  - `dtsave_dat=0.01`
  - `ditregrid=1`
  - `time_max=4.0`

#### 2. 边界条件

- 内边界 `r=min`：
  - `rho` 和 `p` 用多方 Parker 背景重置。
  - `mom(:)=0`，也就是 3 个动量分量都被置零。
  - 磁场第一层 ghost 固定为 OFF 场，后续 ghost 做外推。
- 外边界 `r=max`：
  - `rho`、`p`、`mom(1)` 做零梯度外推。
  - `mom(2:3)=0`。
  - `mag(2:3)=0`。
  - `Br` 保持 `r^2 Br` 守恒外推。
- 相比原始 `off_2270`，边界的大框架没有变，重点变化主要在背景状态和 AMR 诊断配置上。

#### 3. 比起上一个模拟修改的地方

- 这个分支可以看作相对公共基线 `off_2270` 的“多方背景 + 诊断化 AMR”版本。
- 关键变化有：
  - 初始化与内边界背景都切换到 polytropic Parker 风。
  - 活跃运行配置改成 `amr_probe` 诊断输出，而不是普通 `off` 长跑输出。
  - `ditregrid` 从基线常见的 `8` 改到 `1`，便于更频繁追踪 AMR 变化。
  - 自定义 AMR 使用 `J/B` 判据，并对极区引入 `1 deg` 保护带。
- 但它仍保留了“内边界动量全部置零”这一点，因此和后续 `I1_ppbound` 的区别很明确。

#### 4. 模拟结果（跑的时间）

- `bjobs -l` 显示当前任务仍在运行：
  - 作业号：`201093`
  - 作业名：`off_polytropic`
  - 提交时间：`2026-04-09 02:21:12`
  - 核心数：`960`
- 2026-04-17 16:47 左右的 `out` 尾部显示：
  - 迭代步：`it=419869`
  - 物理时间：`t=1.2189`
  - 时间步长：`dt≈2.8471e-6`
- 当前数据文件已到 `amr_probe0121.dat`。
- 这说明该案例不是简单 smoke test，而是已经持续跑了较长时间的在线诊断分支。

### 3.5 `off_2270_cfl08`（`guoyang2`）

远端目录：

- `/share/home/guoyang2/liyihua/outer_corona_relaxation/off_2270_cfl08`

#### 1. 初始条件

- 与 `guoyang` 目录下的 `off_2270_cfl08` 一致，仍然从 `../relaxation_230502_rho1/data/polytropic0070.dat` fresh start。
- 目标仍是保留原始物理模型，只测试更激进的 CFL 与并行参数。

#### 2. 边界条件

- `mod_usr.t` 未改，边界条件与基线 `off_2270` 保持一致。

#### 3. 比起上一个模拟修改的地方

- 相比基线 `off_2270`：
  - `courantpar: 0.5 -> 0.8`
  - 计划核心数提高到 `1280`
  - 仍然是 fresh start，而不是续跑

#### 4. 模拟结果（跑的时间）

- 当前目录中未发现：
  - `out`
  - `data/out.txt`
  - `data/*.dat`
- 因此这一份 `guoyang2` 侧的 `cfl08` 与 `guoyang` 侧一样，仍属于“已准备、未正式运行”的参数分支。

### 3.6 `off_2270_lts_32_linde`

远端目录：

- `/share/home/guoyang2/liyihua/outer_corona_relaxation/off_2270_lts_32_linde`

#### 1. 初始条件

- 当前 `amrvac.par` 配置为：
  - `restart_from_file='../off_2270/data/off0000.dat'`
  - `firstprocess=.true.`
  - `reset_grid=.true.`
- 磁场从 `./initial/OFF_combined_lmax10_q1.008_nr600.bin` 读取。
- 流体背景仍是原始基线里的 Parker 写法，初始化和边界并没有像 `initwind` 那样切换到多方背景。

#### 2. 边界条件

- 内边界 `r=min`：
  - `rho/p` 按 Parker 背景给定。
  - `mom(:)=0`。
  - 磁场仍是第一层 ghost 固定为 OFF 场、后续向内外推。
- 外边界 `r=max`：
  - `rho`、`p`、`mom(1)` 零梯度外推。
  - `mom(2:3)=0`。
  - `mag(2:3)=0`。
  - `Br` 维持 `r^2 Br` 守恒外推。

#### 3. 比起上一个模拟修改的地方

- 这个分支的目标不是改物理边界，而是测试 LTS 是否能给 `off_2270` 提速。
- 相比基线 `off_2270`，主要变化是：
  - `typedivbfix: glm -> linde`
  - `local_timestep=.true.`
  - 保持 `typecourant='maxsum'`
  - 使用 restart from `off0000.dat`
  - `autoconvert=.true.`
  - `dtsave_dat=0.005`
- 代码计划里还明确说明：这一版把自定义 AMR 改成全局 `J/B` 判据，不再走原来死掉的 `inwin` 路径。

#### 4. 模拟结果（跑的时间）

- `data/out.txt` 显示至少有两次 `off_LTS` 正式尝试：
  - `2026-04-03 16:37:00` 到 `2026-04-03 16:37:17`，运行 `17 sec`
  - `2026-04-03 16:37:29` 到 `2026-04-03 16:39:01`，运行 `92 sec`
- `out` 末尾明确报错：
  - `Time step too small`
  - `small value encountered, run crash`
- 崩溃位置大约在：
  - `step=31`
  - `t=1.6846812506139e-4`
- 当前只留下很少量输出：
  - `lts_smoke0000.dat`、`lts_smoke0001.dat`
  - `lts_off0000.dat`、`lts_off0001.dat`
- 因此这个分支的结论很清楚：LTS + linde 在这套 `off_2270` 设置上没有稳定跑起来。

## 4. 当前判断

从目前能确认到的 dell 服务器测试链条看，`off_2270` 系列已经形成了两条比较清楚的演化路线：

1. `guoyang` 这条链上，`off_2270` 是公共基线，`off_2270_cfl08` 是单纯调 CFL 的加速尝试，而 `off_2270_highstream/B2_ghost_antisym_mom1` 代表“改背景风 + 改 ghost + 改 AMR”的 highstream 修正方向。
2. `guoyang2` 这条链上，`off_2270_initwind` 把背景改成多方 Parker 并长期在线诊断；`off_2270_lts_32_linde` 验证了 LTS/linde 不稳定；`off_2270_highstream/I1_ppbound` 则代表“直接改底边界输入”的 highstream 修正方向。
3. 合并来看，`off_2270_highstream` 目前已经形成两种互补策略：
   - `B2`：反对称 `mom(1)` ghost + 极区 AMR 抑制
   - `I1`：map 驱动的 PKSW/PP 内边界注入
4. 当前从推进效率看，`I1_ppbound` 更占优；从改动保守性看，`B2_ghost_antisym_mom1` 更接近基线 `off_2270`。

如果后面要继续扩展这个文档，最自然的下一步是：

- 把后续新出现的远端 `off_2270` 分支继续按同样四项模板追加；
- 为每个测试再补一列“是否解决了上一版的哪个问题”，这样更适合做连续迭代记录。

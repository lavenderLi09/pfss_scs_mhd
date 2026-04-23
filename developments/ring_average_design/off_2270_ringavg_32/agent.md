# `off_2270_ringavg_32` Ring-Average Agent Notes (Consolidated)

## 0) 目的与范围

本文档把本轮对话和当前测试目录下已有文档整合为一份“可执行、可追溯”的总说明，重点覆盖三件事：

1. Ring Average 原始论文方法（含核心公式）
2. 在 AMRVAC 中当前实现到什么程度、改了哪些代码、参数怎么用
3. 目前参数扫描和对照测试结果（含为什么 `dt` 变化不明显）

本文件位置：  
[`agent.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/agent.md)

---

## 1) 已整合文档来源（本目录内）

以下文档内容已并入本总结：

- [`RING_AVERAGE_METHOD_SUMMARY.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/RING_AVERAGE_METHOD_SUMMARY.md)
- [`RING_AVERAGE_FULL_IMPL_SMOKE_2026-04-10.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/RING_AVERAGE_FULL_IMPL_SMOKE_2026-04-10.md)
- [`RING_AVERAGE_POLE_AMR2_TEST_2026-04-10.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/RING_AVERAGE_POLE_AMR2_TEST_2026-04-10.md)
- [`references/ACCESS_NOTE_zhang_2019_ring_average.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/references/ACCESS_NOTE_zhang_2019_ring_average.md)
- [`analysis/bottleneck_hotspots.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/analysis/bottleneck_hotspots.md)
- [`analysis/ring_debug/bottleneck_hotspots.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/analysis/ring_debug/bottleneck_hotspots.md)
- [`analysis/off_2270_summary_2026-04-07.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/analysis/off_2270_summary_2026-04-07.md)
- [`analysis/keyframe_dt_stability_report_2026-04-06.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/analysis/keyframe_dt_stability_report_2026-04-06.md)
- [`analysis/dt_physics_diagnostics.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/analysis/dt_physics_diagnostics.md)

---

## 2) Ring Average 原始论文方法（Zhang et al., JCP 2019）

### 2.1 解决的问题

球坐标 3D 中极区方位向物理尺度很小：

\[
\Delta s_\phi \sim r\sin\theta\,\Delta\phi
\]

当 \(\theta \to 0\) 或 \(\pi\) 时，\(\Delta s_\phi \to 0\)，显式 CFL 步长受最小网格限制：

\[
\Delta t \le C_{\mathrm{CFL}} \min\left(\frac{\Delta s}{\lambda_{\max}}\right)
\]

### 2.2 论文的核心思路（不是改网格，而是每步后处理）

在选定极区 ring 上，把 \(\phi\) 向细胞按 chunk（每块 \(N_C\) 个 cell）分组：

1. 对 chunk 内保守量做平均
2. 在 chunk 尺度上做单调重构（PCM/PLM/PPM）
3. 再把重构结果回填到原始细网格 cell 平均量
4. CFL 用有效网格口径重新估算

这本质是 conservative averaging-reconstruction，不是简单滤波开关。

### 2.3 流体变量的关键公式

对 chunk 平均：

\[
q_{\mathrm{mean}}=\frac{1}{N_C}\sum_{k_l=1}^{N_C}q(k_l)
\]

对 chunk 内局部坐标 \(y\in[0,1]\) 取二次型：

\[
q(y)=Ay^2+By+C
\]

约束：

\[
q(0)=q_L,\quad q(1)=q_R,\quad \int_0^1 q(y)\,dy=q_{\mathrm{mean}}
\]

得到：

\[
A=3(q_L+q_R-2q_{\mathrm{mean}})
\]
\[
B=2(3q_{\mathrm{mean}}-q_R-2q_L)
\]
\[
C=q_L
\]

对应第 \(k_l\) 个细胞平均（\(k_l=1,\dots,N_C\)）：

\[
q_r(k_l)=\frac{A}{3N_C^2}(3k_l^2-3k_l+1)+\frac{B}{2N_C}(2k_l-1)+C
\]

并保持 chunk 守恒：

\[
\sum_{k_l=1}^{N_C}q_r(k_l)=N_C\,q_{\mathrm{mean}}
\]

### 2.4 磁场处理（论文完整版本）

- Cell-centered \(B\)：可与流体同样处理，之后做散度清理
- Staggered/CT \(B\)：论文采用更复杂的 \(\delta E\) 方案更新 \(\Phi_k\)，避免直接三分量独立重构破坏 \(\nabla\cdot B=0\)

当前 AMRVAC 实现尚未做 CT 路径的完整论文版本（见第 6 节）。

### 2.5 论文的 CFL 有效口径

论文强调应基于 effective grid 重新估算步长，常见写法：

\[
\Delta t \sim \min\left(\frac{\Delta L_{\mathrm{eff}}}{V_{\max}}\right),\qquad
\Delta L_{\mathrm{eff}}=\frac{\Delta V_{\mathrm{eff}}}{\min(A_r^\*,A_\theta^\*,A_\phi^\*)}
\]

---

## 3) AMRVAC 当前实现（`amrvac3.2_ring_average`）

### 3.1 工作仓库与分支

- 代码仓库：`/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average`
- 当前分支：`ra/test-github-push`
- 历史基线提交：
  - `67a5027`：CFL-only ring-average
  - `45499da`：reconstruction/evolution path 接入
- 当前尚有未提交改动（含 `polar_ring_recon`、post-ring clean、`ds_eff=dV/min(A*)` 等）

### 3.2 参数接口（新增/扩展）

定义在 [`mod_global_parameters.t`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/mod_global_parameters.t#L37)：

- `polar_ring_average`（逻辑开关）
- `polar_ring_theta_cap_north`（北极帽弧度）
- `polar_ring_theta_cap_south`（南极帽弧度）
- `polar_ring_max_chunk`（chunk 上限）
- `polar_ring_evolve`（是否启用 averaging-reconstruction 子步骤）
- `polar_ring_recon`（`auto|pcm|plm|ppm`）

在 [`mod_input_output.t`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/io/mod_input_output.t#L280) 接入 `paramlist`、默认值和校验：

- 默认值见 [`mod_input_output.t:472-477`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/io/mod_input_output.t#L472)
- 参数合法性校验见 [`mod_input_output.t:749-763`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/io/mod_input_output.t#L749)

### 3.3 CFL 路径实现（`mod_dt.t`）

关键位置：[`mod_dt.t`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/mod_dt.t#L365)

实现要点：

1. `get_ring_effective_scale` 计算 `hit_ring` 与 `ds_eff`
2. `get_courant_ds` 对命中 ring 的单元，用 `ds_eff` 替换 `block%ds`
3. 在 `typecourant=maxsum` 的 spherical 分支里，对命中 ring 单元采用
   \[
   cmaxtot \leftarrow c_{\max,\mathrm{cell}}/ds_{\mathrm{eff}}
   \]
4. 当前 `ds_eff` 口径：
   \[
   ds_{\mathrm{eff}}=\frac{N_{\mathrm{chunk}}\Delta V}{\min(A_r^\*,A_\theta^\*,A_\phi)}
   \]
   其中 \(A_r^\*=N_{\mathrm{chunk}}A_r,\ A_\theta^\*=N_{\mathrm{chunk}}A_\theta,\ A_\phi\) 不放大

辅助调试输出：

- `RING_DEBUG ...`
- `CFL_DEBUG ...`
- `DT_SOURCE_DEBUG ...`
- `DT_CFL_LIMITER ...`

### 3.4 演化路径实现（`mod_finite_volume.t`）

挂接点在 [`reconstruct_LR` 调用后](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/mod_finite_volume.t#L1129)，核心子程序：

- [`apply_polar_ring_average_primitive`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/mod_finite_volume.t#L1235)

实现范围：

- 仅 `spherical_3D`
- 仅 `idims==phi_`
- 仅非 staggered（cell-centered）路径

实现流程：

1. primitive \(\to\) conserved（`wcons`）
2. chunk 平均（`qbar`）
3. chunk 边界重构（PCM/PLM/PPM）
4. chunk 内二次型回填细胞平均（保守）
5. `phys_post_ring_clean`（可选）后再转回 primitive

### 3.5 `polar_ring_recon` 与 `type_limiter` 对齐

在 [`mod_finite_volume.t:1294-1320`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/mod_finite_volume.t#L1294)：

- `polar_ring_recon='pcm'/'plm'/'ppm'`：强制指定
- `polar_ring_recon='auto'`：按 `type_limiter(block%level)` 自动映射
  - `ppm/mp5/weno/teno` 家族 \(\to\) ring `PPM`
  - 其他 \(\to\) ring `PLM`

### 3.6 post-ring divB clean hook（cell-centered）

物理接口定义与缺省绑定：

- [`mod_physics.t:93`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/physics/mod_physics.t#L93)

MHD 绑定：

- [`mod_mhd_phys.t:680`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/mhd/mod_mhd_phys.t#L680)

实现：

- [`mhd_post_ring_clean`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/mhd/mod_mhd_phys.t#L5496)
- 对 `divb_linde*` 类型调用 `add_source_linde`
- 为避免几何算子越界，对操作区间做了 `ixO ± 3` 的内缩
- 首次触发输出 `POST_RING_CLEAN_DEBUG type_divb=...`

### 3.7 一个关键修正：PPM 系数符号

当前代码在 [`mod_finite_volume.t:1437`](/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average/src/mod_finite_volume.t#L1437) 使用：

\[
A=3(q_L+q_R-2q_{\mathrm{mean}})
\]

这是满足 \((q_L,q_R,q_{\mathrm{mean}})\) 约束的一致写法。  
此前若采用 \(q_L-q_R-2q_{\mathrm{mean}}\) 会导致 PPM 子路径不稳定（见第 5.4 节临时测试记录）。

---

## 4) `par` 层如何调用（当前版本）

### 4.1 最小开关组合

可直接用已有覆盖文件：

- 开：[`ringavg_full_on.par`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_full_on.par)
- 关：[`ringavg_full_off.par`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_full_off.par)
- 无 dat 的短跑：[`ringavg_nodat_1step.par`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_nodat_1step.par), [`ringavg_nodat_4steps.par`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_nodat_4steps.par)

推荐起步参数：

```fortran
&paramlist
  polar_ring_average          = .true.
  polar_ring_evolve           = .true.
  polar_ring_theta_cap_north  = 0.10d0
  polar_ring_theta_cap_south  = 0.00d0
  polar_ring_max_chunk        = 32
  polar_ring_recon            = 'auto'
/
```

### 4.2 参数设计建议（与物理案例耦合）

1. `theta_cap`：先小后大扫（`0.05 -> 0.10 -> 0.15`），避免一开始改动过宽极区
2. `max_chunk`：先 `8/16` 再 `32`，过大 chunk 会明显增耗散
3. `recon`：先 `auto`，稳健优先时可固定 `plm`，追求更高阶再试 `ppm`
4. 若目标是只验证 CFL 逻辑，可设 `polar_ring_evolve=.false.`（仅 CFL-only）

---

## 5) 参数测试结果（已保留日志 + 临时日志记录）

> 注：本节区分两类证据  
> A. 已保留在本目录的 `out_*` 文件（可直接复查）  
> B. 当时在 `/tmp` 的临时日志摘要（目前临时目录已清理）

### 5.1 Smoke（on/off）对照

日志：

- [`out_ringavg_smoke_off`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_ringavg_smoke_off)
- [`out_ringavg_smoke_on`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_ringavg_smoke_on)

结果（`it=0,1`）：

- off：`dt0=1.2131E-05`, `dt1=1.2142E-05`
- on ：`dt0=1.2131E-05`, `dt1=1.2142E-05`

结论：该工况下 ring 开关生效不等于立刻抬升全局最小 `dt`。

### 5.2 full-impl 4 步（off vs on）

日志：

- [`out_fullimpl_off_4steps`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_fullimpl_off_4steps)
- [`out_fullimpl_on_4steps`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_fullimpl_on_4steps)

关键信号：

- on 跑出 `RING_DEBUG active=T max_nchunk=32 min_ds_phi=5.3578E-04 min_ds_phi_eff=1.7145E-02`
- on 跑出 `RING_EVOLVE_DEBUG ...`（说明演化子步骤确实执行）
- `CFL_DEBUG` 中 \(\phi\) 方向项显著下降（`1.1501E+04 -> 1.4765E+03`）
- 但 `DT_CFL_LIMITER` 仍在同一 block（`igrid=209`）且主导项是径向项，导致 `dt` 基本不变

### 5.3 极区放开 AMR 到 L2 的对照（1 步）

日志：

- [`out_poleamr2_fulloff_1step`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_poleamr2_fulloff_1step)
- [`out_poleamr2_fullon_1step`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_poleamr2_fullon_1step)
- 说明文档：[`RING_AVERAGE_POLE_AMR2_TEST_2026-04-10.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/RING_AVERAGE_POLE_AMR2_TEST_2026-04-10.md)

关键信号：

- on：`min_ds_phi=1.3763E-04`, `min_ds_phi_eff=4.4042E-03`（极区更细时 ring 仍命中）
- off：`CFL_DEBUG ... dirmax_phi=4.0855E+04`
- on ：`CFL_DEBUG ... dirmax_phi=2.6419E+03`
- 两者 `dt` 仍同为 `6.0612E-06`，因为 `DT_CFL_LIMITER` 仍由同一径向主导项控制

### 5.4 绝对值修复与 fullimpl-v2 单步

日志：

- [`out_fullimpl_on_1step_absfix`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_fullimpl_on_1step_absfix)
- [`out_fullimpl_v2_on_1step_nodat`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/out_fullimpl_v2_on_1step_nodat)

关键信号：

- chunk 比值改成 `abs(ds_theta)/abs(ds_phi)` 后：
  - `RING_EVOLVE_DEBUG nchunk=16 ds_phi=+1.6027E-03 ds_phi_eff=2.5643E-02`
- 运行可稳定完成 1 步，无额外 dat 输出污染

### 5.5 临时测试记录（`/tmp`，已清理）

该阶段用于验证：

1. `polar_ring_recon = auto/pcm/plm/ppm` 是否都能命中执行
2. `auto` 是否与 `type_limiter` 对齐
3. `linde + post-ring clean` hook 是否执行

当时日志摘录结论：

- `auto + minmod limiter`：`RING_EVOLVE_DEBUG recon=PLM`，可稳定完成
- `auto + ppm limiter`：命中 `recon=PPM`；早期曾出现 `SIGILL`，定位到 PPM 二次系数实现错误；修正后可稳定完成
- `pcm/plm/ppm` 显式模式均可触发对应 `recon=...`
- `typedivbfix=linde` 时出现 `POST_RING_CLEAN_DEBUG type_divb=4`，说明 hook 被执行

---

## 6) 目前实现与论文“完整版本”的差距

当前进展（已实现）：

- CFL 有效口径（含 `dV/min(A*)`）
- 非 staggered / cell-centered 路径的 averaging-reconstruction（PCM/PLM/PPM）
- post-ring divB clean hook（Linde 家族）

尚未实现（论文完整 CT 版本核心）：

- staggered/CT 磁通量的 ring-average 保守重构
- 论文中的 \(\delta E\) 驱动 \(\Phi_k\) 更新链路
- 为严格保持 CT 局部离散 \(\nabla\cdot B=0\) 的完整电场闭合流程

结论：当前版本是“cell-centered 简化完整版 + CFL 完整口径”，不是论文 CT 完整版。

---

## 7) 为什么很多测试里 `dt` 不变（当前最清晰判断）

即使 ring 已明显削弱了极区 \(\phi\) 项，`dt` 仍可能不变，根因是：

1. 全局 `dt` 取决于最坏单元的综合约束，不只看某个方向
2. `DT_CFL_LIMITER` 显示最坏单元常被径向项主导（本例多次是 `igrid=209` 或 `2543`）
3. 因此“ring 命中且有效”与“全局 `dt` 立刻增大”不是一回事

这和 [`analysis/off_2270_summary_2026-04-07.md`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/analysis/off_2270_summary_2026-04-07.md) 里“瓶颈长期锁定在极区几何 + AMR 交互”的诊断是一致的。

---

## 8) 可直接复用的运行模板

### 8.1 基线与 ring 对照（4 步无 dat）

```bash
mpirun -np 2 ./amrvac -i amrvac.par ringavg_nodat_4steps.par ringavg_full_off.par > out_fullimpl_off_4steps 2>&1
mpirun -np 2 ./amrvac -i amrvac.par ringavg_nodat_4steps.par ringavg_full_on.par  > out_fullimpl_on_4steps 2>&1
```

### 8.2 极区 AMR L2 条件下对照（1 步）

```bash
mpirun -np 2 ./amrvac -i amrvac.par ringavg_poleamr2_test.par ringavg_nodat_1step.par ringavg_full_off.par > out_poleamr2_fulloff_1step 2>&1
mpirun -np 2 ./amrvac -i amrvac.par ringavg_poleamr2_test.par ringavg_nodat_1step.par ringavg_full_on.par  > out_poleamr2_fullon_1step 2>&1
```

### 8.3 临时输出目录（避免污染 case）

可叠加：

- [`ringavg_nodat_tmpout.par`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/ring_average_design/off_2270_ringavg_32/ringavg_nodat_tmpout.par)

---

## 9) 当前状态一句话

本目录对应的 ring-average 已从“仅 CFL 放宽”推进到“cell-centered averaging-reconstruction + 可选 post-ring clean + recon 模式对齐 limiter”的阶段；功能已可用，证据链完整，但若要达到论文 CT 完整版，下一阶段必须进入 staggered/CT 磁场与电场更新主路径。


# off_2270_initwind Notes

## 定位

- 该分支相对 `off_2270` 的核心变化，是把背景与初始化切换到 polytropic Parker 风。
- 它同时承担了近底边界 dt / AMR bottleneck 的长期诊断角色。

## 初始条件与边界

- 主运行 `amrvac.par`：
  - `base_filename='data/off'`
  - `saveprim=.true.`
  - `autoconvert=.true.`
  - `typedivbfix='glm'`
  - `ditregrid=8`
- 流体背景采用：
  - `parker_solar_wind_polytropic`
  - `set_background_polytropic_state`
- 内边界仍然保留：
  - `rho / p` 按背景重置
  - `mom(:)=0`
  - OFF 场第一层 ghost 固定、其余外推

## 诊断性质最强的变化

- 当前分支把 AMR 诊断做得比基线更强：
  - `J/B` 判据
  - 极区保护带
  - 额外 `init_check`、`convert`、`dt_cfl_debug` 等参数文件
- `off_2270_initwind_hpc/` 进一步把 `base_filename` 切到 `data/amr_probe`，把 `ditregrid` 提高到 1，用于更密集的 dt 追踪。

## 已有测试与结论

- 本地已有：
  - `data/amr_probe*.log`
  - `data/init_check.log`
  - `analysis/check_init_atmosphere.py`
- 已整理的服务器记录显示：
  - 2026-04-09 在 `guoyang2` 提交诊断长跑；
  - 曾推进到 `t=1.2189`，`dt≈2.8471e-6`；
  - 数据已到 `amr_probe0121.dat`。
- 当前结论是：
  - 这个分支更像“多方背景 + 在线 dt/AMR 诊断”的平台；
  - 但仍保留 `mom(:)=0`，因此它也是后续 `I1_ppbound` 继续修改底边界的直接上游。

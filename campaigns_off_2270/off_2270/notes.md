# off_2270 Notes

## 定位

- 这是 `campaigns_off_2270` 的公共基线 case。
- 当前本地 `mod_usr.t` 和 `amrvac.par` 已在 2026-04-11 对齐 HPC 基线版本。

## 初始条件与边界

- 主运行 `amrvac.par` 使用 `base_filename='data/off'`。
- `restart_from_file='../relaxation_230502_rho1/data/polytropic0070.dat'`，属于从既有 relaxed 状态继续发展。
- `mhd_gamma=1.05`，`typedivbfix='glm'`，`refine_max_level=2`，`ditregrid=8`。
- 内边界延续基线 Parker / OFF 场处理，是后续 highstream 与 initwind 分支的共同参照。

## 已有测试与诊断结论

- 历史长程运行日志和输出已经存在于 `data/` 与 `logs/`。
- 原始诊断记录表明，后期近底边界异常高速流更像“边界/数值主导伪加速”，而不是稳健物理结构。
- 现有 `analysis/` 已整理出高流速触发、dt bottleneck、AMR 覆盖等分析脚本，可继续复用。

## 相对后续分支的意义

- `off_2270_initwind`：在本 case 基础上切到 polytropic Parker 背景并强化 AMR 诊断。
- `off_2270_lts_32_linde`：在本 case 基础上测试 `linde + local_timestep`。
- `off_2270_highstream/*`：围绕底边界高速流问题做边界/初始化/AMR 定向修正。

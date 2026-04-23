# off_2270_initwind_hpc Notes

## 定位

- 这是 `off_2270_initwind` 的诊断化子 case，而不是新的物理主线。
- 目标是更高频率地追踪 dt-CFL 控制区、AMR 变化和近底边界残差来源。

## 主配置特征

- `amrvac.par` 使用：
  - `base_filename='data/amr_probe'`
  - `saveprim=.false.`
  - `autoconvert=.false.`
  - `time_max=4.0`
  - `fix_small_values=.true.`
  - `ditregrid=1`
- 额外提供多组 `diag_nodat_*` 与 `diag_speed_regrid8.par`，用于 very short run 或对照诊断。

## 现有结果

- `analysis/` 已经整理了完整的 dt-CFL 与 rho bottleneck 脚本。
- `data/` 下保留了：
  - `amr_probe.log`
  - 多张 `maglines*` / `vr_maglines*` 图
  - `*.mp4` 汇总视频
- 历史说明表明，该分支曾用于 `amr_probe0000~0027` 的局部 dt 控制区分析，并确认瓶颈常锁定在近内边界壳层。

## 当前结论

- 这是一个“诊断资产很完整”的 case。
- 后续如需继续定位 dt collapse 或极区数值问题，应优先复用这里的参数文件和 `analysis/` 脚本。

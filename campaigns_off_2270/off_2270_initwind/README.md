# off_2270_initwind

`off_2270_initwind` 是在 `off_2270` 基线之上，把背景与初始化切到 polytropic Parker 风并强化 AMR 诊断的 case。

- 角色：作为 `I1_ppbound` 之前的重要上游分支
- 当前主运行文件：`amrvac.par`
- 额外参数文件：`init_check.par`、`convert.par`、`dt_cfl_debug*.par`
- 当前提交脚本：`submit.lsf`（由历史 `main.job` 保留而来）
- 主要输出前缀：`data/off`
- 诊断子 case：`off_2270_initwind_hpc/`

更详细的物理和诊断背景见 `notes.md`。

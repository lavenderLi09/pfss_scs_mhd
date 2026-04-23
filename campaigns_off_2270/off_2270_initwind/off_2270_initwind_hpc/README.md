# off_2270_initwind_hpc

`off_2270_initwind_hpc` 是 `off_2270_initwind` 下面的诊断子 case。

- 角色：高频 regrid / 无 autoconvert 输出的 dt-CFL 诊断版本
- 当前主运行文件：`amrvac.par`
- 额外参数文件：`diag_nodat_1step.par`、`diag_nodat_2steps.par`、`diag_nodat_4steps.par`、`diag_speed_regrid8.par`、`ring_on.par`
- 当前提交脚本：`submit.sh`
- 主要输出前缀：`data/amr_probe`
- 分析脚本：`analysis/` 下已经有完整的 dt-CFL 诊断代码

更详细的说明见 `notes.md`。

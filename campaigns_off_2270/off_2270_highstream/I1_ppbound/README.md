# I1_ppbound

`I1_ppbound` 是 `off_2270_highstream` 链条里最明确的 PP / map 驱动底边界 case。

- 角色：从 `off_2270_initwind` 继续演化，改成按 `(theta, phi)` map 给定的 Parker / PP 内边界注入
- 当前主运行文件：`amrvac.par`
- 额外参数文件：`amrvac_p85.par`、`smoke_1step.par`、多组 `dt_cfl_debug*.par`
- 当前提交脚本：`submit.sh`
- 主要输出前缀：`data/off`、`data_p85/off`
- 辅助说明：`pp_strategy_Wang.md`

本 case 的技术背景和已有结果见 `notes.md`。

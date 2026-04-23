# off_2270_lts_32_linde

`off_2270_lts_32_linde` 是基于 AMRVAC 3.2 的 LTS + `linde` 加速试验 case。

- 角色：检验 `off_2270` 是否能通过 `local_timestep=.true.` 提速
- 当前主运行文件：`amrvac.par`
- 额外参数文件：`amrvac_smoke.par`
- 当前提交脚本：`submit.lsf`（由历史 `main.job` 保留而来）
- 主要输出前缀：`data/off`，smoke 对照用 `data/lts_off`
- 当前本地结果：保留了 `smoke.out`、`data/lts_smoke.log`、`data/lts_off.log`

详细结论见 `notes.md`。

# I1_ppbound_B0split_ab

`I1_ppbound_B0split_ab` 是 `I1_ppbound` 的 B0 splitting / short timestep 校验分支。

- 角色：比较 B0 split 与 no-split 两种磁场演化写法在 PP 边界下的短程稳定性
- 当前主运行文件：`amrvac_B0split.par`
- 额外参数文件：`amrvac.par`、`amrvac_nosplit.par`
- 当前提交脚本：`submit.sh`
- 主要输出前缀：`data/b0split_ab`、`data/nosplit_ab`
- 当前本地结果：已有多组 `b0split_*` / `nosplit_*` 日志与 `out_nosplit_ab`

具体测试背景和本地观察见 `notes.md`。

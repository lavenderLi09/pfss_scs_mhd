# off_2270_energy_split

`off_2270_energy_split` 是在 `off_2270` highstream/PP 边界链条上继续做能量方程与经验加热试验的 case。

- 角色：检验在 PP 底边界基础上引入 energy split / heating 后，短程 smoke 是否可稳定推进
- 当前主运行文件：`amrvac.par`
- 当前提交脚本：`submit.sh`
- 主要输出前缀：`data/off`
- 当前本地可见结果：多组 `smoke_1step*` / `mpirun_dbg*` 短程输出

本 case 的推断性说明较多，细节见 `notes.md`。

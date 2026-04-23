# off_2270_highstream

`off_2270_highstream` 是围绕 `off_2270` 近底边界异常高速流所展开的一组定向修正 case 集合。

- `B1_mom_soft/`
  - highstream 早期试探版，主要转向 polytropic 背景与更软的动量处理。
- `B2_ghost_antisym_mom1/`
  - 在 B1 基础上进一步测试径向动量反对称 ghost。
- `I1_ppbound/`
  - 从 `off_2270_initwind` 分出，改为 PP / map 驱动的底边界注入。
- `I1_ppbound_B0split_ab/`
  - `I1_ppbound` 的 B0 split 与短程 timestep 校验版。
- `original_code/`
  - 保留 highstream 开始前的原始参考代码快照。
- `archive/`
  - 保存特定日期的归档试验目录。

每个真正的 case 都已经单独补齐 `README.md`、`manifest.yaml`、`notes.md` 和 `logs/` 结构。

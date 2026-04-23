# B2_ghost_antisym_mom1 Notes

## 定位

- 这是 `off_2270_highstream` 路线里“保守修正”的代表 case。
- 它从 `B1_mom_soft` 继续往前走，核心目标是把内边界径向动量改成反对称 ghost，而不是直接做 PP 注入。

## 相对上一个 case 的关键变化

- 相对 `B1_mom_soft`，`mod_usr.t` 增加了针对 `mom(1)` 的反对称 ghost 构造。
- 代码里明确留下了 `B2-smoke` 样式的边界诊断打印，用于检查 ghost / physical cell 间的奇对称关系。
- 保持了 polytropic Parker 背景、`typedivbfix='glm'`、`ditregrid=1` 和诊断型 `amr_probe` 输出。

## 已有测试记录

- 本地已有：
  - `smoke_1step.out`
  - `data/amr_probe.log`
- 已整理的 dell 服务器记录显示：
  - 2026-04-13 在 `guoyang` 上正式提交；
  - 远端曾推进到 `t=0.84143`，`dt≈3.2141e-6`；
  - 本地 `smoke_1step.out` 曾出现 `MPI_ABORT`，但不代表远端长程作业已失效。

## 当前结论

- 这是一个“本地 smoke 有波动、远端长跑曾持续推进”的边界试验分支。
- 它比 `off_2270` 更快，但推进效率仍低于后续的 `I1_ppbound`。

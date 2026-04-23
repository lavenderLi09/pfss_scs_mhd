# I1_ppbound Notes

## 定位

- 这是 `off_2270_highstream` 路线里的 PP / map 驱动底边界版本。
- 相对 `off_2270_initwind`，它不再使用统一标量背景，而是按 `(theta, phi)` 读取 `rhob / Tiso` map。

## 关键变化

- 内边界 `rho / p` 按 `pp_inner_bc_rho_Tiso*.bin` 的二维表给定。
- 内边界 `mom(1)` 改成 `rho * vr` 的 Parker / PP 注入，而不是 `mom(:)=0`。
- 初始化与边界共享同一套 map，使初值与底边界更一致。
- 主 `amrvac.par` 和 `amrvac_p85.par` 都保留 `typedivbfix='glm'`，输出前缀分别为 `data/off` 与 `data_p85/off`。

## 已有本地与远端记录

- 本地已有：
  - `data/off.log`
  - `data_p85/off.log`
  - 多组 `out_smoke_*` 与 `out_dtdebug_*`
  - `analysis/` 中的 dt limiter 结果图和窗口文本
- `pp_strategy_Wang.md` 已完整记录该方案的 notebook 到 Fortran case 的实现思路。
- 已整理的 dell 服务器记录显示：
  - 2026-04-15 在 `guoyang2` 提交；
  - 曾推进到 `t=1.2388`，`dt≈9.7663e-6`；
  - 相比 `B2_ghost_antisym_mom1`，推进效率更高。

## 当前结论

- `I1_ppbound` 是 highstream 路线里“直接改底边界输入”的代表方案。
- 当前证据支持它在推进效率上优于 `B2`，因此后续派生出了 `I1_ppbound_B0split_ab`。

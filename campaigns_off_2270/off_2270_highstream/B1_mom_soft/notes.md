# B1_mom_soft Notes

## 定位

- 这是 `off_2270_highstream` 链条里的早期 case。
- 本地没有独立说明文档，因此这里的判断主要来自代码差分和后继 `B2` 的上下游关系。

## 相对 `off_2270` 的可确认变化

- `mod_usr.t` 已改用 polytropic Parker 背景生成 `rho / vr / p`。
- 初始 OFF 场读取路径改为 `../../off_2270_lts_32_linde/initial/OFF_combined_lmax10_q1.008_nr600.bin`。
- 主 `amrvac.par` 使用诊断型输出：
  - `base_filename='data/amr_probe'`
  - `saveprim=.false.`
  - `autoconvert=.false.`
  - `fix_small_values=.true.`
  - `ditregrid=1`

## 当前状态

- 本地没有看到独立的 `out` / `err` / `log` 结果。
- 因此它更像 `B2_ghost_antisym_mom1` 之前的准备版或中间态代码快照。

## 当前结论

- 该 case 适合作为 `B2` 的父版本保留。
- 若后续继续启用，建议先做一次最小 smoke，确认它与 `B2` 的差异是否仍然有比较价值。

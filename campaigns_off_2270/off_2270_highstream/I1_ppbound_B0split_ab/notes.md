# I1_ppbound_B0split_ab Notes

## 定位

- 这是 `I1_ppbound` 的直接派生 case。
- 目标是在相同 PP 内边界下，对比：
  - 总磁场直接演化；
  - `B0field=.true.` 的 B0 split 扰动演化。

## 相对 `I1_ppbound` 的关键变化

- `mod_usr.t` 增加了 `usr_set_B0 => specialset_B0`。
- 当 `B0field=.true.` 时：
  - 初始化磁场改成只演化扰动部分；
  - 边界构造会显式减掉背景 `B0`。
- `amrvac_B0split.par` 里：
  - `base_filename='data/b0split_ab'`
  - `autoconvert=.false.`
  - `it_max=2`
  - `B0field=.true.`
- `amrvac_nosplit.par` 提供对应的 no-split 对照。

## 已有本地测试

- 当前本地已有：
  - `data/b0split_ab.log`
  - `data/b0split_ab.out`
  - `data/b0split_ab.err`
  - `data/nosplit_ab.out`
  - `data/nosplit_ab.err`
  - `data/b0split_timestep.*`
  - `out_nosplit_ab`
- 从这些文件命名可以确认，这个分支主要用于 very short run 的稳定性与 timestep 对照，而不是正式长跑。

## 当前结论

- 这是一个明确的局部验证分支，适合保留为后续 B0 split 回归检查模板。
- 若后续重新启用，应继续用 `amrvac_B0split.par` 与 `amrvac_nosplit.par` 成对比较。

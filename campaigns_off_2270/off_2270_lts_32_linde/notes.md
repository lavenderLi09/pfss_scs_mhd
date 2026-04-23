# off_2270_lts_32_linde Notes

## 定位

- 该分支的目的不是改物理边界，而是测试 `off_2270` 在 AMRVAC 3.2 下是否能通过：
  - `typedivbfix='linde'`
  - `local_timestep=.true.`
  获得速度提升。

## 相对 `off_2270` 的关键变化

- `restart_from_file='../off_2270/data/off0000.dat'`
- `typedivbfix='linde'`
- `local_timestep=.true.`
- 仍保留 `mhd_gamma=1.05` 和原始 Parker / OFF 场框架。

## 已有测试结果

- 本地已有：
  - `smoke.out`
  - `data/lts_smoke.log`
  - `data/lts_off.log`
  - `out`
- 已整理的运行记录表明：
  - `amrvac_smoke.par` 的 5 step smoke 可以完成；
  - 正式运行在大约 `step=31`、`t≈1.6847e-4` 附近出现 `Time step too small`；
  - 同时伴随 `small value encountered, run crash`。

## 当前结论

- 这个分支已经给出比较明确的负结果：当前 `linde + LTS` 组合对该 case 并不稳定。
- 它仍然值得保留，因为它回答了“单靠 LTS 是否能给 `off_2270` 提速”这个问题。

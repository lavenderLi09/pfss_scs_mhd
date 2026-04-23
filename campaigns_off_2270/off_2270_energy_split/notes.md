# off_2270_energy_split Notes

## 定位

- 该分支没有找到独立的长文档说明。
- 依据 `mod_usr.t` 差分、PP map 路径和输出文件命名，当前更合理的定位是：
  - 它继承了 `I1_ppbound` 的 PP 内边界框架；
  - 又额外加入了 energy / heating 方向的试验性修改。

## 代码层面的关键变化

- `mod_usr.t` 仍读取 `pp_inner_bc_rho_Tiso_p85.bin`，说明底边界 map 驱动策略被保留。
- 相对 `I1_ppbound`，新增了 `usr_source => empirical_heating_source`。
- 文件内出现 Fan (2017) 风格经验加热参数：
  - `heat_flux_cgs`
  - `heat_Lh_cgs`
  - `heat_rmin`, `heat_rmax`
  - `heat_Bweight_on`, `heat_alpha`
- 主 `amrvac.par` 采用：
  - `mhd_gamma=1.6666667`
  - `typedivbfix='glm'`
  - `time_max=2.0`
  - `saveprim=.true.`

## 当前本地测试

- 仅看到一组本地短程 smoke / debug 输出，没有发现正式长程提交记录。
- 现有输出文件包括：
  - `out_smoke_1step_*`
  - `out_smoke_1step_fast_*`
  - `out_smoke_1step_isothermal_*`
  - `out_smoke_1step_mpirun_*`
- 这些文件名表明，当前主要工作是对 1 step 级别稳定性和参数匹配做快速排查。

## 当前结论

- 这是一个“已有短程测试、尚未形成稳定长跑结论”的方法分支。
- 后续如果继续使用，应先查看 `logs/` 中各次 smoke 输出，再决定是否保留当前 heating 配置。

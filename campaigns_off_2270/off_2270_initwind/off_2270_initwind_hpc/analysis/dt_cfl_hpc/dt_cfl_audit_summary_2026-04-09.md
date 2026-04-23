# DT/CFL 审查总结 (2026-04-09)

输入: `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/off_2270_initwind_hpc/analysis/dt_cfl_hpc/dt_cfl_per_frame.csv`
样本帧数: `2`

## 最小 dt 空间轨迹
- r 范围: `1.000643 -> 1.000643`
- theta 范围(deg): `77.812500 -> 77.812500`
- 起点: `amr_probe0000` it=0 r=1.000643 theta=77.812500deg
- 终点: `amr_probe0001` it=1 r=1.000643 theta=77.812500deg

## CFL 主导项演化
- 平均 `C_ph/C_sum`: `0.0371`
- 平均 `C_r/C_sum`: `0.9266`
- 平均 `C_th/C_sum`: `0.0363`
- `dominant_direction` 计数: `{'radial': 2}`
- `dominant_physics` 计数: `{'geometry': 2}`

## cf_max 与 |v| 关系
- `cf_max/|v|`: p50=350.2990, p90=351.6427, max=351.9786

## 与 AMR(n1/n2) 的对应关系
- n2 范围: `5536 -> 5912`
- 最小dt位置落入 1deg 极区保护带的帧数: `0/2`

## 时间步与一致性
- dt_log: min=6.052866e-06, p50=6.053981e-06, max=6.055096e-06
- dt_pred_min: min=1.211692e-05, p50=1.211721e-05, max=1.211751e-05
- dt_ratio(dt_pred/dt_log): p50=2.0015, p90=2.0018

## 最小 dt 前5帧
- `amr_probe0000` it=0 dt_pred_min=1.211692e-05 (r=1.000643, theta=77.812500deg, dir=radial, phys=geometry)
- `amr_probe0001` it=1 dt_pred_min=1.211751e-05 (r=1.000643, theta=77.812500deg, dir=radial, phys=geometry)


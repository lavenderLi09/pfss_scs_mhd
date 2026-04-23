# amrvac_polytropic Directory Structure

```text
amrvac_polytropic/
├── analysis/
│   包括：工作区级别的综合分析输出与汇总目录。
│
├── archive/
│   包括：历史归档分支与旧实验快照。
│   ├── 2026-03-27_bipolar_hpc_test_jb_50/
│   │   包括：`bipolar_hpc_test` 的 `jb_50` 归档版本。
│   ├── 2026-03-27_bipolar_hpc_test_rgfreeze/
│   │   包括：`bipolar_hpc_test` 的 `rgfreeze` 归档版本。
│   └── 2026-03-30_bipolar_hpc_test_jb_50_clean/
│       包括：`bipolar_hpc_test` 的 `jb_50_clean` 归档版本。
│
├── campaigns_bipo/
│   包括：双极场（bipolar）相关 AMRVAC case 的主工作区。
│   ├── analytic_bipolar_1to20/
│   │   包括：早期解析双极场 `1..20 Rsun` case，以及其 `data/` 输出目录。
│   ├── analytic_bipolar_1to20_stretched/
│   │   包括：带拉伸径向网格的解析双极场 case，以及 `analysis/`、`data/` 子目录。
│   ├── bipolar_hpc_test/
│   │   包括：双极场 HPC 变体 case，以及 `data/`、`data_jb_50/`、`data_jb_50_clean/`、`data_rgfreeze/` 子目录。
│   ├── polytropic_bipolar/
│   │   包括：参考/原始 polytropic bipolar case，以及 `data/` 子目录。
│   └── real_bipo_test/
│       包括：更大规模的 bipolar 变体分支，以及 `data_large/` 子目录。
│
├── campaigns_off_2270/
│   包括：`off_2270` 主线与其派生实验的主工作区。
│   ├── off_2270/
│   │   包括：主基线 case，以及 `analysis/`、`data/`、`hao_code/`、`hpc_baseline_2026-04-11/`、`initial/` 子目录。
│   ├── off_2270_energy_split/
│   │   包括：energy split 相关分支，以及 `archive/`、`data/` 子目录。
│   ├── off_2270_highstream/
│   │   包括：highstream 边界/初始化修正链条。
│   │   ├── B1_mom_soft/
│   │   │   包括：soft momentum 边界变体。
│   │   ├── B2_ghost_antisym_mom1/
│   │   │   包括：反对称 `mom(1)` ghost 变体，以及 `data/` 子目录。
│   │   ├── I1_ppbound/
│   │   │   包括：PP inner boundary 变体，以及 `analysis/`、`data/`、`data_p85/`、`initial/` 子目录。
│   │   ├── I1_ppbound_B0split_ab/
│   │   │   包括：`I1_ppbound` 的 B0 split 扩展变体，以及 `data/`、`initial/` 子目录。
│   │   ├── archive/
│   │   │   包括：highstream 历史归档分支。
│   │   │   └── 2026-04-16_I1_ppbound_B0split_ab_timestep/
│   │   │       包括：`I1_ppbound_B0split_ab` 的特定 timestep 归档版本。
│   │   └── original_code/
│   │       包括：highstream 修改前的原始代码参考。
│   ├── off_2270_initwind/
│   │   包括：initwind 相关分支，以及 `analysis/`、`data/`、`off_2270_initwind_hpc/` 子目录。
│   ├── off_2270_initwind_B0split/
│   │   包括：initwind 的 B0 split 变体。
│   └── off_2270_lts_32_linde/
│       包括：LTS + linde 实验分支，以及 `data/`、`initial/`、`off/` 子目录。
│
├── developments/
│   包括：开发型专题目录与方法实验目录。
│   └── ring_average_design/
│       包括：ring-average 设计与实现工作区。
│       ├── off_2270_ringavg_32/
│       │   包括：ring-average 主 case，以及 `analysis/`、`data/`、`hao_code/`、`initial/`、`references/` 子目录。
│       └── ring_average_site/
│           包括：ring-average 相关展示/文档站点内容。
│
├── docs/
│   包括：工作区级文档、目录说明与整理后的 Markdown 文档。
│
├── knowledge/
│   包括：知识整理、背景说明、长期说明性文档。
│
├── plan/
│   包括：计划、实施步骤、阶段性方案文档。
│
├── rules/
│   包括：工作区内使用的规则、约定和整理规范。
│
└── outline.md
    包括：旧版项目梳理与上下文说明文件；当前目录文档以本文件和实际目录结构共同为参考。
```

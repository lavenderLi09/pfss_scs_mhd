# analysis toolbox

本目录现在只收纳“可复用的分析代码”，不再镜像整个 case 目录。

## 任务分组

- `task_dt_amr/`
  - 定位最小 dt、AMR 覆盖、CFL/flux residual/rho bottleneck 相关分析脚本。
- `task_dat_diagnostics/`
  - 面向 `dat` 数据的 highstream、动量/β、径向线探针等诊断脚本。
- `task_init_checks/`
  - 初始化状态与背景大气检查脚本。
- `common_io/`
  - `dat` 读写和 AMR Morton 辅助模块。
- `tests/`
  - 分析脚本对应的测试文件。

## 来源规则

- 每个复制过来的 `.py` 文件顶部都注明了原始来源 case 和原文件路径。
- 当前这些脚本主要来自 `campaigns_off_2270`，但已经按“任务目标”重组，不按 case 存放。

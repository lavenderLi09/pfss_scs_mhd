# Workspace Gaps

本文件记录扫描整个 `amrvac_polytropic` 后，当前仍值得继续补齐或整理的地方。

更新时间：2026-04-20

## 1. 最高优先级

### 1.1 `campaigns_bipo/` 仍未整体标准化

当前扫描结果显示，以下 case 还没有补齐当前标准结构：

- `campaigns_bipo/analytic_bipolar_1to20`
  - 缺：`manifest.yaml`、`notes.md`、`logs/`、`figs/`、`analysis/`
- `campaigns_bipo/analytic_bipolar_1to20_stretched`
  - 缺：`manifest.yaml`、`notes.md`、`logs/`、`figs/`
- `campaigns_bipo/bipolar_hpc_test`
  - 缺：`manifest.yaml`、`notes.md`、`logs/`、`figs/`、`analysis/`
- `campaigns_bipo/polytropic_bipolar`
  - 缺：`README.md`、`manifest.yaml`、`notes.md`、`logs/`、`figs/`、`analysis/`
- `campaigns_bipo/real_bipo_test`
  - 缺：`README.md`、`manifest.yaml`、`notes.md`、`logs/`、`data/`、`figs/`、`analysis/`

结论：

- 如果下一步要继续把整个仓库统一成当前 `campaigns_off_2270` 的标准，`campaigns_bipo/` 是最需要补的部分。

## 2. 中优先级

### 2.1 reference / snapshot 目录需要明确身份

这些目录不是标准实验 case，但现在仍容易被误当成可直接继续跑的 case：

- `campaigns_off_2270/off_2270/hpc_baseline_2026-04-11`
- `campaigns_off_2270/off_2270_highstream/original_code`
- `archive/2026-03-27_bipolar_hpc_test_jb_50`
- `archive/2026-03-27_bipolar_hpc_test_rgfreeze`
- `archive/2026-03-30_bipolar_hpc_test_jb_50_clean`

建议：

- 至少补一个简短 `README.md`，明确它是：
  - baseline snapshot
  - original reference
  - archive only
- 避免后续被当成“应该补齐标准结构的 active case”。

### 2.2 `campaigns_off_2270/off_2270_highstream/` 作为 case group 可再加一层说明

当前已经有 `README.md`，但如果后续使用更频繁，可以考虑再补一个：

- `manifest.yaml` for case-group

不是必须，但有助于把“目录组”和“真正 case”区分开。

## 3. 低优先级

### 3.1 `analysis/` 工具箱还可以继续补索引

现在总 `analysis/` 已经按任务分组，但还缺：

- 每个任务子目录自己的简短 `README.md`
- 常用入口脚本与依赖关系说明

### 3.2 `docs/` 与 `knowledge/` 可以继续做更强聚合

目前已复制一部分 case 文档上收到总目录，但后续还可以继续整理成：

- 一个 `docs/campaign_summary_bipo.md`
- 一个 `docs/campaign_summary_off_2270.md`
- 一个 `knowledge/numerical_failure_patterns.md`

## 4. 当前不建议马上动的地方

- `campaigns_bipo/real_bipo_test/data_large/*`
- `campaigns_bipo/bipolar_hpc_test/data_*`
- 各类 `.ipynb_checkpoints/`

原因：

- 这些目录更像历史输出、旧实验副本或 legacy 快照。
- 直接原地规范化容易把“输出目录”和“标准 case 目录”混淆。
- 更稳妥的做法是：需要继续工作时，复制成新的标准 case。

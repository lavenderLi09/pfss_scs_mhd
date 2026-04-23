# Case Creation And Layout

## 1. 基本原则

- 最小管理单位是 `case directory`。
- 一个 case 只服务一个明确目标：
  - 基线
  - 初始化变体
  - 边界变体
  - 数值方法变体
  - smoke / diagnostic 对照
- 新工作默认从已有 case 复制，不从零散文件重新拼装。

## 2. 新 case 如何创建

### 2.1 先选父 case

- 默认选择“物理模型最接近、改动最少”的父 case。
- 不要为了方便直接从老基线复制，再一次性塞入很多互不相干的改动。
- 如果目标只是：
  - 改一个边界条件
  - 改一个 `typedivbfix`
  - 改一个初始化 map
  - 改一个输出策略
  应该从最近的相关 case 继续分叉。

### 2.2 创建方式

- 新 case 应创建在对应 `campaigns_*` 下。
- 推荐流程：
  1. 复制父 case 到新目录。
  2. 立即删除不应继承的临时输出或过时日志副本。
  3. 保留有价值的 `initial/`、参考脚本、分析脚本。
  4. 立刻补齐标准文档与目录。

### 2.3 新 case 必须具备的固定结构

- `README.md`
- `manifest.yaml`
- `notes.md`
- `mod_usr.t`
- 主 `amrvac.par`
- 相关 `.par` 文件
- `submit.lsf` 或 `submit.sh`
- `logs/`
- `data/`
- `figs/`
- `analysis/`

## 3. 标准目录含义

- `data/`
  - 原始运行输出：`.dat`、`.vtu`、`.vts`、convert 输出、运行时生成的主数据。
- `logs/`
  - 标准输出、标准错误、运行日志、调试 `out*`、作业日志副本。
- `figs/`
  - 快速图、诊断图、动画、导出 PDF。
- `analysis/`
  - case 自己的分析脚本、局部分析结果、临时分析说明。
- `initial/`
  - 初始场、边界 map、外部输入二进制。
- `archive/`
  - case 内部的历史小分支、一次性保留快照。

## 4. README 要写什么

- 这个 case 是什么。
- 它从哪个父 case 复制而来。
- 这次只改哪一类东西。
- 主运行 `.par` 是哪个。
- 提交脚本入口是哪个。
- 主要输出前缀是什么。

## 5. manifest 至少要写什么

- `name`
- `campaign`
- `parent_case`
- `role`
- `amrvac_dir`
- `main_par`
- `extra_par_files`
- `submit_script`
- `output_prefixes`
- `analysis_dir`
- `figs_dir`
- `logs_dir`
- `notes_file`

## 6. legacy case 的处理方式

仓库扫描显示，旧目录里存在一些不完全符合当前标准结构的情况，例如：

- `campaigns_bipo/real_bipo_test/data_large/*`
- `campaigns_bipo/bipolar_hpc_test/data_*`
- `campaigns_off_2270/off_2270_highstream/original_code/`

对这些目录的规则是：

- 默认不强行原地规范化。
- 若要继续开展新实验，先复制成新的标准 case。
- legacy 目录只作为参考、归档或对照输入，不作为新实验主目录。

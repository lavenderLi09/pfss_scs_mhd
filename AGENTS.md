# amrvac_polytropic Agents Guide

## 1. 这个仓库是什么

`amrvac_polytropic` 是同一个 MHD 模型的演化型项目工作区。

- 它不是一个单一 case，而是一组沿着同一物理问题不断派生、对照、修正的 AMRVAC case 集合。
- 这里的最小管理单位是 `case directory`。
- 每个 case 默认应自带：
  - `mod_usr.t`
  - 一个主 `amrvac.par` 和相关 `.par`
  - 提交脚本
  - 日志
  - 数据
  - 图
  - 总结文档
- agent 在这个项目里默认采用“复制已有 case，再做小改动”的方式工作，而不是从零新建全套物理模型。

## 2. 快速导航

当前工作区按功能分层组织，详细树状结构见 [DIRECTORY_STRUCTURE.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/DIRECTORY_STRUCTURE.md)。

- `campaigns_bipo/`
  - 双极场主线 case 工作区。
- `campaigns_off_2270/`
  - `off_2270` 主线及其派生实验工作区。
- `developments/`
  - 方法开发型目录，如 `ring_average_design`。
- `archive/`
  - 历史测试输入与旧分支归档。
- `docs/`
  - 工作区级说明、命令文档、整理文档。
- `knowledge/`
  - 背景知识和长期说明。
- `plan/`
  - 计划、实施步骤和方案文档。
- `analysis/`
  - 工作区级分析输出与可复用分析代码工具箱。
- `rules/`
  - 规则性材料的存放目录。
- `outline.md`
  - 旧版总纲与上下文说明。

## 3. campaigns 中每个 case 的固定结构

一个新 case 必须包含以下文件或目录：

- `README.md`
- `manifest.yaml`
- `mod_usr.t`
- `amrvac.par`
- 相关 `.par` 文件
- `submit.lsf` 或 `submit.sh`
- `logs/`
- `data/`
- `figs/`
- `analysis/`
- `notes.md`

一个新 case 至少必须更新这些内容：

- `README.md`
  - 写明这个 case 从哪个父 case 复制而来。
- `manifest.yaml`
  - 记录 case 名称、父 case、AMRVAC 外部代码目录、主运行 `.par`、提交脚本、输出前缀。
- `mod_usr.t`
  - 记录当前物理/边界/初始化代码。
- `amrvac.par` 和相关 `.par`
  - 记录当前运行参数、输出前缀、restart 设置。
- `submit.lsf` 或 `submit.sh`
  - 记录当前提交方式与作业入口。
- `notes.md`
  - 记录本次相对父 case 改了什么、为什么改、初步观察到什么。

一个新 case 一定要提前建好的目录：

- `logs/`
- `data/`
- `figs/`
- `analysis/`

## 4. agent 工作规则与 workflow

本项目的详细工作规则、case 复制规则、归档规则、输出保护规则、运行后记录规则，统一见：

- [rules/README.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/README.md)

agent 在开始任何 case 工作前，默认先读 `rules/README.md`。

## 5. 常用命令

常用命令、编译顺序、运行方式、restart、case 对比方式，统一见：

- [CASE_COMMANDS.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/docs/CASE_COMMANDS.md)

这里先记住最核心的几条：

- 每次修改 `mod_usr.t` 之后，必须重新构建。
- 构建顺序固定为：
  - `make clean`
  - `setup.pl -d=3 -arch=debug`
  - `make`
- 运行时使用外部 AMRVAC 代码目录，而不是把源码复制进每个 case。
- restart、提交、case diff 都按 `CASE_COMMANDS.md` 的模板执行。

## 6. 什么算做完一个测试

一个 case 的一次测试，至少满足以下条件才算完成：

- case 可以成功编译
- 作业可以正常提交
- 主要输出文件已生成
- 快速图已生成
- `manifest.yaml` 已更新
- `notes.md` 已记录参数变化和初步结论

如果以上条件没有全部满足，默认视为“未完成测试”，不能直接当作稳定分支。

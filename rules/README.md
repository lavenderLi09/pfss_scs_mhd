# amrvac_polytropic Rules Index

本目录整理的是适用于整个 `amrvac_polytropic` 工作区的共用规则。

这些规则面向：

- `campaigns_bipo/`
- `campaigns_off_2270/`
- 后续新增的任何 `campaigns_*`

## 适用边界

- 这些规则优先约束“新建 case”和“继续演化中的 case”。
- 仓库里存在一批 legacy 目录，尤其是：
  - `campaigns_bipo/real_bipo_test/data_large/*`
  - `campaigns_bipo/bipolar_hpc_test/data_*`
  - 一些历史 `archive/` 快照
- 对 legacy case，默认不要原地重构；需要继续做新实验时，应复制成新的标准 case 再继续。

## 文档导航

- [CASE_CREATION_AND_LAYOUT.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/CASE_CREATION_AND_LAYOUT.md)
  - 新 case 如何创建、标准目录结构、父 case 选择规则。
- [CHANGE_AND_BUILD_RULES.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/CHANGE_AND_BUILD_RULES.md)
  - 什么时候改 `mod_usr.t`，什么时候只改 `.par`，什么时候必须重编译。
- [RUN_LOG_ARCHIVE_RULES.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/RUN_LOG_ARCHIVE_RULES.md)
  - `manifest` 何时更新，日志/数据/图如何归档，finished case 是否冻结。
- [NAMING_STATUS_AND_FAILURES.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/NAMING_STATUS_AND_FAILURES.md)
  - 命名规范、状态管理、失败 case 如何记录、哪些操作禁止做。

## 最短口径

如果只记 6 条，先记下面这些：

1. 新工作默认从最接近目标的父 case 复制，不从零散文件拼装。
2. 只要改了 `mod_usr.t`，就必须 `make clean -> setup.pl -> make`。
3. 所有新运行都必须用新的 `base_filename` 或新的输出目录，禁止覆盖旧输出。
4. 每次测试后至少更新 `notes.md`；值得保留的节点同时更新 `manifest.yaml`。
5. finished case 默认冻结，不在原目录继续做新实验。
6. 删除、清空、覆盖已有 `data/`、日志、图文件，默认都是禁止操作。

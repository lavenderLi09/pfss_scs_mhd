# amrvac_polytropic Rules

详细规则已拆分到 [`rules/`](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules) 目录。

优先阅读顺序：

1. [rules/README.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/README.md)
2. [rules/CASE_CREATION_AND_LAYOUT.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/CASE_CREATION_AND_LAYOUT.md)
3. [rules/CHANGE_AND_BUILD_RULES.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/CHANGE_AND_BUILD_RULES.md)
4. [rules/RUN_LOG_ARCHIVE_RULES.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/RUN_LOG_ARCHIVE_RULES.md)
5. [rules/NAMING_STATUS_AND_FAILURES.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/rules/NAMING_STATUS_AND_FAILURES.md)

最短摘要：

- 新工作默认从最接近目标的父 case 复制。
- 新 case 至少补齐 `README.md`、`manifest.yaml`、`notes.md`、`logs/`、`data/`、`figs/`、`analysis/`。
- 改了 `mod_usr.t` 就必须重新编译。
- 只改现有运行参数时，优先只改 `.par`。
- 每次测试后更新 `notes.md`；主运行入口变化时同步更新 `manifest.yaml`。
- finished case 默认冻结，不在原目录继续做新实验。
- 禁止覆盖旧输出，禁止删除已有数据/日志/图文件，除非用户明确允许。

命令模板见 [CASE_COMMANDS.md](/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/docs/CASE_COMMANDS.md)。

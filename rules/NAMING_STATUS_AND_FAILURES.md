# Naming Status And Failures

## 1. campaign 命名规范

- campaign 目录统一使用：
  - `campaigns_<topic>`
- 现有主线：
  - `campaigns_bipo`
  - `campaigns_off_2270`

## 2. case 命名规范

- 统一使用小写 ASCII。
- 不用空格。
- 用短横线或下划线表达含义，当前仓库以下划线为主。
- 名称要能表达“相对父 case 改了什么”。

推荐模式：

- 基线：
  - `off_2270`
  - `polytropic_bipolar`
- 物理/初始化变体：
  - `off_2270_initwind`
  - `off_2270_initwind_B0split`
- 数值方法变体：
  - `off_2270_lts_32_linde`
  - `off_2270_energy_split`
- 边界链条子 case：
  - `B1_mom_soft`
  - `B2_ghost_antisym_mom1`
  - `I1_ppbound`

## 3. `.par` 命名规范

- 主运行文件优先用 `amrvac.par`。
- 对照/短跑/工具型参数文件使用功能名：
  - `smoke_1step.par`
  - `init_check.par`
  - `convert.par`
  - `resume.par`
  - `amrvac_smoke.par`
- 不要使用看不出用途的临时名字，除非立即写进 `notes.md`。

## 4. 输出前缀命名规范

- `base_filename` 应落在 `data/` 或 case 明确约定的输出子目录下。
- 名字要能看出测试目标，例如：
  - `data/off`
  - `data/amr_probe`
  - `data/lts_off`
  - `data/b0split_ab`
  - `data/nosplit_ab`

规则：

- 新实验必须用新的输出前缀。
- 不允许复用旧前缀去覆盖既有结果。

## 5. 状态建议

一个 case 建议至少有以下状态语义：

- `draft`
  - 已创建，但未完成首次有效测试。
- `running`
  - 当前有活跃测试在跑。
- `failed`
  - 最近一次关键测试失败，但 case 仍值得保留。
- `finished`
  - 已得到明确阶段性结论。
- `frozen`
  - 作为结论节点保存，不再原地做新实验。

这些状态可写进 `manifest.yaml`，或至少写进 `notes.md`。

## 6. failed case 如何记录

failed case 不等于无效 case，必须保留证据链。

至少要记录：

- 失败时使用的 `.par`
- 是否重新编译
- 最后跑到哪个 `it` / `t`
- 错误类型：
  - `Time step too small`
  - `small value encountered`
  - `MPI_ABORT`
  - 其他运行时错误
- 对应日志文件在哪里
- 是否有可复现 smoke
- 下一步建议是什么

## 7. failed case 的处理方式

- 不删除失败输出。
- 不覆盖失败现场。
- 不把失败 case 悄悄改成别的实验继续跑。
- 如果要继续试新思路，复制成新 case。

## 8. 明确禁止的操作

- 删除、清空、覆盖已有 `data/`、日志、图文件，除非用户明确允许。
- 重用旧 `base_filename` 覆盖历史结果。
- 在 finished / frozen case 原地继续做新实验。
- 改了 `mod_usr.t` 却不重编译。
- 把多个不相关目标混进同一个 case。
- 不更新 `notes.md` 就提交或结束测试。
- 不更新 `manifest.yaml` 就切换主 `.par`、主提交脚本或 `AMRVAC_DIR`。
- 把临时调试结果只留在终端，不落到 `notes.md` 或 `logs/`。
- 直接修改仓库根 `archive/` 下的历史快照作为新实验入口。

## 9. 推荐但不强制的做法

- 新建 case 后立刻创建 `README.md`、`manifest.yaml`、`notes.md`。
- 失败测试当天就把失败点写清楚。
- 把可复用分析代码同步到工作区级 `analysis/`。

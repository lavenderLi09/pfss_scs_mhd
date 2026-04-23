# Run Log Archive Rules

## 1. manifest 什么时候更新

以下任一情况发生时，`manifest.yaml` 都应更新：

- 新建 case。
- 父 case 关系确定或改变。
- 主运行 `.par` 改了。
- 新增或删除重要 `.par` 文件。
- 提交脚本入口改了。
- 外部 `AMRVAC_DIR` 改了。
- 主要输出前缀改了。
- case 状态发生变化：
  - draft
  - running
  - failed
  - finished
  - frozen

最简单的规则：

- 只要“别人重新接手这个 case 时会被误导”，就该更新 `manifest.yaml`。

## 2. `notes.md` 什么时候更新

- 每次测试后都要更新。
- 失败测试也必须更新。
- 不要等到做完一大串实验后再一次性回填。

至少要记录：

- 相对父 case 改了什么
- 本次用了哪个 `.par`
- 输出前缀是什么
- 是否重新编译
- 是否成功提交 / 运行
- 成功时看到什么
- 失败时卡在哪里

## 3. 日志、数据、图如何归档

### `data/`

- 保存原始输出：
  - `.dat`
  - `.vtu`
  - `.vts`
  - convert 结果
- 主数据优先留在 `data/`，不要搬到 `logs/`。

### `logs/`

- 保存：
  - `out`
  - `*.out`
  - `*.err`
  - `*.log`
  - 作业调度日志
  - 调试 stdout/stderr
- 若日志最初散落在 case 根目录或 `data/` 里，允许复制一份到 `logs/` 整理；默认不要删除原件。

### `figs/`

- 保存：
  - `.png`
  - `.jpg`
  - `.pdf`
  - `.mp4`
  - 其他可视化产物

### `analysis/`

- 保存 case 自己的分析脚本和局部分析说明。
- 可复用的通用分析代码，应复制到工作区级 `analysis/`，按任务类别整理。

## 4. archive 怎么用

### case 内 `archive/`

- 用于保留：
  - 某次关键 timestep 快照
  - 某轮失败但仍有参考价值的短分支
  - 已被新 case 取代但还不想移到仓库总 `archive/` 的内容

### 仓库根 `archive/`

- 用于保存：
  - 已经退出主线的大分支
  - 历史 HPC 快照
  - 有日期戳的归档版本

## 5. archive 命名规范

- 推荐：`YYYY-MM-DD_case_reason`
- 仓库里现有例子：
  - `2026-03-27_bipolar_hpc_test_jb_50`
  - `2026-04-16_I1_ppbound_B0split_ab_timestep`

## 6. finished case 是否冻结

是。默认规则：

- 一个 case 一旦成为“阶段性结论的承载目录”，就视为 finished。
- finished case 默认冻结，不继续在原目录做新实验。
- 若还想继续探索，应该复制成新的子 case。

## 7. 什么叫冻结

- 不再修改 `mod_usr.t`。
- 不再修改主运行 `.par` 去承载新实验。
- 不再覆盖它的输出目录。
- 只允许：
  - 补 `README.md`
  - 补 `manifest.yaml`
  - 补 `notes.md`
  - 整理 `logs/` / `figs/`
  - 补来源说明

## 8. 什么时候可以解冻

只有两种情况：

- 用户明确要求继续在原目录追改。
- 该 case 实际上还不是 finished，只是误标。

除此之外，默认不要解冻 finished case。

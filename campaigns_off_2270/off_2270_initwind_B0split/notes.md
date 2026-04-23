# off_2270_initwind_B0split Notes

## 定位

- 这是 `off_2270_initwind` 的直接派生分支。
- 它的核心目标很明确：在不改动 initwind 主体物理框架的前提下，测试 `B0field=.true.`。

## 相对 `off_2270_initwind` 的确定变化

- `amrvac.par` 相比父 case 只新增了：
  - `B0field=.true.`
- `mod_usr.t` 增加了：
  - `usr_set_B0 => specialset_B0`
  - 当 `B0field=.true.` 时，初始化和边界构造中的磁场会转成“背景场 + 扰动场”拆分处理。

## 当前状态

- 本地没有找到独立的 `out` / `err` / `log`。
- 因此这里更像“为后续 B0 split 检验准备好的 case 骨架”，而不是已经完成测试的稳定分支。

## 当前结论

- 若后续要研究 B0 split 对 initwind 的影响，这个 case 可以直接作为入口。
- 在没有新运行结果前，建议把它视为“已完成配置、待重新 smoke”的准备分支。

# Change And Build Rules

## 1. 什么时候改 `mod_usr.t`

只要修改涉及“代码逻辑”而不只是运行参数，就应该改 `mod_usr.t`。

典型包括：

- 初始条件的构造方式变了。
- 边界条件公式变了。
- `specialbound_usr`、`initonegrid_usr`、`special_refine_grid` 等逻辑变了。
- 新增或删除 source term。
- 新增对外部 map / bin 文件的读取。
- 新增输出变量、aux 变量、特殊诊断量。
- 新增 B0 splitting、ring average、特殊 limiter、特殊清理步骤。

一句话判断：

- 改的是“物理/数值逻辑”本身，就改 `mod_usr.t`。
- 改的是“当前逻辑的运行参数”，优先改 `.par`。

## 2. 什么时候只改 `.par`

当代码逻辑已经存在，只需要切换现有参数时，只改 `.par`。

典型包括：

- `base_filename`
- `restart_from_file`
- `time_max`
- `dtsave_dat`
- `ditregrid`
- `typedivbfix`
- `local_timestep`
- `B0field`
- `courantpar`
- `saveprim`
- `autoconvert`
- `convert_type`
- `fix_small_values`

只改 `.par` 的前提是：

- 当前 `mod_usr.t` 和外部 `AMRVAC_DIR` 已支持该参数。
- 没有新增代码变量、没有新增未实现的 namelist 项。

## 3. 什么时候既要改 `mod_usr.t` 也要改 `.par`

以下情况通常两者都要改：

- 新增了代码路径，同时要用 `.par` 开关控制。
- 新增外部输入文件，同时 `.par` 要切换输出/运行策略。
- 新边界逻辑需要配套改变 `base_filename`、restart 或 save 频率。
- 从 total-B 改到 `B0field=.true.` 这类代码和参数同时联动的方案。

## 4. 什么时候必须重编译

- 只要 `mod_usr.t` 改了，就必须重编译。
- 只改 `.par`、`README.md`、`manifest.yaml`、`notes.md` 时，不需要重编译。

## 5. 重编译固定顺序

```bash
make clean
$AMRVAC_DIR/setup.pl -d=3 -arch=debug
make
```

规则：

- 不要跳过 `setup.pl`。
- 不要假设旧 `makefile` 或旧 `amrvac` 二进制仍然可用。
- 只改了 `mod_usr.t` 但没重编译，不算有效测试。

## 6. 提交前的最小检查

### 改了 `mod_usr.t`

- 编译成功。
- 本地 smoke 至少做过一次，或确认该 case 只允许集群跑。
- `notes.md` 已记录这次代码逻辑改动。

### 只改了 `.par`

- 核对 `base_filename` 不覆盖旧输出。
- 核对 `restart_from_file`、`firstprocess`、`snapshotnext`。
- 核对提交脚本调用的是正确的 `.par`。

## 7. 不要做的事

- 改了 `mod_usr.t` 但直接沿用旧二进制。
- 同时在多个 `.par` 里改不同方向的实验，却不记录哪一个是主运行文件。
- 只改 `.par` 就声称“方法已修改”。

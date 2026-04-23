# AMRVAC Case Commands

本文档只整理本项目里最常用的 case 命令模板。

## 1. 外部 AMRVAC 代码目录

本项目的 case 默认使用外部 `AMRVAC_DIR`，而不是把完整源码放进每个 case。

当前项目现在常用两套外部代码目录：

```bash
export AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac3.2
```

用于 ring-average 开发或依赖 ring-average 改动的 case：

```bash
export AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average
```

默认情况：

- 普通 case 使用 `amrvac3.2`
- 依赖 ring-average 修改的 case 使用 `amrvac3.2_ring_average`

如果某个 case 改用别的外部 AMRVAC 目录，也必须写进该 case 的 `manifest.yaml` 和 `notes.md`。

## 2. 编译

只要修改过 `mod_usr.t`，就必须重新构建。

标准顺序：

```bash
cd /path/to/case
export AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac3.2
make clean
$AMRVAC_DIR/setup.pl -d=3 -arch=debug
make
```

不要跳过 `setup.pl`。

如果当前 case 是 ring-average 版本，则先切到对应代码树：

```bash
export AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average
```

## 3. 本地运行

最小本地运行模板：

```bash
cd /path/to/case
OMP_NUM_THREADS=1 ./amrvac -i amrvac.par
```

如果是某个专用参数文件：

```bash
cd /path/to/case
OMP_NUM_THREADS=1 ./amrvac -i smoke.par
OMP_NUM_THREADS=1 ./amrvac -i init_check.par
OMP_NUM_THREADS=1 ./amrvac -i custom_test.par
```

## 4. 集群提交

如果 case 使用 LSF：

```bash
cd /path/to/case
bsub < submit.lsf
```

如果仓库里仍沿用旧文件名，例如 `main.job`：

```bash
cd /path/to/case
bsub < main.job
```

如果 case 使用 shell 提交脚本：

```bash
cd /path/to/case
bash submit.sh
```

## 5. restart

restart 前先检查以下参数：

- `restart_from_file`
- `firstprocess`
- `snapshotnext`
- `base_filename`

典型 restart 命令：

```bash
cd /path/to/case
./amrvac -i resume.par -resume
```

或在提交脚本里：

```bash
mpirun -np 960 ./amrvac -i resume.par -resume > out
```

如果主 `.par` 本身就是 restart 版本，也可以直接：

```bash
cd /path/to/case
./amrvac -i amrvac.par
```

## 6. convert

如果需要从已有 `.dat` 做 convert：

```bash
cd /path/to/case
./amrvac -i convert.par -convert
```

使用前先确认：

- `restart_from_file`
- `convert_type`
- 输出前缀

## 7. 比较两个 case

最常见的比较分三类。

### 7.1 比输入

比较 `mod_usr.t`：

```bash
diff -u /path/to/caseA/mod_usr.t /path/to/caseB/mod_usr.t
```

比较主参数文件：

```bash
diff -u /path/to/caseA/amrvac.par /path/to/caseB/amrvac.par
```

比较提交脚本：

```bash
diff -u /path/to/caseA/submit.lsf /path/to/caseB/submit.lsf
```

### 7.2 比输出目录

```bash
find /path/to/caseA/data -maxdepth 1 -type f | sort
find /path/to/caseB/data -maxdepth 1 -type f | sort
```

### 7.3 比分析结果

优先比较：

- `notes.md`
- `analysis/` 下的摘要
- `figs/` 下的快速图

例如：

```bash
diff -u /path/to/caseA/notes.md /path/to/caseB/notes.md
```

## 8. 新建 case 的最小命令流程

```bash
cp -R /path/to/parent_case /path/to/new_case
cd /path/to/new_case
mkdir -p logs data figs analysis
touch README.md manifest.yaml notes.md
export AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac3.2
make clean
$AMRVAC_DIR/setup.pl -d=3 -arch=debug
make
```

之后再修改：

- `mod_usr.t`
- `amrvac.par` 和相关 `.par`
- `submit.lsf` 或 `submit.sh`
- `manifest.yaml`
- `notes.md`

如果新 case 依赖 ring-average 代码改动，把上面的 `AMRVAC_DIR` 改成：

```bash
export AMRVAC_DIR=/Users/zhaoyan/Documents/codes/amrvac3.2_ring_average
```

## 9. 什么算一次测试完成

至少满足：

- case 成功编译
- 作业可正常运行或提交
- 主要输出文件生成
- 快速图生成
- `manifest.yaml` 更新
- `notes.md` 更新

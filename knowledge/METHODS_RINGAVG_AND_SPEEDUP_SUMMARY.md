# off_2270_initwind：Ring Average 与提速方法汇总（2026-04-10）

## 1) 当前测试结论（基于本地诊断）

在当前 `off_2270_initwind_hpc` 的诊断中，`it=0,1,2` 均表现为：

- 时间步由 `CFL` 限制，而非 `phys_get_dt`
- 最小 `dt` 的 limiter 位置稳定在同一 block（`igrid=1997`）
- 方向贡献中 `r` 主导，`phi` 不是主导项

对应 CFL 结构：

\[
dt \approx \frac{C_{\text{courant}}}{\max\left(\frac{c_r}{ds_r}+\frac{c_\theta}{ds_\theta}+\frac{c_\phi}{ds_\phi}\right)}
\]

由于当前最小步长主导项是 `c_r/ds_r`，只放松极区 `ds_phi` 的 ring-average-CFL 改动不会显著抬升全局 `dt`。

## 2) 在该测试下，是否还有必要实现 Ring Average？

### 短期（当前 case）

- **优先级不高**。当前瓶颈不在 `phi`，而在 `r` 主导单元。
- ring average 仍可保留为能力储备，但不应作为本 case 的首要提速手段。

### 中长期（更高分辨率/更强极区限制）

- **仍有价值**：当 limiter 回到极区 `phi` 主导时，ring average 才会直接转化为更大 `dt`。

## 3) 更值得优先推进的提速路线（按收益/改动比）

1. AMR 开销治理（优先）
- 增大 `ditregrid`、必要时用 `itfixgrid/tfixgrid` 冻结网格
- 这通常比改核心求解器更快见效

2. 面向 `r` limiter 区的网格策略
- 针对内边界附近的最限步区域调整 refinement 策略（阈值/限级）
- 目标是先降低 `c_r/ds_r` 的峰值

3. 性能评测使用优化编译
- 避免用 debug 编译来评估 wall-time

4. 物理与边界平滑
- 若内边界局部波速尖峰过高，可优化初值/边界过渡，降低局部刚性

5. 中长期网格架构方案
- 若长期受极区结构限制，可评估 pole-free 网格（Yin-Yang / cubed-sphere）

## 4) 文献与参考链接

- Ring Average 原始论文（JCP 2019）  
  https://www.sciencedirect.com/science/article/abs/pii/S0021999118305436

- GMD 2021（ring/滤波等多策略比较）  
  https://gmd.copernicus.org/articles/14/859/2021/

- 双 FFT polar filter（示例）  
  https://arxiv.org/abs/2305.01537

- Yin-Yang overset grid  
  https://arxiv.org/abs/physics/0403123

- MPI-AMRVAC 参数文档（`ditregrid` 等）  
  https://amrvac.org/md_doc_2par.html


# PFSS and SCS model grid discussion
_Exported on 2026/3/25 at GMT+8 10:37:59 from Cursor (2.6.20)_

---

**User**

我的需求是根据一个综合磁图计算pfss，然后用pfss计算SCS模型，一直计算到10-20个太阳半径。查看 pfss_code @gychen/yihua/solar_wind/pfss_code/ 中计算pfss和scs的部分，告诉我是否支持拉伸网格，以及你认为是否需要拉伸网格。

---

**Cursor**

Let me first explore the code structure and then read the relevant source files.

Let me also check the PFSS script that actually performs the computation, to see how the radial grid is constructed.

我已经仔细阅读了 `pfss_code` 中的所有核心代码。以下是我的分析：

---

## 一、当前代码是否支持拉伸网格？

**不支持。** 当前代码在 r、θ、φ 三个方向上都使用**均匀网格**。具体证据如下：

### PFSS 部分

径向网格在 `pfss_script.py` 和 `pfss_module.py` 中都是均匀的：

```425:427:gychen/yihua/solar_wind/pfss_code/codes/pfss/pfss_script.py
            nr = n_r
            r_list = np.linspace(1, Rs, nr)
            rr, tt, pp = np.meshgrid(r_list, T_resample, P_resample, indexing='ij')
```

`pfss_module.py` 的 `get_rtp()` 同样如此：

```275:277:gychen/yihua/solar_wind/pfss_code/codes/pfss/pfss_module.py
        r_list   = np.linspace( 1     , self.Rs, self.n_r)
        t_list   = np.linspace(-90*deg, 90*deg , self.n_t+1)[:-1]+dth/2
        p_list   = np.linspace( 0     , 2*np.pi, self.n_p+1)[:-1]+dph/2
```

### SCS 部分

SCS 的径向网格也是均匀的：

```200:201:gychen/yihua/solar_wind/pfss_code/codes/pfss/scs_module.py
        r_list = np.linspace(self.Rcp,self.Rtp,Nr)
        rr,tt,pp = np.meshgrid(r_list,t_list,p_list,indexing='ij')
```

### 插值假设了均匀网格

更关键的是，`get_Brtp()` 中的三线性插值**硬编码了均匀网格的索引映射**。例如 PFSS 中：

```435:437:gychen/yihua/solar_wind/pfss_code/codes/pfss/pfss_module.py
            ir       = (rr   - 1)/(self.Rs-1)*(self.n_r-1)
            it       = (np.pi-tt)/np.pi*(self.n_t-1)
            ip       = (pp   - 0)/2/np.pi*(self.n_p-1)
```

SCS 中也一样：

```272:273:gychen/yihua/solar_wind/pfss_code/codes/pfss/scs_module.py
            ir          = (rr-self.Rcp)/(self.Rtp-self.Rcp)*(Nr-1)
            it          = (np.pi-tt)/np.pi*(Nt-1)
```

这些公式 `(r - r_min) / (r_max - r_min) * (N-1)` 只对均匀网格成立。如果换成拉伸网格，这些索引映射会给出错误的插值结果。

此外，`geometry.py` 中的 Lame 系数计算也只取了第一个网格间距作为 dr：

```47:49:gychen/yihua/solar_wind/pfss_code/codes/pfss/geometry.py
        dr    = r[1,0,0]-r[0,0,0]
        dt    = t[0,1,0]-t[0,0,0]
        dp    = p[0,0,1]-p[0,0,0]
```

这同样隐含了均匀网格假设。

> 注意：`pfss_script.py` 有一个 `--rtp` 参数可以加载自定义网格用于**场的求值**，但后续的插值和磁力线追踪仍然假设均匀网格，所以这并不是真正的拉伸网格支持。

---

## 二、是否需要拉伸网格？

**对于你的需求（PFSS + SCS 计算到 10-20 R_sun），我认为拉伸网格是有价值的，但并非绝对必要。** 具体分析：

### 有利的一面

1. **径向范围巨大**：从 1 R_sun 到 20 R_sun，径向跨度为 19 R_sun。磁场强度按 \( \sim r^{-(l+2)} \) 衰减——在 cusp surface（~2.5 R_sun）附近变化剧烈，而在 10-20 R_sun 处已经非常平滑。均匀网格在远处浪费大量网格点，近处分辨率却不够。

2. **内存效率**：如果要在 1-20 R_sun 上保持近表面的分辨率（如 PFSS 域中 400 个点/1.5 R_sun ≈ dr=0.00375），均匀网格需要约 5000 个径向点。配合 200×400 的角向网格，单个场分量就需要 ~5000×200×400 × 8B ≈ **3 GB**，三个分量就是 ~9 GB，这对内存和计算都是巨大负担。

3. **磁力线追踪精度**：在 cusp surface 附近（PFSS→SCS 过渡区），场梯度最大，拉伸网格可以集中分辨率在这里，提高追踪精度。

### 不太必要的一面

1. **PFSS 和 SCS 的解是解析的**：两个模型都是基于球谐函数展开，解析地计算 Alm/Blm 和 glm/hlm 系数，然后在网格点上求值。网格只是用于**存储和插值**，不是用于解 PDE。理论上你可以在任意点直接从球谐系数计算场值（代码中已有 `method='harmonics'` 选项），完全绕过网格。

2. **SCS 的 lmax 很低**（默认只有 10）：低阶球谐意味着场本身就很光滑，即使用中等分辨率的均匀网格也能较好表示。

3. **替代方案**：你可以不用拉伸网格，而是在 PFSS 和 SCS 域分别使用不同分辨率的均匀网格（代码已经支持这一点——`Nrtp_scs` 参数），在 SCS 域使用较少的径向点即可。

### 我的建议

**如果你的主要目标是：**
- 计算 PFSS+SCS 并导出磁场 → 现有代码基本够用，只需调整 `Nrtp_scs` 参数
- 进行高精度的磁力线追踪 → 建议添加拉伸网格支持（或直接使用 `harmonics` 方法逐点计算，避免插值误差）
- 需要高分辨率且省内存 → 拉伸网格价值最大

如果你决定要添加拉伸网格支持，需要修改的地方包括：
1. `get_rtp()` 方法：生成非均匀 r 数组（如对数间距或双曲正切间距）
2. `get_Brtp()` 方法：将均匀索引映射替换为二分查找（`np.searchsorted`）
3. `geometry.py`：`Lame_coefficient` 中用逐点间距替代全局常数 dr

需要我帮你实现拉伸网格的支持吗？

我的上一条回复实际上已经完整了，最后是问你是否需要我帮你实现拉伸网格的支持。让我简要总结一下结论：

---

## 总结

| 问题 | 回答 |
|------|------|
| **是否支持拉伸网格？** | **不支持。** 三个方向全部使用 `np.linspace` 均匀网格，且插值索引映射和 Lame 系数计算都硬编码了均匀间距假设。 |
| **是否需要拉伸网格？** | **有价值但非绝对必要。** 对于 10-20 R_sun 的大径向范围，拉伸网格能显著节省内存并提升 cusp surface 附近的分辨率；但由于 PFSS/SCS 都是解析解，也可以通过其他方式规避（如分域不同分辨率、直接球谐求值）。 |

## 三种可行方案（按实现难度排序）

1. **最简单**：保持现有代码不变，利用 PFSS 和 SCS 分域各自设置分辨率（`n_r` 和 `Nrtp_scs` 独立控制），SCS 域用较少径向点即可覆盖到 10-20 R_sun。
2. **中等难度**：在磁力线追踪等需要高精度的场景，使用 `method='harmonics'` 直接从球谐系数逐点求值，完全绕过网格和插值。
3. **最完整**：给代码添加拉伸网格支持，修改 `get_rtp()`、`get_Brtp()` 和 `Lame_coefficient()` 三处。

需要我帮你实现其中某个方案吗？

---

**User**

我后需要用这个pfss+scs作为磁场初始条件在AMRVAC里进行MHD松弛

---

**Cursor**

这个信息很关键，改变了对拉伸网格必要性的判断。让我看一下现有的 example notebook 来了解当前的使用流程。

我看到你的 `initial_magnetic.ipynb` 还在搭建初期。既然最终目标是为 AMRVAC MHD 松弛提供初始磁场，这让情况变得更清晰了。以下是我的综合分析：

---

## 对于 AMRVAC MHD 松弛，拉伸网格的问题其实不重要

原因是 **PFSS 和 SCS 都是解析解**，工作流程应当是：

```
综合磁图 → PFSS (Alm/Blm) → SCS (glm/hlm) → 在 AMRVAC 网格点上直接求值 → MHD 松弛
```

而**不是**：

```
综合磁图 → PFSS 存到均匀网格 → SCS 存到均匀网格 → 插值到 AMRVAC 网格 → MHD 松弛
```

第二种方式多了一层插值，既损失精度又引入 ∇·B ≠ 0 的数值误差。

## 关键注意事项

### 1. 直接在 AMRVAC 网格上求值

AMRVAC 通常使用**径向拉伸网格**（对数或幂律分布），你应该：
- 先确定 AMRVAC 的网格坐标 (r, θ, φ)
- 用 `get_Brtp(rtp, method='harmonics')` 直接在这些点上从球谐系数计算场值
- 这样完全绕过了 PFSS/SCS 代码自身网格的限制

现有代码已经支持这一点——`pfss_module.py` 的 `get_Brtp()` 有 `method='harmonics'` 选项，可以在任意 (r,θ,φ) 点上从 Alm/Blm 直接计算。但 **SCS 的 `get_Brtp()` 目前只支持 `interpolation` 方法**，没有 harmonics 方法。

### 2. ∇·B = 0 保证

这是 MHD 松弛中最关键的问题之一：
- PFSS/SCS 的解析解本身严格满足 ∇·B = 0
- 但**插值到另一个网格上后**，∇·B ≠ 0 是不可避免的
- AMRVAC 通常有散度清理机制（如 GLM 方法），但初始条件的散度越小越好
- 最佳做法是通过**矢量势 A** 来给定初始条件，然后让 AMRVAC 自己算 B = ∇×A，这样可以保证机器精度的 ∇·B = 0

### 3. PFSS→SCS 过渡区处理

在 cusp surface（~2.5 R_sun）附近，PFSS 和 SCS 的场需要平滑衔接。当前代码在 `scs_module.py` 的 `get_Brtp()` 中通过判断 r < Rcp 或 r >= Rcp 来分区，但在边界处可能存在不连续。MHD 松弛对这种不连续性很敏感。

## 推荐的工作流程

```
步骤1: PFSS 计算
   综合磁图 → 球谐展开 → 得到 Alm, Blm 系数

步骤2: SCS 计算  
   在 cusp surface 上取 PFSS 场 → 重新展开 → 得到 glm, hlm 系数

步骤3: 准备 AMRVAC 初始条件
   方案A（推荐）: 从解析解推导矢量势 A(r,θ,φ)，在 AMRVAC 网格面上求值，
                  B = curl(A) 保证 divB = 0
   方案B（简单）: 在 AMRVAC 网格点上用球谐系数直接求值 B(r,θ,φ)，
                  依赖 AMRVAC 的散度清理

步骤4: AMRVAC MHD 松弛
   从 PFSS+SCS 初始条件出发，松弛到稳态太阳风解
```

## 你需要我做什么？

基于你的场景，我建议优先做以下工作：

1. **给 SCS 添加 `method='harmonics'` 直接求值功能**——让你能在任意点（包括 AMRVAC 网格点）上直接从 glm/hlm 计算磁场，避免先存网格再插值
2. **编写 AMRVAC 初始条件生成脚本**——读取 AMRVAC 网格，调用 PFSS+SCS 求值，输出为 AMRVAC 可读格式（`.dat` 或通过 `usr_init_one_grid` 子程序）
3. **（可选）实现矢量势 A 的输出**——从 PFSS/SCS 球谐系数推导 A，保证 ∇·B = 0

你想先从哪个方向开始？

---

**User**

在当前工作目录2026_scs_relaxation下拷贝一份我们使用的pfss_code代码，并在这个基础上修改，确认SCS可以直接用glm和hlm计算磁场，并给SCS添加利用harmonics求值的功能

---

**Cursor**

代码已拷贝到 `/home/aiia/newDisk/liyihua/2026_scs_relaxation/codes/` 下。现在让我仔细分析 SCS 的 `build_SCS_Brtp` 函数，确认它可以直接用 glm/hlm 在任意点计算磁场。

现在让我分析 `build_SCS_Brtp` 函数并实现修改。

**确认 `build_SCS_Brtp` 的能力：**

`build_SCS_Brtp` 函数（第99-121行）是一个纯解析计算函数，它直接从 glm/hlm 系数和 Pnm/DPnm 球谐函数在任意 (r, θ, φ) 点上计算 SCS 磁场。公式为：

```99:121:liyihua/2026_scs_relaxation/codes/pfss/scs_module.py
def build_SCS_Brtp(rr,tt,pp,glm,hlm,lmax=10, **kwargs):
    # ... 在任意 (rr, tt, pp) 点上解析计算
    for l in range(lmax+1):
        for m in range(l+1):
            br+=(l+1)*(Rcp/rr)**(l+2)*plm*(glm[l][m]*torch.cos(m*pp)+hlm[l][m]*torch.sin(m*pp))
            bt+=-(Rcp/rr)**(l+2)*dplm*(glm[l][m]*torch.cos(m*pp)+hlm[l][m]*torch.sin(m*pp))
            bp+=(Rcp/rr)**(l+2)*plm*m/torch.sin(tt)*(glm[l][m]*torch.sin(m*pp)-hlm[l][m]*torch.cos(m*pp))
    # ...
```

这个函数完全可以在任意坐标上求值。但现有的 `scs_solver.get_Brtp()` 却**强制使用插值**，并且在第250-251行直接拒绝了 harmonics 方法：

```250:251:liyihua/2026_scs_relaxation/codes/pfss/scs_module.py
        if method!='interpolation':
            warnings.warn("`method` in SCS module supports only 'interpolation'.", UserWarning)
```

现在我来添加 harmonics 支持。需要注意 SCS 有一个极性反转的 mask——glm/hlm 是在 Br 取绝对值后拟合的，所以求值后需要根据 mask 恢复原始极性。

让我再通读一遍修改后的文件，确认逻辑正确。

全部完成。以下是所做工作的总结：

---

## 完成的修改

### 1. 拷贝代码

已将 `gychen/yihua/solar_wind/pfss_code/codes/` 完整复制到 `/home/aiia/newDisk/liyihua/2026_scs_relaxation/codes/`。

### 2. 确认 `build_SCS_Brtp` 的能力

`build_SCS_Brtp` 函数是纯解析求值，直接从 glm/hlm 系数 + 归一化连带 Legendre 函数在**任意** (r, θ, φ) 点上计算 SCS 磁场，无需预计算网格。其核心公式：

\[
B_r = \sum_{l,m} (l+1)\left(\frac{R_{cp}}{r}\right)^{l+2} \tilde{P}_l^m(\theta) \left[g_l^m \cos(m\phi) + h_l^m \sin(m\phi)\right]
\]

### 3. 给 `scs_solver.get_Brtp()` 添加 `method='harmonics'`

在 `scs_module.py` 中新增了两个功能：

**`_get_mask_at(tt, pp)`** — 将极性 mask 从 cusp surface 网格映射到任意 (θ, φ)，使用最近邻查找。因为 SCS 的 glm/hlm 是在对 Br 取绝对值后拟合的，所以需要通过 mask 恢复原始极性。

**`get_Brtp(..., method='harmonics')`** — 新的 harmonics 分支：
- **r < Rcp** (PFSS 区域)：调用父类 `pfss_solver.get_Brtp(method='harmonics')` 从 Alm/Blm 直接求值
- **r >= Rcp** (SCS 区域)：调用 `build_SCS_Brtp()` 从 glm/hlm 直接求值，再通过 `_get_mask_at` 恢复极性
- 支持 numpy 数组和标量两种输入
- 原有的 `method='interpolation'` 完全保留不变

### 使用方式

```python
from pfss.scs_module import scs_solver

ss = scs_solver(fits_file, Rtp=20.0, ...)

# 在任意网格上直接从球谐系数求值（无需预计算网格文件）
rtp = np.stack([rr, tt, pp])
Brtp = ss.get_Brtp(rtp, method='harmonics', device='cuda')
```

上一条回复实际上已经完整了。所有修改都已完成，代码无 lint 错误。

简要回顾一下，你现在可以在 `initial_magnetic.ipynb` 中这样使用：

```python
from pfss.scs_module import scs_solver

# 初始化（会自动计算 PFSS 的 Alm/Blm 和 SCS 的 glm/hlm）
ss = scs_solver(fits_file, lmax=120, Rtp=20.0, ...)

# 在 AMRVAC 网格上直接从球谐系数求值，无需预先存网格文件
rtp = np.stack([rr, tt, pp])  # AMRVAC 的网格坐标
Brtp = ss.get_Brtp(rtp, method='harmonics', device='cuda')

Br, Bt, Bp = Brtp
```

`method='harmonics'` 会自动根据 r 值分区：r < Rcp 走 PFSS（Alm/Blm），r >= Rcp 走 SCS（glm/hlm），并正确恢复极性 mask。

---

**User**

你修改了代码之后我还可以直接from pfss.scs_module import scs_solver？

---

**Cursor**

好问题。让我检查一下包的结构和安装情况。

`pfss` 包目前没有通过 pip 安装。让我看看原来的 notebook 是怎么导入的。

`pfss` 包目前**不在** Python 的搜索路径中，直接 `from pfss.scs_module import scs_solver` 会报 `ModuleNotFoundError`。你有两个选择：

**方案 A（推荐）：在本地 `codes/` 目录中以开发模式安装**

```bash
cd /home/aiia/newDisk/liyihua/2026_scs_relaxation/codes && pip install -e .
```

这样安装后，Python 全局都能 `from pfss.scs_module import scs_solver`，而且修改代码后立即生效（`-e` 是 editable 模式）。

**方案 B：在 notebook 中手动添加路径**

在 notebook 开头加：
```python
import sys
sys.path.insert(0, './codes')
```

---

要我帮你用方案 A 安装吗？上一条回复已经完整了。我给了你两个方案来解决 `from pfss.scs_module import scs_solver` 的导入问题。

你想用哪个方案？

- **方案 A**：`pip install -e ./codes` 开发模式安装（全局可用，改代码即时生效）
- **方案 B**：notebook 开头加 `sys.path.insert(0, './codes')`（不需要安装，局部生效）

---

**User**

在，pfss在zzsolar_rht这个conda环境里。但是我现在需要使用的是你修改之后的@liyihua/2026_scs_relaxation/codes/pfss/scs_module.py ，对吧？告诉我应该如何。

---

**Cursor**

对，你需要用修改后的版本。让我先确认一下现有安装的情况。

现在清楚了。当前 `zzsolar_rht` 环境中的 `pfss` 指向的是**原始版本** `/home/aiia/newDisk/gychen/yihua/solar_wind/pfss_code/codes/pfss`，不是你修改后的。

最简单的做法：在 notebook 开头用 `sys.path.insert` 让修改后的版本优先于已安装的版本。我直接帮你写到 notebook 里。

这样就行了。原理很简单：

- `sys.path.insert(0, './codes')` 把你的本地 `codes/` 目录插到 Python 搜索路径的**最前面**
- 后续 `from pfss.scs_module import scs_solver` 会**优先**找到 `./codes/pfss/`（你修改后的版本），而不是 conda 环境里安装的原始版本
- `%autoreload 2` 保证你修改 `scs_module.py` 后不用重启 kernel，自动重新加载

---

**User**

@/home/aiia/newDisk/liyihua/2026_scs_relaxation/initial_magnetic.ipynb:1-4 利用这个构造一个拉伸网格，并在pfss范围内计算pfss，利用R_hat高度的pfss磁场作为scs的顶边界，在R_hat和R_scs之间计算scs，注意要在同一个拉伸网格上计算。先询问我细节然后再执行

---

**Cursor**

让我先确认几个关键细节。

**拉伸网格与计算细节**

> **径向拉伸网格的类型？**
> - 对数拉伸 r = R_min * (R_max/R_min)^(i/N)，近表面密、远处疏
> - 双曲正弦拉伸（类似 AMRVAC 常用方案），可控制集中程度
> - 幂律拉伸 r = R_min + (R_max-R_min)*(i/N)^alpha
> - 其他方案，我会说明
>
> Answer: 告诉我几种不同拉伸网格的优劣性和数值稳定性，并从a m r v a camrvac的角度分析判断更适合哪一个

> **网格分辨率大致需要多少？（径向 × θ × φ）**
> - 低分辨率：Nr=200, Nt=180, Np=360
> - 中等分辨率：Nr=400, Nt=180, Np=360
> - 高分辨率：Nr=600, Nt=360, Np=720
> - 自定义，我会说明
>
> Answer: r方向高分辨率（1-15个太阳半径至少需要600，但不确定，我们来计算不同网格对应的最大和最小的r分辨率决定），t和p分别180和360暂时就可以

> **PFSS 球谐展开的最大阶数 lmax？**
> - 80（默认值，适合一般精度） ✔
> - 120（较高精度）
> - 150 或更高
> - 自定义，我会说明

> **在拉伸网格上求值用什么方法？**
> - harmonics 直接从球谐系数求值（精确，不需要预存网格，但较慢） ✔
> - 先在均匀网格上计算，再插值到拉伸网格（快，但有插值误差）
> - 你建议一个方案

好问题，让我详细分析三种拉伸网格在你的场景下的表现。

## 三种拉伸网格对比分析

### 参数设定

R_min = 1 R_sun，R_max = 15 R_sun，总径向跨度 14 R_sun。

### 1. 对数拉伸（Logarithmic）

$$r_i = R_{\min} \cdot \left(\frac{R_{\max}}{R_{\min}}\right)^{i/N}$$

**核心特性：相邻格子的比值恒定**，即 $r_{i+1}/r_i = \text{const} = (R_{\max}/R_{\min})^{1/N}$

| N | dr_min (r=1) | dr_max (r=15) | 膨胀比 (恒定) |
|---|---|---|---|
| 400 | 0.00677 | 0.1016 | 1.00678 |
| 600 | 0.00451 | 0.0677 | 1.00452 |
| 800 | 0.00339 | 0.0508 | 1.00339 |

**优点：**
- 膨胀比恒定，数值最稳定
- **与 AMRVAC 的 `qstretch_baselevel` 几何拉伸完全对应**（AMRVAC 的标准拉伸方式就是让每个格子比前一个大固定倍数）
- 实现最简单

**缺点：**
- 不能独立调节集中程度，完全由 R_max/R_min 决定

### 2. 双曲正弦拉伸（Sinh）

$$r = R_{\min} + (R_{\max} - R_{\min}) \cdot \frac{\sinh(q_s \cdot \xi)}{\sinh(q_s)}, \quad \xi = i/N$$

$q_s$ 控制集中程度，$q_s \to 0$ 退化为均匀网格。

| N | q_s | dr_min | dr_max | 最大膨胀比 |
|---|---|---|---|---|
| 600 | 3 | 0.00699 | 0.0703 | ~1.005 |
| 600 | 5 | 0.00157 | 0.1167 | ~1.0084 |
| 800 | 3 | 0.00524 | 0.0527 | ~1.004 |

**优点：**
- 可通过 $q_s$ 灵活调节集中程度
- 膨胀比平滑变化
- 天体物理模拟中常用

**缺点：**
- AMRVAC 没有内置 sinh 拉伸，需要自定义坐标变换
- 多一个参数需要调节

### 3. 幂律拉伸（Power-law）

$$r_i = R_{\min} + (R_{\max} - R_{\min}) \cdot (i/N)^\alpha$$

| N | α | dr_min (i=1) | dr_max | 第一格膨胀比 |
|---|---|---|---|---|
| 600 | 2.0 | 0.0000778 | 0.0467 | **2.0** |
| 600 | 1.5 | 0.00143 | 0.035 | **1.414** |

**致命缺点：** 内边界附近膨胀比 $(i+1/i)^{\alpha-1}$ 在 i=1 时极大（α=2 时为 2.0）。**完全不适合 MHD 模拟**，会导致数值不稳定和 CFL 时间步极小。

---

## 从 AMRVAC 角度的判断

**推荐：对数拉伸。** 理由：

1. **AMRVAC 原生兼容**：AMRVAC 的标准径向拉伸就是几何级数（`qstretch_baselevel`），与对数拉伸数学等价。你在这里生成的初始条件网格可以**精确匹配** AMRVAC 的计算网格，零插值误差。

2. **恒定膨胀比**：保证 CFL 条件在径向上变化平滑，不会出现局部时间步骤急剧变小的问题。

3. **N=600 时，dr_min ≈ 0.0045 R_sun**：太阳表面附近分辨率约 3100 km，对于 lmax=80 的 PFSS 解来说已经足够。

4. **幂律拉伸应排除**：内边界膨胀比过大，不适合 MHD。

5. **Sinh 可以作为备选**：如果你需要近表面更高的分辨率，可以用 sinh 拉伸（q_s=3~5），但需要在 AMRVAC 中手动实现对应的坐标变换。

---

你觉得**对数拉伸 + N=600**是否合适？还是你希望调整 N 或者考虑 sinh 方案？


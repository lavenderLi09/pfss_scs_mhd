> Copied to project knowledge from case: `off_2270_highstream`
> Original file: `/Users/zhaoyan/Documents/05_heliosphere/pfss_scs_mhd/amrvac_polytropic/campaigns_off_2270/off_2270_highstream/high_speed_streams_coronal_mhd_notes.md`
> Copied on: 2026-04-20# 强磁场区异常高速流问题整理：全球日冕 MHD 模型中的文献、机制、处理方法与对当前模拟的判断

## 0. 这份笔记在回答什么问题

这份笔记整理的是这样一个问题：在全球日冕到内日球层范围的 MHD 模拟中，尤其是计算域覆盖约 $1\sim 10$--$20\,R_\odot$、并且近太阳强磁场区存在活动区或低阶球谐重建出的强场结构时，模型有时会在强场区域上方出现非常强、很窄、且看起来不太物理的高速流。这里整理四部分内容：

1. 文献检索：哪些模型和文章明确讨论过这类问题。
2. 过往文章中的讨论与解决方法：作者们把问题归因为什么，采用了哪些数值处理。
3. 这类高速流从数值上为什么容易出现，以及哪些方法通常有效。
4. 结合我当前模拟的已知信息，对我自己的模型做出的判断，以及下一步应该检查什么。

下文中所有公式都写成可直接显示的 markdown 数学格式。

---

## 1. 文献检索：哪些文章明确讨论过强磁场区附近的异常高速流或相关数值稳定性问题

### 1.1 COCONUT 2025：太阳活动峰年 time-evolving 全球日冕模拟

这是目前最直接点名“活动区上方异常高速流”的文献之一。

文献：

- Wang et al., *Time-evolving coronal modelling of the solar maximum around the solar storms in May 2024 by COCONUT*, 2025, A&A / arXiv:2505.11990.
- 链接：<https://arxiv.org/pdf/2505.11990>

原文位置：

> “its numerical stability is still poor for solar maximum coronal simulations that involve strong complex magnetic fields with low-β values ..., and it risks developing abnormally high-speed streams above active regions (ARs) when the magnetograms are only poorly preprocessed”  
> 见 p.1, lines 81–87.

这段话非常重要，因为它直接说明：

- 问题出现的背景是 **solar maximum、强而复杂的磁场、低 $\beta$ 区域**；
- **magnetogram 预处理不充分** 会明显增加这种风险；
- 作者讨论的现象就是 **活动区上方的异常高速流**。

同文还给出了它们采用的处理：

- 对磁图做滤波，限制近表面最大磁场强度；
- 对密度和热压使用 positivity-preserving 处理；
- 对近边界速度使用 mass-flux / velocity limitation；
- 在边界直接引入与 $B^2$ 和最小 $\beta$ 相关的压力和密度修正。

相关原文位置：

> “where a filter strength ξ = 1 × 10−4 was adopted. This approach effectively limited the maximum magnetic field near the solar surface to approximately 30 Gauss during the simulation.”  
> 见 p.2, lines 246–248.

> “Under the presence of low-β regions ( β < 10−3) and rapid variation in the magnetic field at the inner boundary during the time-evolving coronal simulations, an unphysical negative thermal pressure can easily occur during the calculation when deriving the thermal pressure from energy density.”  
> 见 p.2, lines 268–270.

> “we used the 2nd-order accurate backward Euler (BDF2) method ... When the updated plasma density or thermal pressure fell beyond a threshold ...”  
> 见 p.3, lines 325–336.

边界修正的形式写成：

$$
\rho_{\mathrm{BC}} = \Upsilon_\rho \frac{B_{\mathrm{BC}}^2}{V_{A,\mathrm{BCmax}}^2} + (1-\Upsilon_\rho)\rho_s
$$

$$
p_{\mathrm{BC}} = \Upsilon_p \frac{B_{\mathrm{BC}}^2}{2}\beta_{\min} + (1-\Upsilon_p)p_s
$$

其中同文明确给出：

$$
\beta_{\min}=10^{-3}, \qquad V_{A,\mathrm{BCmax}} = 3000\ \mathrm{km\ s^{-1}}
$$

### 1.2 COCONUT 2024：time-evolving MHD corona 的 positivity-preserving 改进

文献：

- Wang et al., *Efficient magnetohydrodynamic modelling of the time-evolving solar corona...*, 2024, arXiv:2409.02043.
- 链接：<https://arxiv.org/pdf/2409.02043>

这篇文章把“异常高速流”和“增密以限制 Alfvén 速度”的关系说得非常清楚。

原文位置：

> “the active region density was carefully increased to avoid an abnormally high Alfvénic speed ... This method can prevent some non-physical negative thermal pressures and very high-speed streams developed in the domain above the active regions”  
> 见 p.2, lines 186–191.

并且进一步说明：

> “we further extended this method to all the computational domains and applied a smooth hyperbolic tangent function to gradually increase the plasma density when the local Alfvénic speed exceeds 2000 km s−1”  
> 见 p.2, lines 192–194.

这篇文章给出的逻辑很明确：

1. 强场区会带来过高的 $v_A$；
2. 过高的 $v_A$ 会破坏数值稳定性；
3. 通过在强场区增密来压低 $v_A$；
4. 这样可以同时减少非物理负热压和活动区上方异常高速流。

### 1.3 COCONUT 2026：energy decomposition 版本进一步处理低-$\beta$ 问题

文献：

- Wang et al., *COCONUT: A coronal model with an energy decomposition strategy*, 2025/2026, arXiv:2508.20423.
- 链接：<https://arxiv.org/pdf/2508.20423>

这篇文章更系统地总结了低-$\beta$、强场、负热压、磁场分解与能量分解之间的关系。

原文位置：

> “limiting its strength in the low corona to not significantly exceed 30 Gauss”  
> 见 p.8, lines 620–621.

> “allowing magnetic field strengths to exceed 100 Gauss near the active regions in the low corona during the time-evolving simulation”  
> 见 p.8, lines 622–623.

> “successfully calculating the solar maximum simulation with magnetic field strength exceeding 50 Gauss, in which COCONUT, adopting the traditional full energy equation, crashed”  
> 见 p.19, lines 930–932.

同文还直接总结了为什么低-$\beta$ 容易出问题：

> “this approach helps prevent the occurrence of non-physical negative thermal pressures in low-β regions resulting from catastrophic cancellation”  
> 见 p.2, lines 113–116.

以及它们为什么采用 decomposed energy：

> “to balance numerical stability with the conservative property of the energy equation”  
> 见 p.2, lines 118–130.

这篇文章给出的信息很有价值：

- 以前的做法常常需要把低日冕强场压到 **约 30 G** 左右；
- 当低日冕强场超过 **50 G** 时，传统 full-energy 版本会崩；
- 引入新的 energy decomposition 后，可以把低日冕活动区附近的磁场强度推进到 **超过 100 G**。

### 1.4 SIP-IFVM 2025：观测约束 CME 模型中的 extended magnetic field decomposition

文献：

- Wang et al., *SIP-IFVM: An observation-based magnetohydrodynamic model of coronal mass ejection*, 2025, arXiv:2506.19711.
- 链接：<https://arxiv.org/pdf/2506.19711>

这篇文章没有把“活动区上方异常高速流”作为主标题写出来，但它在数值稳定性背景介绍里把相关措施写得很集中。

原文位置：

> “solar maximum with a magnetic field as strong as about 30 Gauss. Some specific measures, such as switching to a magneto-frictional model when driving coronal evolution below 1.15 Rs ..., and employing a solution-preserving method that artificially increases plasma density to broaden the transition region ...”  
> 见 p.2, lines 94–95.

同文对 extended magnetic field decomposition 的描述是：

> “we adjust B01 to accommodate B1 and then reset B1 to zero at the end of each physical time step”  
> 见 p.3, lines 139–140.

这说明对于 time-evolving、低-$\beta$、强场问题，他们采取的是：

- 下层区域用更“温和”的演化模型，例如 magneto-friction；
- 提高密度、加宽 transition region；
- 通过 extended magnetic field decomposition 控制数值上真正被推进的扰动磁场 $B_1$；
- 每个物理步后重置 $B_1$，避免它积累得太大。

### 1.5 AWSoM-R 2025：强磁图下的能量输入问题

文献：

- Liu et al., *Forecasting the 8 April 2024 Total Solar Eclipse with Multiple Solar Photospheric Magnetograms*, 2025, arXiv:2503.10974.
- 链接：<https://arxiv.org/pdf/2503.10974>

AWSoM-R 的公开描述里没有像 COCONUT 那样直接写“active region 上方异常高速流”，但它明确点名了强磁图条件下边界能量输入太大带来的数值问题。

原文位置：

> “using a relatively large (SA/B)⊙ for AWSoM-R near the solar maximum can result in unstable chromospheric evaporation”  
> 见 p.3, lines 150–153.

这说明 AWSoM-R 的风险更多表述为：

- 强磁通磁图下，边界波能输入标度过大；
- 低层蒸发响应不稳定；
- 需要把 $(S_A/B)_\odot$ 随太阳活动水平和磁图总磁通重新调节。

也就是说，AWSoM 这一路更强调“加热/驱动过强导致的不稳定”，而 COCONUT/SIP-IFVM 这一路更直接强调“低-$\beta$ 强场区的数值脆弱性和异常高速流”。

### 1.6 MAS：过渡区加宽与半隐式推进

文献：

- MAS User Guide / documentation.
- 链接：<https://www.predsci.com/mas/doc/User-Guide.pdf>

MAS 文档没有像 COCONUT 一样直接点名“活动区上方异常高速流”，但它明确把 transition region 的极端陡梯度视为主要数值难点，并采用人工加宽和半隐式处理。

原文位置：

> “extremely steep temperature and density gradients arise in the transition region ... To make these calculations feasible, we have developed a way of artificially broadening the transition region.”  
> 见 p.71, lines 1980–1987.

> “Implicit and semi-implicit time differencing”  
> 见 p.10, lines 98–101.

> “Switch between fully explicit and semi-implicit algorithm”  
> 见 p.30–31, lines 710–712.

这说明 MAS 对这类问题的思路是：

- 不硬碰低层极端刚性结构；
- 通过热传导/辐射函数修改，把 transition region 人工加宽；
- 数值推进上使用隐式或半隐式来提高稳定性。

---

## 2. 过往文章中的讨论与解决方法

这一部分把上面的文献按“问题—讨论—处理”整理成统一图景。

### 2.1 文献中对问题的共同判断

不同代码的表述不一样，但它们反复指向的核心条件基本一致：

$$
\text{强磁场} + \text{低 }\beta + \text{低密度} + \text{强边界变化} 
$$

常见表现包括：

1. 活动区上方长出窄而快的高速流；
2. 热压接近 0，甚至出现非物理负热压；
3. 近边界速度过大；
4. 时间步被极端高的 Alfvén 速度压得很小；
5. Newton 迭代或整体松弛变得困难，甚至直接崩溃。

### 2.2 各文献中提出的主要解决方法

#### 2.2.1 磁图预处理 / 磁场滤波

目标：不要把过强、过尖锐的小尺度场直接送进低日冕。

典型做法：

- 降低球谐阶数；
- 对高阶模加滤波；
- 将低日冕最大磁场控制在某个更稳的量级，比如 30 G 左右。

文献依据：

- COCONUT 2025：通过滤波把近表面最大场压到约 30 G。
- COCONUT 2026：回顾旧方案中“not significantly exceed 30 Gauss”的处理。

这一类方法的思想很简单：**先把数值上最危险的场梯度削平，再逐步提高复杂度。**

#### 2.2.2 在强场区人为增密，限制 Alfvén 速度

这是文献中最常见、也最直接针对“异常高速流”的办法。

Alfvén 速度是

$$
v_A = \frac{B}{\sqrt{\mu_0\rho}}.
$$

当 $B$ 大、$\rho$ 小时，$v_A$ 会非常高。于是很多模型在强场区直接增大密度，等价于压低 $v_A$。

典型做法：

- 当局地 $v_A$ 超过某个阈值，就平滑提高密度；
- 在边界密度里显式加入 $B^2/V_A^2$ 型项；
- 对活动区单独增密。

文献依据：

- COCONUT 2024：当 $v_A$ 超过 2000 km/s 时增密；并明确说这样可以防止活动区上方出现 very high-speed streams。
- COCONUT 2025：边界密度显式写成 $\rho_{BC}\sim B^2/V_A^2$ 的混合形式。

#### 2.2.3 给热压设置最低支撑 / 保证最小 $\beta$

目标：避免强场区热压被抽得过小，导致

$$
\beta = \frac{2p}{B^2} \ll 1
$$

进入极端数值脆弱区。

典型做法：

- 在边界热压中显式加入 $\frac{B^2}{2}\beta_{\min}$ 项；
- 在 Newton 迭代或更新过程中对热压设置 floor；
- 当热压超出允许区间时回退到最近的有效值。

文献依据：

- COCONUT 2025：直接设置 $\beta_{\min}=10^{-3}$，并在更新超阈值时回退。

#### 2.2.4 近边界质量通量和速度限幅

目标：防止强边界场与内边界插值一起推出过大的局地速度。

文献依据：

- COCONUT 2025：采用 mass-flux limitation；限制 inner boundary plasma velocity 不显著超过 characteristic speeds of supergranulation or sunspots。

这一步的核心作用是：**即使强场区附近数值驱动力偏大，也不让速度在最底层一下子失控。**

#### 2.2.5 人工加宽 transition region

目标：缓和近底边界极端陡峭的温度和密度梯度，让低层不至于过于 stiff。

典型做法：

- 提高低层密度；
- 修改导热系数和辐射函数；
- 通过等效方法让 transition region 的数值尺度变宽。

文献依据：

- MAS：明确采用 artificially broadened transition region；
- SIP-IFVM 背景介绍：人工增密以 broaden the transition region；
- COCONUT 2025 背景讨论：把 increased plasma density broadened transition region 作为现有稳定化方法之一。

#### 2.2.6 使用磁场分解或扩展磁场分解

写成

$$
\mathbf{B} = \mathbf{B}_0 + \mathbf{B}_1
$$

其中 $\mathbf{B}_0$ 是大背景场，$\mathbf{B}_1$ 是真正由代码推进的小扰动场。

目标：

- 让数值求解器不要直接推进一个非常大的总磁场；
- 减小磁场离散误差；
- 在低-$\beta$ 强场区保持数值量级更温和。

文献依据：

- SIP-IFVM 2025：extended magnetic field decomposition，并在每个物理步后重置 $B_1$；
- COCONUT 2026：明确把 extended magnetic decomposition 视为 time-evolving low-$\beta$ 问题的重要方案。

#### 2.2.7 用能量分解替代传统 full-energy 形式

这是 COCONUT 2026 的核心改进之一。

目标：避免在强磁场低-$\beta$ 区域中，热压由大数相减得到而导致 catastrophic cancellation。

文献依据：

- COCONUT 2026：指出传统 full-energy 版本在 $>50$ G 情况下会崩，而 decomposed-energy 可以稳定得多。

#### 2.2.8 用更“温和”的低层模型

典型做法：

- $r<1.1\,R_\odot$ 用 1D field-line / threaded field-line；
- $r<1.15\,R_\odot$ 用 magneto-friction；
- 低层与上层 3D MHD 分区耦合。

文献依据：

- SIP-IFVM 2025：总结了 $1.15\,R_\odot$ 以下 magneto-friction 的方案；
- COCONUT 2025：背景讨论中提到 solving 1D equations below $1.1\,R_\odot$ 和 magneto-friction below $1.15\,R_\odot$。

#### 2.2.9 调整强磁图下的能量输入参数

这是 AWSoM-R 特别强调的一点。

当边界加热或波能输入与 $B$ 成正比时，强磁图本身会把能量源项推得很强。此时需要重新调节注入参数，而不是简单照搬太阳活动低时的默认值。

文献依据：

- AWSoM-R 2025：$(S_A/B)_\odot$ 过大可导致 unstable chromospheric evaporation。

---

## 3. 出现高速流的数值原因，以及有效解决方法

这一部分不再按文献逐条摘录，而是从数值机制本身来解释。

### 3.1 总体图景

这类问题最常出现在下面这个组合里：

$$
|B|\ \text{大}, \qquad \rho\ \text{小}, \qquad p_{\mathrm{th}}\ \text{小}, \qquad \beta \ll 1.
$$

此时：

$$
v_A = \frac{B}{\sqrt{\mu_0\rho}}
$$

很大，而

$$
c_s = \sqrt{\gamma p/\rho}
$$

可能很小。系统就进入“磁主导、热压几乎失去缓冲能力”的状态。

### 3.2 为什么强磁场区尤其容易在数值上长出高速流

#### 3.2.1 动量方程里，大力相消后的小残差决定运动

动量方程写成

$$
\rho \frac{d\mathbf v}{dt}
= -\nabla p + \frac{(\nabla\times\mathbf B)\times\mathbf B}{\mu_0} + \rho\mathbf g.
$$

在强磁场低日冕区域里，经常是这样一种情况：

- 压强梯度很大；
- 洛伦兹力也很大；
- 重力也不小；
- 真正决定加速度的是它们相加后剩下的那个小残差。

也就是说，真实情况可能类似于

$$
120 - 100 - 20 = 0.
$$

这个地方本来应该几乎不加速。

但如果数值离散把其中一项算偏了一点，比如磁力误差从 $-100$ 变成 $-92$，那总力就变成

$$
120 - 92 - 20 = 8.
$$

如果此时密度又很低，比如 $\rho=0.1$，那么加速度就是

$$
a = \frac{8}{0.1}=80.
$$

于是这里就会被推成一股明显的高速流。

所以第一个核心原因是：

**强场区里，真实运动由几项大力之间的微小差别决定；数值误差只要改动其中一项，这个小差别就会被明显改写。**

#### 3.2.2 磁压和磁张力对 $B$ 的误差很敏感

洛伦兹力可以理解成磁压梯度和磁张力的组合：

$$
\frac{(\nabla\times\mathbf B)\times\mathbf B}{\mu_0}
\sim
-\nabla\left(\frac{B^2}{2\mu_0}\right)
+
\frac{(\mathbf B\cdot\nabla)\mathbf B}{\mu_0}.
$$

这里直接包含 $B^2$ 和 $\nabla B$。

所以当 $B$ 很大时：

- 磁压本身大；
- 磁张力本身大；
- 同样百分比的离散误差会变成更大的绝对误差；
- 最后会直接污染净力残差。

这就是为什么活动区、强场区、球谐截断后仍然很尖的局地强场，是最容易出问题的地方。

#### 3.2.3 低密度把误差对应的加速度进一步放大

同样的净力误差，密度越小，加速度越大：

$$
a = \frac{F_{\mathrm{net}}}{\rho}.
$$

所以 $|B|$ 大和 $\rho$ 小这两件事叠在一起时，最危险。

#### 3.2.4 过高的 Alfvén 速度让显式推进和松弛都变得很困难

如果是显式推进，时间步受 CFL 限制：

$$
\Delta t \lesssim \frac{\Delta x}{|\mathbf v| + c_f},
$$

而快磁声速 $c_f$ 在低-$\beta$ 强场区会非常高。于是：

- 全局时间步被局地强场区卡死；
- 松弛需要很多很多步；
- 局地尖峰即使很窄，也会被长期保存并逐步发展。

如果是隐式方法，虽然可以跨过 CFL 限制，但方程组会更 stiff，Newton 迭代和局地修正也更容易在强场区成为瓶颈。

#### 3.2.5 如果解的是总能量方程，还会有“大数相减”的灾难性消减

总能量形式里

$$
E = \frac{p}{\gamma-1} + \frac{\rho v^2}{2} + \frac{B^2}{2\mu_0}.
$$

反推出热压时需要计算

$$
p = (\gamma-1)\left(E - \frac{\rho v^2}{2} - \frac{B^2}{2\mu_0}\right).
$$

低-$\beta$ 区域里，$\frac{B^2}{2\mu_0}$ 很大，而 $p$ 很小，于是热压是由两个很大的量相减得到的。磁能只要有一点点离散误差，热压就可能被减成接近零甚至负值。COCONUT 2026 里把这一点直接称为 **catastrophic cancellation**。

#### 3.2.6 如果解的是内能方程，问题会以另一种形式出现

如果求解的是内能方程，那么

$$
p_{\mathrm{th}} = (\gamma-1)e.
$$

此时不再是“由总能量减磁能得到热压”的问题，而是：

- 内能 $e$ 自己就被推进得很小；
- 或者边界和源项没有把强场区热压托起来；
- 或者膨胀项持续把内能抽走。

一旦出现

$$
p_{\mathrm{th}} \to 0,
$$

就有

$$
\beta = \frac{2p}{B^2} \to 0,
$$

同时

$$
c_s = \sqrt{\gamma p / \rho} \to 0.
$$

这样一来，气体几乎失去缓冲和支撑作用，流动对磁力残差、边界插值误差和任何局地数值修正都会变得非常敏感。这个时候，即使没有“总能量大数相减”的问题，也还是会长出很尖的高速流。

### 3.3 对应的有效解决方法：按优先级排序

下面按“最直接、最对症”的顺序整理。

#### 方法 A：在强场区增密，主动限制 $v_A$

这是最直接的。

目标是把

$$
v_A = \frac{B}{\sqrt{\mu_0\rho}}
$$

压下来。它通常同时改善：

- CFL 限制；
- 低密度导致的高加速度放大；
- 低层 transition region 的 stiff 程度；
- 活动区上方异常高速流。

#### 方法 B：给热压设置最低支撑，避免极低 $\beta$

例如在边界或迭代更新中引入

$$
p \gtrsim \frac{B^2}{2}\beta_{\min},
$$

其中常用量级是

$$
\beta_{\min} \sim 10^{-3}
$$

左右。这样可以让热压不会在强场区被抽得太低。

#### 方法 C：对 magnetogram 或球谐结果做适度滤波

如果最危险的高速流总是在最强磁图像素上方出现，先把这些地方的峰值削平，是最省事也最有效的做法之一。文献里控制低日冕最大场到 30 G 左右，就是这一类思路。

#### 方法 D：近边界速度或质量通量限幅

这一步主要是防止“边界一上来就把速度推出去”。尤其当你的边界磁场每步都在更新，或者采用插值磁图驱动时，很有必要。

#### 方法 E：磁场分解 / 扩展磁场分解

适合需要保留较强磁场而又不想过度滤波的情况。它的价值在于：

- 不直接推进完整的大场；
- 让求解器主要处理相对小的扰动量；
- 保持对强背景场更好的数值控制。

#### 方法 F：若用总能量方程，可改用能量分解或内能/热压更新

这是为了解决 catastrophic cancellation 的根本问题。

#### 方法 G：加宽 transition region，或把低层单独处理

适合底层特别 stiff 的情况。MAS、AWSoM、SIP-IFVM 都有类似思路。

#### 方法 H：重新审视强场区的加热/驱动源项

当加热项和 $B$ 成正比，强场区会天然被“过热”。这一步对 AWSoM 类模型尤其关键，但对任何加热强依赖 $B$ 的模型都值得检查。

---

## 4. 对我当前模拟的判断，以及下一步应该检查的地方

这一部分只基于当前已经明确的信息来判断，不引入额外推测。

### 4.1 当前已经明确的信息

根据前面的讨论，我现在已经知道：

1. 我的问题区域位于相对较强磁场的位置；
2. 那里出现了很强的高速流；
3. 我检查到这个区域的热压 $p_{\mathrm{th}}$ 很小，已经接近 0；
4. 我的 MHD 模型求解的是 **内能方程**，不是总能量方程。

这四条已经足够做出一个比较明确的数值判断。

### 4.2 对当前问题的直接判断

因为我解的是内能方程，所以这里最主要的问题 **不是** “总能量减磁能得到热压”的 catastrophic cancellation。

对我现在的模型，更合理的判断是：

$$
p_{\mathrm{th}} = (\gamma-1)e
$$

中的内能 $e$ 自己被推进得过低，导致这个位置进入了一个非常极端的低-$\beta$、低声速、磁主导状态。这样一来：

- 热压几乎失去支撑作用；
- 动量方程中的压强梯度项变得很弱；
- 运动主要由磁力和重力的微小残差决定；
- 任何局地离散误差、边界插值误差、source term 不平衡或密度修正，都会更容易被积累成明显的高速流。

因此，我现在最核心的判断是：

**当前高速流问题，很可能是“强场区内能/热压被抽得过低”之后，在极低 $\beta$ 状态下被数值误差放大的结果。**

### 4.3 最应该先检查的物理量

在出问题的位置，建议沿时间连续输出下面这些量：

$$
p_{\mathrm{th}},\ e,\ \rho,\ \beta,\ v_A,\ c_s,\ |B|.
$$

如果可能，最好还要输出：

$$
\nabla\cdot\mathbf v,
$$

因为它直接告诉我这块区域是不是在持续膨胀。如果

$$
\nabla\cdot\mathbf v > 0,
$$

那内能方程里的压缩功项会持续降低内能。

### 4.4 最应该检查的方程项

#### 4.4.1 内能方程里的各项贡献

最应该拆开的，是内能演化各项。例如如果我的方程大致形如

$$
\frac{\partial e}{\partial t} + \nabla\cdot(e\mathbf v)
= -p\,\nabla\cdot\mathbf v + Q_{\mathrm{heat}} + Q_{\mathrm{other}},
$$

那么我要检查：

1. 是不是 $-p\nabla\cdot\mathbf v$ 项持续把 $e$ 往下拉；
2. 加热项 $Q_{\mathrm{heat}}$ 在强场区是否过弱，托不住热压；
3. 是否有其它 cooling / relaxation / artificial terms 在这里持续抽走内能。

#### 4.4.2 动量方程中的受力分解

建议直接把出问题区域的动量方程三大项画出来：

$$
-\nabla p,
\qquad
\frac{(\nabla\times\mathbf B)\times\mathbf B}{\mu_0},
\qquad
\rho\mathbf g.
$$

如果能沿一条磁力线或沿径向切线输出它们的分量，会更有用。

这里要看的不是单项有多大，而是：

- 哪一项突然出现了不连续或尖峰；
- 哪一项开始和其它两项失衡；
- 高速流长出来之前，哪个残差最先偏离了近似平衡。

### 4.5 最应该检查的边界与数值设置

#### 4.5.1 边界热压和密度是否与强磁场匹配

如果边界磁场很强，但边界给定的密度和热压还按“普通背景区”的量级设定，那么强场区一开始就会进入：

$$
|B|\ \text{大},\quad p\ \text{偏小},\quad \rho\ \text{偏小}
$$

这正是最危险的组合。

所以要检查：

- 强场像素处边界 $\rho$ 是否偏小；
- 强场像素处边界 $p$ 是否偏小；
- 边界温度是否一刀切，导致强场区没有额外热支撑。

#### 4.5.2 初始速度和磁拓扑是否匹配

如果闭合场区一开始也被赋予明显向外流速，那么这些地方会先经历一轮很不自然的再调整。这个问题未必是当前高速流的唯一主因，但它会明显增加前期松弛的混乱程度。

建议检查：

- 高速流根部是位于开放场区还是闭合场区；
- 如果是在闭合场附近，初始速度是否过大；
- 是否值得让闭合场区初始速度更接近静止。

#### 4.5.3 是否存在显式或隐式修正把热压/密度推到边界值

如果代码中存在 pressure floor、density floor、回退到上一步有效值、局地 limiter 等处理，建议记录这些修正在问题区域触发的时间和次数。因为：

- 如果高速流总是在修正器频繁触发之后长出来；
- 或者问题区域的 $p$ 与 $\rho$ 总被频繁夹到允许区间端点；

那就说明该位置已经进入数值脆弱区，后面的高速流多半是症状，不是起点。

### 4.6 最值得做的几组对照实验

#### 实验 A：只改强场区密度

保持磁场不变，只把强场区或高 $v_A$ 区域的密度提高一倍、几倍，或者采用平滑的 $v_A$-triggered 增密。

如果高速流显著减弱，那么主因就是：

$$
\rho\ \text{太低} \Rightarrow v_A\ \text{过高} \Rightarrow \text{低-} \beta\ \text{数值敏感性过强}.
$$

#### 实验 B：只改强场区热压 / 最小 $\beta$

例如给边界或内部更新加一个较弱的最小热压支撑，或者人为保证

$$
\beta \gtrsim 10^{-3}
$$

这一类量级。

如果高速流明显减弱，说明主因更偏向：

$$
p_{\mathrm{th}}\ \text{过小} \Rightarrow \text{压强支撑消失}.
$$

#### 实验 C：只对 magnetogram 或球谐高阶部分做更强滤波

若高速流对磁图平滑非常敏感，说明其根源很大程度上来自：

- 强场峰值；
- 尖锐场梯度；
- 或边界驱动过于激烈。

#### 实验 D：比较开放场区和闭合场区上的行为

如果高速流只在开放场线根部或准开放区域上方出现，那么它更可能由低-$\beta$ 强场区的受力失衡沿开放方向放大而来；如果它反而从闭合场顶端或分界附近出现，则更值得检查拓扑、边界速度和重构误差。

### 4.7 我现在对这个问题的工作判断

基于当前信息，我会把这件事按下面的优先顺序来判断：

#### 第一判断

**问题已经明确和强磁场区的极低热压有关。**

因为我已经看到 $p_{\mathrm{th}}$ 接近 0，这不是普通背景风速偏高，而是强场区热力学支撑已经接近塌陷。

#### 第二判断

**对于我的代码，主因更可能是内能自己掉得太低，而不是 total energy 的 catastrophic cancellation。**

因为我解的是内能方程。

#### 第三判断

**下一步最有希望直接改善问题的手段，是强场区增密和最低热压/最小 $\beta$ 支撑。**

因为这两步最直接对应当前暴露出来的两个症状：

- $v_A$ 过高；
- $p_{\mathrm{th}}$ 过低。

#### 第四判断

**如果问题对 magnetogram 滤波非常敏感，那就说明边界强场峰值和场梯度本身也是重要来源。**

这种情况下，先把峰值压到 20–30 G 量级测试，是合理而且和文献一致的工程做法。

---

## 5. 可以直接执行的检查清单

为了方便实际操作，这里把前面的判断压缩成一个最实用的检查清单。

### 5.1 先做输出

在问题区域输出随时间变化的：

$$
p_{\mathrm{th}},\ e,\ \rho,\ \beta,\ v_A,\ c_s,\ |B|,\ \nabla\cdot\mathbf v.
$$

### 5.2 再拆项

沿问题流束根部，输出：

$$
-\nabla p,
\qquad
\mathbf J\times \mathbf B,
\qquad
\rho\mathbf g.
$$

以及内能方程中的各 source term。

### 5.3 然后做最小对照实验

1. 只增密；
2. 只抬高最低热压或最小 $\beta$；
3. 只平滑 magnetogram / 降球谐高阶；
4. 分开测试开放场区和闭合场区初始速度。

### 5.4 最后决定长期方案

如果 A 最有效：优先走 **强场区增密 / $v_A$ 限幅**。  
如果 B 最有效：优先走 **最小热压 / 最小 $\beta$ 支撑**。  
如果 C 最有效：优先走 **磁图预处理 / 低日冕强场控制**。  
如果只有 A+B+C 一起才有效：说明问题确实是典型的 **强场 + 低密度 + 极低热压 + 边界场梯度** 共同作用。

---

## 6. 参考文献与原文位置索引

### COCONUT / time-evolving / high-speed streams / PP

- Wang et al. 2025, *Time-evolving coronal modelling of the solar maximum around the solar storms in May 2024 by COCONUT*  
  链接：<https://arxiv.org/pdf/2505.11990>
  - p.1, lines 81–87：strong complex magnetic fields with low-β values；abnormally high-speed streams above active regions
  - p.2, lines 246–248：filtered magnetogram，near-surface field limited to about 30 G
  - p.2, lines 268–270：low-β and rapid boundary-field variation lead to unphysical negative thermal pressure
  - p.2–3, lines 273–324：boundary density and pressure formulas with $V_A$ limit and $\beta_{\min}=10^{-3}$
  - p.3, lines 325–336：threshold/revert style PP treatment
  - p.1, lines 107–113：mass-flux limitation and thermal-pressure PP

- Wang et al. 2024, *Efficient magnetohydrodynamic modelling of the time-evolving solar corona...*  
  链接：<https://arxiv.org/pdf/2409.02043>
  - p.2, lines 186–194：increase active-region density to reduce Alfvén speed; prevents negative thermal pressure and very high-speed streams above active regions; smooth density increase when local $v_A>2000$ km/s
  - p.2, lines 205–209：magnetic-pressure discretisation error comparable to thermal pressure in low-$\beta$ regions

- Wang et al. 2025/2026, *COCONUT: A coronal model with an energy decomposition strategy*  
  链接：<https://arxiv.org/pdf/2508.20423>
  - p.2, lines 113–116：non-physical negative thermal pressures in low-$\beta$ regions resulting from catastrophic cancellation
  - p.8, lines 620–623：old setup limited low-corona field to about 30 G; new setup can exceed 100 G near active regions
  - p.18, lines 909–911：older PP algorithm stable when low-corona field not obviously exceeding 30 G
  - p.19, lines 930–932：traditional full-energy COCONUT crashed when field exceeded 50 G

### SIP-IFVM / extended decomposition / transition-region broadening

- Wang et al. 2025, *SIP-IFVM: An observation-based magnetohydrodynamic model of coronal mass ejection*  
  链接：<https://arxiv.org/pdf/2506.19711>
  - p.2, lines 94–95：magneto-friction below $1.15\,R_\odot$; artificially increased density to broaden transition region
  - p.3, lines 139–140：adjust $B_{01}$ to accommodate $B_1$ and reset $B_1$ to zero each physical step

### AWSoM-R / boundary energy input

- Liu et al. 2025, *Forecasting the 8 April 2024 Total Solar Eclipse with Multiple Solar Photospheric Magnetograms*  
  链接：<https://arxiv.org/pdf/2503.10974>
  - p.3, lines 145–153：$(S_A/B)_\odot$ proportional to local radial field strength; too large near solar maximum may cause unstable chromospheric evaporation

### MAS / broadened transition region / semi-implicit

- MAS User Guide  
  链接：<https://www.predsci.com/mas/doc/User-Guide.pdf>
  - p.71, lines 1980–1989：transition region gradients extremely steep; artificially broadened transition region
  - p.10, lines 98–101：implicit and semi-implicit time differencing
  - p.30–31, lines 710–712：switch between fully explicit and semi-implicit algorithm

---

## 7. 一句话总结

从现有文献和我当前模型的表现来看，这类“强磁场区上方异常高速流”最常见的数值根源是：

$$
\text{强场} + \text{低密度} + p_{\mathrm{th}}\to 0 + \text{低 }\beta
$$

使得系统进入极端敏感状态；此时动量方程由大力相消后的微小残差决定，任何边界不匹配、离散误差或热力学支撑不足，都会被放大成明显的局地高速流。

对我当前这套**求解内能方程**的模型来说，最值得优先检查和测试的是：

1. 强场区内能为什么掉得这么低；
2. 强场区是否需要增密以限制 $v_A$；
3. 是否需要给热压设置最低支撑或最小 $\beta$；
4. magnetogram / 球谐重建的强场峰值是否需要进一步平滑。

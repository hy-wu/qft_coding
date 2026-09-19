# 从 Dyson 级数到费曼图：$T\{\phi^4\cdots\phi^4\}$ 的展开 (From the Dyson Series to Diagrams)

> **主要目标**：把 $\langle f|S|i\rangle$ 完整归约到"处理 $T\{\phi^4(x_1)\cdots\phi^4(x_n)\}$"这一个
> 对象上，并给出处理它的全部规则：三类收缩、组合计数、对称因子的来源，以及 $n=1,2$ 的
> 完全显式展开。最后给出汇总的大公式链条。
>
> 配套阅读：[`feynman_diagrams.md`](feynman_diagrams.md)（规则与观测量）、
> [`asymptotic_series.md`](asymptotic_series.md)（为什么可以这样展开）。

**约定**（与 [`feynman_diagrams.md`](feynman_diagrams.md) 一致，Peskin 型）：

$$
\mathcal L=\tfrac12(\partial_\mu\phi)^2-\tfrac12m^2\phi^2-\frac{\lambda}{4!}\phi^4
\quad\Longrightarrow\quad
\mathcal H_{\rm int}(x)=\frac{\lambda}{4!}\phi^4(x)\quad(\lambda>0),
$$

$$
D_F(x-y)=\int\!\frac{d^4k}{(2\pi)^4}\frac{i}{k^2-m^2+i\epsilon}e^{-ik\cdot(x-y)},
\qquad
\langle\mathbf p|\mathbf q\rangle=(2\pi)^3\,2\omega_p\,\delta^3(\mathbf p-\mathbf q).
$$

---

## 一、 归约：为什么全部困难都落在 $T\{\phi^4\cdots\}$ 上

$|i\rangle$、$|f\rangle$ 是**渐近态**（$t\to\mp\infty$ 时演化到自由哈密顿量 $H_0$ 的本征态）。
用相对论归一化 $|i\rangle=\sqrt{2\omega_{p_1}}\sqrt{2\omega_{p_2}}\;a^\dagger_{\mathbf p_1}a^\dagger_{\mathbf p_2}|0\rangle$，
则

$$
\langle f|S|i\rangle=\langle f|\,T\exp\Big[-i\!\int dt\,H_{\rm int}^I(t)\Big]|i\rangle
=\sum_{n=0}^\infty\frac{(-i\lambda)^n}{n!\,(4!)^n}
\int\prod_{a=1}^{n}d^4x_a\;
\underbrace{\big\langle f\big|\,T\big\{\phi^4(x_1)\cdots\phi^4(x_n)\big\}\big|i\big\rangle}_{\textbf{全部困难所在}}
$$

**问题被完整归约到**：对每个 $n$ 计算这个由 $4n$ 个场算符组成的编时矩阵元。其中 $\phi$ 已是
相互作用绘景下的**自由场** —— 这点很重要，唯有如此才能做模式展开并用 Wick 定理。

### 1.1 一个容易被跳过的逻辑环节：外态也必须进入 T-乘积

$a^\dagger,a$ 不是"旁观者"，它们**可以和 $\phi$ 收缩**。理由是 $a^\dagger_{\mathbf p}$ 可以写成
某个表面对 $\phi$ 的积分，

$$
a^\dagger_{\mathbf p}=i\int d^3x\;e^{ip\cdot x}\overleftrightarrow{\partial_0}\,\phi(x),
$$

而渐近态定义在 $t=\mp\infty$，位于整个相互作用区域**之外**，所以把它移进 T-乘积内部是合法的。

> 这一点严格来说是 **LSZ 约化公式**的内容，不是"显然"的（见
> [`feynman_diagrams.md`](feynman_diagrams.md) §8.1）。实践中先用下面的规则，再用 $n=1$
> 的例子验证它确实给出正确答案。

---

## 二、 处理它的完整规则：三类收缩

在 $\langle f|T\{\phi(x_1)\cdots\phi(x_m)\}|i\rangle$ 中，**先把外态摊开成场算符**，
然后对所有场（内部的 $4n$ 个 + 外态带来的那些）做 Wick 收缩。

| 收缩类型 | 结果 | 图上 |
|---|---|---|
| 场 — 场 | $\displaystyle D_F(x_i-x_j)=\int\!\frac{d^4k}{(2\pi)^4}\frac{i}{k^2-m^2+i\epsilon}e^{-ik(x_i-x_j)}$ | 内线 |
| 场 — 外态 | $e^{-ip\cdot x_i}$（入射）/ $e^{+ip\cdot x_i}$（出射） | 外线 |
| 未收缩的场 | 进入 $:\cdots:$，作用在 $\langle 0\vert$ 或 $\vert 0\rangle$ 上 | 不存在 |

**第三行是关键的选择定则**：$:\cdots:$ 的真空期望值为零。所以**每一个场都必须被收缩掉**，
否则整个项为零。正是这条定则把"求和于**所有配对**"约化为"求和于**所有图**"。

一旦外腿收缩掉，$e^{\pm ip\cdot x}$ 里每个 $x_a$ 的指数就用来做 $\int d^4x_a$，逐个生成
$(2\pi)^4\delta^4(\sum_{\rm in}p-\sum_{\rm out}p)$；最后整体剩一个守恒 $\delta$，
其系数定义为 $i\mathcal M$。

---

## 三、 组合爆炸与对称因子的来源

$4n$ 个场，完全配对的数目是

$$
(4n-1)!!=\frac{(4n)!}{2^{2n}(2n)!}:\qquad
n=1:\;3,\quad n=2:\;105,\quad n=3:\;10\,395,\quad n=4:\;2\,027\,025 .
$$

对 $2\to2$ 散射还要再乘上"哪些场分配给外腿"的选择数。**逐项手算是不可能的。**
但分母里的两个因子会大幅回收：

- $1/n!$ 抵消"把 $n$ 个顶点标上 $x_1,\dots,x_n$"的 $n!$ 种标号方式；
- 每个顶点的 $1/4!$ 抵消"4 个场分配到 4 条线"的 $4!$ 种方式。

剩下的净过度计数恰好是图的**对称因子**：

$$
(\text{给出拓扑 }G\text{ 的收缩方式数})=\frac{(4!)^n\,n!}{S(G)}
\quad\Longrightarrow\quad
\frac{1}{n!\,(4!)^n}\times\frac{(4!)^n n!}{S}=\frac{1}{S}
$$

**这就是费曼图这套记法的全部价值**：把 $(4n-1)!!$ 量级的组合求和压缩成
"少数几个拓扑 + 一个对称因子"。

**校验**：

| 图 | $n$ | 收缩方式数 | $S$ |
|---|---|---|---|
| $2\to2$ 树图（4 场全给外腿） | 1 | $4!=24$ | $1$ |
| tadpole 自能（2 外腿 + 1 自圈） | 1 | $4\cdot3=12$ | $2$ |
| 八字形真空图（两两自配） | 1 | $3$ | $8$ |
| sunset 自能（两顶点、三平行线） | 2 | $192$ | $6$ |

八字形：$\frac{3}{4!}=\frac18=\frac1S$ ✓（与 [`feynman_diagrams.md`](feynman_diagrams.md) §7 及
`qft_basics.tex` 里 $\frac{\lambda}{4!}\times3=\frac{\lambda}{8}$ 一致）。

---

## 四、 常见误解澄清

| 容易误以为 | 实际 |
|---|---|
| 对 $\mathcal L$ 做 $n$ 重四维积分 | $n$ 重积分来自"相互作用可以发生在 $n$ 个**不同时空点**"；$\mathcal L$ 是**局域**的，只依赖同一点的场，作用量只积一重 |
| 某点的贡献 $=$ 一个传播子 | 某点的贡献 $=$ **一串传播子相乘**（每条内线一个） |
| "高阶可以拆成二阶" | 是"**多点**拆成**两点收缩的乘积**"（Wick 定理），不是"高阶拆成低阶" |
| $n$ 是场的个数 | $n$ 是**顶点个数**；场的个数是 $4n$ |
| $4$ 个场的 $3$ 种配对是全部 | $4$ 个场时是 $3$ 种；$8$ 个场时是 $105$ 种；且配对数还取决于有几个场分给外腿 |

**"阶"这个字有两种含义，不要混**：

| 说法 | 含义 | 决定什么 |
|---|---|---|
| **耦合阶** $n$ | Dyson 展开到第 $n$ 阶 $=$ 顶点个数 | 有几重积分、几个顶点 |
| **几点**（2 点 / 4 点） | 场的个数 | 有几条线、几个传播子 |

$n=2$ 是"2 个顶点、8 个场"，不是"2 个场"。

---

## 五、 $n=1$ 完全显式地做一遍

$$
\langle f|S^{(1)}|i\rangle=\frac{-i\lambda}{4!}\int d^4x\;\big\langle f\big|T\{\phi^4(x)\}\big|i\big\rangle
$$

**同一个 $\phi^4(x)$**，四个场怎么分配决定了是哪个过程：

| 场的分配 | 收缩方式数 | $x$ 处的被积函数 | 结果 |
|---|---|---|---|
| 4 个全给外腿（$2\to2$） | $4!=24$ | $24\,e^{-ip_1x}e^{-ip_2x}e^{+ip_3x}e^{+ip_4x}$ | $i\mathcal M=-i\lambda$ |
| 2 个给外腿，2 个自成一对（tadpole） | $4\cdot3=12$ | $12\,D_F(y-x)D_F(z-x)D_F(0)$ | $-i\Sigma=\frac{-i\lambda}{2}D_F(0)$ |
| 0 个给外腿，两两自配（真空八字） | $3$ | $3\,D_F(0)^2$ | 真空能 $\frac{\lambda}{8}\int D_F(0)^2$ |

### 5.1 $2\to2$ 散射

$$
\begin{aligned}
\langle p_3p_4|S|p_1p_2\rangle
&=\frac{-i\lambda}{4!}\int d^4x\;\langle p_3p_4|\,T\{\phi^4(x)\}|p_1p_2\rangle\\[4pt]
&=\frac{-i\lambda}{4!}\int d^4x\;\underbrace{4!\,e^{-ip_1x}e^{-ip_2x}e^{+ip_3x}e^{+ip_4x}}_{24\ \text{种收缩：四个场各接一条外腿}}\\[4pt]
&=\frac{-i\lambda}{4!}\cdot 24\int d^4x\;e^{-i(p_1+p_2-p_3-p_4)x}\\[4pt]
&=-i\lambda\,(2\pi)^4\delta^4(p_1+p_2-p_3-p_4)
\qquad\Longrightarrow\qquad \boxed{\mathcal M=-\lambda},\quad S=1
\end{aligned}
$$

**这一项只有 $n=1$**：没有内线，位置积分只产出一个 $\delta$ 函数。

### 5.2 tadpole 自能（两点函数，不是散射）

$$
\begin{aligned}
G^{(2)}(y,z)\Big|_{\mathcal O(\lambda)}
&=\frac{-i\lambda}{4!}\int d^4x\;\langle 0|\,T\{\phi(y)\phi(z)\phi^4(x)\}|0\rangle\\[4pt]
&=\frac{-i\lambda}{4!}\int d^4x\;\underbrace{12\;D_F(y-x)D_F(z-x)D_F(0)}_{\text{注意：三个传播子相乘}}
=\frac{-i\lambda}{2}D_F(0)\int d^4x\,D_F(y-x)D_F(z-x)
\end{aligned}
$$

**$12$ 怎么来的**：把顶点的四个场记作 $A,B,C,D$，外场 $\phi(y),\phi(z)$。$y$ 的伙伴 $4$ 选 $1$，
$z$ 的伙伴剩 $3$ 选 $1$，剩下两个**只能**互相配 → $4\cdot3=12$。

做掉 $x$ 积分（卷积 → 动量空间）：

$$
\int d^4x\,D_F(y-x)D_F(z-x)
=\int\!\frac{d^4k}{(2\pi)^4}\,e^{-ik(y-z)}
\frac{i}{k^2-m^2+i\epsilon}\frac{i}{k^2-m^2+i\epsilon}
$$

与 $\frac{i}{k^2-m^2}(-i\Sigma)\frac{i}{k^2-m^2}$ 对比，读出

$$
\boxed{\;-i\Sigma_{\rm tad}=\frac{-i\lambda}{2}\,D_F(0)
=\frac{-i\lambda}{2}\int\!\frac{d^4k}{(2\pi)^4}\frac{i}{k^2-m^2+i\epsilon}\;}
$$

> 这正是 [`feynman_diagrams.md`](feynman_diagrams.md) 附录 A 表格里那一行的来源，
> 也解释了为什么 `qft_basics.tex` 写成 $\frac\lambda2\int\frac{d^3p}{(2\pi)^32\omega_p}$
> —— 同一张图、Kapusta 型记账。

---

## 六、 $n=2$ 一瞥：sunset 自能

$$
\frac{1}{2!\,(4!)^2}\int d^4x_1d^4x_2\;
\big\langle 0\big|T\{\phi(y)\phi(z)\underbrace{\phi^4(x_1)\phi^4(x_2)}_{8\ \text{个场}}\}\big|0\big\rangle
$$

$8$ 个场完全配对共 $7!!=105$ 种。其中给出 sunset 拓扑的：

$$
N_{\rm 收缩}=\frac{(4!)^2\cdot 2!}{S}=\frac{576\cdot 2}{6}=192,
\qquad
\frac{1}{2!\,(4!)^2}\times192=\frac{192}{1152}=\frac16=\frac1S
$$

$S=3!=6$ 的来源：两顶点间的三条平行内线可任意排列。动量空间：

$$
i\mathcal M_{\rm sunset}=\frac16\,(-i\lambda)^2
\int\!\frac{d^4k}{(2\pi)^4}\frac{d^4l}{(2\pi)^4}\;
\frac{i}{k^2-m^2+i\epsilon}\frac{i}{l^2-m^2+i\epsilon}\frac{i}{(k+l-p)^2-m^2+i\epsilon},
\qquad L=3-2+1=2
$$

---

## 七、 三件真正费劲的事

处理 $T\{\phi^4\cdots\phi^4\}$ 的困难可以拆成三层，**性质完全不同**：

1. **组合记账** —— $1/n!$、$1/4!$ 与配对数的抵消 → 对称因子。纯组合问题，靠图解决。
   **这是唯一被费曼图"彻底解决"的部分。**
2. **圈积分发散** —— 收缩给出 $L=I-V+1$ 个未定四维动量积分，UV 发散，必须正则化 + 重整化
   （见 [`feynman_diagrams.md`](feynman_diagrams.md) §11）。**图只负责写出积分，不负责让它收敛。**
3. **外腿处理** —— 未截腿的 Green 函数含外腿传播子极点，入/出的相位因子 $e^{\pm ipx}$ 也在这里。
   靠 LSZ 削掉。

---

## 八、 汇总：大公式链条

$$
\begin{aligned}
\langle f|S|i\rangle
&=\Big\langle f\Big|\,T\exp\Big[-i\!\int\! d^4x\,\mathcal H_{\rm int}(x)\Big]\Big|i\Big\rangle
&&\text{① 定义}\\[4pt]
&=\sum_{n=0}^{\infty}\frac{(-i)^n}{n!}\int\! d^4x_1\cdots d^4x_n\;
\Big\langle f\Big|\,T\big\{\mathcal H_{\rm int}(x_1)\cdots\mathcal H_{\rm int}(x_n)\big\}\Big|i\Big\rangle
&&\text{② 展开指数}\\[4pt]
&=\sum_{n=0}^{\infty}\frac{(-i\lambda)^n}{n!\,(4!)^n}\int\! d^4x_1\cdots d^4x_n\;
\Big\langle f\Big|\,T\big\{\underbrace{\phi^4(x_1)\cdots\phi^4(x_n)}_{4n\ \text{个场}}\big\}\Big|i\Big\rangle
&&\text{③ 代入 }\mathcal H_{\rm int}=\tfrac{\lambda}{4!}\phi^4\\[4pt]
&=\sum_{n=0}^{\infty}\frac{(-i\lambda)^n}{n!\,(4!)^n}\int\! d^4x_1\cdots d^4x_n\;
\sum_{C}\ \prod_{\text{内线 }(a,b)}D_F(x_a-x_b)\prod_{\text{外线 }j}e^{i\eta_j p_j\cdot x_{v(j)}}
&&\text{④ Wick 定理}\\[4pt]
&=\sum_{n=0}^{\infty}\frac{(-i\lambda)^n}{n!\,(4!)^n}\sum_{G}\ \frac{(4!)^n\,n!}{S(G)}\;\mathcal J_G
&&\text{⑤ 把 }C\text{ 按图 }G\text{ 归并}\\[4pt]
&=\sum_{n=0}^{\infty}\sum_{G}\frac{(-i\lambda)^n}{S(G)}\;\mathcal J_G
&&\text{⑥ 约掉 }n!\,(4!)^n\\[4pt]
&=\sum_{n=0}^{\infty}\sum_{G}\frac{(-i\lambda)^n}{S(G)}\;
\int\! d^4x_1\cdots d^4x_n\prod_{\text{内线}}D_F(x_a-x_b)\prod_{\text{外线}}e^{i\eta_jp_j\cdot x}
&&\text{⑦ 写开 }\mathcal J_G\\[4pt]
&=(2\pi)^4\delta^4\Big(\sum_{\rm in}p-\sum_{\rm out}p\Big)\;i\mathcal M
&&\text{⑧ 位置积分 }\to\delta
\end{aligned}
$$

$$
\boxed{\;
i\mathcal M=\sum_{n}\sum_{G}\frac{1}{S(G)}\;(-i\lambda)^n
\int\prod_{l=1}^{L}\frac{d^4k_l}{(2\pi)^4}\prod_{\text{内线 }e}\frac{i}{q_e^2-m^2+i\epsilon},
\qquad L=I-V+1
\;}
$$

**符号说明**：

- $C$ 跑遍**所有完全收缩方案**：$4n$ 个场每个恰好配一次，配成"场–场"（内线）或"场–外线"（外线）。
- $\eta_j=-1$（入射）/ $+1$（出射）；$v(j)$ 是外线 $j$ 所接的顶点。
- 同一对顶点之间有几条线就写几个 $D_F$；自圈写 $D_F(x_a-x_a)=D_F(0)$。
- $S(G)$ 只在第 ⑤ 步出现一次，来源是 $\frac{(4!)^nn!/S}{n!(4!)^n}=\frac1S$。

---

## 九、 三条绕开"逐项手算"的路

同一件事有三种等价表述，各有各的省力之处：

1. **LSZ 约化公式** —— 对象仍是 $T\{\phi^4\cdots\}$，但把外腿处理从"收缩外态"换成"对截腿
   Green 函数削极点"，外腿部分变干净、与内线结构解耦。
2. **路径积分** —— 最优雅。因为

   $$
   Z[J]=\exp\Big[-\frac{i\lambda}{4!}\int d^4x\Big(\frac1i\frac{\delta}{\delta J(x)}\Big)^4\Big]Z_0[J],
   \qquad Z_0[J]=\exp\Big[-\frac i2\int d^4x\,d^4y\,J(x)D_F(x-y)J(y)\Big]
   $$

   里的 $\exp$ 其泰勒级数**自带 $1/n!$**，正好抵消 Dyson 的 $1/n!$；而 $Z_0[J]$ 是高斯积分，
   $n$ 点函数自动就是全部配对之和。$T\{\ \}$ 这个符号与 Wick 定理的组合结构被指数形式完全吸收。
3. **Wick 定理的生成泛函形式** ——

   $$
   T\Big\{e^{\int J\phi}\Big\}=:e^{\int J\phi}:\;\exp\Big[\tfrac12\int J D_F J\Big],
   $$

   一行包含**所有**配对方式；对 $J$ 求 $4n$ 次泛函导数即回到原式。
   这正是路径积分版本能奏效的代数原因。

---

## 一句话总结

$\langle f|S|i\rangle$ 的**外态部分**（谁入谁出、归一化）与**内线部分**（哪些场互相收缩）
完全解耦；后者只需按"**所有场必须被收缩**"这一条定则枚举 —— 枚举的结果就是图，
而图与积分之间只差一个对称因子。所谓"必须处理 $T\{\phi^4\cdots\}$"，就是必须完成这次枚举
与随后的圈动量积分。

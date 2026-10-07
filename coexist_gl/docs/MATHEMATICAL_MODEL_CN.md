# coexist-gl 数学模型

合同版本：`coexist_gl_v2_no_ale_unit_edge`

## 1. 科学问题与状态空间

给定一个已经限定为“含有该gene family”的host集合，检验同一host内的同源物共存状态是否表现出超过普通gain/resolution过程的跨species-tree边持续性。

\[
S:n_h=1,\qquad C:n_h\ge2.
\]

`n_h`是host `h`中的正式family成员数。零拷贝host不在v2观测空间中，因此这是presence-conditioned模型。

## 2. 参数化

\[
G=D+T,\qquad A=G/L.
\]

species tree的每条边统一设为unit length，整体时间尺度不可单独识别。实现固定`L=1,G=A`。`G`只是`single -> coexist`的collapsed acquisition强度，不区分duplication和transfer。

基础生成矩阵为：

\[
Q(A)=\begin{pmatrix}-A&A\\2&-2\end{pmatrix}.
\]

第二行系数2表示二拷贝状态可由任一copy的resolution回到single，是当前collapsed模型的定义，不是独立估计参数。

## 3. Unit-edge基础转移矩阵

定义：

\[
\lambda=A+2,\quad \pi_S=\frac{2}{A+2},\quad
\pi_C=\frac{A}{A+2},\quad e=\exp[-(A+2)].
\]

则`P_1=exp(Q)`为：

\[
P_1=\begin{pmatrix}
\pi_S+\pi_Ce&\pi_C(1-e)\\
\pi_S(1-e)&\pi_C+\pi_Se
\end{pmatrix}.
\]

行表示parent状态，列表示child状态。

## 4. Coexistence persistence参数

唯一新增参数为`r>=1`。它只改变`C`行，将基础`C -> C`概率的odds乘以`r`：

\[
P_r(C\to C)=\frac{rP_1(C\to C)}{P_1(C\to S)+rP_1(C\to C)},
\]

\[
P_r(C\to S)=\frac{P_1(C\to S)}{P_1(C\to S)+rP_1(C\to C)}.
\]

`S`行不变，且`r=1`逐项恢复`P_1`。所以`r`是每条unit edge上的coexistence-continuation odds multiplier，不是概率倍数、copy数倍数或祖先来源概率。

## 5. Root prior

根节点使用基础生成矩阵的平稳分布：

\[
\Pr(root=S)=\frac{2}{A+2},\qquad
\Pr(root=C)=\frac{A}{A+2}.
\]

该prior不施加`r`倾斜；`r`只描述从parent `C`向child的持续。

## 6. 全树likelihood

对tip `v`：

\[
L_v(s)=\mathbf1(s=y_v).
\]

对内部节点：

\[
L_v(s)=\prod_{c\in children(v)}
\sum_{t\in\{S,C\}}P_r(s,t)L_c(t).
\]

整树likelihood为：

\[
L(A,r)=\sum_s\pi_s(A)L_{root}(s).
\]

程序在log space计算。固定`A,r`时复杂度为`O(E)`，`E`为species-tree边数。

## 7. Branch posterior与宏观事件期望

Inside/outside递推给出每条边`e=(u,v)`的联合posterior：

\[
\gamma_e(i,j)=\Pr(X_u=i,X_v=j\mid Y,A,r).
\]

四类宏观转移期望为：

\[
E[N_{ij}]=\sum_e\gamma_e(i,j).
\]

对应`S -> S`、`S -> C`、`C -> S`和`C -> C`。这些是species-state macro-edge expectations，不是reconciliation D/T/L事件数，也不必为整数。

## 8. 三层coexistence证据

### 8.1 观测存在

`n_C>0`只是数据事实，不计算p值。只有S和C同时存在时才拟合`r`。

### 8.2 Species-tree聚集

固定`n_S,n_C`并随机置换tip labels。统计量为unit-edge Fitch最小变化数`M`：

\[
p_{cluster}=\frac{1+\sum_{b=1}^{B}\mathbf1(M_b\le M_{obs})}{B+1}.
\]

较小`M`表示相同数量的C集中于更少树上区域。这是pattern test，不拟合`A,r`，也不确定真实转移历史。

### 8.3 `r>1` persistence检验

\[
H_0:r=1,\qquad H_1:r>1.
\]

Null固定`r=1`并优化`A`；free在`r>=1`下联合优化`A,r`：

\[
T=2\{\ell(\hat A,\hat r)-\ell(\hat A_0,1)\}.
\]

由于`r=1`是边界，报告：

\[
p=\tfrac12\Pr(\chi_1^2\ge T),\quad T>0.
\]

它检验普通gain/resolution比例重新优化后，是否仍需额外`C -> C`倾斜。

## 9. Profile与surface

Profile在每个固定`r`下重新优化`A`：

\[
\ell_p(r)=\max_A\ell(A,r).
\]

当前95%区间由离散profile网格中满足`2[ell(hat r)-ell_p(r)] <= chi-square(1,0.95)`的最小和最大点定义，是profile-grid interval，而不是插值后的精确端点。

二维surface直接计算预设`r x A`点的likelihood和宏观事件期望，用于诊断参数混淆及事件负担变化。

## 10. 可选SpeciesCut

完整树是主模型。可选SpeciesCut将hosts压缩为拓扑blocks，以beta-binomial emission描述block内C比例，并拟合nuisance parameter `kappa`。`kappa`不是进化率。不同cut改变观测空间，因此absolute likelihood不得跨cut比较。actin当前示例未启用非平凡SpeciesCut。

## 11. 解释边界

模型支持检验“coexistence状态是否较少形成、较长持续”。它不支持单独声称：

- 某copy来自duplication或transfer；
- `S -> C`等于一次D或T；
- `C -> S`等于一次真实gene loss；
- 所有coexistence来自同一祖先事件；
- 新模型logL可与ALE logL直接比较；
- `r>1`本身证明selection、功能分化或retention机制。
## 理论边界补充（2026-08-18）

普通 ALE 可以表示 duplication 后两个 daughter lineage 沿 species tree 长期共同存在；其“independent lineage propagation”不等于不能产生 ancestral multicopy，而是没有额外的 host-level joint coexistence state。coexist-gl 在 S=1 与 C=2+ 的状态层检验这种状态是否聚集并具有额外持续性。

正式参数空间固定为 `r>=1`：`r=1` 是 null，`r>1` 是 C→C persistence odds multiplier。`r<1` 已取消，不作为 accelerated resolution 分支。S→C、C→S、C→C 是宏观状态转移，不直接等同于逐事件 D/T/L；ALE 仅作可匹配的外部 benchmark。

边界参数说明：若未来在同一 coexist-gl 生成模型内检验 `H0:r=1` 对 `H1:r>1`，应考虑边界参数的非标准 LRT 分布或使用 parametric bootstrap；不能把常规无界一参数卡方近似直接套用。`r=1` 与普通 ALE 仅作概念上的无额外 persistence baseline，不代表两种 likelihood 数值等价。coexist-gl likelihood 不得与普通 ALE likelihood 直接相减或做跨观察模型 LRT。

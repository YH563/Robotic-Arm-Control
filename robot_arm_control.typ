#import "@preview/zh-kit:0.1.0" as zh-kit: setup-base-fonts,zhnumber-lower, zhnumber-upper 
#import "@preview/ergo:0.1.1": *
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

// 设置页面样式
#show : doc => setup-base-fonts(
  doc,
  cjk-mono-family: ("霞鹜文楷 Mono", "JetBrains Mono"),
)
#set math.mat(delim: "[")
#set par(leading: 1.2em, justify: true)
#set text(12pt)
#show heading.where(level:1):set align(center)
#set heading(
  numbering: "1.",
)
#set page(
  footer: context[
    #set align(center)
    #counter(page).display()
  ]
)
#set list(indent: 1.0em)
#set enum(indent: 1.0em)

#show heading.where(level:2): it =>{
  set text(fill: blue,size: 16pt, weight: "bold")
  it
  v(0.5em)
}

#show heading.where(level: 3):it =>{
  set text(fill: eastern,size: 14pt, weight: "bold")
  it
  v(0.5em)
}

#show heading.where(level: 1): it =>{
  set text(size: 20pt, weight: "bold")
  it
}

#let d(name, content) = defn[
  *#name*
][#v(0.5em)#content]

#let t(name, content, proof) = theorem[
  *#name*
][#content][#proof]

#let e(name, content) = ex[*#name*][#v(0.5em)#content]

#align(center, text(size: 25pt)[机器人理论笔记])

#align(center, text(size: 18pt)[引言])

此笔记可谓有生之年系列，在生物医学工程专业中承受巨量苦难后，笔者有幸参与到机械臂相关项目的开发，之前所学的李理论也算是有了用武之地了。

同时为了便于队友的机器人控制理论学习，特记录此笔记，作为个人对相关理论的总结。

为保证通俗易懂，本笔记尽量避开抽象的数学概念以及定理证明，主要针对在机器人控制中可能出现的数学公式以及引理进行论述，同时配合图片说明(数学理论部分图片较少，尽管李群同时也有几何上作为流形的特性，但这部分相关性较低，且直观的几何图片会误导读者，因此数学理论部分不会出现过多的几何图片)，以保证直观性。

不过需要说明的是，本笔记默认读者有较好的线性代数基础，因此不会重复叙述线性代数基本概念，若部分内容令读者感到生分，还望能自行查找资料查漏补缺。

笔者以李理论和螺旋运动作为笔记开篇，叙述在描述空间旋转运动时的常规方法的局限性，进而引出李理论，并基于李理论搭建螺旋理论的框架。螺旋理论是描述机器人运动的基础，同时也为后续机器人运动学提供了统一的数学理论。

#pagebreak()

#outline()

#pagebreak()

= 李理论与旋量理论

== 刚体运动

=== 旋转矩阵和齐次坐标

对于刚体运动，我们知道，可以将其运动拆解为平移运动与旋转运动两种，纯粹的平移运动是最为简单的，但是实际上的刚体运动是平移运动与旋转运动混合的结果，所以我们需要一种理论框架来描述这种混合运动。

我们首先讨论最简单的一种旋转，二维平面上绕原点的逆时针旋转。

假设一个二维平面中的三角形，绕原点进行旋转，假设其旋转矩阵为 $A$，二维平面的正交基为 $e_1=[1, 0]^T,e_2=[0, 1]^T$，显然有，旋转操作不改变三角形的面积，则对于旋转后的基 $e_1^prime, e_2^prime$ 有，
$ A [e_1, e_2]=[e_1^prime, e_2^prime] $
满足正交关系，
$ e_1^prime dot e_2^prime=0 $
其中 $dot$ 表示内积。且矩阵 $A$ 满足，
$ det(A)=1 $

因此我们可以发现，二维平面上的旋转矩阵满足正交关系(正交基转完还是正交基)以及其行列式必须为1。那么我们就得到了旋转矩阵的定义。
#d[旋转矩阵][
  在 $n$ 维空间 $RR^n$ 上的旋转矩阵定义为，
  $ R in {R in RR^(n times n)|R R^T=I, det(R)=1} $
]

那么我们知道在二维平面上的旋转矩阵可以写为，
$ R=mat(cos theta, -sin theta;sin theta, cos theta) $
或者将二维平面转为复平面，那么根据欧拉公式，旋转矩阵转为一维复数形式，即
$ R=e^(i theta) $

对于平移运动，设空间中点 $x in RR^n$，其进行平移运动，平移量为 $t in RR^n$，那么我们可以知道，平移后的点为，$y=x+t$，然而，平移并非线性变换，我们更加希望能够使用矩阵语言对平移进行描述，因此提出了齐次坐标的概念。

#d[齐次坐标][
  对 $n$ 维空间中点的坐标 $x=(x_1, x_2, dots, x_n)^T$，有齐次坐标为，
  $ x^prime = (x_1, x_2, dots, x_n, 1)^T $
  那么 $n$ 维空间中的平移变换 $t$ 便可转为 $n+1$ 维空间中的线性变换，写作矩阵为，
  $ M=mat(I,t;0,1) $
  其中 $I$ 为 $RR^(n times n)$ 中的单位矩阵。
]

以三维空间中的平移为例。
#e[][
  设点 $x=(1, 2, 3)^T$，平移量为 $t=(2, 2, 2)^T$，那么齐次坐标矩阵为，
  $ M=mat(1, 0, 0, 2;0, 1, 0, 2;0 ,0, 1, 2;0,0,0,1) $
  带入可得，
  $ M x^prime = mat(I, t;0, 1)mat(x;1)=mat(1, 0, 0, 2;0, 1, 0, 2;0 ,0, 1, 2;0,0,0,1)mat(1;2;3;1)=mat(3;4;5;1)=mat(x+t;1) $
]

齐次坐标的提出，将平移变换与矩阵表达进行了统一，也就是，我们可以用一个矩阵同时描述平移与旋转，即有矩阵 $M$，满足
$ M = mat(R,t;0,1) $
其中 $R$ 为旋转矩阵，$t$ 为平移量。而矩阵 $M$ 也被称为刚体变换矩阵，是机器人运动学中最为重要的矩阵。

=== 四元数

在处理三维空间旋转问题中，四元数一直都是一个无法忽略的数学工具，其具有比旋转矩阵更少的参数，且可以避免万向节锁问题，因此四元数在机器人控制中有着广泛的应用。



== 基础李理论

在这一章节我们将引入一些新颖的数学概念作为我们描述刚体运动的工具，同时也解决了不少在传统表示方法下不易解决的问题。

=== 李群与指数映射

首先需要引入群的概念。
#d[群][
  对于集合 $G$ 和二元运算 $dot$，构成*群* $(G,dot)$ 满足以下四条性质。
  - 封闭性：$forall a , b in G$ 有 $c = a dot b in G$
  - 幺元：$forall a in G$ 有 $exists e in G$，满足 $e dot a = a dot e = a$
  - 逆元：$forall a in G$ 有 $exists a^(-1) in G$，满足 $a^(-1) dot a = a dot a^(-1) = e$
  - 结合律：$forall a, b, c in G$ 满足 $a dot (b dot c) = (a dot b) dot c$

  倘若同时满足交换律，即
  $ forall a , b in G, a dot b = b dot a $
  成立，则该群也称为*阿贝尔群*。
]

为方便书写，后续提到群的概念则不强调二元运算 $dot$ 的存在，默认乘法(或矩阵乘法)，同时省略不写。根据群的定义，我们可以很轻松地证明前面讲到的旋转矩阵和刚体变换矩阵对矩阵乘法构成群。
#e[旋转群和特殊欧式群][
  对于 $n$ 维空间中的*旋转群*，记作 $S O(n)$，描述了旋转操作
  $ S O(n)={R in RR^(n times n)|R R^T=I, det(R)=1} $
  而对于 $n$ 维空间中的*特殊欧式群*，记作 $S E(n)$，描述了刚体运动，其中包括平移，旋转以及平移和旋转的耦合
  $ S E(n) = {mat(R,t;0,1)|R in S O(n), t in RR^n} $
  而这两种群，同时也被称为*李群*(Lie Group，一种既为群也为流形的数学结构)。
]

重新思考在三维空间中的旋转操作，确定一个物体的旋转操作只需要确定旋转轴以及绕轴旋转的角度，而这个实际上可以使用一个三维向量来表示，向量的模长表示旋转角度，而该向量的方向表示旋转轴。

不过，若使用旋转矩阵进行表达的话，实际上需要9个参数以及一个约束方程，在实际工程中会极大地引入不必要的麻烦。那么我们是否存在一种方式，可以从一个三维向量生成旋转矩阵并且天然地满足旋转矩阵的约束方程呢？答案是存在的，这种方式为指数映射。

#figure(
  image("figure/fig1.png", width: 45%),
  caption: [三维旋转示意图]
)

在引入指数映射前，我们需要先研究反对称矩阵，其与指数映射息息相关。
#d[反对称矩阵][
  对于 $n$ 维空间中的反对称矩阵 $A$，满足 $A^T=-A$，即 $A_(i j)=-A_(j i)$，其中 $A_(i j)$ 表示矩阵 $A$ 中第 $i$ 行第 $j$ 列的元素。

  因此，我们知道一个反对称矩阵只需要 $(n (n-1))/2$ 个参数即可确定，且对角线元素必定为0。

  在三维空间中，反对称矩阵由三个参数确定，我们可以将其对应为一个三维向量，设三维向量$bold(u) in RR^3$，其对应的反对称矩阵记作 $u^and$，而其逆运算记作 $u^or$，称 $and$ 为反对称矩阵算子。
  $ u^and = mat(0, -u_3, u_2;u_3, 0, -u_1;-u_2, u_1, 0) $
]

根据反对称矩阵的定义，我们可以得到以下几条重要性质，作为简化运算的依据。

#t[反对称矩阵重要性质][
  - 在三维空间中，叉乘与反对称矩阵乘法等价，即
  $ (bold(a) times bold(b))^and = a^and bold(b) $
  - 反对称矩阵的平方满足如下关系，
  $ (a^and)^2 = -abs(bold(a))^2 I + bold(a)bold(a)^T $
  进而可以得到，若反对称矩阵由单位向量生成，对于反对称矩阵的偶数次幂有，
  $ (a^and)^(2k) = (-1)^(k-1) (a^and)^2 $
  - 反对称矩阵的立方满足如下关系，
  $ (a^and)^3=-abs(bold(a))^2 a^and $
  进而可以得到，若反对称矩阵由单位向量生成，对于反对称矩阵的奇数次幂有，
  $ (a^and)^(2k-1) = (-1)^(k-1) a^and $
][]

接下来我们仅对 $S O(3)$ 进行讨论，其结论可以自然地推广到 $S E(3)$ 上。

#v(0.5em)
#line(length: 100%)
#v(0.5em)

对于旋转矩阵 $R(t) in S O(3)$，其中 $t$ 为参数，可以由唯一的参数确定唯一的 $R(t)$，同时有 $R(0)=I$，$R(t)$ 满足李群性质，即
$ R^T R = I $

则两边同时对 $t$ 求导，有
$ dot(R)^T R+R^T dot(R)=O $

设 $omega^and=R^T dot(R)$，则根据上式可知，$omega^and$ 为一个反对称矩阵(即矩阵对角位相加为0)，则有，
$ dot(R) = R R^T dot(R)=R omega^and $

假设 $omega^and$ 为不随时间变化的常量，将上式看作常微分方程，发现其解为指数映射，对于矩阵也同样成立，则解为(也可写作 $e^(omega^and t)$ 的形式)，
$ R(t) = exp(omega^and t)=sum_(k=0)^infinity (omega^and t)^k/(k!) $

为了便于化简，我们将 $omega^and t$ 重写为 $omega^and theta$，其中 $omega$ 为单位向量，代表了旋转轴方向，而 $theta$ 为旋转角度，再根据反对称矩阵的性质，我们带入即可整理得到，
$ R(theta)&=I + omega^and (theta-1/(3!)theta^3+1/(5!)theta^5-dots) + (omega^and) ^2(1/(2!)theta^2-1/(4!)theta^4+1/(6!)theta^6-dots)\
&=I + (omega^and) sin theta + (omega^and)^2 (1-cos theta) $

因此，我们得到了在* $S O(3)$ 群上的指数映射表达式*，即同时也推广得到了三维空间中的旋转公式，即*Rodrigues旋转公式*。

而推广到 $S E(3)$ 上，三维向量 $omega$ 则需改写为六维向量 $xi = (rho^T,phi.alt^T)^T$，其中 $rho$ 为平移量，$phi.alt$ 为旋转量，则 $xi$ 对应的反对称矩阵(并非严格意义上的反对称矩阵)为，
$ xi^and = mat(phi.alt^and, rho;0,0) $

指数映射的结果为，
$ exp(xi^and) = mat(exp(phi.alt^and),bold(J)rho;0,1) $
其中 $bold(J)$ 为雅可比矩阵，满足如下关系式，
$ bold(J) = (sin abs(phi.alt))/abs(phi.alt) I + (1-(sin abs(phi.alt))/abs(phi.alt)) (phi.alt phi.alt^T)/abs(phi.alt)^2 + (1-cos abs(phi.alt))/abs(phi.alt)^2 phi.alt^and $

因此我们得到了在* $S E(3)$群上的指数映射表达式*，这也是研究机器人运动的重要工具。(证明留作练习，不是很难)

最后我们再介绍几个指数映射的性质，证明同样留作练习。

#t[指数映射性质][
  + $"d"(exp(A t)) \/"d"t = A e^(A t) = e^(A t) A$
  + 若 $A = P D P^(-1)$，则有 $exp(A t) = P exp(D t) P^(-1)$
  + 若有 $A B = B A$，则有 $exp(A)exp(B) = exp(A + B)$ (由于矩阵乘法不具备对称性，所以常规指数映射的性质不再成立)
  + $(exp(A))^(-1) = exp(-A)$
][]

#v(0.5em)
#line(length: 100%)
#v(0.5em)

根据上述推导我们可以知道，在三维空间中，如果我们已知旋转轴以及旋转角度，则可以用一个三维向量表示旋转信息，再基于三维向量生成反对称矩阵，通过指数映射自然生成满足约束的旋转矩阵。

同样的，我们也能融合平移量，共同构成刚体变换矩阵，不过从 $S E(3)$ 指数映射的结果可以显然推出，刚体变换运动的复杂性在于平移量与旋转量的耦合，因此我们无法直接通过一个三维向量生成满足约束的刚体变换矩阵。

在几何上，指数映射构建了一个从线性空间到流形空间的映射关系，这使得我们可以利用在线性空间上构建的加法关系定义在流形空间的求导操作，使其在涉及旋转矩阵的优化问题上提供了便利性。(这段话看不懂建议跳过喵，这里只是笔者基于数学角度给出的理解)

=== 李代数

在这一节中，我们要引入*李代数*的概念，其描述了上一节中所提到的反对称矩阵。

#d[李代数][
  李代数定义为李群 $G$ 在幺元 $e$ 处的切空间 $T_e G$ ，记作 $g$。同时在李代数上定义了*李括号*运算 $[X,Y]=X Y - Y X$，且满足如下性质
  - 封闭性：$forall X,Y in g$ 有 $[X,Y] in g$
  - 双线性：$[a X_1 + b X_2, Y]=a[X_1, Y]+b[X_2,Y]$
  - 自反性：$[X,X]=O$
  - 雅可比等价：$[X,[Y,Z]]+[Y,[Z,X]]+[Z,[X,Y]]=0$
]

这个定义看不懂也没关系，只需要知道李代数为切空间，所以允许我们像处理一般线性空间那样处理李代数，而且，切空间的定义，使得我们可以将李代数看作速度空间，包含了线速度和角速度两个物理量。

同时李括号运算作为一种特殊的运算，其衡量了李代数空间中两个元素的可交换性。

为了和李群进行区分，$S O(3)$ 对应的李代数为 $frak(s o(3))$，$S E(3)$ 对应的李代数为 $frak(s e(3))$。

在 $frak(s o(3))$ 中，李括号同时也等价于叉乘运算，即对于 $phi.alt_1, phi.alt_2 in RR^3$，满足如下关系，
$ [phi.alt_1^and, phi.alt_2^and] = (phi.alt_1 times phi.alt_2)^and $

#v(0.5em)
#line(length: 100%)
#v(0.5em)

既然存在指数映射，那么也就同时存在对数映射，允许我们从李群元素映射到李代数元素，在 $S O(3)$ 群中，即为从旋转矩阵，得到旋转轴和旋转角度。

根据Rodrigues旋转公式，对于三维旋转矩阵 $R$，有三维单位向量 $omega in RR^3$ 满足如下关系，
$ R = I + (omega^and) sin theta + (omega^and)^2 (1-cos theta) $

计算 $R$ 的迹，
$ tr(R) &= tr(I) + tr(omega^and) sin theta + tr((omega^and)^2)(1 - cos theta)\
&=3 + 0 + tr((omega^and)^2)(1 - cos theta) $

由于 $(omega^and)^2 = omega omega^T - I$，则有
$ tr((omega^and)^2) = tr(omega omega^T) - tr(I) = 2 abs(omega)^2 - 3 = -2 $

带入可得，
$ tr(R) = 1 + 2cos theta $

因此，旋转角度为
$ theta = arccos((tr(R)-1)/2) $

而考虑 $R - R^T$，有
$ R - R^T &= I + (omega^and) sin theta + (omega^and)^2 (1-cos theta) - (I - (omega^and) sin theta + (omega^and)^2 (1-cos theta))\
&=2 (omega^and) sin theta $

根据反对称矩阵性质，可以得到旋转轴，
$ omega = 1/(2 sin theta)mat(R_(32)-R_(23);R_(13)-R_(31);R_(21)-R_(12)) $

因此，我们得到了 $S O(3)$ 群到 $frak(s o)(3)$ 的*对数映射*，即
$ log(R) = omega^and theta = theta/(2 sin theta)(R - R^T) $

同理，有$S E(3)$ 群到 $frak(s e)(3)$ 的*对数映射*，对于 $S E(3)$ 群中元素，
$ T = mat(R,t;0,1) $

则有 $frak(s e)(3)$ 中元素 $xi^and$ 满足，
$ xi^and = mat(log(R),bold(J)^(-1) t;0,0) $

#v(0.5em)
#line(length: 100%)
#v(0.5em)

对数映射建立了李群到李代数的联系，使得我们可以将李群上的优化问题转化至李代数空间上，从而利用李代数上的性质进行求解。

=== 伴随表示

在处理机器人运动时，我们需要频繁地进行坐标系变换，而伴随表示就是有效处理坐标系变换的工具。

首先是对于坐标变换问题的研究。这种问题广泛存在于机器人运动学中，我们希望从空间坐标系(Space Frame)的变换得到在本体坐标系(Body Frame)的变换，从而得到关节角对末端执行器末端位置的影响。

假设我们现在有两个相同维度的坐标系 $A,B$，两个坐标系的标准正交基分别为 $bold(alpha),bold(beta)$，且已知坐标变换矩阵 $T$，满足
$ bold(beta) = T bold(alpha) $

假设在坐标系 $A$ 中存在一个线性变换矩阵 $M_A$，我们希望得到坐标系 $B$ 中对应的线性变换矩阵 $M_B$。

在 $A$ 中应用线性变换，然后再变换到 $B$ 中，等价于直接在 $B$ 中进行线性变换，因此可以得到等式，
$ M_B bold(beta) = T M_A bold(alpha) $

由于 $bold(beta) = T bold(alpha)$，因此有
$ M_B T bold(alpha) = T M_A bold(alpha) $

因此，
$ M_B = T M_A T^(-1) $

进而通过相似变换，我们得到了在新坐标系下的线性变换矩阵，也就是李群上的*伴随表示*的定义。

#d[李群伴随表示][
  设李群 $G$，对应的李代数为 $frak(g)$，对于李代数元素 $xi in frak(g)$，有李群元素 $g in G$，定义*伴随表示*为映射，满足
  $ "Ad"_g (xi) = g xi g^(-1) $
  伴随表示的定义，将李群上的线性变换，转化为了李代数上的线性变换。
]

接下来我们看一下在 $S O(3)$ 群和 $S E(3)$ 群上的伴随表示。

#e[$S O(3)$ 群上的伴随表示][
  设三维向量 $omega in RR^3$，三维旋转矩阵 $R$，则对应的伴随表示为，
  $ "Ad"_R (omega^and) = R omega^and R^T = (R omega)^and $
  
  对于任意向量 $x in RR^3$ 有，
  $ (R omega^and R^T)x = R (omega times (R^T x)) $

  根据旋转不改变叉乘，我们有，
  $ R(omega times (R^T x)) &= (R omega) times (R R^T x)\
  &= (R omega) times x $

  因此，
  $ (R omega^and R^T)x = (R omega) times x = (R omega)^and x $

  所以有，
  $ "Ad"_R (omega^and) = (R omega)^and $

  或使用矢量形式，$"Ad"_R (omega) = R omega$

  为方便，我们记作伴随表示矩阵 $"Ad"_R = R$。
]

#e[$S E(3)$ 群上的伴随表示][
  设刚体变换矩阵 $T = mat(R,t;0,1) in S E(3)$，李代数元素 $xi^and in frak(s e)(3)$，则对应的伴随表示为，
  $ "Ad"_T (xi^and) = T xi^and T^(-1) $

  在矢量形式下，即
  $ xi = mat(v;omega) in RR^6 $
  对应的伴随表示为，
  $ "Ad"_T (xi) = T xi T^(-1) = mat(R v + t times (R omega);R omega) = mat(R, -R t^and;0, R) xi $

  可以化简为伴随表示矩阵
  $ "Ad"_T = mat(R, -R t^and;0, R) in RR^(6 times 6) $

  证明留作练习。
]

从上面的两个例子，我们知道了，如何计算伴随表示矩阵，利用伴随表示矩阵就可以很轻易的计算坐标变换，也就是说，在空间坐标系下的速度量，可以通过伴随表示矩阵，映射到本体坐标系下的速度量。

除了在李群上有伴随表示，在李代数上也有伴随表示，其被广泛应用于机器人动力学中。在此仅做概念阐述，不过度展开。

#d[李代数伴随表示][
  对于李代数 $frak(g)$ 中两个元素 $xi,eta in frak(g)$，定义*伴随表示*为，其中 $[,]$ 为李括号。
  $ "ad"_xi (eta) = [xi,eta] $
  
  其与李群的伴随表示为导数关系，即
  $ "d"/("d"t) "Ad"_(exp(t xi)) (eta) |_(t=0) = "ad"_xi (eta) $
]

== Screw Theory(旋量理论)

旋量理论是处理机器人运动学中，坐标变换问题的工具，其将空间坐标系和本体坐标系之间的变换，统一为旋量的形式，从而简化了计算。

旋量理论主要包括运动旋量和力旋量两个部分，共同构成了机器人动力学的理论基础。

=== Twist(运动旋量)

在上一节中我们提到了空间坐标系和本体坐标系之间的变换，所谓空间坐标系，可以简单理解为固定的世界坐标系，也就是“上帝视角”。而本体坐标系，可以理解为固定在机械臂末端的坐标系，也就是“机器人视角”。我们需要研究在不同坐标系之间速度量的变换。

*PS*：本节数学符号上与前几节略有差异。

// 先留着
#pagebreak()

#v(0.5em)
#line(length: 100%)
#v(0.5em)

在此我们将空间坐标系记为 ${s}$，本体坐标系记为 ${b}$。
#figure(
  image("figure/fig2.png", width: 55%),
  caption: [坐标变换示意图]
)

根据上图，我们可以记从空间坐标系到本体坐标系之间的变换矩阵为，
$ T_(s b)(t)=T(t)=mat(R(t), p(t);0,1) $

其中 $R(t)$ 为旋转矩阵，$p(t)$ 为平移向量，并暂时忽略角标。

我们对时间求导，计算 $T^(-1)dot(T)$ 的结果
$ T^(-1)dot(T) &= mat(R^T, -R^T p; 0,1)mat(dot(R), dot(p); 0,1)\
&= mat(R^T dot(R), R^T dot(p) ; 0,1)\
&= mat(omega_b^and, v_b; 0,1) $

其中 $omega_b$ 为本体坐标系下的角速度，$v_b$ 为本体坐标系下的线速度。

因此，我们定义在本体坐标系下的*本体运动旋量(Body Twist)* $cal(V)_b$ 为(注意在此对旋量表述为先角速度后线速度，这与前文 $frak(s e)(3)$ 中的顺序相反)，
$ cal(V)_b = mat(omega_b; v_b) in RR^6\
cal(V)^and_b=T^(-1)dot(T)=mat(omega_b^and, v_b; 0,1)\ $

对于这里的旋量定义，忽略掉传统的线速度和角速度等物理量，我们需要明确的一点是，我们研究运动的核心在于其对于速度求导的结果，以及我们希望该运动在对应坐标系下是固定的。

而本体旋量，就是描述了在本体坐标系下，运动模型相对于刚体固定的运动，将运动模式“拉回(右变换)”至本体坐标系，而在“拉回”下保持不变的即为本体旋量。(实在不理解就记住)

那么同理，我们可以计算 $dot(T)T^(-1)$ 的结果。
$ dot(T)T^(-1) &= mat(dot(R), dot(p); 0,1)mat(R^T, -R^T p; 0,1)\
&=mat(dot(R)R^T, dot(p)-dot(R)R^T p;0,0)\
&=mat(omega_s^and, v_s; 0,0) $

其中 $omega_s$ 为空间坐标系下的角速度，$v_s$ 为空间坐标系下的线速度，满足
$ v_s = dot(p) - omega_s times p $

因此，我们定义在空间坐标系下的*空间运动旋量(Space Twist)* $cal(V)_s$ 为，
$ cal(V)_s = mat(omega_s; v_s) in RR^6\
cal(V)_s^and = dot(T)T^(-1) =mat(omega_s^and, v_s; 0,0) $

而空间旋量也是同理，将运动模式“推前(左变换)”到空间坐标系，而在“推前”下保持不变的即为空间旋量。

那么我们现在就得到了在空间坐标系描述和本体坐标系描述下的两个旋量，然后我们希望能够将二者联系起来，也就是说，我们希望得到两者的变换关系。

根据定义，我们很容易能得到，
$ cal(V)_s^and &= dot(T) T^(-1)\
&=T cal(V)_b^and T^(-1) $

所以，我们发现，通过伴随表示可以实现两个坐标系旋量之间的变换，即
$ cal(V)_s^and = "Ad"_T cal(V)_b^and\
mat(omega_s^and; v_s) = mat(R, 0;p^and R,R) mat(omega_b^and;v_b) $

#v(0.5em)
#line(length: 100%)
#v(0.5em)

根据上述推导，得到了在不同坐标系中旋量的定义，以及如何在空间坐标系和本体坐标系之间进行旋量的变换。

在一些时候，我们并不需要关注真实的角速度和线速度，而是关注角速度和线速度的方向以及在对应方向上旋转和平移的结果，为方便，我们定义了旋量轴(Screw Axis)。

=== Wrench(力旋量)


= 机器人运动学

== Forward Kinematics

== Inverse Kinematics

= 控制理论

= 轨迹生成

= 姿态估计

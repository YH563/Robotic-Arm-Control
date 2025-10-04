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

为保证通俗易懂，本笔记尽量避开抽象的数学概念以及定理证明，主要针对在机器人控制中可能出现的数学公式以及引理进行论述，同时配合图片说明，以保证直观性。不过需要说明的是，本笔记默认读者有较好的线性代数基础，因此不会重复叙述线性代数基本概念，若部分内容令读者感到生分，还望能自行查找资料查漏补缺。

笔者以李理论和螺旋运动作为笔记开篇，叙述在描述空间旋转运动时的常规方法的局限性，进而引出李理论，并基于李理论搭建螺旋理论的框架。螺旋理论是描述机器人运动的基础，同时也为后续机器人运动学提供了统一的数学理论。

#pagebreak()

#outline()

#pagebreak()

= 李理论与螺旋运动

== 刚体运动

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

齐次坐标的提出，首先将平移变换与矩阵表达进行了统一，同时，我们也可以用一个矩阵同时描述平移与旋转，即有矩阵 $M$，满足
$ M = mat(R,t;0,1) $
其中 $R$ 为旋转矩阵，$t$ 为平移量。而矩阵 $M$ 也被称为刚体变换矩阵，是机器人运动学中最为重要的矩阵。

#pagebreak()

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
  而这两种群，同时也被称为*李群*(Lie Group，一种即为群也为流形的数学结构)。
]

重新思考在三维空间中的旋转操作，确定一个物体的旋转操作只需要确定旋转轴以及绕轴旋转的角度，而这个实际上可以使用一个三维向量来表示，向量的模长表示旋转角度，而该向量的方向表示旋转轴。

不过，若使用旋转矩阵进行表达的话，实际上需要9个参数以及一个约束方程，在实际工程中会极大地引入不必要的麻烦。那么我们是否存在一种方式，可以从一个三维向量生成旋转矩阵并且天然地满足旋转矩阵的约束方程呢？答案是存在的，这种方式为指数映射。

#figure(
  image("figure/image.png", width: 45%),
  caption: [三维旋转示意图]
)

在引入指数映射前，我们需要先研究反对称矩阵，其与指数映射息息相关。
#d[反对称矩阵][
  对于 $n$ 维空间中的反对称矩阵 $A$，满足 $A^T=-A$，即 $A_(i j)=-A_(j i)$，其中 $A_(i j)$ 表示矩阵 $A$ 中第 $i$ 行第 $j$ 列的元素。

  因此，我们知道一个反对称矩阵只需要 $(n (n-1))/2$ 个参数即可确定，且对角线元素必定为0。

  在三维空间中，反对称矩阵由三个参数确定，我们可以将其对应为一个三维向量，设三维向量 $hat(omega) in RR^3$，我们将其确定的反对称矩阵记作 $[hat(omega)]$，即
  $ [hat(omega)] = mat(0, -omega_z, omega_y;omega_z, 0, -omega_x;-omega_y, omega_x, 0) $

  对于三维向量 $bold(u) in RR^3$，其对应的反对称矩阵记作 $u^and$，而其逆运算记作 $u^or$，则有 $(u^and)^or = bold(u)$ (此种符号表达仅在需要区分正逆运算时使用，一般使用前一种记法)。
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

对于旋转矩阵 $R(t) in S O(3)$，其中 $t$ 为参数，可以由唯一的参数确定唯一的 $R(t)$，同时有 $R(0)=I$，$R(t)$ 满足李群性质，即
$ R^T R = I $

则两边同时对 $t$ 求导，有
$ dot(R)^T R+R^T dot(R)=O $

设 $[hat(omega)]=R^T dot(R)$，则根据上式可知，$[hat(omega)]$ 为一个反对称矩阵(即矩阵对角位相加为0)，则有，
$ dot(R) = R R^T dot(R)=R[hat(omega)] $

假设 $[hat(omega)]$ 为不随时间变化的常量，将上式看作常微分方程，发现其解为指数映射，对于矩阵也同样成立，则解为，
$ R(t) = exp([hat(omega)] t)=sum_(k=0)^infinity ([hat(omega)] t)^k/(k!) $

为了便于化简，我们将 $[hat(omega)]t$ 重写为 $[hat(omega)]theta$，其中 $hat(omega)$ 为单位向量，代表了旋转轴方向，而 $theta$ 为旋转角度，再根据反对称矩阵的性质，我们带入即可整理得到，
$ R(theta)&=I + [hat(omega)](theta-1/(3!)theta^3+1/(5!)theta^5-dots) + [hat(omega)]^2(1/(2!)theta^2-1/(4!)theta^4+1/(6!)theta^6-dots)\
&=I + [hat(omega)] sin theta + [hat(omega)]^2 (1-cos theta) $

因此，我们得到了指数映射的定义，同时也推广得到了三维空间中的旋转公式，即*Rodrigues旋转公式*。

根据上述推导我们可以知道，在三维空间中，如果我们已知旋转轴以及旋转角度，则可以用一个三维向量表示旋转信息，再基于三维向量生成反对称矩阵，通过指数映射自然生成满足约束的旋转矩阵。

在几何上，指数映射构建了一个从线性空间到弯曲光滑空间的映射关系，这使得我们可以利用在线性空间上构建的加法关系定义在弯曲光滑空间的求导操作，使其在涉及旋转矩阵的优化问题上提供了便利性。(这段话看不懂建议跳过喵，这里只是笔者基于数学角度给出的理解)



=== 李代数

== 螺旋理论

= 机器人运动学

== Forward Kinematics

== Inverse Kinematics

= 控制理论

= 轨迹生成

= 姿态估计

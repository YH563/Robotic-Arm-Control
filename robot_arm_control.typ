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

笔者以李理论和旋量理论作为笔记开篇，叙述在描述空间旋转运动时的常规方法的局限性，进而引出李理论，并基于李理论搭建旋量理论的框架。旋量理论是描述机器人运动的基础，同时也为后续机器人运动学提供了统一的数学理论。

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

一般使用 $T$ 来描述刚体变换，设空间中存在点 $p_a,p_b$，分别为{a},{b}两个坐标系，设 $T_(a b)$ 为从坐标系{a}到坐标系{b}的变换，那么有，
$ p_b = T_(a b) p_a $

而反过来变换时，则有，
$ p_a = T_(b a)p_b = T_(a b)^(-1)p_b $

=== 四元数

在处理三维空间旋转问题中，四元数一直都是一个无法忽略的数学工具，其具有比旋转矩阵更少的参数，且可以避免万向节死锁问题，因此四元数在机器人控制中有着广泛的应用，不过一般仅应用于纯旋转问题。

在描述三维物体的旋转时，我们一般有三种方式，分别为旋转矩阵，欧拉角，以及四元数。

旋转矩阵我们在前一节中已经讲过，就不再重复，我们先简单介绍一下欧拉角。

*欧拉角*是指通过绕固定轴旋转来表示任意方向的旋转的办法，一般固定轴为物体自身的参考系，以下图为例，其中蓝色坐标系为世界坐标系(固定的)，红色坐标系为物体坐标系通过旋转后的新坐标系，而初始物体坐标系与蓝色坐标系重合。
#figure(
image("figure/fig3.png", width: 30%),
caption: [欧拉角示意图])

从图片上可以看出，物体首先绕z轴旋转 $alpha$，将x轴转到 $N$ 方向上，然后绕x轴旋转 $beta$，最后再绕z轴旋转 $gamma$，旋转至新坐标系。(这里的x轴，z轴都是物体坐标系)。

也就是说，实际上的旋转矩阵 $R$ 可以写作为，
$ R &= R o t(z,gamma) R o t(x,beta) R o t(z,alpha)\
&=mat(cos gamma, -sin gamma, 0;sin gamma, cos gamma, 0;0, 0, 1)mat(1, 0, 0;0, cos beta, -sin beta;0, sin beta, cos beta)mat(cos alpha, -sin alpha, 0;sin alpha, cos alpha, 0;0, 0, 1) $

欧拉角十分直观，描述旋转也非常方便，但是欧拉角的旋转顺序中，倘若绕第二个轴的旋转 $plus.minus 90^circle.small$，那么会导致其直接丢失一个自由度，进而导致万向节死锁问题的发生(在此不做详细解释，读者可以自行查阅资料)。

而四元数有效解决了欧拉角的问题，在描述纯旋转运动中，四元数往往是一个更好的选择。

#d[四元数][
  四元数类似于复数的扩展，不过四元数并不满足交换律，其形式为，
  $ q = w + x i+ y j + z k quad (w,x,y,z in RR) $
  其中 $i, j, k$ 为虚部，满足，
  $ i^2 = j^2 = k^2 = i j k = -1 $

  并定义四元数共轭为，
  $ q^* = w - (x i+y j+z k) $

  有四元数模长公式，其中模长为1的四元数称为*单位四元数*
  $ abs(q)^2 = q q^*=w^2+x^2+y^2+z^2 $

  四元数的逆为，
  $ q^(-1) = q^* / abs(q)^2 $
]

根据四元数的定义，我们可以得到四元数的乘法公式，设四元数 $q_1 = w_1 + x_1i + y_1j + z_1k,q_2 = w_2 + x_2i + y_2j + z_2k$，那么 $q_1 q_2$为，
$ q_1 q_2 &= (w_1 + x_1i + y_1j + z_1k)(w_2 + x_2i + y_2j + z_2k)\
&=(w_1w_2 - arrow(v_1)dot arrow(v_2))+(w_1 arrow(v_2)+w_2 arrow(v_1)+arrow(v_1)times arrow(v_2)) $

其中，$arrow(v_1), arrow(v_2)$ 分别表示四元数的虚部，$dot,times$ 分别表示内积与外积。

同理可以得到纯虚四元数(实部为0)的乘法结果为，
$ q_1q_2 = - arrow(v_1)dot arrow(v_2) + arrow(v_1)times arrow(v_2) $

#v(0.5em)
#line(length: 100%)
#v(0.5em)

#figure(
image("figure/fig4.png", width: 30%),
caption: [旋转示意图])

接下来我们将推导如何利用四元数计算三维空间中的旋转。为便利，我们将四元数写作 $q=[w, arrow(v)]$ 的形式，使用 $w$ 表示实部，使用 $arrow(v)$ 表示虚部。

以上图为例，我们假设旋转轴 $A$ 对应的单位向量为 $arrow(u)$，而 $P$ 为待旋转的向量，假设其为 $arrow(v)$，并将其分解为平行于旋转轴和垂直于旋转轴的两个向量 $arrow(v)_parallel,arrow(v)_ tack.t$，显然有，旋转不改变平行于旋转轴的向量，只改变垂直于旋转轴的向量。

因此，我们可以设旋转后的向量为 $arrow(v)^prime$，那么有
$ arrow(v)^prime_parallel = arrow(v)_parallel $

而对于垂直于旋转轴的向量，我们可以通过四元数乘法实现。设 $v_tack.t = [0, arrow(v)_tack.t]$，$q = [cos theta, (sin theta) arrow(u)]$，$v^prime_tack.t = [0, arrow(v)^prime_tack.t]$ (不带箭头默认为纯虚四元数)，那么有，
$ v^prime_tack.t = q v_tack.t $

同时，根据四元数乘法的定义，不难证明以下公式成立，
$ q^2 = [cos 2theta, (sin 2theta) arrow(u)]\
q v_parallel = v_parallel q\
q v_tack.t = v_tack.t q^* $

因此，我们就可以得到四元数旋转公式，
$ v^prime &= v^prime_parallel + v^prime_tack.t\
&= v_parallel + q v_tack.t $

由于 $abs(q)=1$，则有
$ q^(-1) = q^* $

设 $p = [cos theta/2,(sin theta/2)arrow(u)]$，进而有，
$ q = p^2 $ 

带入四元数旋转公式，
$ v^prime &=p p^(-1) v_parallel + p^2 v_tack.t\
&= p v_parallel p^* + p v_tack.t p^*\
&=p v p^* $

因此，我们得到了四元数旋转公式。对向量 $arrow(v)$ 绕旋转轴 $arrow(u)$ 旋转 $theta$，那么设四元数$p = [cos theta/2,(sin theta/2)arrow(u)]$，则旋转后的结果为，
$ v^prime = p v p^* $

需要注意的是，旋转 $theta$ 时，参与计算使用的是 $theta/2$。

这样我们就顺利得到了四元数旋转公式，能够发现，四元数只使用了四个参数便可以描述旋转操作，远少于旋转矩阵的参数，同时又避免了万向节死锁问题，所以在描述纯旋转操作时，四元数是一个很好的选择，不仅在描述旋转操作中，在设计旋转运动插值算法时，也会使用四元数进行球面线性插值(Slerp)，实现了平滑插值的效果。


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

#t[反对称矩阵重要性质(三维空间下)][
  - 在三维空间中，叉乘与反对称矩阵乘法等价，即
  $ bold(a) times bold(b) = a^and bold(b) $
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
$ exp(xi^and) = mat(exp(phi.alt^and),G rho;0,1) $
其中矩阵 $G$ 满足如下关系式，
$ G = (sin abs(phi.alt))/abs(phi.alt) I + (1-(sin abs(phi.alt))/abs(phi.alt)) (phi.alt phi.alt^T)/abs(phi.alt)^2 + (1-cos abs(phi.alt))/abs(phi.alt)^2 phi.alt^and $

因此我们得到了在* $S E(3)$群上的指数映射表达式*，在描述机器人时，我们需要使用 $S E(3)$ 群来描述，并将其描述的结果称为 *位姿*，即位置(处在什么位置，平移量)与姿态(朝向什么方向，旋转量)。

最后我们再介绍几个指数映射的性质，证明留作练习。

#t[指数映射性质][
  + $"d"(exp(A t)) \/"d"t = A e^(A t) = e^(A t) A$
  + 若 $A = P D P^(-1)$，则有 $exp(A t) = P exp(D t) P^(-1)$
  + 若有 $A B = B A$，则有 $exp(A)exp(B) = exp(A + B)$ (由于矩阵乘法不具备交换性，所以常规指数映射的性质不再成立)
  + $(exp(A))^(-1) = exp(-A)$
][]

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

接下来我们将介绍对数映射，从李群空间映射回到李代数空间。

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
$ xi^and = mat(log(R),G^(-1) t;0,0) $

对数映射建立了李群到李代数的联系，使得我们可以将李群上的优化问题转化至李代数空间上，从而利用李代数上的性质进行求解。

#pagebreak()

=== 伴随表示

在处理机器人运动时，我们需要频繁地进行坐标系变换，而伴随表示就是有效处理坐标系变换的工具。

首先是对于坐标变换问题的研究。这种问题广泛存在于机器人运动学中，我们希望从空间坐标系(Space Frame，有时也叫基坐标系)的变换得到在本体坐标系(Body Frame，有时也叫末端坐标系)的变换，从而得到关节角对末端执行器末端位置的影响。

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

  为方便，我们记作伴随表示矩阵 $["Ad"_R] = R$。
]

#e[$S E(3)$ 群上的伴随表示][
  设刚体变换矩阵 $T = mat(R,t;0,1) in S E(3)$，李代数元素 $xi^and in frak(s e)(3)$，则对应的伴随表示为，
  $ "Ad"_T (xi^and) = T xi^and T^(-1) $

  在矢量形式下，即
  $ xi = mat(v;omega) in RR^6 $
  对应的伴随表示为，
  $ "Ad"_T (xi) = T xi T^(-1) = mat(R v + t times (R omega);R omega) = mat(R, -R t^and;0, R) xi $

  可以化简为伴随表示矩阵
  $ ["Ad"_T] = mat(R, -R t^and;0, R) in RR^(6 times 6) $

  为将伴随表示与伴随表示矩阵做区分，我们对伴随表示矩阵添加了中括号。

  证明留作练习。
]

从上面的两个例子，我们知道了，如何计算伴随表示矩阵，利用伴随表示矩阵就可以很轻易的计算坐标变换，也就是说，通过伴随表示，我们可以计算在两个李代数之间的转换，也就是描述刚体运动的速度量之间的转换。

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

在上一节中我们提到了空间坐标系和本体坐标系之间的变换，所谓空间坐标系，可以简单理解为固定的世界坐标系，也就是“上帝视角”。而本体坐标系，可以理解为固定在机械臂末端的坐标系，也就是“机器人视角”。我们需要研究在不同坐标系之间速度量的变换。(*PS*：本节数学符号上与前几节略有差异。)

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

变换矩阵满足约束为，
$ T T^(-1) = I $

对时间求导有，
$ dot(T)T^(-1) + T dot(T)^(-1) = O $

取 $dot(T)T^(-1)$ 描述空间速度，并计算 $dot(T) T^(-1)$ 的结果
$ dot(T)T^(-1) &= mat(dot(R), dot(p); 0,1)mat(R^T, -R^T p; 0,1)\
&=mat(dot(R)R^T, dot(p)-dot(R)R^T p;0,0)\
&=mat(omega_s^and, v_s; 0,0) $

其中 $omega_s$ 为空间坐标系下的角速度，$v_s$ 为空间坐标系下的线速度，满足
$ v_s = dot(p) - omega_s times p $

因此，我们定义在空间坐标系下的*空间运动旋量(Space Twist)* $cal(V)_s$ 为，
$ cal(V)_s = mat(omega_s; v_s) in RR^6\
cal(V)_s^and = dot(T)T^(-1) =mat(omega_s^and, v_s; 0,0) $

此处的角速度与线速度排列同前文相反。

在得到空间运动旋量后，我们希望计算在末端坐标系下的运动旋量，我们现在有末端位于本体坐标系下的运动旋量，我们希望得到末端位于末端坐标系下的运动旋量，那么对应的速度转换就要将空间运动旋量“逆”到末端坐标系下，因此根据伴随表示有，
$ cal(V)_b^and &=  "Ad"_(T_(b s))cal(V)_s^and = "Ad"_(T_(s b)^(-1))cal(V)_s^and = T^(-1)dot(T) $

从而得到了*本体运动旋量(Body Twist)*的定义，以及本体运动旋量与空间运动旋量之间的相互转换公式。

为方便，省略下标，写成如下形式，
$ cal(V)_s^and = "Ad"_(T)cal(V)_b^and $

或向量形式，
$ cal(V)_s &= ["Ad"_T]cal(V)_b\
mat(omega_s;v_s) &= mat(R,0;p^and R,R)mat(omega_b;v_b) $

旋量的定义，允许我们使用一个李代数元素来描述机器人的速度量，这在后续的速度运动学以及逆运动学中，都有广泛的应用。

在一些时候，我们并不需要关注真实的角速度和线速度，而是关注角速度和线速度的方向以及在对应方向上旋转和平移的结果，为方便，我们定义了*旋量轴(Screw Axis)*。

#d[旋量轴(Screw Axis)][
  对于运动旋量 $cal(V)=mat(omega;v) in RR^6$，为方便，我们将对其进行归一化处理。其中分为两种情况。
  - 当 $omega$ 不为 $0$ 时，我们定义旋量轴为 
  $ cal(S)=cal(V)/abs(omega)=mat(omega/abs(omega);v/abs(omega)) in RR^6 $
  - 当 $omega$ 为 $0$ 时，我们定义旋量轴为
  $ cal(S)=cal(V)/abs(v)=mat(0;v/abs(v)) in RR^6 $

  最终都可以写作如下形式，$dot(theta) = abs(omega)$ 或 $dot(theta) = abs(v)$
  $ cal(V) = cal(S)dot(theta) $

  有时也会写作，
  $ cal(V) = cal(S)theta $
]

根据旋量轴的定义，我们考虑纯旋转机械臂的旋量轴，一般的计算方式如下，设单位旋转轴为 $omega$，在旋转轴上选取任意一点 $p$，那么旋量轴表示为，
$ cal(S) = mat(omega^and;-omega times p) $

同时可以证明，点 $p$ 只要位于旋转轴线上(空间中的旋转轴线)，则无论点 $p$ 的选取如何，最终旋量轴的计算结果不变。

设在旋转轴线上选取新的一点 $p^prime = p + r$，那么新旋量轴中的线速度结果为( $r$ 与 $omega$ 共线，所以叉乘结果为0)，
$ v^prime = - omega times p^prime = - omega times (p + r) = -omega times p $

因此，只要点 $p$ 位于旋转轴线上，并不影响其旋量轴的计算结果。

接下来推导一下旋量轴指数映射的结果。作为对 $S E(3)$ 群指数映射的回顾。

设运动旋量为 $cal(V)=cal(S)theta$，那么其对应的李代数元素为，
$ cal(S)^and theta = mat(omega^and, v; 0,0)theta $

其中 $cal(S)$ 为旋量轴，所以有 $abs(omega)=1$。

那么根据指数映射定义，我们可以得到其对应的 $S E(3)$ 群元素为，
$ T &= exp(cal(S)^and theta)\
 &= sum_(k=0)^infinity (cal(S)^and theta)^k/(k!)\
 &= mat(exp(omega^and theta), G(theta) v;0,1) $

根据 $(omega^and)^3=-omega^and$，可以对矩阵 $G(theta)$ 进行化简，
$ G(theta) &= sum_(k=0)^infinity ((omega^and)^k theta^(k+1))/(k+1)!\
&=I theta + (omega^and theta^2)/2! + ((omega^and)^2 theta^3)/3! + dots\
&=I theta + omega^and (theta^2/2!-theta^3/3!+dots) + (omega^and)^2(theta^3/3!-theta^5/5!+dots)\
&=I theta + (1 - cos theta) omega^and + (theta - sin theta) (omega^and)^2 $

总的来说，指数映射的结果为，
$ T = exp(cal(S)^and theta) = mat(exp(omega^and theta), (I theta + (1 - cos theta) omega^and + (theta - sin theta) (omega^and)^2) v;0,1) $

若 $abs(omega) = 0$，则退化为，
$ T = mat(I, v theta; 0, 1) $

而对数映射在此则不再重复。

=== Wrench(力旋量)

#figure(
  image("figure/fig5.png", width: 45%),
  caption: [力旋量示意图]
)

在{a}坐标系下，设力矩为 $m_a = r_a times f_a$，其中 $f_a$ 为力，则我们可以得到在{a}坐标系下的力旋量为，
$ cal(F)_a = mat(m_a;f_a) in RR^6 $

在不同的坐标系下，功率(功率计算使用内积)应保持不变，因此我们有，
$ cal(V)_a^T cal(F)_a = cal(V)_b^T cal(F)_b $

根据速度旋量变换公式，$cal(V)_a = ["Ad"_(T_(a b))] cal(V)_b$，我们可以得到，
$ cal(V)_a^T cal(F)_b &= (["Ad"_T_(a b)] cal(V)_b)^T cal(F)_a\
&=cal(V)_b^T ["Ad"_T_(a b)]^T cal(F)_a\
&=cal(V)_b^T cal(F)_b $

因此有，*力旋量变换公式*为，
$ cal(F)_b = ["Ad"_T_(a b)] cal(F)_a $

到此，我们的数学基础理论部分就基本搭建完毕，从下一章开始，我们就要开始逐步讲述机器人理论相关的内容了。

#pagebreak()

= 机器人动力学

== Forward Kinematics(前向运动学)

在前向运动学中，我们主要使用PoE(Product of Exponentials)公式进行描述，不过在此之前，我们先以一个示例讲述一下如何推导从基坐标系到末端坐标系的变换矩阵。

#figure(
  image("figure/fig6.png", width: 60%),
  caption: [开链前向运动学示例]
)<openchain>

此时的x轴与y轴构成了二维平面，而z轴垂直于纸面。我们接下来要推导从坐标系{0}到坐标系{4}的变换矩阵，根据@openchain ，可以很轻松得到相邻坐标系间的变换矩阵分别为，

#figure(
  table(
    columns: 2,
    stroke: none,
    inset: (y:5pt),
    [$ T_(01) = mat(cos theta_1, -sin theta_1,0,0;sin theta_1, cos theta_1,0,0;0,0,1,0;0,0,0,1) $], [$ T_(12) = mat(cos theta_2, -sin theta_2,0,L_1;sin theta_2, cos theta_2,0,0;0,0,1,0;0,0,0,1) $],
    [$ T_(23) &= mat(cos theta_3, -sin theta_3,0,L_2;sin theta_3, cos theta_3,0,0;0,0,1,0;0,0,0,1) $],[$ T_(34) = mat(1, 0,0,L_3;0, 1,0,0;0,0,1,0;0,0,0,1) $]
  )
)

根据上述变换矩阵，我们可以得到，从坐标系{0}到坐标系{4}的变换矩阵为，
$ T_(04) = T_(01)T_(12)T_(23)T_(34) $

也就是，我们可以通过关节的角度，确定末端坐标系的姿态。

=== 基于基坐标系的PoE

接下来我们尝试将@openchain 的变换矩阵改写为指数乘积的形式。

想象(这里真的懒的画图了，虽然前边的也没画过)将机械臂完全放平，此时定义为零时，可以得到其变换矩阵为，
$ M = mat(1, 0, 0, L_1 + L_2 +L_3;0,1,0,0;0,0,1,0;0,0,0,1) $

然后，关节3先进行运动，设其旋转 $theta_3$，我们可以写出其运动旋量轴为，
$ cal(S)_3^and = mat(omega_3;v_3) = mat(0;0;1;0;-(L_1+L_2);0) $

其中，$v_3 = -omega_3 times p_3 = (0, -(L_1+L_2), 0)^T$

然后，我们就可以得到此时的变换矩阵为，
$ exp(cal(S)_3^and theta_3)M $

然后，我们再让关节2和关节1进行运动，使用同样的计算方式，那么最终变换矩阵为，
$ T_(04) = exp(cal(S)_1^and theta_1)exp(cal(S)_2^and theta_2)exp(cal(S)_3^and theta_3)M $

由此，我们就得到了在*基坐标系下的PoE公式*。

#figure(
  image("figure/fig7.png", width: 80%),
  caption: [n链前向运动学]
)

公式表述为，
$ T = product_(i=1)^n exp(cal(S)_i^and theta_i)M $

=== 基于末端坐标系的PoE

除了在基坐标系下的前向运动学，我们还希望得到在末端坐标系下的PoE公式。

针对这个，我们只需要将基坐标系下的运动旋量转换到末端坐标系下，便可以从基于基坐标系下的PoE公式转变为基于末端坐标系下的PoE公式。

根据运动旋量变换公式，我们可以知道，每个关节在末端坐标系下的运动旋量为，
$ cal(B)^and_i = "Ad"_(M^(-1))cal(S)_i^and = M^(-1)cal(S)_i^and M $

将其带入基坐标系下的PoE公式中，
$ T &= product_(i=1)^n exp(cal(S)_i^and theta_i)M\
&= product_(i=1)^n exp(M cal(B)_i^and M^(-1) theta_i)M\
&= product_(i=1)^n (M exp(cal(B)_i^and theta_i) M^(-1)) M\
&= M product_(i=1)^n exp(cal(B)_i^and theta_i) $

由此，我们便得到了在*末端坐标系下的PoE公式*。

到此，前向运动学就基本讲述完毕，还有一种广泛使用的用于描述前向运动学的公式为*Denavit–Hartenberg parameters (D–H parameters)*。这个公式使用的参数量相较于PoE公式要少，但是在理论计算中，我们更倾向于使用PoE公式，所以在此不对此公式进行展开，感兴趣的读者可以自行查阅相关资料。

== Velocity Kinematics and Statics(速度运动学与静力学)

在上一节中，我们了解了前向运动学中广泛使用的PoE公式，并推导了其在基坐标系与末端坐标系下的形式，这主要研究的是位置与方向的问题，而在这一节中，我们将研究开链机械臂的速度量也就是运动旋量(twist)以及静力学的相关问题。

=== Jacobian(雅可比矩阵)

设自由度为 $n$ ($n$ 为关节数)的机械臂的末端姿态为 $x in RR^m$ ($m$ 通常为6，包含位置分量与方向分量)，关节角度用 $theta in RR^n$ 表示，设 $x(theta) = f(theta(t))$，根据复合函数求导法则，可以得到，
$ dot(x)(theta) &= (partial f)/(partial theta) ("d"theta)/("d"t)\
&=J(theta)dot(theta) $

而其中，$J(theta) in RR^(m times n)$ 为*雅可比矩阵*，展开为，
$ J(theta) = mat((partial f_1)/(partial theta_1), dots, (partial f_1)/(partial theta_n);fence.dotted, dots.down, fence.dotted;(partial f_m)/(partial theta_1), dots, (partial f_m)/(partial theta_n)) $

雅可比矩阵包含的是速度的方向信息，而关节速度大小信息由 $dot(theta)$ 所定义，通过雅可比矩阵，我们找到了*关节速度与末端速度之间的联系*。

#e[雅可比矩阵示例][
  设 $x in RR^2$，且满足如下公式，
  $ x = mat(x_1;x_2) = mat(L_1 cos theta_1 + L_2 cos(theta_1 + theta_2);L_1 sin theta_1 + L_2 sin(theta_1 + theta_2)) $

  其中 $theta_1, theta_2$ 都是对时间 $t$ 的函数，则 $x$ 对 $t$ 求导，有
  $ dot(x)_1 = -L_1 dot(theta)_1 sin theta_1 - L_2(dot(theta)_1 + dot(theta)_2) sin(theta_1 + theta_2)\
  dot(x)_2 = L_1 dot(theta)_1 cos theta_1 + L_2(dot(theta)_1 + dot(theta)_2) cos(theta_1 + theta_2) $

  写成矩阵形式便得到了雅可比矩阵 $dot(x) = J(theta)dot(theta)$ 为
  $ mat(dot(x)_1;dot(x)_2) = mat(-L_1 sin theta_1 - L_2 sin (theta_1 + theta_2), - L_2 sin (theta_1 + theta_2); L_1 cos theta_1 + L_2 cos(theta_1 + theta_2), L_2 cos(theta_1+theta_2))mat(dot(theta)_1;dot(theta)_2) $

  也可以写作，
  $ v = J_1(theta)dot(theta)_1 + J_2(theta)dot(theta)_2 $
]

#v(0.5em)
#line(length: 100%)
#v(0.5em)

根据前向运动学的PoE公式，我们知道，基坐标系下的变换矩阵为，
$ T = product_(i=1)^n exp(cal(S)_i^and theta_i)M $

两边同时对时间求导，我们可以得到，
$ dot(T) &= (("d")/("d"t)exp(cal(S)_1^and theta_1))dots exp(cal(S)^and_n theta_n) M + exp(cal(S)^and_1 theta_1)(("d")/("d"t)exp(cal(S)_2^and theta_2))dots exp(cal(S)^and_n theta_n) M + dots\
&=(cal(S)_1^and dot(theta)_1exp(cal(S)_1^and theta_1)) dots exp(cal(S)^and_n theta_n) M + exp(cal(S)^and_1 theta_1)(cal(S)_2^and dot(theta)_2exp(cal(S)_2^and theta_2)) dots exp(cal(S)^and_n theta_n) M + dots $

同时有，
$ T^(-1) = M^(-1) product_(i=0)^(n-1) exp(-cal(S)_(n-i)^and theta_(n-i)) $

因此可以得到基坐标系下向量形式的运动旋量 $cal(V)_s$ 为，
$ cal(V)_s = (dot(T)T^(-1))^or = cal(S)_1 dot(theta)_1 + ["Ad"_(exp(cal(S)_1^and theta_1)) ]cal(S)_2 dot(theta)_2 + ["Ad"_(exp(cal(S)_1^and theta_1)exp(cal(S)_2^and theta_2)) ]cal(S)_3 dot(theta)_3 + dots $

同时，对于基坐标系下的运动旋量，我们有雅可比矩阵表示，
$ cal(V)_s = J_s (theta) dot(theta) = mat(J_(s 1)(theta), J_(s 2)(theta), dots, J_(s n)(theta))mat(dot(theta)_1;dot(theta)_2; fence.dotted; dot(theta)_n) $

因此我们得到了在基坐标系下的雅可比矩阵为，
$ J_s (theta) = mat(J_(s 1)(theta), J_(s 2)(theta), dots, J_(s n)(theta))\
J_(s 1) = cal(S)_1\
J_(s i) = ["Ad"_(H_(s i))]cal(S)_i quad i = 2, 3, dots, n $

其中，
$ H_(s i) = product_(j=1)^(i-1) exp(cal(S)_j^and theta_j) $

这样我们便得到了*空间雅可比*。

而根据末端坐标系下的PoE公式，
$ T = M product_(i=1)^n exp(cal(B)_i^and theta_i) $

可以计算得到，末端运动旋量 $cal(V)_b = (T^(-1)dot(T))^or$ 的计算公式，推导同上，证明留作练习。

$ cal(V)_b = (T^(-1)dot(T))^or = cal(B)_n dot(theta)_n + ["Ad"_(exp-(cal(B)_n^and theta_n)) ]cal(B)_(n-1) dot(theta)_(n-1) + dots + ["Ad"_(exp(-cal(B)_n^and theta_n) dots exp(-cal(B)_2^and theta_2)) ]cal(B)_1 dot(theta)_1  $

同时，对于末端坐标系下的运动旋量，我们有雅可比矩阵表示，
$ cal(V)_b = J_b (theta) dot(theta) = mat(J_(b 1)(theta), J_(b 2)(theta), dots, J_(b n)(theta))mat(dot(theta)_1;dot(theta)_2; fence.dotted; dot(theta)_n) $

因此我们得到了在末端坐标系下的雅可比矩阵为，
$ J_b (theta) = mat(J_(b 1)(theta), J_(b 2)(theta), dots, J_(b n)(theta))\
J_(b n) = cal(B)_n\
J_(s i) = ["Ad"_(H_(b i))]cal(B)_i quad i = 1, 2, dots, n-1 $

其中，
$ H_(b i) = product_(j=0)^(n-i-1) exp(-cal(B)_(n-j)^and theta_(n-j)) $

这样我们便得到了*本体雅可比*。

#figure(
  image("figure/fig8.png", width: 90%),
  caption: [空间雅可比与本体雅可比对比图]
)

在得到本体雅可比与空间雅可比后，我们可以很自然地通过旋量的相互转换得到雅可比矩阵间的相互转换公式。根据旋量转换公式，
$ cal(V)_s = ["Ad"_(T_(s b))] cal(V)_b $

将其带入雅可比公式中有，
$ cal(V)_s = J_s (theta) dot(theta) = ["Ad"_(T_(s b))] cal(V)_b = ["Ad"_(T_(s b))] J_b (theta) dot(theta) $

因此，得到了*雅可比矩阵转换公式*，
$ J_s (theta) = ["Ad"_(T_(s b))] J_b (theta) $

根据雅可比矩阵，我们可以轻松地根据关节的运动速度 $dot(theta)$ 来得到机械臂末端的运动速度与方向，同样的，我们可以使用雅可比矩阵应用在轨迹规划与逆运动学中。

=== Statics(静力学)

当机器人在进行装配、打磨等与环境有接触的任务时，我们就需要去实时计算每个关节应当输出的扭矩的大小，而根据虚功原理以及雅可比矩阵，我们就可以很轻松地计算得到实际的扭矩输出值。

假设一个自由度为 $n$ 的机械臂，则关节变量为，
$ theta = (theta_1, theta_2, dots, theta_n)^T $

其末端的位姿为 $x in RR^6$，设末端受到广义力
$ cal(F) = (f_x, f_y, f_z, tau_x, tau_y, tau_z)^T $

关节力矩为
$ tau = (tau_1, tau_2, dots, tau_n)^T $

我们考虑在惯性系下，且不考虑重力的情况，若机械臂系统要保持平衡，那么根据虚功原理，所有作用在系统上的力在任意虚位移上做的总虚功为0。

对于关节虚位移有，$delta theta$，对于末端虚位移有，
$ delta x = J(theta) delta theta $

那么根据虚功原理可得，
$ tau^T delta theta + cal(F)^T delta x = 0\
[tau^T + cal(F)^T J(theta)] delta theta = 0 $

化简可得，关节力矩为，
$ tau = - J^T (theta) cal(F) $

因此我们得到了*机械臂静力学关系*如下(通常省略负号)，

$ tau = J^T (theta) cal(F) $

=== Singularity Analysis(奇点分析)

根据之前的讨论，我们知道了可以利用雅可比矩阵从关节速度计算得到末端速度。在这一节中，我们将对雅可比矩阵 $J in RR^(m times n)$ 进行讨论，其性质关系到了机械臂能够运动的最大限度。

对于自由度为 $n$ 的机械臂，其雅可比矩阵的尺寸为 $n times 6$，我们在此考虑 $n >= 6$ 的情况，因此，雅可比矩阵的秩满足，
$ "rank"(J)<=6 $

同时我们知道，雅可比矩阵联系了关节速度与末端速度有，
$ dot(x) = J(theta) dot(theta) $

而末端速度 $dot(x) in RR^6$，也就意味着，如果要保持上式解的存在性，那么需要有 $"rank"(J)=6$，而当 $"rank"(J)<6$ 时，这意味着末端执行器一定的自由度，导致其会存在“方向的缺失”，在对应方向上的运动是无法产生的。

而这种秩的亏损而导致失去运动自由度的现象就被称为*Singularity(奇点)*，或者说奇点构型，指机器人处于某些特定的关节构型时，其雅可比矩阵的秩会降低。

对于常见的六自由度串联机器人，其奇点可以根据物理意义和位置分为三大类。

+ *边界奇点*
  - 定义：当机器人完全伸展开或完全收缩回来，使得末端执行器处于工作空间的边界时发生。
  - 物理意义：这类似于人的手臂完全伸直去够一个远处的物体。在这一点，你无法再沿着手臂的方向继续延伸。
+ *内部奇点*
  - 定义：发生在机器人工作空间内部，通常是由于多个关节轴共线或共面导致的。
  - 物理意义：机器人“折叠”起来，失去了一个或多个方向的移动能力。
+ *姿态奇点*
  - 定义：这类奇点不影响末端执行器的位置，但影响其姿态（方向）。
  - 物理意义：机器人无法绕操作空间的某个特定轴旋转。

对于奇点问题，我们可以通过一些方式来检测是否为奇点。

对于自由度为6的机械臂，雅可比矩阵为方阵，可以通过计算行列式是否为0，得到其是否满秩的结果。但对于机械臂冗余，也就是自由度大于6的情况下，就需要使用数值方式来分析其是否处于奇点位置，而这里就需要使用到*奇异值分解*的内容。对此我们简单补充一下奇异值分解的内容，读者可自行阅读或跳过。

#t[奇异值分解][
  对于矩阵 $A in RR^(m times n)$，存在如下分解，
  $ A = U Sigma V^T $
  其中 $U in RR^(m times m),Sigma in RR^(m times n),V in RR^(n times n)$
][
  对于矩阵 $A in RR^(m times n)$，有 $A^T A in RR^(n times n)$ 为半正定矩阵，则其存在 $n$ 个特征值满足，
  $ lambda_1 >= lambda_2 >= lambda_n >= dots>=0 $
  设其中非零特征值个数为 $r$，则矩阵 $A^T A$ 的秩为 $r$，同时有 $"rank"(A) = "rank"(A^T A)$，则矩阵 $A$ 的秩也为 $r$。由于 $A^T A$ 是对称矩阵，可以产生 $n$ 个标准正交的特征向量 $bold(v)_i$ (使用施密特正交化)，
  $ lambda_1, lambda_2, dots, lambda_r eq.not 0->bold(v)_1, bold(v)_2, dots, bold(v)_r\
  lambda_(r+1), lambda_(r+2), dots, lambda_n = 0->bold(v)_(r+1), bold(v)_(r+2), dots, bold(v)_n\ $
  
  进而有，
  $ A^T A bold(v)_i = lambda_i bold(v)_i quad i = 1,2,dots, n $

  设列向量 $bold(u)_i$ 为
  $ bold(u)_i = 1/sqrt(lambda_i)A bold(v)_i quad i = 1, 2, dots, r $
  可以证明列向量组 $bold(u)_i$ 构成了 $RR^m$ 中的一个标准正交向量组，然后对其进行扩充为 $m$ 维空间的标准正交基，进而可以写成矩阵形式，
  $ A(bold(v)_1, dots, bold(v)_r, bold(v)_(r+1), dots, bold(v)_n)&=(sqrt(lambda_1)bold(u)_1, dots, sqrt(lambda_2)bold(u)_r, 0, dots, 0)\
  &=(bold(u)_1, dots, bold(u)_r, bold(u)_(r+1), dots, bold(u)_m)mat("diag"(sqrt(lambda_1), dots, sqrt(lambda_r)), O;O,O) $

  从而得到了矩阵的奇异值分解为，
  $ A = U Sigma V^T $

  并记 $sigma_i = sqrt(lambda_i)$ 为*奇异值*。
]

对雅可比矩阵进行奇异值分解，
$ J = U Sigma V^T $

若其非0奇异值的个数小于 $m=6$，那么雅可比矩阵处于奇点构型。我们可以用*可操作度*来判断其是否接近或处于奇点，即计算
$ w = sqrt(det(J J^T)) = product_(i=1)^m sigma_i $

当 $w$ 越小，机器人越接近奇点。

还可以使用*条件数*衡量矩阵求逆的灵敏度，定义为，
$ kappa(J) = sigma_("max")/sigma_("min") $

条件数越大，距离奇点越近。

在具体工程中，奇点问题无法被完全规避，但仍有应对的策略。可以在轨迹规划之前，计算一些可能的奇点位置，达到事先避免奇点的目的。也可以使用自由度冗余的机器人，使用多余的自由度来避开奇点。在逆运动学中，会使用阻尼最小二乘法，防止因小奇异值导致的巨大关节速度，但代价是引入了跟踪误差，这是一种在速度和稳定性之间的权衡。

== Inverse Kinematics(逆向运动学)

在正向运动学中，我们希望通过关节的角度确定末端执行器的位姿，很显然，这个映射关系是唯一的。但是在一些场景下，我们希望得到在已知末端位姿的情况下，求解可能的关节角度，同样显然的是，这个问题并不是一定唯一的，可能存在多组解达到我们期望的位姿，所以逆运动学是复杂且难以求解的，同时也并不存在十分普适的求解方法，在此我们主要围绕几种主流方式进行讨论，在分类上，我们将其分为了解析法与数值法两种。

=== 解析法

解析法一般来并非对所有机器人都适用，我们以2连杆平面臂为例。

#figure(
  image("figure/fig9.png", width: 45%),
  caption: [2连杆平面臂]
)

假设末端坐标为 $(x,y)$，由余弦定理得，
$ x^2 + y^2 = L_1^2 + L_2^2 - 2 L_1 L_2 cos beta\
L_2^2 = L_1^2 + x^2 + y^2 - 2 L_1 sqrt(x^2+y^2)cos alpha $

则有，
$ beta = arccos((L_1^2 + L_2^2 - x^2 - y^2)/(2 L_1 L_2))\
alpha = arccos((L_1^2 + x^2 + y^2-L_2^2)/(2 L_1 sqrt(x^2+y^2))) $

因此可以得到关节角度为，
$ theta_1 = arctan(y/x) - alpha quad theta_2 = pi - beta $

不过这种计算方式泛用性并不好，以及我们并不清楚到底在什么情况下才能是用解析法求得解析解，不过，在6自由度机械臂的情况下，存在一个重要的判断依据：

- *Pieper准则*：如果一个*6自由度机械臂*的后三个关节轴*相交于一点*（即构成一个球形腕），或者后三个关节轴*相互平行*，那么该机械臂的逆运动学具有封闭解析解。

#v(0.5em)
#line(length: 100%)
#v(0.5em)

接下来就后三个关节轴相交于一点的情况进行分析，

因为后三个关节轴交于一点，在计算线速度时点 $q in RR^3$ 选取为三个轴相交的点，称为*腕部中心点*。那么有后三个旋量轴满足，
$ cal(S)_i = mat(omega_i; - omega_i times q) quad i=4,5,6 $

接下来我们可以证明，后三个关节的运动不会改变腕部中心点的位置。对于空间中的任意一点，其变换公式可以由PoE公式给出，
$ X(theta) = product_(i=1)^n exp(cal(S)^and_i theta_i) X(0) $
其中，$X(0)$ 为初始坐标(在计算时需转换为齐次坐标)。

那么对于6自由度机械臂的腕部中心点的坐标变换公式为,
$ q(theta) = product_(i=1)^6 exp(cal(S)^and_i theta_i) q $

我们考虑 $exp(cal(S)_4^and theta_4)q$ 的计算结果，
$ exp(cal(S)_4^and theta_4)q &= mat(exp(omega_4^and theta_4), - G(theta_4) (omega_4 times q);0,1)mat(q;1)\
&=mat(exp(omega_4^and theta_4) q - G(theta_4) (omega_4 times q);1)\ $

对上式进行化简，使用罗德里格斯旋转公式，带入可得，
$ exp(omega_4^and theta_4) q - G(theta_4) (omega_4 times q) &= I q + sin theta (omega^and_4) q + (1-cos theta)(omega^and_4)^2 q \
&- (I theta + (1-cos theta) omega^and_4 + (theta - sin theta)(omega^and_4)^2) (omega_4 times q)\
&= I q + sin theta (omega^and_4) q + (1-cos theta)(omega^and_4)^2 q \
&- (I theta omega_4^and q + (1 - cos theta) (omega_4^and)^2 q + (1-cos theta)(omega^and_4)^3 q)\
&= q quad("使用三维反对称算子平方和立方的性质化简") $

因此，后三个关节的运动并不会改变腕部中心点坐标，其坐标仅依赖于前三个关节，即，
$ product_(i=4)^6 exp(cal(S)_i^and theta_i)q = q $

在零位时，我们可以根据目标位姿 $T_d$ 满足，
$ product_(i=1)^6 exp(cal(S)^and_i theta_i) = T_d $

求出目标位姿下腕部中心点的位置 $q_t$，那么根据目标腕部中心点，带入有，
$ q_t = product_(i=1)^3 exp(cal(S)_i^and theta_i)q $

这样可以通过几何法(与具体机械臂构型有关)对前三个关节的角度进行求解，并且，此时相当于确定了目标位姿的位置逆解(确定下来在哪了，但是没有确定面向哪)，然后带入位姿的PoE公式中，我们就有，
$ product_(i=4)^6 exp(cal(S)_i^and theta_i) = [product_(i=1)^3 exp(cal(S)_i^and theta_i)]^(-1) T_d M^(-1) $

因为我们知道后三个关节的运动不改变腕部中心点坐标，因此，后三个关节的运动决定的是姿态，也就是得的是一个纯旋转矩阵，这个问题一般也存在解析解。

因此，我们知道了，当后三个关节轴线交于一点时，能够将逆运动学问题分解为*位置逆解*与*姿态逆解*的两个子问题，进而简化了求解过程，允许我们求得解析解。

一般来说，机械臂厂商会设计满足Pieper准则的6自由度机械臂，并针对自家产品开发针对性的逆运动学解析解求解算法，因此我们在使用中往往会选择调包。

=== 数值法

除了解析法，我们还会使用一些数值方法对逆运动学进行求解，我们将介绍几个比较常见的数值求解方法。

- *Newton–Raphson法*：将待求解的逆运动学方程写为一种通用格式，$g(theta)=0,theta in RR^n$，那么对函数进行一阶泰勒展开可以得到，
$ g(theta) = g(theta_0) + (partial g)/(partial theta)|_(theta=theta_0) (theta - theta_0) + omicron(theta^2) $

令其等于0，有，
$ theta = theta_0 - ((partial g)/(partial theta)|_(theta=theta_0))^(-1) g(theta_0) $

因而有迭代式，
$ theta_(i+1) = theta_i - ((partial g)/(partial theta)|_(theta=theta_i))^(-1) g(theta_i) $

可以不断计算到约定误差下，但是在奇点构型或非方阵的情况下，无法使用。

- *雅可比转置法*：雅可比转置法，其本质为梯度下降法的应用。在介绍这个方法前，首先讲述一下何为梯度下降法。

#figure(
  image("figure/fig10.png",width: 80%),
  caption: [梯度下降法]
)

梯度下降法被广泛应用于人工智能当中。所谓导数，即为斜率或者梯度方向，在地理上我们知道，等高线越密的地方，越陡峭，而用数学描述，等高线越密即为方向导数越大，而等高线最密的方向就是方向导数最大的方向，也就是梯度方向。假设我们要向山上走去，或者我们要下山，很显然，最近的路线一定是沿梯度的路线，因此，就有了*梯度下降法*。

假设我们现在有一个损失函数 $J(theta)$ (衡量模型拟合程度的，越小越好，此处符号并不代表雅可比矩阵)，那么我们希望找到一个特定的 $theta$ 满足使损失函数最小的条件，而这个过程，就是我们前面描述的下山的过程，用数学公式表达即为，
$ theta_(i+1) = theta_i - eta (partial J)/(partial theta) $

其中 $eta$ 表示步长，称作*学习率*，也就是说，沿着梯度方向走 $eta$ 距离，因为越到最值点，其梯度一定是越小(越平缓)，所以不断迭代下去，就会到达梯度为0的点，也就是我们要到达的损失函数最小值点。

现在理解了梯度下降法的原理，但是距离实现还有一定距离。因为在通常情况下，参与计算的是矩阵，而求导的自变量也有可能是矩阵，所以我们需要对矩阵求导有一定的了解，在此仅作简单介绍。

先考虑一个从向量映射到标量的函数，
$ f:RR^n->RR\
X in RR^n, f(X) in RR $
那么我们知道其微分为，
$ "d"f = sum_(i=1)^n (partial f)/(partial x_i) "d"x_i $

改写为矩阵形式，即，
$ "d"f = ((partial f)/(partial X))^T "d"X $

其中，$(partial f)/(partial X)$ 为我们要求的梯度，而 $"d"X$ 为，
$ "d"X = ("d"x_1, dots, "d"x_n)^T $

那么我们可以推广到自变量为矩阵时，微分为，
$ "d"f = tr(((partial f)/(partial X))^T "d" X) $

同时结合迹的循环性质，
$ tr(A B C) = tr(B C A) = tr(C A B) $

那么我们就可以很简单地求解矩阵的梯度了，接下来给出一个简单的例子。

#e[矩阵导数示例][
  $ alpha = bold(a)^T X bold(b) $
  对其两边同时求微分有，
  $ "d"alpha = tr(bold(a)^T "d" X bold(b)) = tr(bold(b)bold(a)^T"d" X ) $
  因此有，
  $ ((partial alpha)/(partial X))^T = bold(b)bold(a)^T $
  所以，我们可以得到矩阵导数为，
  $ (partial alpha)/(partial X) = bold(a) bold(b)^T $
]

接下来，我们将根据上述内容，推导雅可比转置法求解逆运动学的公式。设前向运动学函数为 $f:theta in RR^n|->x_d in RR^6 $，目标位姿为 $x_d$ (使用六维矢量进行表示)，定义需要最小化的误差函数为，
$ E(e) = 1/2 ||e(theta)||^2 = 1/2 e^T (theta) e(theta)\
e(theta) = x_d - f(theta) $

对 $E(e)$ 求微分有，
$ "d"E(e) = tr(1/2 "d"e^T (theta) e(theta) + 1/2 e^T (theta) "d"e(theta)) = tr(e^T (theta) "d"e(theta)) $

同时，$e(theta)$ 对 $theta$ 求微分有，
$ "d"e(theta) = - J(theta) "d"theta $

带进去可得，
$ "d"E(e) = tr(- e^T (theta) J(theta) "d"theta) $

同时有，
$ "d"E = tr(((partial E)/(partial theta))^T "d" theta) $

因此，
$ ((partial E)/(partial theta))^T = - J^T (theta) e(theta) $

由此，根据梯度下降法，我们得到了关节角度的迭代公式为，
$ theta_(i+1) = theta_i + eta J^T (theta_i) e(theta_i) $

=== 逆速度运动学

除了希望根据目标位姿求解关节角度，我们也希望通过目标速度求解关节速度，在很多情况下，我们需要对机械臂运动进行实时微调，这是与机械臂的关节速度实时相关的，因此我们需要研究逆速度运动学，让机器人的运动变得动态、平滑、响应迅速且智能。

因为我们知道末端速度 $dot(x)$ 满足如下关系式，
$ dot(x) = J(theta) dot(theta) $

那么目标即为，
$ dot(theta) = J^(-1) (theta) dot(x) $

不过对于矩阵求逆的问题，多半会有问题，所以接下来介绍几种更加方便通用的求解方法。



== Dynamics(动力学)

= 控制理论

= 轨迹生成

= 姿态估计

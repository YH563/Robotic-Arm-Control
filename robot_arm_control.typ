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

笔者以李理论和旋量理论作为笔记开篇，叙述在描述空间旋转运动时的常规方法的局限性，进而引出李理论，并基于李理论搭建旋量理论的框架。螺旋理论是描述机器人运动的基础，同时也为后续机器人运动学提供了统一的数学理论。

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
$ exp(xi^and) = mat(exp(phi.alt^and),G rho;0,1) $
其中矩阵 $G$ 满足如下关系式，
$ G = (sin abs(phi.alt))/abs(phi.alt) I + (1-(sin abs(phi.alt))/abs(phi.alt)) (phi.alt phi.alt^T)/abs(phi.alt)^2 + (1-cos abs(phi.alt))/abs(phi.alt)^2 phi.alt^and $

因此我们得到了在* $S E(3)$群上的指数映射表达式*，这也是研究机器人运动的重要工具。

最后我们再介绍几个指数映射的性质，证明留作练习。

#t[指数映射性质][
  + $"d"(exp(A t)) \/"d"t = A e^(A t) = e^(A t) A$
  + 若 $A = P D P^(-1)$，则有 $exp(A t) = P exp(D t) P^(-1)$
  + 若有 $A B = B A$，则有 $exp(A)exp(B) = exp(A + B)$ (由于矩阵乘法不具备对称性，所以常规指数映射的性质不再成立)
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
  $ ["Ad"_T] = mat(R, -R t^and;0, R) in RR^(6 times 6) $

  为将伴随表示与伴随表示矩阵做区分，我们对伴随表示矩阵添加了中括号。

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
$ cal(V)_s^and = "Ad"_T (cal(V)_b^and)\
mat(omega_s; v_s) = mat(R, 0;p^and R,R) mat(omega_b;v_b) $

根据上述推导，得到了在不同坐标系中旋量的定义，以及如何在空间坐标系和本体坐标系之间进行旋量的变换。

在这需要强调的是，这里的变换矩阵 $T=T_(s b)$，可以看到其与坐标系变换矩阵是“相反”的，因为我们是从本体坐标系变换到空间坐标系，因此变换矩阵是“反”的。

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

在旋量轴的定义下，我们可以推导一下其指数映射的结果。作为对 $S E(3)$ 群指数映射的回顾。

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

= 机器人运动学与静力学

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

设一个物体在 $n$ 维空间中的位置可以表示为 $x(theta) = f(theta(t))$，根据复合函数求导法则，可以得到，
$ dot(x)(theta) &= (partial f)/(partial theta) ("d"theta)/("d"t)\
&=J(theta)dot(theta) $

而其中，$J(theta) in RR^(m times n)$ 为*雅可比矩阵*，展开为，
$ J(theta) = mat((partial f_1)/(partial theta_1), dots, (partial f_1)/(partial theta_n);fence.dotted, dots.down, fence.dotted;(partial f_m)/(partial theta_1), dots, (partial f_m)/(partial theta_n)) $

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

=== Singularity Analysis(奇点分析)

== Inverse Kinematics(逆向运动学)

= 控制理论

= 轨迹生成

= 姿态估计

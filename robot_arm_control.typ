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
  $ M in {M in RR^(n times n)|M M^T=I, det(M)=1} $
]

那么我们知道在二维平面上的旋转矩阵可以写为，
$ M=mat(cos theta, -sin theta;sin theta, cos theta) $
或者将二维平面转为复平面，那么根据欧拉公式，旋转矩阵转为一维复数形式，即
$ M=e^(i theta) $

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

== 基础李理论

= 机器人运动学

== Forward Kinematics

== Inverse Kinematics

= 控制理论

= 轨迹生成

= 姿态估计

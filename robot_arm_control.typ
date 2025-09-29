#import "@preview/zh-kit:0.1.0" as zh-kit: setup-base-fonts,zhnumber-lower, zhnumber-upper 
#import "@preview/ergo:0.1.1": *
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

// 设置页面样式
#show : doc => setup-base-fonts(
  doc,
  cjk-mono-family: ("霞鹜文楷 Mono", "JetBrains Mono"),
)
#set math.mat(delim: "[")
#set par(leading: 1.2em, justify: true, first-line-indent: 2em)
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

为保证通俗易懂，本笔记尽量避开抽象的数学概念以及定理证明，主要针对在机器人控制中可能出现的数学公式以及引理进行论述，同时配合图片说明，以保证直观性。

笔者以李理论和螺旋运动作为笔记开篇，叙述在描述空间旋转运动时的常规方法的局限性，进而引出李理论，并基于李理论搭建螺旋理论的框架。螺旋理论是描述机器人运动的基础，同时也为后续机器人运动学提供了统一的数学理论。

#pagebreak()

#outline()

#pagebreak()

= 李理论与螺旋运动

== 空间旋转问题

== 基础李理论

= 机器人运动学

== Forward Kinematics

== Inverse Kinematics

= 控制理论

= 轨迹生成

= 姿态估计

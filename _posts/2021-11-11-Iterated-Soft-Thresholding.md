---
title: Iterated Soft Thresholding
author: Elsa Riachi
categories: [Notes]
tags: [sparsity]
date: 2021-11-08
math: true
---
<div style="display:none">
$$
\newcommand\testmacro[2]{\mathbf{F\alpha}(#1)^{#2}}
\def\norm#1{\left\|{#1}\right\|} % A norm with 1 argument
\newcommand\zeronorm[1]{\norm{#1}_0} % L0 norm
\newcommand\onenorm[1]{\norm{#1}_1} % L1 norm
\newcommand\twonorm[1]{\norm{#1}_2} % L2 norm
\def\<{\left\langle} % Angle brackets
\def\>{\right\rangle}
\newcommand\inner[1]{\langle #1 \rangle} % inner product
\newcommand\argmax{\mathop\mathrm{arg max}} % Defining math symbols
\newcommand\argmin{\mathop\mathrm{arg min}}
\newcommand{\pd}[2]{\frac{\partial{#1}}{\partial{#2}}}
\newcommand{\pdd}[2]{\frac{\partial^2{#1}}{\partial{#2}^2}}
$$

### Introduction

$$
\begin{equation}
\min_x \twonorm{Dx - y}^2 \quad \text { subject to } \zeronorm{x} \leq k
\end{equation}
$$

$$
\begin{equation}
\min_x \frac{1}{2}\twonorm{Dx - y}^2 \quad \text { subject to } \onenorm{x} \leq \tau
\end{equation}
$$

$$\begin{equation}
L(x, \lambda) = \frac{1}{2} \twonorm{Dx - y}^2  - \lambda (\tau - \onenorm{x})
\end{equation}$$

At the optimal solution (which is unique if $$D^TD$$ is full rank):

$$\begin{equation}
D^T(Dx - y)  + \lambda sign(x) = 0
\end{equation}$$


Proximal Method:

objective $$F(x) = f(x) + g(x)$$, with $$f(x)$$ continuously differentiable and $$g(x)$$ not differentiable.


Recall Newton's method for finding roots:
x_{k+1} = x_k - (D^T D)^{-1}(D^T(Dx - y)  + \lambda sign(x_k))
x_{k+1} =  (D^T D)^{-1}D^Ty  + \lambda sign(x_k)


$$F(x + p) \sim f(x) + \nabla f(x)^T p + \frac{1}{2}p^T H p + g(x + p)$$

Proximal operator with respect to g(x + p).

---
title: Spectral Clustering from the Min Cut Problem
author: Elsa Riachi
categories: [Notes]
tags: [graphs]
date: 2021-10-31
math: true
---

<div style="display:none">
Spectral clustering, often used as a community detection method, has its origins in the min-cut problem in graph theory. While the min-cut problem is NP-hard in the general case, a relaxation of the problem can be obtained from the graph Laplacian.
</div>

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

</div>

Spectral clustering, often used as a community detection method, has its origins in the min-cut problem in graph theory. While the min-cut problem is NP-hard in the general case, a relaxation of the problem can be obtained from the graph Laplacian.

We will start by introducing the objective for the min-cut problem which intuitively captures some notion of clusters or communities on graphs. Since this objective is NP-hard in the general case, we look at the derivation of an alternate objective based on an ideal case where the graph consists of $k$ connected components. In the following discussion we describe a cut $S$ as the set of edges that connect nodes in different partitions.

### Graph Cut Based on Conductance

The conductance objective for a cut S is as follows:

$$
\begin{equation}
\phi(A_1, A_2, ... A_k) = \sum_{l=1}^{k}\frac{\sum_{(i, j) \in A_l \times \bar{A}_l} w_{ij}}{vol(A_l)}
\end{equation}
$$

where $$vol(A_l)$$ is the total weighted degree (sum of weights of edges connected to nodes in $$A_l$$) of the nodes in $$A_l$$.

Note that for a fixed cut cost $$\sum_{(i, j) \in S} w_{ij}$$,  the conductance between the two partitions is higher if either $vol(A)$ or $$vol(B)$$ is too small. That is, we want to keep a substantial edge mass in both partitions.

Intuitively, computing the optimal cut in the general case is NP-hard, since there are exponentially many candidate graph cuts for a graph with $$N$$ edges. One can compare candidate solutions in polynomial time but there is no easy way to compute an improvement since the problem does not have a greedy or recursive sub-problem structure.

### Spectral Graph Partitioning

In the special case where the graph already consists of 2 connected, d-regular components, the optimal conductance is $0$, and we will show how to find the partitions by looking at the graph's adjacency matrix and graph Laplacian. We will then compare other cases to this ideal case and derive an approximate solution to the minimum conductance criterion.

#### The Eigenvectors of the Ideal Case:

For a d-regular graph, the most obvious eigenvector of the adjacency matrix $$A$$ is the all ones vector $$\mathbb{1}$$, with eigenvalue $d$ which is the largest eigenvalue of $$A$$.

For a graph with k disjoint d-regular components, the largest eigenvalue is still $$d$$, however this eigenvalue now has multiplicity $$k$$. The eigenspace of $$d$$ is spanned by eigenvectors consisting of all 1's at the indices corresponding to nodes in a single connected component. In this ideal case, if we know the eigenvectors corresponding to the largest eigenvalue of the adjacency matrix, we can obtain the $$k$$ connected components corresponding to the min-cut partitions.

#### What if the graph is connected?

What if we were to move a few edges to connect the separate partitions, while maintaining d-regularity? Now our graph is connected, and its largest eigenvalue is still $$d$$ corresponding to the eigenvector of all ones.

Now the eigenvalue $$d$$ has multiplicity $$1$$ but intuitively we can **partition** or permute the rows of $$A$$ into an almost block-diagonal matrix. Since by construction the off-block diagonal entries are limited compared to $$d$$  (remember conductance?), the second largest eigenvalue is close to $d$.

Also, since eigenvectors are orthogonal, the entries of the second largest eigenvector must sum to 0. So we have entries that are non-negative and entries that are negative. The sign of the entries can be used to determine the partition of the graph into two components.

#### What if we require more than 2 partitions?

Understanding the k-partition problem requires us to understand the graph Laplacian $$L = D - W$$ where $$D$$ is the weighted degree matrix and $$W$$ is the adjacency matrix.
Note that:

$$
\begin{align}
x^T L x &= x^T(D - W)x \\
&= \sum_{i}d_i x_i^2  - 2\sum_{(i, j)\ in E} w_{ij} x_i x_j \\
& = \sum_{(i, j)\ in E}\left( \sum_j w_{ij} x_i^2 + \sum_i w_{ji} x_j^2 - 2 w_{ij} x_i x_j \right)\\
& = \sum_{(i, j)\ in E} w_{ij}(x_i  - x_j)^2
\end{align}
$$

We can therefore think of the Laplacian as an operator on $$x \in \mathbb{R}^n$$ which is a vector of values assigned to nodes. If adjacent nodes are assigned similar values, then $$x^T L x$$ takes on a small value. Note that the all ones vector $$\mathbb{1}$$ is an eigenvector of the Laplacian matrix, with eigenvalue $$0$$.

Let \
$$cut(A_l, \bar{A}_l) := \sum_{(i, j) \in A_l \times \bar{A}_l} w_{ij}$$.

Note the similarity between the expressions for $$cut(A_l, \bar{A}_l)$$ and $$x^T L x$$.

In fact we can write $$cut(A_l, \bar{A}_l)$$ in terms of the Laplacian matrix, if we know how to define $$x$$.

Specifically if we set $$x$$ such that $$x_i = x_j$$ for any $$i, j$$ in $$A_l$$ or any $$i, j$$ in $$\bar{A}_l$$, then $$x^T L x$$ reduces to $$\sum_{(i, j) \in A_l \times \bar{A}_l} w_{ij} (x_i - x_j)^2$$.

If we also set
$$x_i = \frac{1}{\sqrt{vol(A_l)}}$$ if $$i \in A_l$$, and $$x_i = 0$$ otherwise, we obtain

$$
\begin{equation}
N \cdot vol(A) \cdot cut(A_l, \bar{A}_l) = x^T L x
\end{equation}
$$

The vector $$x$$ acts like an indicator vector for partition $$A_l$$, we can define such vectors $$x_l$$ for each partition obtaining a set of $$k$$ vectors in $$\mathbb{R}^n$$, $$\{x_l\}_{l=1}^{k}$$.

We obtain a different formulation for the min-cut problem:

$$
\begin{equation}
min_{(A_1, ... A_l)} \sum_{l=1}^{k} cut(A_l, \bar{A}_l) = min_{(A_1, ... A_l)} \sum_{l=1}^{k} x_l^T L x_l
\end{equation}
$$

Let $$H \in \mathbb{R}^{n \times k}$$ be a matrix whose columns are $$\{x_l\}_{l=1}^{k}$$. Then we can write:

$$
\begin{equation}
 \sum_{l=1}^{k} x_l^T L x_l = Tr(H^T L H)
\end{equation}
$$

Note that we still have a discrete constraint on $$H$$, where $$H_{ij} = \frac{1}{\sqrt{vol(A_j)}}$$ or $$H_{ij} = 0$$. Also $$vol(A_j)$$ depends on $$x_j$$.
We can relax this constraint by allowing the entries of $$H$$ to take any real values, as long as $$H^TH = I$$. The resulting relaxed objective for the min-cut problem becomes:

$$
\begin{equation}
min_{H: H^T D H = I} Tr(H^T L H)
\end{equation}
$$

Let $$T = D^{\frac{1}{2}}H$$

$$
\begin{equation}
min_{H: T^T T = I} Tr(T^T D^{-1\frac{1}{2}} L D^{-1\frac{1}{2}} T)
\end{equation}
$$

Using the Rayleigh quotient, we can show that the solution to the above problem consists of the top $$k$$ eigenvectors of $$D^{-1\frac{1}{2}} L D^{=1\frac{1}{2}}$$ which is the normalized Laplacian matrix. But the eigenvectors do not have the discrete structure corresponding to partitions which we defined above. Here is where we rely on heuristics.

Since $$T = D^{\frac{1}{2}}H$$, $$T$$ is obtained from $$H$$ by re-weighting the rows of $$H$$. Recall that the rows of $$H$$ (and in turn those of $$T$$) are by design expected to indicate the membership of each node. Once the rows of $$T$$ are normalized, we can cluster them using $$k$$-means and obtain a partition.

### References 
{% bibliography --file spectral_clustering %}

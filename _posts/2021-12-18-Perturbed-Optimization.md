---
title: Perturbed Optimization
author: Elsa Riachi
categories: [Notes]
tags: [optimization]
date: 2021-12-18
math: true
---

<div style="display:none">
Many forward computations whose parameters we would like to optimize may include discrete, or non-differentiable elements. Often, these forward procedures can be expressed as or approximated by the solution of a linear program. While linear programs are not differentiable with respect to their parameters either, their smooth approximations can be used to optimize over the underlying parameters.
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

### Introduction

Many forward computations whose parameters we would like to optimize may include discrete, or non-differentiable elements. Often, these forward procedures can be expressed as or approximated by the solution of a linear program. While linear programs are not differentiable either, their smooth approximations can be used to optimize over the underlying parameters.

Consider a Linear Program of the form:

$$
\begin{align}
F(\theta) &= \max _{y \in \mathcal{C}}\langle y, \theta\rangle \\
y^{*}(\theta) &= \underset{y \in \mathcal{C}}{\arg \max }\langle y, \theta\rangle
\end{align}
$$

Where $$C$$ is a convex polytope. The solution to the above LP is *almost* always uniquely determined by $$\theta$$ and corresponds to a vertex of $$C$$. To see why, consider the following inequality-form LP:

$$
\begin{equation}
\max _{y \in \mathcal{C}}\langle y, \theta\rangle,\ \text{ such that }  Q y \geq b
\end{equation}
$$

for $$y \in \mathbb{R}^n$$. And recall the following KKT condition, where $$\mathcal{A}$$ denotes the set of active inequality constraints:

$$
\begin{align}
\theta &= - Q^T \lambda^* \\  
\lambda_i^* &\neq 0 \ \ \forall i \in \mathcal{A}
\end{align}
$$

In general, $$n$$ columns of $$Q^T$$ are required to construct the n-dimnesional $$\theta$$. So $$\lambda^*$$ should have $$n$$ non-zero entries. This means that $$y^*$$ is the solution of a system of $$n$$ linear equations $$\bar{Q}y = \bar{b}$$, where $$\bar{Q}$$ and $$\bar{b}$$ denote a subset of the rows of $$Q$$ and $$b$$ respectively. Effectively, $$y^*$$ is the intersections of $$n$$ sides of the polytope, corresponding to a vertex point.

However, the maximizer $$y^{*}$$ is not a differentiable function of $$\theta$$. This is because, small perturbations in $$\theta$$ do not change the active constraints at $$y^{*}$$.  While the entries $$\lambda_i$$ might be different, the effective system of $$n$$ linear equations remains the same, unless a large enough perturbation to $$\theta$$ is applied. This means that the Jacobian of $$y^*$$ with respect to $$\theta$$ is zero for almost all $$\theta$$, and undefined at the boundaries of normal cones. This is inconvenient when we wish to optimize the parameter $$\theta$$ for some objective, as is often the case in non-differentiable learning or inverse problems.

<figure>
    <img src="/assets/img/optimization/normal_fan.png" width="650">
    <figcaption>Figure 1: The normal fan of the polytope consists of the regions in which the cost parameters result in the same maximizing vertex.</figcaption>
</figure>

Berthet et al. \cite{} propose to *smooth out* the above LP. The intuition behind this approach is as follows. Since we want $$y^{*}$$ to vary smoothly in $$\mathcal{C}$$ with respect to $$\theta$$, we can consider an alternative solution which is a convex combination of the $$v$$ vertices of $$\mathcal{C}$$,  


$$\begin{equation}
y^{*}(\alpha) = \sum_{i=1}^{v} \alpha_i y_i \text{ where  } \alpha_i \in [0, 1]
\end{equation}
$$

 where $$\alpha \in \mathbb{R}^{v}$$ is a function of $$\theta$$ which will be revealed shortly. Intuitively, looking at the above figure, as $$\theta$$ is moved clock-wise, we expect $$y^*(\alpha)$$ to move closer to $$y_2$$. In a way, we want all vertices to contribute to the resulting estimate, while keeping the estimate $$y^{*}(\alpha)$$ close to the solution. Note that $$\alpha$$ effectively defines a probability mass function $$p_{\theta}$$ over the possible solutions $$\{y_i\}_{i=1}^{v}$$.  

The authors propose to achieve a smooth approximation to $$y^{*}(\theta)$$ by taking the expectation of $$y^{*}(\theta + \epsilon Z)$$ where $$\epsilon$$ controls the size of the added noise, and $$Z$$ is a random variable sampled from a Gibbs distribution $$\mu(Z) \propto e^{-\nu(Z)}$$.

$$
\begin{equation}
y^{*}_\epsilon (\theta) = \mathbb{E}_z \mathrm{arg max}_{y \in C} \langle y,\theta + \epsilon Z \rangle
\end{equation}
$$


$$
\begin{equation}
F_\epsilon (\theta) = \mathbb{E}_z \mathrm{max}_{y \in C} \langle y,\theta + \epsilon Z \rangle
\end{equation}
$$


The gradient of $$F(\theta) = \max_{y \in C}\langle y, \theta \rangle$$ with respect to $$\theta$$ is $$y^*$$, as long as $$\theta$$ remains in the normal cone where $$y^*$$ is optimal. It is undefined on the boundaries of normal cones. There is a similar relationship between $$\nabla_\theta F_{\epsilon}(\theta)$$ and $$y^*_{\epsilon}(\theta)$$:

$$
\begin{align}
\nabla_{\theta}F_{\epsilon}(\theta) &= \nabla_{\theta} \mathbb{E}_z \mathrm{max}_{y \in C} \langle y,\theta + \epsilon Z \rangle = \mathbb{E}_z \nabla_{\theta} \mathrm{max}_{y \in C} \langle y,\theta + \epsilon Z \rangle \\
\\
&= \mathbb{E}_z \nabla_{\theta} \langle y^*(\theta + \epsilon Z), \theta \rangle = \mathbb{E}_z y^*(\theta + \epsilon Z) = y_{\epsilon}^*(\theta)
\end{align}
$$

This leads to the fact that:

$$\begin{equation}
\nabla_{\theta} y_{\epsilon}^*(\theta) = \nabla_{\theta}^2 F_\epsilon(\theta)
\end{equation}
$$

### Derivatives of Soft Maximizers
To see why the smooth approximations are differentiable with respect to $$\theta$$, consider the integral which evaluates $$F_\epsilon (\theta)$$.

$$
\begin{align}
F_\epsilon (\theta) &= \int_{Z} \max_{y \in C}\langle y, \theta + \epsilon z \rangle \mu(z) \, dz \\
\\
\nabla_{\theta}  F_\epsilon (\theta) &= \int_{Z} \nabla_{\theta} \max_{y \in C}\langle y, \theta + \epsilon z \rangle \mu(z) \, dz \\
\end{align}
$$


Since the integral boundaries don't depend on $$\theta$$ all we have to do is differentiate the integrand and integrate the gradient. Since $$\max_{y \in C}\langle y, \theta + \epsilon z \rangle$$ is not everywhere differentiable, we consider $$z' = \theta + \epsilon z$$ to be a fixed parameter vector, so now $$z$$ depends on $$z'$$ and $$\theta$$. A simple change of variable results in:

$$
\begin{align}
 \nabla_{\theta}  F_\epsilon (\theta) &= \int_{Z'} \max_{y \in C}\langle y, z' \rangle \nabla_{\theta}  \mu\left(\frac{z' - \theta}{\epsilon}\right) \, \frac{dz'}{\epsilon} \\
 \\
 &= \int_{Z'} \max_{y \in C}\langle y, z' \rangle \mu\left(\frac{z' - \theta}{\epsilon}\right)  \nabla_{\theta} \nu\left(\frac{z' - \theta}{\epsilon}\right) \, \frac{dz'}{\epsilon} \\
 \\
 &= \mathbb{E} \left[ F(\theta + \epsilon Z) \nabla_{Z} \nu\left(Z \right) / \epsilon \right]
\end{align}
$$

It follows that:

$$
\begin{align}
\nabla_{\theta} y_{\epsilon}^{*}(\theta) &= \nabla_{\theta}^{2} F_{\epsilon}(\theta) \\
\\
\nabla_{\theta} y_{\epsilon}^{*}(\theta) &= \nabla_{\theta}^{2} \int_{Z} \max_{y \in C}\langle y, \theta + \epsilon z \rangle \mu(z) \, dz \\
\\
&=  - \int_{Z'} \max_{y \in C}\langle y, z' \rangle \nabla^{2}_{\theta}  \mu\left(\frac{z' - \theta}{\epsilon}\right) \, \frac{dz'}{\epsilon} \\
\\
&= - \int_{Z'} \max_{y \in C}\langle y, z' \rangle \nabla_{\theta} \left[ \mu\left(\frac{z' - \theta}{\epsilon}\right) \nabla_{\theta} \nu\left(\frac{z' - \theta}{\epsilon}\right) \right] \, \frac{dz'}{\epsilon} \\
\\
&= \mathbf{E} \left[F(\theta+\varepsilon Z)\left(\nabla_{z} \nu(Z) \nabla_{z} \nu(Z)^{\top}-\nabla_{z}^{2} \nu(Z)\right) / \varepsilon^{2} \right]
\end{align}
$$

### Monte-Carlo Estimates of the Gradients

Using the following derived expression:

$$
\begin{equation}
\nabla_{\theta} y_{\epsilon}^{*}(\theta) = \mathbf{E} \left[F(\theta+\varepsilon Z)\left(\nabla_{z} \nu(Z) \nabla_{z} \nu(Z)^{\top}-\nabla_{z}^{2} \nu(Z)\right) / \varepsilon^{2} \right]
\end{equation}
$$

we see that we can estimate the derivative of $$\nabla_{\theta} y_{\epsilon}^{*}(\theta)$$ by:

- sampling $$Z_i \sim \mu(Z)$$ for $$ i = 1, ... ,m$$
- solving each LP $$\max_{y \in C} \langle y, \theta + \epsilon Z_i \rangle$$
- computing $$ \left(\nabla_{z} \nu(Z) \nabla_{z} \nu(Z)^{\top}-\nabla_{z}^{2} \nu(Z)\right) / \varepsilon^{2}$$  at $$Z = Z_i$$.


### Fenchel-Young Losses (??)

First, we write $$F_\epsilon (\theta)$$ with respect to $$\theta$$ as a convex conjugate. Let $$\Omega(y)$$ be the convex conjugate of $$F_1 (\theta)$$. A simple change of variable, where $$\theta' := \epsilon \theta$$ shows that the convex conjugate of $$F_{\epsilon}(\theta)$$ is $$\epsilon \Omega(y)$$.

$$
\begin{align}
\Omega(y) &= \sup_{\theta}\{\langle \theta, y \rangle - F_1(\theta)\} \\
\\
\epsilon \Omega(y) &= \sup_{\theta} \{\langle \epsilon \theta, y \rangle - \mathbb{E}_Z \mathrm{max}_{y \in C} \langle y, \epsilon \theta + \epsilon Z \rangle\} \\
\\
\epsilon \Omega(y) &= \sup_{\theta'} \{\langle \theta', y \rangle - \mathbb{E}_Z \mathrm{max}_{y \in C} \langle y, \theta' + \epsilon Z \rangle\} \\
\\
\epsilon \Omega(y) &= F^{*}_{\epsilon}(y)
\end{align}
$$

Then using the fact that $$F_{\epsilon}(\theta)$$ is strictly convex:

$$
\begin{align}
F_{\epsilon}(\theta)= \sup_{y} \{\langle y, \theta \rangle - \epsilon \Omega(y)\}
\end{align}
$$

The supremum is achieved at $$y_{\epsilon}^{*}(\theta)$$ where $$\epsilon \nabla_{y}\Omega(y_{\epsilon}^{*}) = \theta$$. Using the property $$\nabla f^{*} = (\nabla f)^{-1}$$,

$$
\begin{align}
\left( \epsilon \nabla_{y} \Omega(y_{\epsilon}^{*}) \right)^{-1} &= \nabla_{\theta} F_{\epsilon}(\theta) \\
\\
y_{\epsilon}^{*}(\theta) &= \nabla_{\theta} F_{\epsilon}(\theta)
\end{align}
$$


Consider the problem of maximizing the likelihood of observations $$\{y_i\}$$, with respect to a parameterized model $$p_{\theta}(y)$$, where $$p_{\theta}$$ is a Gibbs distribution. The empirical loss is shown below.

$$
\begin{align}
\bar{\ell}_{n}(\theta) &= \frac{1}{n} \sum_{i=1}^{n} \log p_{\text {Gibbs }, \theta}\left(y_{i}\right) \\
\\
 &= \frac{1}{n} \sum_{i=1}^{n}\left\langle y_{i}, \theta\right\rangle-\log Z(\theta)\\
  \text { with } \nabla_{\theta} \bar{\ell}_{n}(\theta) &= \frac{1}{n} \sum_{i=1}^{n} y_{i}-\mathbf{E}_{\text {Gibbs }, \theta}[Y]
\end{align}
$$

Notice that setting $$\nabla_{\theta} \bar{\ell}_{n}(\theta) = 0$$ points to a moment matching procedure to fit $$\theta$$. However,

$$
\begin{align}
F_{\epsilon}(\theta) &= \mathbb{E}_Z \max_{y} \langle y, \theta + \epsilon Z \rangle \\
                     &\geq \mathbb{E}_Z \log{ \mathbb{E}_Y exp\{\langle Y, \theta + \epsilon Z \rangle\} } \\
                     &\geq Z_{\epsilon}(\theta) \\
\end{align}
$$

So $$\left\langle y_{i}, \theta\right\rangle- F_{\epsilon}(\theta)$$ is a lower bound on $$\ell}_{n}(\theta)$$.


### How Good are MC Estimates of the Gradients?

To answer this question we need to develop asymptotic results about the MC estimates of $$\nabla_{\theta} y_{\epsilon}^{*}(\theta)$$ and $$\nabla_{\theta} F_\epsilon (\theta)$$. But before we do, since asymptotic results depend on properties of function of the random variable, we should first characterize $$F_\epsilon (\theta)$$ and $$y_{\epsilon}^{*}(\theta)$$.

**Properties of $$F_{\epsilon}$$**:

1. $$F_{\epsilon}$$ is strictly convex (with respect to $$\theta$$).

*Proof*:

Since $$F$$ is the supremum of a set of linear functions, it is convex. Let $$\lambda \in \[0, 1\]$$,  $$\theta_{\lambda}=\lambda \theta+(1-\lambda) \theta^{\prime}$$. So:

$$
\begin{equation}
\lambda F_{\varepsilon}(\theta)+(1-\lambda) F_{\varepsilon}\left(\theta^{\prime}\right)=\mathbf{E}\left[\lambda F(\theta+\varepsilon Z)+(1-\lambda) F\left(\theta^{\prime}+\varepsilon Z\right)\right] \geqslant \mathbf{E}\left[F\left(\lambda \theta+(1-\lambda) \theta^{\prime}+\varepsilon Z\right)\right]=F_{\varepsilon}\left(\theta_{\lambda}\right)
\end{equation}
$$

To show that $$F_{\epsilon}$$ is strictly convex, we need to show that the inequality above is strict for all values of $$\theta, \ \theta^{\prime} \in \mathbb{R}^{d}$$. Note that the inequality applies to the arguments inside the expectation. Therefore, for the inequality to hold with equality, it should hold for *almost* all values of $$Z$$.

If  $$\lambda F(\theta+\varepsilon z)+(1-\lambda) F\left(\theta^{\prime}+\varepsilon z\right) \geqslant F\left(\lambda \theta+(1-\lambda) \theta^{\prime}+\varepsilon z\right)$$, for almost all $$z$$, then $$F$$ is linear over the segment $$[\theta + \epsilon z, \theta^{\prime} + \epsilon z]$$ for almost all $$z \in \mathbb{R}^{d}$$. Since $$F$$ is the result of a linear program, $$F$$ is the value obtained by $$\langle \theta + \epsilon z, y^{*} \rangle$$ at $$y^{*}$$ a vertex point in $$\mathcal{C}$$. So if $$F$$ is linear in $$[\theta, \theta^{\prime}]$$ then $$y^{*}(\theta + \epsilon z)$$ remains the same as $$\theta$$ is varied on the line $$[\theta, \theta^{\prime}]$$, for almost all $$z$$. For this to hold, $$\theta + \epsilon z$$ and $$\theta^{\prime} + \epsilon z$$ need to be in the same normal cone of $$\mathcal{C}$$ for almost all $$z$$. This can only happen if $$\theta - \theta^{\prime} = 0$$. So $$F_{\epsilon}$$ is strictly convex.

2. $$F_{\epsilon}$$ is $$R_C-Lipschitz$$.

*Proof*:
Since $$F_{\epsilon}$$ is the supremum of a set of linear functions, and $$\twonorm{y} \leq R_C$$, then $$| F_{epsilon}(\theta) - F_{\epsilon}(\theta^{\prime}) | \leq R_C (\theta - \theta^{\prime})$$.

3. $$F_{\epsilon}$$ is $$R_C M_{\mu} / \epsilon$$-gradient Lipschitz.

*Proof*:

Recall that $$ \nabla_{\theta}  F_\epsilon (\theta) = \mathbb{E} F(\theta + \epsilon Z) \frac{\nabla_{z}\mu(z)}{\epsilon}$$. Then,

$$
\begin{align}
\twonorm{F_\epsilon (\theta) - F_\epsilon (\theta^{\prime})} &= \twonorm{ \mathbb{E} \left[ F(\theta + \epsilon Z) -  F(\theta^{\prime} + \epsilon Z) \frac{\nabla_{z}\mu(z)}{\epsilon} \right]}\\
\end{align}
$$

Using Jensen's inequality since the $l2$ norm is convex, and the Cauchy-Schwarz inequality we obtain:

$$
\begin{align}
& \leq \mathbb{E} \twonorm{ F(\theta + \epsilon Z) - F(\theta^{\prime})} \mathbb{E} \twonorm{\frac{\nabla_{z}\mu(z)}{\epsilon}} \\
& \leq R_C (\theta - \theta^{\prime}) M_{\mu} / {\epsilon} \\
\end{align}
$$

**Asymptotic Convergence of $$\

# Tractography algorithms

```@contents
Pages = ["algos.md"]
Depth = 3
```

The library provides several algorithms to compute the streamlines. More precisely, assume that we are given a field of fibre orientation distribution functions (FODF) 

$$\bm u\to ODF(\bm x,\bm u)$$

where $ODF(\bm x,\cdot)$ is the distribution of directions at position $\bm x$.

In practice, to prevent sharp direction changes, the full FODF is not considered at a given step. Typically, a new direction is obtained by sampling the FODF in a cone around the previous direction.  This can represented as a modification of the incoming direction $\bm{u} \in \mathbb{S}^2$

$$g(\bm{x}, \bm{u}, \bm{u}') = \frac{1}{N(\bm{x}, \bm{u})}f(\bm{x}, \bm{u}')c(\bm{u}, \bm{u}')$$

where $c(\bm{u}, \bm{u}') \in \mathbb{R}^+$ is non--zero in a cone around $\bm{u}$ and $N(\bm{x}, \bm{u}) \in \mathbb{R}^+$ is a normalization factor that ensures

$$\int_{\mathbb{S}^2}g(\bm{x}, \bm{u}, \bm{u}')d\bm{u}' = 1.$$

!!! tip "Cone"
    The cone function is passed to a `TMC` via [`Tractography.TMC`](@ref). The `TMC` also determine how the FODF are computed, see [SH evaluation](@ref sheval).

## 1. Deterministic

The algorithm `alg = Deterministic()` (see [`Tractography.Deterministic`](@ref)) implements the following situation. We compute a sequence $(\bm x_i, \bm u_i)_i$ such that

$$\bm x_{i+1} = \bm x_i + \Delta s \bm u_i$$

$$\bm u_{i+1} = \argmax g(\bm x_i, \cdot)$$

## 2. Cumulative sum distribution

The algorithm `alg = CSD()` (see [`Tractography.CSD`](@ref)) implements the following situation. We compute a sequence $(\bm x_i, \bm u_i)_i$ such that

$$\bm x_{i+1} = \bm x_i + \Delta s \bm u_i$$

$$\bm u_{i+1} \sim \argmax g(\bm x_i, \cdot)$$

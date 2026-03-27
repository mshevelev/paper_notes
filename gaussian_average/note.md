# Gaussian Variables

Let us introduce $n$ Gaussian variables $x_1, \dots, x_n$ with covariance matrix $\Sigma$. We denote the variances of these variables as $\sigma_i^2$, which correspond to the diagonal elements of the covariance matrix (i.e., $\Sigma_{ii} = \sigma_i^2$).

Next, we introduce the standardized variables $$y_i = \frac{x_i}{\sigma_i}$$


Finally, we introduce a variable $b$ defined as a weighted sum of all $x_i$, and normalized to have a unit variance:

$$
b = \frac{\sum_{i=1}^n w_i x_i}{\sqrt{\mathbf{w}^T \Sigma \mathbf{w}}}
$$

Next, we introduce variables $z_i$ as follows:
$$
z_i = y_i -  \underbrace{\text{corr}(y_i, b)}_{\rho_i} \cdot b
$$

$\textcolor{red}{Th1}$. $z_i \perp b$ and $\text{Var}(z_i) = 1 - \rho_i^2$

$$
\text{Cov}(z_i, b) = \text{Cov}(y_i - \rho_i b, b) = \text{Cov}(y_i, b) - \rho_i \underbrace{\text{Cov}(b, b)}_{\text{Var}(b)=1} = \rho_i - \rho_i \cdot 1 = 0
$$

$\textcolor{red}{Th2}$. $\sum_{i=1}^n w_i \sigma_i z_i = 0$. 
$$
\begin{align*}
\sum_{i=1}^n w_i \sigma_i z_i &= \sum_{i=1}^n w_i \sigma_i y_i - b \sum_{i=1}^n w_i \sigma_i \rho_i \\
&= \sum_{i=1}^n w_i x_i - b \sum_{i=1}^n w_i \sigma_i \text{Cov}\left(\frac{x_i}{\sigma_i}, b\right) \\
&= b \sqrt{\mathbf{w}^T \Sigma \mathbf{w}} - b \sum_{i=1}^n w_i \text{Cov}\left(x_i, \frac{\sum_k w_k x_k}{\sqrt{\mathbf{w}^T \Sigma \mathbf{w}}}\right) \\
&= b \sqrt{\mathbf{w}^T \Sigma \mathbf{w}} - b \sum_{i=1}^n w_i \frac{\sum_k w_k \text{Cov}(x_i, x_k)}{\sqrt{\mathbf{w}^T \Sigma \mathbf{w}}} \\
&= b \sqrt{\mathbf{w}^T \Sigma \mathbf{w}} - b \frac{\sum_{i=1}^n \sum_{k=1}^n w_i w_k \text{Cov}(x_i, x_k)}{\sqrt{\mathbf{w}^T \Sigma \mathbf{w}}} \\
&= b \sqrt{\mathbf{w}^T \Sigma \mathbf{w}} - b \frac{\mathbf{w}^T \Sigma \mathbf{w}}{\sqrt{\mathbf{w}^T \Sigma \mathbf{w}}} \\
&= b \sqrt{\mathbf{w}^T \Sigma \mathbf{w}} - b \sqrt{\mathbf{w}^T \Sigma \mathbf{w}} \\
&= 0
\end{align*}
$$

**Geometric Intuition for $\text{Th2}$**

We can understand $\text{Th2}$ elegantly through the lens of linear algebra and geometry:
- Think of the random variables $y_i$ as vectors in a Hilbert space where the inner product corresponds to covariance: $\langle A, B \rangle = \text{Cov}(A, B)$.
- In this space, the variable $b \propto \sum w_i \sigma_i y_i$ is just a specific vector (representing the aggregate portfolio).
- The residual $z_i = y_i - \text{proj}_b(y_i)$ represents the exact component of $y_i$ that is orthogonal to $b$. 
- Therefore, each $z_i$ vector lies perfectly in the subspace perpendicular to $b$.

When we compute $\sum w_i \sigma_i z_i$, we are evaluating the exact same linear combination that constructed $b$ in the first place, but we are applying it to the orthogonal residuals instead. 
Since projection is a linear operator, we can swap the sum and the projection:

$$ \sum_{i=1}^n w_i \sigma_i z_i = \sum_{i=1}^n w_i \sigma_i y_i - \text{proj}_b\left(\sum_{i=1}^n w_i \sigma_i y_i\right) $$

Notice that the term $\sum w_i \sigma_i y_i$ inside the projection is actually just a scaled version of $b$ itself! 

The residual of any vector projected onto itself is exactly zero—it has no orthogonal component to itself. Thus, when we take the orthogonal residuals ($z_i$) and recombine them using the precise weights ($w_i \sigma_i$) that form the original vector $b$, the resulting sum must perfectly collapse to zero.
### Simultaneous Decomposition of Residuals

Suppose the aggregate portfolio $b$ is constructed as a linear combination of two sub-portfolios $b_1$ and $b_2$. We can formalize this decomposition as:
$$ b = a_1 b_1 + a_2 b_2 $$
where we carefully define our notations:
- **$b_1 \perp b_2$**: The two sub-portfolios are strictly orthogonal in our Gaussian space, meaning their covariance is exactly zero ($\text{Cov}(b_1, b_2) = 0$).
- **Unit Variances**: We normalize the variables such that they all have a variance of 1 ($\text{Var}(b) = \text{Var}(b_1) = \text{Var}(b_2) = 1$). 
- **Convex Allocation**: To maintain $\text{Var}(b) = 1$ given the orthogonality of $b_1$ and $b_2$, the allocation weights must naturally satisfy the geometric identity $a_1^2 + a_2^2 = 1$.

Recall that for any variable $y_i$, its residual with respect to $b$ (the component strictly orthogonal to the aggregate portfolio) is defined as:
$$ z_i = y_i - \text{proj}_b(y_i) = y_i - \rho_i b $$
where $\rho_i = \text{Cov}(y_i, b)$ is the standard projection weight.

Similarly, we can define the local residuals of $y_i$ with respect to the individual sub-portfolios:
- **$z_i^{(b_1)} = y_i - \rho_{i1} b_1$**: The residual from projecting $y_i$ strictly onto $b_1$, where $\rho_{i1} = \text{Cov}(y_i, b_1)$.
- **$z_i^{(b_2)} = y_i - \rho_{i2} b_2$**: The residual from projecting $y_i$ strictly onto $b_2$, where $\rho_{i2} = \text{Cov}(y_i, b_2)$.

We want to express the total aggregate residual $z_i$ as a simultaneous combination of the local, partial residuals $z_i^{(b_1)}$ and $z_i^{(b_2)}$.

#### A Simple Proof

Let us construct a convex combination of the two local residuals using the squared allocation weights $a_1^2$ and $a_2^2$. Since $a_1^2 + a_2^2 = 1$, we can expand this combination as follows:

$$ a_1^2 z_i^{(b_1)} + a_2^2 z_i^{(b_2)} = a_1^2 (y_i - \rho_{i1} b_1) + a_2^2 (y_i - \rho_{i2} b_2) $$
$$ = (a_1^2 + a_2^2) y_i - \Big( a_1^2 \rho_{i1} b_1 + a_2^2 \rho_{i2} b_2 \Big) $$
$$ = y_i - \Big( a_1^2 \rho_{i1} b_1 + a_2^2 \rho_{i2} b_2 \Big) $$

Now, let's analyze the true residual $z_i$. By definition, $z_i = y_i - \text{proj}_b(y_i)$. Because covariance itself is a linear operator, we know that $\rho_i = \text{Cov}(y_i, a_1 b_1 + a_2 b_2) = a_1 \rho_{i1} + a_2 \rho_{i2}$. We can fully expand the projection onto the aggregate portfolio $b$:

$$ \text{proj}_b(y_i) = \rho_i b = (a_1 \rho_{i1} + a_2 \rho_{i2}) (a_1 b_1 + a_2 b_2) $$
$$ = a_1^2 \rho_{i1} b_1 + a_2^2 \rho_{i2} b_2 + a_1 a_2 (\rho_{i1} b_2 + \rho_{i2} b_1) $$

Notice a beautiful structural alignment: the first two terms in this true projection exactly match the projection components hidden inside our convex combination derived earlier! 

If we substitute this expanded projection back into our definition $z_i = y_i - \text{proj}_b(y_i)$, we get:
$$ z_i = \underbrace{y_i - \Big( a_1^2 \rho_{i1} b_1 + a_2^2 \rho_{i2} b_2 \Big)}_{\text{convex combination result}} - \underbrace{a_1 a_2 (\rho_{i1} b_2 + \rho_{i2} b_1)}_{\text{missing cross-terms}} $$

Substituting the convex combination from our first step into this equation yields the final, elegant simultaneous decomposition:
$$ \textcolor{blue}{ z_i = a_1^2 z_i^{(b_1)} + a_2^2 z_i^{(b_2)} - a_1 a_2 (\rho_{i1} b_2 + \rho_{i2} b_1) } $$

#### Geometric Intuition

Geometrically, how do we build the total residual $z_i$ from local information? We start by simply taking a variance-weighted average of the local sub-portfolio residuals ($z_i^{(b_1)}$ and $z_i^{(b_2)}$). 

However, because the full portfolio $b$ sits diagonally between $b_1$ and $b_2$, projecting directly onto $b$ inherently generates mixed terms. Because the local projections happen completely isolated along the strictly orthogonal axes of $b_1$ and $b_2$, simply averaging them trivially misses this cross-axial "leakage". 

The strict correction factor $-a_1 a_2 (\rho_{i1} b_2 + \rho_{i2} b_1)$ exactly subtracts out these overlapping, mixed projections, guaranteeing that $z_i$ mathematically preserves its zero-covariance orthogonality to the aggregate portfolio vector $b$.

#### Generalization to $K$ Orthogonal Sub-Portfolios

This result flawlessly and elegantly extends to any number of orthogonal sub-portfolios! 

Suppose we aggregate $K$ strictly orthogonal components:
$$ b = \sum_{k=1}^K a_k b_k $$
where $b_j \perp b_k$ for all $j \ne k$, all components have unit variance, and the allocation weights naturally satisfy $\sum_{k=1}^K a_k^2 = 1$. Let the local residual for the $k$-th sub-portfolio be $z_i^{(b_k)} = y_i - \rho_{ik} b_k$, where $\rho_{ik} = \text{Cov}(y_i, b_k)$.

If we construct the same $a^2$-squared convex combination across all $K$ residuals, we find:
$$ \sum_{k=1}^K a_k^2 z_i^{(b_k)} = \sum_{k=1}^K a_k^2 (y_i - \rho_{ik} b_k) = \left( \sum_{k=1}^K a_k^2 \right) y_i - \sum_{k=1}^K a_k^2 \rho_{ik} b_k = y_i - \sum_{k=1}^K a_k^2 \rho_{ik} b_k $$

When substituting this into the true aggregate residual $z_i$, we expand the total projection $\rho_i b$:
$$ \text{proj}_b(y_i) = \left( \sum_{j=1}^K a_j \rho_{ij} \right) \left( \sum_{k=1}^K a_k b_k \right) = \underbrace{ \sum_{k=1}^K a_k^2 \rho_{ik} b_k }_{\text{diagonal terms}} + \underbrace{ \sum_{j \ne k} a_j a_k \rho_{ij} b_k }_{\text{off-diagonal cross-terms}} $$

Just as before, the diagonal unmixed terms exactly match the partial projection extracted by our convex combination! By grouping the off-diagonal terms symmetrically into unique $(j, k)$ pairs where $j < k$, we derive the fully generalized simultaneous decomposition:

$$ \textcolor{blue}{ z_i = \sum_{k=1}^K a_k^2 z_i^{(b_k)} - \sum_{j < k} a_j a_k (\rho_{ij} b_k + \rho_{ik} b_j) } $$

**Geometric intuition extends perfectly:** The aggregate residual is entirely the variance-weighted convex combination of the $K$ distinct local residuals, strictly minus all unique pairwise "cross-leakages" between every orthogonal subspace pair.

#### Using Strict Orthogonalized Components ($b_k^{ort}$)

We can alternatively reformulate the entire $K$-dimensional decomposition by explicitly utilizing the orthogonalization of each individual sub-portfolio $b_k$ with respect to the aggregate portfolio $b$.

Let $b_k^\perp = b_k - a_k b$ be the unnormalized component of $b_k$ that is strictly orthogonal to $b$. Because $b$ naturally scales to $1$ and $b_k$ sits perpendicularly to the rest of the mix, the variance of this perpendicular component resolves beautifully to $1 - a_k^2$. The fully normalized unit-variance component is thus $b_k^{ort} = \frac{b_k - a_k b}{\sqrt{1 - a_k^2}}$.

Starting again from our residual alignment $z_i - \sum_{k=1}^K a_k^2 z_i^{(b_k)} = \sum_{k=1}^K a_k^2 \rho_{ik} b_k - \rho_i b$, we carefully separate $b_k$ inside the sum into its scaled trace along the aggregate portfolio ($a_k b$) and its completely orthogonal complement ($b_k^\perp$):

$$ \sum_{k=1}^K a_k^2 \rho_{ik} b_k = \sum_{k=1}^K a_k^2 \rho_{ik} (b_k^\perp + a_k b) = \sum_{k=1}^K a_k^2 \rho_{ik} b_k^\perp + \left( \sum_{k=1}^K a_k^3 \rho_{ik} \right) b $$

Substituting this back bundles all of the strictly non-orthogonal terms beautifully onto the exact $b$ axis:
$$ z_i = \sum_{k=1}^K a_k^2 z_i^{(b_k)} + \sum_{k=1}^K a_k^2 \rho_{ik} b_k^\perp - \left( \rho_i - \sum_{k=1}^K a_k^3 \rho_{ik} \right) b $$

Since covariance is completely linear ($\rho_i = \sum a_k \rho_{ik}$), the correction acting exclusively along the $b$ axis mathematically isolates:
$$ \rho_i - \sum_{k=1}^K a_k^3 \rho_{ik} = \sum_{k=1}^K a_k (1 - a_k^2) \rho_{ik} $$

Substituting this back and scaling our unnormalized $b_k^\perp$ vectors smoothly to up the strict normalized vectors $b_k^\perp = \sqrt{1 - a_k^2} \, b_k^{ort}$ guarantees the final magnificent shape:

$$ \textcolor{blue}{ z_i = \sum_{k=1}^K a_k^2 z_i^{(b_k)} + \sum_{k=1}^K a_k^2 \sqrt{1-a_k^2} \, \rho_{ik} b_k^{ort} - \left( \sum_{k=1}^K a_k (1 - a_k^2) \rho_{ik} \right) b } $$

This reconstructs the total orthogonal residual exactly across three uniquely segregated dimensions:
1. **Local Baselines**: The variance-weighted convex combination directly driving off the $K$ isolated local residuals ($\sum a_k^2 z_i^{(b_k)}$).
2. **Orthogonal Shifts**: The sum acting purely along normalized, orthogonal individual spaces $b_k^{ort}$ that accounts for orthogonal dispersion against $b$.
3. **Line-of-Sight Anchorage**: The single consolidated correction pulling back strictly along linearly $b$ to secure the non-negotiable definition that $z_i \perp b$.

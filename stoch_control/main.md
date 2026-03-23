# Stochastic Control Optimization for Autocorrelated Returns

## Model Setup

We consider a process of returns $R_t$ composed of a forecast signal $F_t$ and white noise $\epsilon_t$:
$$
R_t = F_t + \epsilon_t,
$$
where $E[F_t] = 0$ and $E[\epsilon_t] = 0$, $\text{cov}(\epsilon_t, \epsilon_{t-k}) = \sigma_\epsilon^2 \delta_k$.

The forecast $F_t$ is modeled as an AR(1) process:
$$
F_t = \nu F_{t-1} + \eta_t, \quad |\nu| < 1,
$$
where $\nu = \text{corr}(F_t, F_{t-1})$ represents the implicit decay in the forecast. 
The autocovariance function of $F_t$ is given by:
$$
\rho(k) = \text{cov}(F_t, F_{t-k}) = \rho(0) \nu^{|k|}.
$$

We define a generating function $\phi(\lambda)$ for the autocovariance:
$$
\phi(\lambda) \coloneqq \sum_{k=0}^\infty \lambda^k \rho(k) = \frac{\rho(0)}{1 - \lambda \nu}.
$$

## Optimization Problem Formulation

The objective is to maximize the expected utility over a long trading horizon of $T$ steps, balancing the expected Ideal P&L (IPNL) against trading costs (Self-Self Impact, or SSI) and holding costs (Self-Self Risk, or SSR):
$$
U = E[\text{IPNL}] - E[\text{SSI}] - E[\text{SSR}].
$$

Instead of solving the full time-varying dynamic programming problem to maximize this utility step-by-step, we optimize the steady-state performance over a long horizon. This yields a time-invariant trading strategy with constant coefficients, which is highly practical to execute. 

*Remark on Optimality:* One might wonder if a more complex strategy could achieve a higher utility. However, a fundamental result in stochastic control theory for Linear-Quadratic (LQ) problems (linear AR(1) dynamics, quadratic impact and risk costs) establishes that the optimal unconstrained policy is always a linear function of the state variables. Furthermore, as the trading horizon $T \to \infty$, the optimal dynamic policy asymptotically converges to a steady-state class with constant coefficients. Therefore, in the stationary limit, no better solution exists outside this family.

## Optimal Strategy Class

Guided by the steady-state approach, we search for an optimal position process $P_t$ in the class of processes defined by a scale parameter $A$ and an exponential decay parameter $\lambda$:
$$
\begin{aligned}
    P_t &= A G_t, \\
    G_t &= F_t + \lambda G_{t-1} = \sum_{k=0}^\infty \lambda^k F_{t-k}.
\end{aligned}
$$
where $G_t$ is the exponentially weighted sum of past forecasts.

Equivalently, 
$$
P_t = A F_t + \lambda P_{t-1}
$$

Note: $P_t$ is an $AR(2)$ process: 
$$
P_t =  A \eta_t + (\lambda + \nu) P_{t-1} - \lambda \nu P_{t-2}.
$$

## Derived Quantities and Costs

We compute the unconditional (stationary) expectations of the utility components over a long trading horizon of $T$ steps. By stationarity, the expected sum over $T$ steps is simply $T$ times the unconditional expectation of a single step.

### Ideal P&L (IPNL)
The expected ideal P&L is:
$$
E[\text{IPNL}] = E\left[ \sum_{t=1}^T P_t R_t \right] = E\left[ \sum_{t=1}^T P_t F_t \right] = A T \sum_{k=0}^\infty \lambda^k E[F_t F_{t-k}] = A T \phi(\lambda) = A T \frac{\rho(0)}{1 - \lambda \nu}.
$$

### Self-Self Impact Cost (SSI)
SSI refers to the quadratic trade cost. The trades are $P_t - P_{t-1} = A(G_t - G_{t-1}) = A(F_t + (\lambda-1)G_{t-1})$.
$$
E[\text{SSI}] = \mu \frac{\sigma}{\text{ADV}} E\left[ \sum_{t=1}^T (P_t - P_{t-1})^2 \right]
$$
For the AR(1) case, this simplifies to:
$$
E[\text{SSI}] = A^2 T \frac{2 \mu \sigma}{\text{ADV}} \rho(0) \frac{1-\nu}{(1+\lambda)(1-\lambda\nu)}.
$$

### Self-Self Risk Cost (SSR)
SSR is the quadratic holding cost (variance cost):
$$
E[\text{SSR}] = \gamma \sigma^2 E\left[ \sum_{t=1}^T P_t^2 \right] = A^2 T \gamma \sigma^2 E[G_t^2].
$$
For AR(1) forecasts:
$$
E[\text{SSR}] = A^2 T \gamma \sigma^2 \frac{\rho(0)}{1-\lambda^2} \frac{1+\lambda \nu}{1 - \lambda \nu}.
$$

## Maximizing Utility

Having derived the components of the expected utility $U(A, \lambda) = E[\text{IPNL}] - E[\text{SSI}] - E[\text{SSR}]$, suppressing the overall factor $T \rho(0)$ and substituting $b(\lambda) = \frac{1}{1-\lambda\nu}$, we can write the scaled utility as:
$$
\frac{U(A, \lambda)}{T \rho(0)} = A b(\lambda) - A^2 a(\lambda),
$$
where
$$
a(\lambda) = \frac{1}{1 - \lambda \nu} \left( \frac{2 \mu \sigma}{\text{ADV}} \frac{1-\nu}{1+\lambda} + \gamma \sigma^2 \frac{1+\lambda\nu}{1-\lambda^2} \right).
$$

### Solving for $A$
For a fixed $\lambda$, the optimal $A_*$ is:
$$
A_*(\lambda) = \frac{b(\lambda)}{2 a(\lambda)}.
$$

### Solving for $\lambda$
Substituting $A_*$ back into the utility:
$$
U^*(\lambda) = \frac{b^2(\lambda)}{4 a(\lambda)}.
$$
Maximizing $U^*(\lambda)$ is equivalent to minimizing $V(\lambda) \coloneqq \frac{4 a(\lambda)}{b^2(\lambda)}$:
$$
V(\lambda) = 4(1 - \lambda \nu) \left[ \frac{2 \mu \sigma}{\text{ADV}} \frac{1-\nu}{1+\lambda} + \gamma \sigma^2 \frac{1+\lambda\nu}{1-\lambda^2} \right].
$$
The optimal $\lambda_*$ is found by setting the derivative to zero, which recaps to the quadratic solution shown in the notes:
$$
\lambda_* = \text{argmin}_\lambda V(\lambda).
$$

# Plasticity

Here we describe the update procedure and tangent operator for non-associative plasticity.

!include application_development/identities.md

!include application_development/yield_function.md

## Flow rules

For plastic deformation, the yield function follows the following rule:

\begin{equation}
  \mathcal{F} \leqslant 0
\end{equation}

The flow rule for the plastic strain increment is written as follow:

\begin{equation}
  \Delta \varepsilon_{ij}^{p} = \Delta \lambda\frac{\partial \mathcal{G}}{\partial \sigma_{ij}}
\end{equation}

This expression can be used to express the invariants of the plastic strain increment:

\begin{equation}
  \begin{aligned}
    &\Delta \varepsilon_{v}^{p} = - \beta \Delta \lambda \\
    &\Delta e_{d}^{p} = \Delta \lambda
  \end{aligned}
\end{equation}

where $\Delta \lambda$ is the plastic increment. It follows the following rule:

\begin{equation}
  \mathcal{F}\Delta \lambda = 0
\end{equation}

which can be expressed as:
\begin{equation}
  \begin{aligned}
    \Delta \lambda &= 0 \quad \text{if} \quad \mathcal{F} < 0 \\
    \Delta \lambda &> 0 \quad \text{if} \quad \mathcal{F} = 0
  \end{aligned}
\end{equation}

Using Eq. 3, the plastic increment can be written:

\begin{equation}
  \Delta \lambda = \frac{<\mathcal{F}^{tr}>}{3G + \alpha \beta K}
\end{equation}

## Tangent operator modulus

Here we extend the tangent operator modulus given in [cite:dunne2005] to account for volumetric deformation and non-associative models.
The tangent operator modulus is defined as:

\begin{equation}
  \mathbb{D}_{ijkl} = \frac{\delta \sigma_{ij}}{\delta \varepsilon_{kl}}
\end{equation}

!bibtex bibliography
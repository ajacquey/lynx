# Viscoplasticity

Here we describe the update procedure and tangent operator for non-associative viscoplasticity.

!include application_development/identities.md

!include application_development/yield_function.md

## Flow rules

The flow rules for the viscoplastic strain increment is defined as:

\begin{equation}
  \Delta \varepsilon_{ij}^{vp} = \frac{<\mathcal{F}> \Delta t}{\eta}\frac{\partial \mathcal{G}}{\partial \sigma_{ij}}
\end{equation}

This expression can be rearranged by considering the magnitude of the viscoplastic strain rate $\Delta\lambda$:

\begin{equation}
  \begin{aligned}
    &\Delta \varepsilon_{v}^{vp} = - \beta \Delta \lambda \\
    &\Delta e_{d}^{vp} = \Delta \lambda \\
    &\Delta \lambda = \frac{<\mathcal{F}> \Delta t}{\eta}
  \end{aligned}
\end{equation}

Using [eq:identity_stress_strain], the increment of viscoplastic strain can be written:

\begin{equation}
  \Delta \lambda = \frac{1}{\left(1 + \frac{\left(3G + \alpha\beta K\right)\Delta t}{\eta}\right)}\frac{<\mathcal{F}^{tr}>\Delta t}{\eta}
\end{equation}

Furthermore, using the expression of the plastic increment $\Delta \lambda^{p}$ described in [Plasticity](application_development/plasticity.md), the previous expression can be expressed as:

\begin{equation}
  \Delta \lambda = \frac{\frac{\left(3G + \alpha\beta K\right)\Delta t}{\eta}}{\left(1 + \frac{\left(3G + \alpha\beta K\right)\Delta t}{\eta}\right)} \Delta \lambda^{p}
\end{equation}

## Tangent operator modulus


!bibtex bibliography
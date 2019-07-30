# Damage Rheology Development

Here we describe the base for the damage rheology implementation in LYNX.

## Stress update

We rely of the elastic energy and stress formulation of [!cite](lyakhovsky1997):

\begin{equation}
  \sigma_{ij} = \left(\lambda - \frac{\alpha \gamma}{\xi}\right) \varepsilon_{kk} \delta_{ij} + \left(2G - \alpha \gamma \left(\xi - 2\xi_{0}\right)\right) \varepsilon_{ij}.
\end{equation}

LYNX uses an incremental formulation for the stress update. Based on a Taylor expansion, the stress at a new time step $n+1$ is written:

\begin{equation}
  \sigma_{ij}^{n+1} = \sigma_{ij}^{n} + \frac{\partial \sigma_{ij}}{\partial \varepsilon_{kl}}:\left(\varepsilon_{kl}^{n+1} - \varepsilon_{kl}^{n}\right) + \frac{\partial \sigma_{ij}}{\partial \alpha} \left(\alpha^{n+1} - \alpha^{n}\right).
\end{equation}

The previous expression can be rearranged as:

\begin{equation}
  \sigma_{ij}^{n+1} = \sigma_{ij}^{n} + \mathbb{C}^{0}_{ijkl} : \left(\varepsilon_{kl}^{n+1} - \varepsilon_{kl}^{n}\right) - \mathbb{C}^{\xi}_{ijkl} : \left(\varepsilon_{kl}^{n+1} - \varepsilon_{kl}^{n}\right) - \sigma_{ij}^{\alpha} \left(\alpha^{n+1} - \alpha^{n}\right),
\end{equation}

where:

\begin{equation}
\label{eq:stiffness}
  \mathbb{C}^{0}_{ijkl} = \lambda \delta_{ij}\delta_{kl} + \left[G - \alpha^{n} \gamma\left(\frac{\xi^{n}}{2} - \xi_{0}\right)\right]\left(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}\right),
\end{equation}

\begin{equation}
  \mathbb{C}^{\xi}_{ijkl} = \alpha^{n} \gamma \left(\frac{\varepsilon_{ij}^{n}}{\lVert\varepsilon_{ij}^{n}\rVert}\delta_{kl} + \frac{\varepsilon_{kl}^{n}}{\lVert\varepsilon_{ij}^{n}\rVert}\delta_{ij} - \xi^{n} \frac{\varepsilon_{ij}^{n} \varepsilon_{kl}^{n}}{{\lVert \varepsilon_{ij}^{n}\rVert}^{2}}\right),
\end{equation}

and 

\begin{equation}
  \sigma_{ij}^{\alpha} = \gamma\left[\left(\xi^{n} - 2\xi_{0}\right)\varepsilon_{ij}^{n} + \lVert\varepsilon_{ij}^{n}\rVert \delta_{ij} \right]
\end{equation}

The parameter $\xi_{0}$ represents the modified internal friction of the material. With some algebra, one can express this coefficient based on the friction angle and the elastic moduli (intact) of the material as:

\begin{equation}
  \xi_{0} = -\sqrt{\frac{3}{1 + \frac{3}{2}\frac{K^{2}}{G^{2}}\sin^{2}\left(\varphi\right)}}
\end{equation}

## Inelastic model

In LYNX, we rely on a different formulation as the one presented in [!cite](lyakhovsky1997,lyakhovsky2015) for the inelastic update. The yield function presented in [!cite](lyakhovsky2015) reads:

\begin{equation}
  f = D\varepsilon_{v}^{3} + {\lVert \varepsilon_{ij} \rVert}^{2} \left(\xi - \xi_{0}\right)
\end{equation}

One can express this yield in the stress space in several ways. Here, we chose to use the following formulation because it includes as less as possible undefined regions for the yield. This formulation reads:

\begin{equation}
\label{eq:critical_strain_ratio}
 f = \xi - \xi_{cr}, \quad \text{with} \quad \xi_{cr} = \xi_{0} - \left(\sqrt{3} + \xi_{0}\right)\frac{\varepsilon_{v}}{\varepsilon_{vc}} \frac{\xi^{2}}{3}, \quad \text{and} \quad \varepsilon_{vc} = \frac{\left(\sqrt{3} + \xi_{0}\right)}{3D}.
\end{equation}

where the values of $\varepsilon_{v}$ and $\xi$ are calculated based on the previous timestep. This expression can therefore be used to express a yield in the form:

\begin{equation}
  f = e_{d} - \frac{\sqrt{2}}{3}\sqrt{\frac{3 - \xi_{cr}^{2}}{\xi_{cr}^{2}}} \varepsilon_{v},
\end{equation}

which can be expressed in stress space using non-relevant elastic moduli. Here we choose to use the ones forming the stiffness defined in [eq:stiffness]:

\begin{equation}
\label{eq:damage_yield}
  \mathcal{F} = \sigma_{e} - \frac{\sqrt{2}G_{e}}{K_{e}} \sqrt{\frac{3 - \xi_{cr}^{2}}{\xi_{cr}^{2}}} p.
\end{equation}

The inelastic strain update is conducted following a non-associative viscoplastic model using [eq:damage_yield] as a yield function. 
The update therefore follows the procedure describe in [Viscoplasticity](application_development/viscoplasticity.md). The general update procedure is therefore conducted in the following order:

1. Trial State: $p^{tr} = \tilde{p}^{tr} + K_{e} \Delta \varepsilon_{v}^{tr}$ and $\sigma_{e}^{tr} = \tilde{\sigma}_{e}^{tr} + G_{e} \Delta e_{d}^{tr}$

2. Inelastic strain update:

3. Inelastic correction: 

3. Damage and strain correction: $\sigma_{ij} = \bar{\sigma}_{ij} - \mathbb{C}^{\xi}_{ijkl} : \Delta\varepsilon_{kl} - \sigma_{ij}^{\alpha} \Delta \alpha$

!alert note prefix=False
The friction coefficient in [eq:damage_yield] is not always defined due to the square roots in its expression and in [eq:critical_strain_ratio]. 

## Damage evolution

After updating the stress following the previously described procedure, we compute the rate of damage accumulation as:

\begin{equation}
  \frac{\partial \alpha}{\partial t} = \frac{\mathcal{F}}{\eta_{\alpha}}
\end{equation}

The expression of the damage accumulation shows that the damage accumulation is controlled by the amount of overstress after the inelastic update. If the deformation is mainly brittle, the accumulation of inelastic strain is quite small and therefore leading to a significant damage accumulation whereas if the deformation is mainly ductile, the inelastic strain will be significant and therefore reduce the amount of accumulated damage.

!bibtex bibliography
# AMP_StateEvolution

Run plost_SE.m file.

we assume the signal is Bernouli-Gaussian $x \sim (1-P) \delta_0 + P \mathcal N(s; \mu_g, \nu_g)$, 

and the denoiser is thre related MMSE denoiser.

We implement the state evolution by numberical integration.

Specially, we characterize the state evolution as

$$
\begin{aligned}
\tau^{(t)}_r &= \frac{1}{\delta} \mathcal E^{(t)} + \tau_w \\
\mathcal E^{(t+1)} &= \mathbb E \left [ \eta^{(t)} (X + \mathcal N (0, \tau^{(t)}_r)) - X \right]^2
\end{aligned}
$$
where the $\mathcal E^{(t)}$ is the denoiser output error, and $\tau^{(t)}_r$ is the denoiser input variance.

We note that as for the calculation of $\mathcal E^{(t+1)} &= \mathbb E \left [ \eta^{(t)} (X + \mathcal N (0, \tau^{(t)}_r)) - X \right]^2$, the expectation is taken over both variables $X$ and $Z$, and the aligned distribution is $p(X, Z)$.

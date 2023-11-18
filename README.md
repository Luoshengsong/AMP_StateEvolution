# AMP_StateEvolution

Run plost_SE.m file.

we assume the signal is Bernouli-Gaussian, and the denoiser is thre related MMSE denoiser.

We implement the state evolution by numberical integration.

Specially, we characterize the state evolution as

$$
\begin{aligned}
\tau^{(t)}_r &= \frac{1}{\delta} \mathcal E^{(t)} + \tau_w \\
\mathcal E^{(t+1)} &= \mathbb E \left [ \eta^{(t)} (X + \mathcal N (0, \tau^{(t)}_r)) - X \right]^2
\end{aligned}
$$

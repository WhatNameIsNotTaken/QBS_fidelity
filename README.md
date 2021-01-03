# Introduction
<big>This repository is used for describing the features of a certain quantum blind signature scheme under different noise, including *phase flipping noise*, *depolarizing noise*, *amplitude damping noise* and *bit flipping noise*.</big>

[TOC]

# 1. The outcome, when $\varphi_{in}=|00\rangle$
## 1.1 Phase flipping noise
The outcome of QBS_Bell_GHZ_phaseflipping.m was as follows:

$$\rho_{out}^{pf}=
\left[\begin{array}{}
4p^3 - 6p^2 + 3p & 0 & 0 & 0 \\
0 & - 4p^3 + 6p^2 - 3p + 1 & 0 & 0  \\ 
0 & 0 & 0 & 0  \\  
0 & 0 & 0 & 0 
\end{array}\right]
$$

$$f_{out}^{pf}=4p^3 - 6p^2 + 3p$$

## 1.2 Amplitude damping noise
 The outcome of QBS_Bell_GHZ_ampdamping.m was as follows:
 
 $$\rho_{out}^{ad}=
\left[\begin{array}{}
\frac 1 2(p^2-p+(1-p)^{ 3/2}) + \frac{\sqrt2} {16} + \sqrt2/16) & p/2 & 0 \\
p/2 & \frac 1 2(p^2-p+(1-p)^{ 3/2}) - \frac{\sqrt2} {16} & 0 & 0  \\ 
0 & 0 & -\frac 1 2(p^2-p) & 0  \\  
0 & 0 & 0 &-\frac 1 2(p^2-p)
\end{array}\right]
$$

$$f_{out}^{ad}=\frac 1 2(p^2-p+(1-p)^{ 3/2}) + \frac{\sqrt2} {16}$$

## 1.3 Bit flipping noise
 The outcome of QBS_Bell_GHZ_bitflipping.m was as follows:

$$\rho_{out}^{bf}=
\left[\begin{array}{}
p^4 + 2p(p - 1)^2 - p^2(p - 1) - p(p - 1)^3 & 0 & 0 & 0 \\
0 & 0 & 0 & 0  \\ 
0 & 0 & (p - 1)^4 + p*(p - 1)^2 - 2p^2(p - 1) - p^3(p - 1) & 0  \\  
0 & 0 & 0 & 0 
\end{array}\right]
$$

$$f_{out}^{bf}=4p^3 - 6p^2 + 3p$$

## 1.4 Depolarizing noise
The outcome of QBS_Bell_GHZ_depolarizing.m was as follows:
 
$$\rho_{out}^{dp}=
\left[\begin{array}{}
(64p^4)/81 - (32p^3)/9 + (16p^2)/3 - (10p)/3 + 1 & 0 & 0 & 0 \\
0 & - (64p^4)/81 + (64p^3)/27 - (8p^2)/3 + (4p)/3 & 0 & 0  \\ 
0 & 0 & - (64p^4)/81 + (64p^3)/27 - (8*p^2)/3 + (4p)/3 & 0  \\  
0 & 0 & 0 &(64p^4)/81 - (32p^3)/27 + (2p)/3
\end{array}\right]
$$

$$f_{out}^{dp}=(64p^4)/81 - (32p^3)/9 + (16p^2)/3 - (10p)/3 + 1 $$

## 1.5 The figure of the parameter $p$ and fidelity


# Introduction
<big>This repository is used for describing the features of a certain quantum blind signature scheme under different noise, including *phase flipping noise*, *depolarizing noise*, *amplitude damping noise* and *bit flipping noise*.</big>

[TOC]

# 1. The outcome, when $\varphi_{in}=|00\rangle$
## 1.1 Phase flipping noise
The outcome of QBS_Bell_GHZ_phaseflipping.m was as follows:

$$\rho_{out_1}^{pf}=
\left[\begin{array}{}
-4p^3 + 6p^2 - 3p+1 & 0 & 0 & 0 \\
0 &  4p^3 - 6p^2 + 3p  & 0 & 0  \\ 
0 & 0 & 0 & 0  \\  
0 & 0 & 0 & 0 
\end{array}\right]
$$

$$f_{out_1}^{pf}=-4p^3 + 6p^2 - 3p+1$$

## 1.2 Amplitude damping noise
 The outcome of QBS_Bell_GHZ_ampdamping.m was as follows:
 
 $$\rho_{out_1}^{ad}=
\left[\begin{array}{}
(p^2-p+(1-p)^{ 3/2}+1)/2  & p/2 & 0 & 0\\
p/2 & (p^2-p-(1-p)^{ 3/2}+1)/2 & 0 & 0  \\ 
0 & 0 & -(p^2-p)/2 & 0  \\  
0 & 0 & 0 &-(p^2-p)/2
\end{array}\right]
$$

$$f_{out_1}^{ad}=\frac 1 2(p^2-p+(1-p)^{ 3/2}+1) $$

## 1.3 Bit flipping noise
 The outcome of QBS_Bell_GHZ_bitflipping.m was as follows:

$$\rho_{out_1}^{bf}=
\left[\begin{array}{}
-4p^3 + 6p^2 - 3p + 1 & 0 & 0 & 0 \\
0 & 0 & 0 & 0  \\ 
0 & 0 &  4p^3 - 6p^2 + 3p & 0  \\  
0 & 0 & 0 & 0 
\end{array}\right]
$$

$$f_{out_1}^{bf}= - 4p^3 + 6p^2 - 3p + 1$$

## 1.4 Depolarizing noise
The outcome of QBS_Bell_GHZ_depolarizing.m was as follows:
 
$$\rho_{out_1}^{dp}=
\left[\begin{array}{}
64p^4/81 - 32p^3/9 + 16p^2/3 - 10p/3 + 1 & 0 & 0 & 0 \\
0 & - 64p^4/81 + 64p^3/27 - 8p^2/3 + 4p/3 & 0 & 0  \\ 
0 & 0 & - 64p^4/81 + 64p^3/27 - 8p^2/3 + 4p/3 & 0  \\  
0 & 0 & 0 &64p^4/81 - 32p^3/27 + 2p/3
\end{array}\right]
$$

$$f_{out_1}^{dp}=64p^4/81 - 32p^3/9 + 16p^2/3 - 10p/3 + 1 $$

## 1.5 The figure of the parameter $p$ and fidelity

![](https://cdn.mathpix.com/snip/images/09ltHE8Jvq2T3rXdxbeONpqLwhjk459qMtCRwMYT6sM.original.fullsize.png)

# 2. The outcome, when $\varphi_{in}=|01\rangle$
## 2.1 Phase flipping noise
The outcome of QBS_Bell_GHZ_phaseflipping.m was as follows:

$$\rho_{out_2}^{pf}=
\left[\begin{array}{}
4p^3 - 6p^2 + 3p & 0 & 0 & 0 \\
0 & - 4p^3 + 6p^2 - 3p + 1 & 0 & 0  \\ 
0 & 0 & 0 & 0  \\  
0 & 0 & 0 & 0 
\end{array}\right]
$$

$$f_{out_2}^{pf}=-4p^3 + 6p^2 - 3p + 1$$

## 2.2 Amplitude damping noise
 The outcome of QBS_Bell_GHZ_ampdamping.m was as follows:
 
 $$\rho_{out_2}^{ad}=
\left[\begin{array}{}
(p^2-p-(1-p)^{ 3/2}-1)/2  & p/2 & 0 \\
p/2 & (p^2-p+(1-p)^{ 3/2}+1)/2 & 0 & 0  \\ 
0 & 0 & -(p^2-p)/2 & 0  \\  
0 & 0 & 0 &-(p^2-p)/2
\end{array}\right]
$$

$$f_{out_2}^{ad}=\frac 1 2(p^2-p+(1-p)^{ 3/2}+1) $$

## 2.3 Bit flipping noise
 The outcome of QBS_Bell_GHZ_bitflipping.m was as follows:

$$\rho_{out_2}^{bf}=
\left[\begin{array}{}
0 & 0 & 0 & 0  \\ 
0 & - 4p^3 + 6p^2 - 3p + 1 &  0 & 0 \\
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 4p^3 - 6p^2 + 3p 
\end{array}\right]
$$

$$f_{out_2}^{bf}= - 4p^3 + 6p^2 - 3p + 1$$

## 2.4 Depolarizing noise
The outcome of QBS_Bell_GHZ_depolarizing.m was as follows:
 
$$\rho_{out_2}^{dp}=
\left[\begin{array}{}
-64p^4/81 + 64p^3/27 - 8p^2/3 + 4p/3 & 0 & 0 & 0 \\
0 & 64p^4/81 - 32p^3/9 + 16p^2/3 - 10p/3 + 1 & 0 & 0  \\ 
0 & 0 &  64p^4/81 - 32p^3/27 + 2p/3 & 0  \\  
0 & 0 & 0 & -64p^4/81 + 64p^3/27 - 8p^2/3 + 4p/3
\end{array}\right]
$$

$$f_{out_2}^{dp}=64p^4/81 - 32p^3/9 + 16p^2/3 - 10p/3 + 1 $$

## 2.5 The figure of the parameter $p$ and fidelity

![](https://cdn.mathpix.com/snip/images/09ltHE8Jvq2T3rXdxbeONpqLwhjk459qMtCRwMYT6sM.original.fullsize.png)

# 3. The outcome, when $\varphi_{in}=|10\rangle$
## 3.1 Phase flipping noise
The outcome of QBS_Bell_GHZ_phaseflipping.m was as follows:

 $$\rho_{out_{3}}^{pf}=
\left[\begin{array}{}
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0  \\ 
0 & 0 & -4p^3+6p^2-3p+1  & 0  \\  
0 & 0 & 0 &  4p^3-6p^2+3p
\end{array}\right]
$$

$$f_{out_3}^{pf}=-4p^3+6p^2-3p+1$$

## 3.2 Amplitude damping noise
 The outcome of QBS_Bell_GHZ_ampdamping.m was as follows:
 
$$\rho_{out_3}^{ad}=
\left[\begin{array}{}
p^3 - 3p^2/2 + ((1 - p)^{ 3/2}/2 + 1)p & p^2/2 & 0 & 0 \\
p^2/2 & p^3 - 3p^2/2 + (1-((1 - p)^{3/2})/2)p & 0 & 0  \\ 
0 & 0 & 3p^2/2 - p^3 - p + ((1 - p)^{5/2}) + 1)/2 & p^2/2 - p/2  \\  
0 & 0 & p^2/2 - p/2 & 3p^2/2 - p^3 - p - ((1 - p)^{5/2} - 1)/2
\end{array}\right]
$$

$$f_{out_3}^{ad}=3p^2/2 - p^3 - p +((1 - p)^{5/2} + 1)/2$$

## 3.3 Bit flipping noise
 The outcome of QBS_Bell_GHZ_bitflipping.m was as follows:

$$\rho_{out_3}^{bf}=
\left[\begin{array}{}
4p^3 - 6p^2 + 3p& 0 & 0 & 0 \\
0 & 0 & 0 & 0  \\ 
0 & 0 & - 4p^3 + 6p^2 - 3p + 1 & 0  \\  
0 & 0 & 0 & 0 
\end{array}\right]
$$

$$f_{out_3}^{bf}=-4p^3 + 6p^2 - 3p+1$$

## 3.4 Depolarizing noise
The outcome of QBS_Bell_GHZ_depolarizing.m was as follows:
 
$$\rho_{out_3}^{dp}=
\left[\begin{array}{}
-64p^4/81 + 64p^3/27 - 8p^2/3 + 4p/3  & 0 & 0 & 0 \\
0 &64p^4/81-32p^3/27 + 2p/3 & 0 & 0  \\ 
0 & 0 & 64p^4/81 - 32p^3/9 + 16p^2/3 - 10p/3+1 & 0  \\  
0 & 0 & 0 &- 64p^4/81 + 64p^3/27 - 8p^2/3 + 4p/3
\end{array}\right]
$$

$$f_{out_3}^{dp}=64p^4/81 - 32p^3/9 + 16p^2/3 - 10p/3 + 1 $$

## 3.5 The figure of the parameter $p$ and fidelity
![](https://cdn.mathpix.com/snip/images/rFTRgXuxlY6wxrSWY82ioYikGqGPRiLscdmLpp9BV-o.original.fullsize.png)

# 4. The outcome, when $\varphi_{in}=|11\rangle$
## 4.1 Phase flipping noise

$$\rho_{out_{4}}^{pf}=
\left[\begin{array}{}
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0  \\ 
0 & 0 & 4p^3 - 6p^2 + 3p  & 0  \\  
0 & 0 & 0 &  -4p^3 + 6p^2 - 3p+1 
\end{array}\right]
$$

$$f_{out_4}^{pf}=-4p^3+6p^2-3p+1$$

## 4.2 Amplitude damping noise
 The outcome of QBS_Bell_GHZ_ampdamping.m was as follows:
 
 $$\rho_{out_4}^{ad}=
\left[\begin{array}{}
p^3 - 3p^2/2 + p - p(1 - p)^{ 3/2}/2  & p^2/2 & 0 & 0 \\
p^2/2 & p^3 - 3p^2/2 + p + p(1 - p)^{ 3/2}/2 & 0 & 0  \\ 
0 & 0 & 3p^2/2 - p^3 - p - ((1 - p)^{5/2} - 1)/2 & p^2/2 - p/2  \\  
0 & 0 & p^2/2 - p/2 & 3p^2/2 - p^3 - p + ((1 - p)^(5/2) + 1)/2
\end{array}\right]
$$

$$f_{out_4}^{ad}=3p^2/2 - p^3 - p + (\sqrt2(1 - p)^{5/2} + 1)/2$$

## 4.3 Bit flipping noise
 The outcome of QBS_Bell_GHZ_bitflipping.m was as follows:

$$\rho_{out_{4}}^{bf}=
\left[\begin{array}{}
0 & 0 & 0 & 0 \\
0 & 4p^3 - 6p^2 + 3p & 0 & 0  \\ 
0 & 0 & 0  & 0  \\  
0 & 0 & 0 &  - 4p^3 + 6p^2 - 3p + 1
\end{array}\right]
$$

$$f_{out_4}^{bf}=-4p^3+6p^2-3p+1$$


## 4.4 Depolarizing noise
The outcome of QBS_Bell_GHZ_depolarizing.m was as follows:
 
$$\rho_{out_4}^{dp}=
\left[\begin{array}{}
64p^4/81 - 32p^3/27 + 2p/3 & 0 & 0 & 0 \\
0 & - 64p^4/81 + 64p^3/27 - 8p^2/3 + 4p/3 & 0 & 0  \\ 
0 & 0 & - 64p^4/81 + 64p^3/27 - 8p^2/3 + 4p/3 & 0  \\  
0 & 0 & 0 &64p^4/81 - 32p^3/9 + 16p^2/3-10p/3+1
\end{array}\right]
$$

$$f_{out_4}^{dp}=64p^4/81 - 32p^3/9 + 16p^2/3 - 10p/3 + 1 $$

## 4.5 The figure of the parameter $p$ and fidelity
![](https://cdn.mathpix.com/snip/images/9luhe7XI1ywvCR82qr80jv_9hHbqFZWFsNN-oAoxReQ.original.fullsize.png)
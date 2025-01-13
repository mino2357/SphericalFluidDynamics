# SphericalFluidDynamics

実験的な流体の方程式。計算メモ。

## 計算メモ

緯度、経度を $(\alpha, \beta)$ 座標とする。

速度ベクトル。緯度方向、経度方向。

$$
u = \begin{pmatrix} u^\alpha \\ u^\beta \end{pmatrix}
$$

計量テンソル。

$$
\begin{align}
g =
\begin{pmatrix}
R^2 \cos^2\beta & 0 \\
0 & R^2
\end{pmatrix}
\end{align}
$$

クリストッフェル記号。

$$
\begin{align}
    \Gamma^\alpha_{\alpha\beta} &= \Gamma^\alpha_{\beta\alpha} = -\tan\beta, \\
    \Gamma^\beta_{\alpha\alpha} &= -\cos\beta \sin\beta.
\end{align}
$$


## 偏微分方程式

実装した式。実験的。ミスもあると思われる。

質量保存則。

$$
\begin{align}
    \frac{\partial \rho}{\partial t} 
    + \frac{1}{R} \left( \frac{\partial (\rho u^\alpha)}{\partial \alpha} 
    + \frac{\partial (\rho u^\beta)}{\partial \beta} \right) = 0
\end{align}
$$

経度方向。

$$
\begin{align}
    \frac{\partial u^\alpha}{\partial t} 
    + u^\alpha \frac{\partial u^\alpha}{\partial \alpha} 
    + u^\beta \frac{\partial u^\alpha}{\partial \beta} 
    - 2 \tan\beta \cdot u^\alpha u^\beta 
    = -\frac{1}{\rho R^2 \cos^2 \beta} \frac{\partial p}{\partial \alpha} + \frac{\nu}{R^2} \left( \frac{1}{\cos^2\beta} \frac{\partial^2 u^\alpha}{\partial \alpha^2} 
    - \tan\beta \frac{\partial u^\alpha}{\partial \beta} 
    + \frac{\partial^2 u^\alpha}{\partial \beta^2} \right)
\end{align}
$$

緯度方向。

$$
\begin{align}
    \frac{\partial u^\beta}{\partial t} 
    + u^\alpha \frac{\partial u^\beta}{\partial \alpha} 
    + u^\beta \frac{\partial u^\beta}{\partial \beta} 
    - \cos\beta \sin\beta \cdot u^\alpha u^\alpha 
    = -\frac{1}{\rho R^2} \frac{\partial p}{\partial \beta} + \frac{\nu}{R^2} \left( \frac{1}{\cos^2\beta} \frac{\partial^2 u^\beta}{\partial \alpha^2} 
    - \tan\beta \frac{\partial u^\beta}{\partial \beta} 
    + \frac{\partial^2 u^\beta}{\partial \beta^2} \right)
\end{align}
$$

状態方程式。

$$
\begin{align}
    p(\rho) = k \rho
\end{align}
$$

ただし極付近は取り除いた。
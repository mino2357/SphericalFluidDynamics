# SphericalFluidDynamics

実験的な球面上の流体シミュレータです。経度 $\alpha$、緯度 $\beta$ で定義された格子を用いて圧縮性流体の運動方程式を単純な陽解法で積分します。極付近は特異性を避けるため除外しています。

## モデル

速度ベクトル $u = (u^\alpha, u^\beta)$ と密度 $\rho$ を考えます。球半径を $R$ とすると計量テンソルは

$$
 g = \begin{pmatrix}
 R^2 \cos^2\beta & 0 \\
 0 & R^2
 \end{pmatrix}
$$

非ゼロクリストッフェル記号は

$$
\Gamma^\alpha_{\alpha\beta} = \Gamma^\alpha_{\beta\alpha} = -\tan\beta,\qquad
\Gamma^\beta_{\alpha\alpha} = -\cos\beta\sin\beta.
$$

これを用いた保存方程式は

$$
\frac{\partial \rho}{\partial t}
+ \frac{1}{R}\left(
 \frac{\partial (\rho u^\alpha)}{\partial \alpha}
 + \frac{\partial (\rho u^\beta)}{\partial \beta}
\right) = 0,
$$

$$
\frac{\partial u^\alpha}{\partial t}
+ u^\alpha\frac{\partial u^\alpha}{\partial \alpha}
+ u^\beta\frac{\partial u^\alpha}{\partial \beta}
- 2\tan\beta\, u^\alpha u^\beta
= -\frac{1}{\rho R^2\cos^2\beta}\frac{\partial p}{\partial \alpha}
  + \frac{\nu}{R^2}\left(
     \frac{1}{\cos^2\beta}\frac{\partial^2 u^\alpha}{\partial \alpha^2}
     - \tan\beta\frac{\partial u^\alpha}{\partial \beta}
     + \frac{\partial^2 u^\alpha}{\partial \beta^2}
  \right),
$$

$$
\frac{\partial u^\beta}{\partial t}
+ u^\alpha\frac{\partial u^\beta}{\partial \alpha}
+ u^\beta\frac{\partial u^\beta}{\partial \beta}
- \cos\beta\sin\beta\, u^\alpha u^\alpha
= -\frac{1}{\rho R^2}\frac{\partial p}{\partial \beta}
  + \frac{\nu}{R^2}\left(
     \frac{1}{\cos^2\beta}\frac{\partial^2 u^\beta}{\partial \alpha^2}
     - \tan\beta\frac{\partial u^\beta}{\partial \beta}
     + \frac{\partial^2 u^\beta}{\partial \beta^2}
  \right).
$$

状態方程式は $p(\rho) = k\rho$ を用います。

## ビルドと実行

```bash
make
./main > field.dat
```

標準エラー出力にはステップ数と最大速度が出力され、標準出力には $(\alpha,\beta)$ 座標と流速・密度が出力されます。


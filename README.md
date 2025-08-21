# SphericalFluidDynamics

実験的な球面上の流体シミュレータです。経度 \( \lambda \) と緯度 \( \varphi \) を用いた球面座標で圧縮性流体の運動方程式を単純な陽解法で積分します。極付近は特異性を避けるため除外しています。

## モデル

球半径を \(R\) とすると，線素は
```math
ds^2 = R^2\cos^2\varphi\,d\lambda^2 + R^2\,d\varphi^2
```
となり，この座標系における非ゼロのクリストッフェル記号は
```math
\Gamma^\lambda_{\lambda\varphi} = \Gamma^\lambda_{\varphi\lambda} = -\tan\varphi,\qquad
\Gamma^\varphi_{\lambda\lambda} = \sin\varphi\cos\varphi.
```

速度ベクトルを \( \boldsymbol{u} = u^\lambda\,\partial_\lambda + u^\varphi\,\partial_\varphi \) とし，密度を \( \rho \) とすると，保存方程式は以下の通りである。

連続の式：
```math
\partial_t \rho + \frac{1}{R\cos\varphi}\left[\partial_\lambda(\rho u^\lambda \cos\varphi) + \partial_\varphi(\rho u^\varphi)\right] = 0
```

経度方向の運動方程式：
```math
\partial_t u^\lambda + u^\lambda\partial_\lambda u^\lambda + u^\varphi\partial_\varphi u^\lambda - 2\tan\varphi\,u^\lambda u^\varphi = -\frac{1}{\rho R^2\cos^2\varphi}\partial_\lambda p + \frac{\nu}{R^2}\left[\frac{1}{\cos^2\varphi}\partial_\lambda^2 u^\lambda - \tan\varphi\,\partial_\varphi u^\lambda + \partial_\varphi^2 u^\lambda\right]
```

緯度方向の運動方程式：
```math
\partial_t u^\varphi + u^\lambda\partial_\lambda u^\varphi + u^\varphi\partial_\varphi u^\varphi + \cos\varphi\sin\varphi\,(u^\lambda)^2 = -\frac{1}{\rho R^2}\partial_\varphi p + \frac{\nu}{R^2}\left[\frac{1}{\cos^2\varphi}\partial_\lambda^2 u^\varphi - \tan\varphi\,\partial_\varphi u^\varphi + \partial_\varphi^2 u^\varphi\right]
```

状態方程式には \( p(\rho) = k\rho \) を用いる。

## ビルドと実行

```bash
make
./main > field.dat
```

標準エラーにはステップ数と最大速度が出力され，標準出力には \((\lambda,\varphi)\) と流速・密度が出力される。

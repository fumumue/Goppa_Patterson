# Goppa符号とPattersonアルゴリズム

### Goppa符号の概要

Goppa符号は有限体上の代数的な誤り訂正符号で、McEliece暗号やNiederreiter暗号などの符号ベース暗号(耐量子計算機暗号の候補)で使われています。

#### 定義(二元の場合)

- $m, t$ を正整数とし、を各成分 $\mathbb{F}_2$ 上の $m$ 次元ベクトルに展開すると、
  $mt \times n$の二元検査行列が得られます。
-  上の次数 $t$ の既約多項式 $g(x)$ を選びます(Goppa多項式)。
- $L = (\alpha_1, \dots, \alpha_n)$ を $\mathbb{F}_{2^m}$ の相異なる元の列とします($n \le 2^m$、すべての $i$ で $g(\alpha_i) \neq 0$)。

このとき、Goppa符号 $\Gamma(L, g)$ は次の集合です。

$$
\Gamma(L, g) = \lbrace c \in \mathbb{F}_2^n \mid \sum_{i=1}^{n} \frac{c_i}{x - \alpha_i} \equiv 0 \pmod{g(x)} \rbrace
$$

#### パラメータ

| 項目 | 値 |
|------|----|
| 符号長 | $n$ |
| 次元 | $k \ge n - mt$ |
| 最小距離 | $d \ge 2t + 1$ |
| 訂正能力 | 最大 $t$ 個の誤り |


#### 検査行列

$g(x) = \sum_{k=0}^{t} g_k x^k$ とします。$\mathbb{F}_{2^m}$ 上の $t \times n$ 行列

$$
H = T V D
$$

を考えます。ここで $T$、$V$、$D$ は次のとおりです。

$$
T = \begin{pmatrix}
g_1 & g_2 & \cdots & g_{t-1} & g_t \\
g_2 & g_3 & \cdots & g_t & 0 \\
\vdots & \vdots & & \vdots & \vdots \\
g_{t-1} & g_t & \cdots & 0 & 0 \\
g_t & 0 & \cdots & 0 & 0
\end{pmatrix}
$$

$$
V = \begin{pmatrix}
1 & \cdots & 1 \\
\alpha_1 & \cdots & \alpha_n \\
\vdots & & \vdots \\
\alpha_1^{t-1} & \cdots & \alpha_n^{t-1}
\end{pmatrix},
\quad
D = \begin{pmatrix}
g(\alpha_1)^{-1} & & \\
& \ddots & \\
& & g(\alpha_n)^{-1}
\end{pmatrix}
$$

$T$ は $g(x)$ の係数からなる正則行列です。$T$ を省いた $VD$ も同じ符号の検査行列になります。

$$
VD = \begin{pmatrix}
g(\alpha_1)^{-1} & \cdots & g(\alpha_n)^{-1} \\
\alpha_1 g(\alpha_1)^{-1} & \cdots & \alpha_n g(\alpha_n)^{-1} \\
\vdots & & \vdots \\
\alpha_1^{t-1} g(\alpha_1)^{-1} & \cdots & \alpha_n^{t-1} g(\alpha_n)^{-1}
\end{pmatrix}
$$

$H$ の各成分を $\mathbb{F}_2$ 上の $m$ 次元ベクトルに展開すると、
$mt \times n$
の二元検査行列が得られます。


## Patterson復号における検査行列とゴッパ多項式の係数

Pattersonアルゴリズムは、シンドロームを多項式として扱います。

$$
S(x) \equiv \sum_{i=1}^{n} \frac{c_i}{x - \alpha_i} \pmod{g(x)}
$$

$S(x)$ の逆元、平方根、拡張ユークリッド互除法は、すべて $\bmod\ g(x)$ の世界で計算します。そのため、シンドロームは $1, x, \dots, x^{t-1}$ を基底とする係数ベクトルの形になっている必要があります。

#### 係数が現れる場所

$g(x) = \sum_{k=0}^{t} g_k x^k$ とします。$\bmod\ g(x)$ では次が成り立ちます。

$$
\frac{1}{x - \alpha} \equiv \frac{g(x) - g(\alpha)}{x - \alpha} \cdot \frac{1}{g(\alpha)} \pmod{g(x)}
$$

右辺の分子を展開すると、

$$
\frac{g(x) - g(\alpha)}{x - \alpha} = \sum_{j=0}^{t-1} x^j \sum_{k=j+1}^{t} g_k \alpha^{k-1-j}
$$

となります。したがって、$S(x)$ の $x^j$ の係数は次のとおりです。

$$
s_j = \sum_{i=1}^{n} \frac{c_i}{g(\alpha_i)} \sum_{l=0}^{t-1-j} g_{j+l+1} \alpha_i^{l}
$$

これを行列でまとめると、

$$
H = T V D
$$

です。

- $T$ は $g(x)$ の係数からなる $t \times t$ 行列で、$T_{j,l} = g_{j+l+1}$($j+l+1 > t$ のとき 0)です。
- $V$ はヴァンデルモンド行列で、$V_{l,i} = \alpha_i^{l}$ です。
- $D = \mathrm{diag}\bigl(g(\alpha_1)^{-1}, \dots, g(\alpha_n)^{-1}\bigr)$ です。

#### 係数の役割

| 項目 | 役割 |
|------|------|
| $VD$ | $\sum_i c_i \alpha_i^{j} / g(\alpha_i)$ という、評価値にもとづく量を与える |
| $T$ | $VD$ で得た量を、$S(x)$ の係数に変換する基底変換 |
| $g_t$ | $T$ の反対角の成分。$g_t \neq 0$ なので $T$ は正則 |
| $g_0$ | $T$ には現れず、$D$ の $g(\alpha_i)^{-1}$ の中に入る |

$T$ は正則なので、$H = TVD$ と $VD$ は同じ符号を定めます。違いは、シンドロームを評価値の形で持つか、$S(x)$ の係数として持つかだけです。

##### 例: $t = 2$ の場合

$g(x) = x^2 + g_1 x + g_0$ とすると、

```math
T = \begin{pmatrix} g_1 & 1 \\ 1 & 0 \end{pmatrix}
```

であり、シンドロームの係数は次のようになります。

$$
s_0 = \sum_{i=1}^{n} \frac{c_i}{g(\alpha_i)} (g_1 + \alpha_i), \qquad
s_1 = \sum_{i=1}^{n} \frac{c_i}{g(\alpha_i)}
$$

##### 暗号での位置づけ

Niederreiter暗号の公開鍵は、この行列の各成分を $\mathbb{F}_2$ 上の $m$ 次元ベクトルに展開し、標準形に変形したものです。公開鍵だけからは $T$ や $g(x)$ の構造は見えません。秘密鍵の $g(x)$ を知っている復号者だけが、シンドロームを $S(x)$ の係数に戻して Patterson を適用できます。

---

## Patterson アルゴリズム

二元Goppa符号に対する効率的な復号アルゴリズムです(Patterson, 1975)。受信語 $r = c + e$($c \in \Gamma(L,g)$、
$\mathrm{wt}(e) \le t$)から誤り $e$ を復元します。

#### 準備

- シンドローム:
  $$S(x) \equiv \sum_{i=1}^{n} \frac{r_i}{x - \alpha_i} \equiv \sum_{i \in E} \frac{1}{x - \alpha_i} \pmod{g(x)}$$
  ($E$ は誤り位置の集合)
- 誤り位置多項式:
  $$\sigma(x) = \prod_{i \in E} (x - \alpha_i)$$
- 二元の場合、形式微分 $\sigma'(x)$
  は $\sigma$ を偶数部と奇数部に分けて
  $$\sigma(x) = a(x)^2 + x\, b(x)^2,\quad \sigma'(x) = b(x)^2$$
  と書けます($\deg a \le \lfloor t/2 \rfloor$、
  $\deg b \le \lfloor (t-1)/2 \rfloor$)。
- 鍵方程式:
  $$S(x)\,\sigma(x) \equiv \sigma'(x) \pmod{g(x)}$$

#### 手順

1. **シンドローム計算**: $S(x)$ を求める。$S(x) = 0$ なら誤りなし。
2. **逆元**: $T(x) = S(x)^{-1} \bmod g(x)$ を求める。
3. **特殊ケース**: $T(x) = x$ なら $\sigma(x) = x$ として手順6へ。
4. **平方根**: $R(x) = \sqrt{T(x) + x} \bmod g(x)$ を求める。
   - $T + x = T_0(x)^2 + x\,T_1(x)^2$と偶奇分解し、事前計算した $\sqrt{x} \bmod g$ を用いて $R = T_0 + \sqrt{x}\,T_1$ とする。
     （ここで、$T_0$は$R$の偶数部分、$T_1=xT_2$は$R$の奇数部分とし、$R=T_0+xT_2$）
5. **拡張ユークリッド互除法**: $a(x) \equiv b(x)\,R(x) \pmod{g(x)}$ を満たす $a, b$ を、次数条件
   $\deg a \le \lfloor t/2 \rfloor$、
   $\deg b \le \lfloor (t-1)/2 \rfloor$を満たした時点で打ち切って求め、

   $$\sigma(x) = a(x)^2 + x\,b(x)^2$$

   とする。
7. **根の探索**: すべての $\alpha_i \in L$ について $\sigma(\alpha_i) = 0$ を判定し、根となる位置 $i$ で $e_i = 1$ とする。
8. **訂正**: $c = r + e$ を出力する。

#### 導出の要点

鍵方程式に $\sigma = a^2 + xb^2$、
$\sigma' = b^2$
を代入し、両辺に $T = S^{-1}$ を掛けると

$$
a^2 \equiv (T + x)\, b^2 \pmod{g}
$$

となります。

$\mathbb{F}_{2^m}[x]/(g)$

は標数2の体なので平方根が一意に存在し、$a \equiv R\,b$($R^2 = T + x$)と変形できます。これにより鍵方程式が $a, b$ の線形合同式に帰着され、拡張ユークリッド互除法で解けます。

#### 計算量

シンドローム計算・逆元・平方根・根探索を合わせて、主要項は根探索の $O(nt)$ 程度の体演算です。

## ベンチマーク結果

5回実行(`./build/bench 5`)の結果です。全回で復号の検証に PASS しています。

| 工程 | 最小 (ms) | 中央値 (ms) | 平均 (ms) | 最大 (ms) |
|------|----------:|------------:|----------:|----------:|
| 鍵生成 (keygen) | 489.3 | 538.7 | 528.9 | 576.3 |
| 暗号化 (encrypt) | 4.6 | 7.1 | 9.2 | 14.5 |
| 復号 (decode) | 313.7 | 341.0 | 342.7 | 376.7 |
| 合計 (total) | 807.7 | 876.5 | 880.8 | 942.1 |

- パラメーター: `N=8192, K=256, E=13`
- 計測環境: CPU=Athlon220G、コンパイラ=gcc、スレッド数(`OMP_NUM_THREADS=4`)

### 20261030:新記録
|iter  4: keygen 0.200s  encrypt 0.000s  decode 0.015s  [PASS]
|------|--------:|------:|------:|-------:|
|stage        | min(ms) |median(ms)  | mean(ms)   | max(ms)|
|keygen         | 119.8     | 183.4     | 165.4     | 203.4|
|encrypt          | 0.1       | 0.1       | 0.1       | 0.2|
|decode          | 12.9      | 15.4      | 19.3      | 36.3|
|total          | 133.3     | 200.2     | 184.8     | 219.0|

verify: PASS 5 / FAIL 0

=== libmceliece 8192128f / 1コア固定 (CPU 0) ===
|------:|------:|------:|------:|------:|
|mceliece8192128f  (runs=20)
|stage      | min(ms) |median(ms)  | mean(ms)   | max(ms)|
|keygen       | 157.9     | 197.7     | 203.4     | 255.3|
|encrypt        | 0.2       | 0.2       | 0.2       | 0.3|
|decode         | 0.2       | 0.3       | 0.2       | 0.4|
|total        | 158.3     | 198.1     | 203.9     | 255.8|

verify: PASS 20 / FAIL 0

=== libmceliece 8192128 / 1コア固定 (CPU 0) ===
|------:|------:|------:|------:|------:|
mceliece8192128  (runs=20)
|stage      | min(ms) |median(ms)   |mean(ms)    |max(ms)|
|keygen       | 171.1     | 242.9      |261.3     | 551.6|
|encrypt        | 0.2       | 0.2       | 0.2       | 0.3|
|decode         | 0.2       | 0.2       | 0.2       | 0.5|
|total        | 171.5     | 243.4     | 261.7    | 552.1|

verify: PASS 20 / FAIL 0


#### 注意

- Goppa多項式 $g(x)$ は既約でなければなりません(Pattersonアルゴリズムの前提)。
- 誤りが $t$ 個を超える場合は復号に失敗します。得られた $\sigma$ の根の数が $\deg \sigma$ と一致するかを確認すると、失敗検出に使えます。

### 参考文献

- N. J. Patterson, "The algebraic decoding of Goppa codes," *IEEE Trans. Inf. Theory*, 21(2), 1975.
- V. D. Goppa, "A new class of linear error-correcting codes," *Probl. Inf. Transm.*, 6(3), 1970.
- R. J. McEliece, "A public-key cryptosystem based on algebraic coding theory," *DSN Progress Report*, 1978.
  
# ビルド

`make` は各ソースを個別にコンパイルし、`build/patterson` にリンクします。
`make test` で固定入力の多項式演算テストを実行し、`make clean` で生成物を削除できます。

共通の型・パラメーター・モジュール間宣言は `common.h`、符号パラメーターと有限体テーブルの宣言は `global.h` にあります。
各 `.c` ファイルを直接 `#include` せず、Makefile からオブジェクトファイルとしてコンパイルしてください。

現在の設定は `N=8192, K=256, E=13` です。大きな行列は `MTX` の連続したヒープ領域として確保し、`matrix_init` / `matrix_free` で初期化・解放します。鍵生成時の有限体行列やシンドローム行列もヒープ上に置き、スクランブラーの行列生成・反転にはビットパックしたワーク領域を使います。

多項式の乗算、評価、誤り位置探索、およびシンドロームの入出力には、係数を次数添字で保持する密表現 `vec` を使います。
復号の逆元・剰余・平方根・有界GCDは、既存の動作を保つため引き続き `OP` 表現です。これらを `vec` に移す場合は、復号位置の一致を確認しながら個別に変更してください。


# 20261028
timeで１秒を切る０．８秒を記録。（N=8192,K=256のとき）

AIの力を借りてボトルネックを解消！コードサイズを半分にした。
実に２年後の進化を遂げました。

更に、エラーの数がK以下の任意の個数でも正しく訂正することに成功。

２年ぶりなんですが、これは公開すべきことだと思うなり。

バイナリGoppa符号を復号するパターソンアルゴリズムと拡張パターソンアルゴリズムのテスト実装です。



# 20240221
参考文献を追加しました。

global.h　のパラメーターを設定するだけで、GF(2^E)の有限体を計算するたびに生成するように方針を変更した。

今迄、有限体で演算テーブルで参照していたのはパラメータが固定の場合にのみ早くするためでしたが、
セキュリティパラメータを任意にすることで可変式の鍵生成ができるようになりました。

# 20240220
電子署名もいいけど、やっぱり古巣の誤り訂正符号は良い感じ。
コードの誤りも訂正しないとｗ

# 20230330

旧来の開発GoppaDecorderを非公開に移行し、リポジトリを分離しました。  
これからはこれをメンテナンスします。

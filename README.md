# Goppa_Patterson

encode:Niederreiter

decode:Patterson

## Goppa符号とPattersonアルゴリズム

### Goppa符号の概要

Goppa符号は有限体上の代数的な誤り訂正符号で、McEliece暗号やNiederreiter暗号などの符号ベース暗号(耐量子計算機暗号の候補)で使われています。

#### 定義(二元の場合)

- $m, t$ を正整数とし、$\mathbb{F}_{2^m}$ 上の次数 $t$ の既約多項式 $g(x)$ を選びます(Goppa多項式)。
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

$\mathbb{F}_{2^m}$ 上の $'t \times n'$ 行列

$$
H = \begin{pmatrix}
g(\alpha_1)^{-1} & \cdots & g(\alpha_n)^{-1} \\
\alpha_1 g(\alpha_1)^{-1} & \cdots & \alpha_n g(\alpha_n)^{-1} \\
\vdots & & \vdots \\
\alpha_1^{t-1} g(\alpha_1)^{-1} & \cdots & \alpha_n^{t-1} g(\alpha_n)^{-1}
\end{pmatrix}
$$

を各成分 $\mathbb{F}_2$ 上の $m$ 次元ベクトルに展開すると、
$mt \times n$
の二元検査行列が得られます。

---

### Patterson アルゴリズム

二元Goppa符号に対する効率的な復号アルゴリズムです(Patterson, 1975)。受信語 $r = c + e$($c \in \Gamma(L,g)$、
$\mathrm{wt}(e) \le t$)から誤り $e$ を復元します。

#### 準備

- シンドローム:
  $$S(x) \equiv \sum_{i=1}^{n} \frac{r_i}{x - \alpha_i} \equiv \sum_{i \in E} \frac{1}{x - \alpha_i} \pmod{g(x)}$$
  ($E$ は誤り位置の集合)
- 誤り位置多項式:
  $$\sigma(x) = \prod_{i \in E} (x - \alpha_i)$$
- 二元の場合、形式微分 $\sigma'(x)$ は $\sigma$ を偶数部と奇数部に分けて
  $$\sigma(x) = a(x)^2 + x\, b(x)^2,\quad \sigma'(x) = b(x)^2$$
  と書けます($\deg a \le \lfloor t/2 \rfloor$、$\deg b \le \lfloor (t-1)/2 \rfloor$)。
- 鍵方程式:
  $$S(x)\,\sigma(x) \equiv \sigma'(x) \pmod{g(x)}$$

#### 手順

1. **シンドローム計算**: $S(x)$ を求める。$S(x) = 0$ なら誤りなし。
2. **逆元**: $T(x) = S(x)^{-1} \bmod g(x)$ を求める。
3. **特殊ケース**: $T(x) = x$ なら $\sigma(x) = x$ として手順6へ。
4. **平方根**: $R(x) = \sqrt{T(x) + x} \bmod g(x)$ を求める。
   - $T + x = T_0(x)^2 + x\,T_1(x)^2$ と偶奇分解し、事前計算した $\sqrt{x} \bmod g$ を用いて $R = T_0 + \sqrt{x}\,T_1$ とする。
5. **拡張ユークリッド互除法**: $a(x) \equiv b(x)\,R(x) \pmod{g(x)}$ を満たす $a, b$ を、次数条件
   $\deg a \le \lfloor t/2 \rfloor$、$\deg b \le \lfloor (t-1)/2 \rfloor$
   を満たした時点で打ち切って求め、

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

#### 注意

- Goppa多項式 $g(x)$ は既約でなければなりません(Pattersonアルゴリズムの前提)。
- 誤りが $t$ 個を超える場合は復号に失敗します。得られた $\sigma$ の根の数が $\deg \sigma$ と一致するかを確認すると、失敗検出に使えます。

### 参考文献

- N. J. Patterson, "The algebraic decoding of Goppa codes," *IEEE Trans. Inf. Theory*, 21(2), 1975.
- V. D. Goppa, "A new class of linear error-correcting codes," *Probl. Inf. Transm.*, 6(3), 1970.
- R. J. McEliece, "A public-key cryptosystem based on algebraic coding theory," *DSN Progress Report*, 1978.
- 
## ビルド

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
早すぎてごめん。



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

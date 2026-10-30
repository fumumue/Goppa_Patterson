#!/usr/bin/env python3
"""
add_prof.py  --  Patterson.c に decode の計測用マーカー(PROF_BEGIN/PROF_END)を差し込む。

使い方(Goppa_Patterson のフォルダで):
    cp prof.h .                       # prof.h を同じフォルダに置く
    python3 add_prof.py Patterson.c   # Patterson_prof.c ができる(元のファイルは変えない)

ビルドの例:
    # Patterson.c の代わりに Patterson_prof.c をコンパイルして -DPROF を付ける
    #  (Makefile の SOURCES の Patterson.c を Patterson_prof.c に差し替えるか、
    #   ファイルを入れ替えて make CPPFLAGS=-DPROF bench)

差し込む場所は「行がちょうど1回だけ現れる」ことを確認してから決める。
1つでも見つからない/複数見つかる場合は、何も書かずにエラーを出して終わる。
"""
import sys

src = sys.argv[1] if len(sys.argv) > 1 else "Patterson.c"
dst = sys.argv[2] if len(sys.argv) > 2 else "Patterson_prof.c"

# 改行コードの違い(CRLF)を吸収して読む
with open(src, encoding="utf-8", newline=None) as f:
    text = f.read()


def B(i):
    return f"  PROF_BEGIN({i});\n"


def E(i, label):
    return f'  PROF_END({i}, "{label}");\n'


L_DEC_EXPAND = "dec: 整数→ビット展開"
L_DEC_MUL = "dec: S^-1 乗算 (GF(2), 1664x1664)"
L_DEC_PACK = "dec: ビット→整数 + Setvec"
L_INV = "patterson: inv(s, w)"
L_SQRT = "patterson: sqrt(h+x) mod w"
L_CHECK = "patterson: 検証 (g1^2, gcd(s,w))"
L_ZGCD = "patterson: zgcd + 位置多項式"
L_CHIEN = "patterson: Chien探索 (chen)"
L_ERO2 = "ero2: 位置の並べ替え"
L_SYND = "zin: シンドローム (公開鍵の行XOR)"

# (アンカー行, 差し込み方, 差し込む文字列)
edits = [
    ("#include <time.h>\n", "after", '#include "prof.h"\n'),
    # --- dec() ---
    ("  // 13ビット整数 -> ビット列\n", "before", B(0)),
    ("  // GF(2) の行列・ベクトル積 (各行の AND の XOR 和 = パリティ)\n", "before",
     E(0, L_DEC_EXPAND) + B(1)),
    ("  // ビット列 -> 13ビット整数\n", "before", E(1, L_DEC_MUL) + B(2)),
    ("  return Setvec(uk, K);\n", "before", E(2, L_DEC_PACK)),
    # --- patterson() ---
    ("  OP ff = inv(f, w);\n", "before", B(3)),
    ("  // g1 = sqrt(ff + x) mod w\n", "before", E(3, L_INV)),
    ("  OP r2 = oadd(ff, tt);\n", "before", B(4)),
    ('  dbg_pol("sqrt(h+x)=======", o2v(g1));\n', "before", E(4, L_SQRT)),
    ("  b2 = omod(omul(g1, g1), w);\n", "before", B(5)),
    ("  // 位置多項式 l = beta^2 + x*alpha^2 (zgcd は1回だけ)\n", "before",
     E(5, L_CHECK)),
    ("  OP ppo = zgcd(w, g1, T);\n", "before", B(6)),
    ("  vec v = chen(o2v(ll));\n", "replace",
     E(6, L_ZGCD) + B(7) + "  vec v = chen(o2v(ll));\n" + E(7, L_CHIEN)),
    # --- ero2() ---
    ("  unsigned short xa[N] = {0};\n", "after", B(8)),
    ("    ya[i] = xa[permutation[i]];\n", "after", E(8, L_ERO2)),
    # --- zin() ---
    ("  unsigned short ss[K] = {0};\n", "after", B(9)),
    ("  // 暗号化されたシンドロームを復元する\n", "before", E(9, L_SYND)),
]

errors = []
for anchor, mode, ins in edits:
    n = text.count(anchor)
    if n != 1:
        errors.append(f"  アンカーが {n} 回見つかりました(1回であるべき): {anchor.strip()!r}")

if errors:
    print("差し込みに失敗しました。ファイルは書き出していません。", file=sys.stderr)
    print("\n".join(errors), file=sys.stderr)
    print("→ このメッセージをそのまま貼ってください。アンカーを直します。", file=sys.stderr)
    sys.exit(1)

for anchor, mode, ins in edits:
    if mode == "before":
        text = text.replace(anchor, ins + anchor, 1)
    elif mode == "after":
        text = text.replace(anchor, anchor + ins, 1)
    else:  # replace
        text = text.replace(anchor, ins, 1)

with open(dst, "w", encoding="utf-8", newline="\n") as f:
    f.write(text)

print(f"{len(edits)} か所に計測用の行を差し込みました → {dst}")

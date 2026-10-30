/*
 * sinv_pack.h  --  S^{-1} (K*E x K*E, GF(2)) をビットパックして掛ける
 *
 * 今の MTX は 1 ビットを unsigned short (2バイト) で持っているため、
 * S^{-1} は 1664 x 1664 x 2B = 約5.5MB。decode のたびにこれを全部なめている。
 * 鍵生成のときにすでにビットパック済みの S^{-1} (346KB) が手元にあるので、
 * それをコピーして持っておき、decode ではそちらを使う。
 *
 * - zin()/dec() の引数や戻り値は変えない。
 * - パック版が使えない(登録されていない/中身が変わった)ときは、
 *   従来と同じ計算に自動で切り替わる。結果は常に同じ。
 * - decode を複数スレッドから同時に呼ぶ場合は、登録と呼び出しを同じスレッドで行うこと。
 */
#ifndef SINV_PACK_H
#define SINV_PACK_H

#include <stddef.h>
#include <stdint.h>
#include "common.h" /* MTX, K, E */

/* random_full_rank_binary_matrix_inv の終わり(S^{-1} の展開のあと)から呼ぶ。
 * inv: n 行 x words 語のビットパック済み S^{-1} (行 r の列 c は inv[r*words + c/64] の c%64 ビット目)
 * n は K*E、words は (K*E+63)/64 であること(違えば何もしない)。 */
void sinv_pack_register(const MTX *Sinv, const uint64_t *inv, size_t n, size_t words);

/* 登録を捨てる。鍵を別の方法(ファイルなど)で差し替えたときに呼ぶと確実。 */
void sinv_pack_invalidate(void);

/* uk = S^{-1} * ss  (ss, uk はどちらも K 個の E ビット整数)
 * 戻り値: パック版を使えたら 1、従来計算に切り替わったら 0 */
int sinv_mul(const MTX *Sinv, const unsigned short ss[K], unsigned short uk[K]);

#endif /* SINV_PACK_H */

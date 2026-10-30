/*
 * prof.h  --  decode の内訳を測る簡易プロファイラ
 *
 * -DPROF を付けてビルドしたときだけ有効。付けなければ何も起きない(空マクロ)。
 * プログラム終了時に、段階ごとの「1回あたりの平均時間」を stderr に表示する。
 *
 * 注意:
 *  - 集計用の変数はスレッドセーフではない。decode を複数スレッドから
 *    同時に呼ぶベンチでは、数字が崩れる。
 *  - 1回の計測に clock_gettime を2回呼ぶ(数十ns)。段階が数ms単位なら誤差は無視できる。
 */
#ifndef PROF_H
#define PROF_H

#ifdef PROF

#include <stdio.h>
#include <time.h>

#define PROF_MAX 16

typedef struct
{
  const char *name;
  double sec;
  long calls;
} prof_slot;

static prof_slot prof_slots[PROF_MAX];

static inline double prof_now(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + t.tv_nsec * 1e-9;
}

#define PROF_BEGIN(id) double prof_t0_##id = prof_now()

#define PROF_END(id, label)                       \
  do                                              \
  {                                               \
    prof_slots[id].name = (label);                \
    prof_slots[id].sec += prof_now() - prof_t0_##id; \
    prof_slots[id].calls++;                       \
  } while (0)

// 終了時に自動で表示する(ベンチ側のコードを触らなくてよい)
static void __attribute__((destructor)) prof_report(void)
{
  double total = 0;
  int any = 0;

  for (int i = 0; i < PROF_MAX; i++)
    if (prof_slots[i].calls > 0)
      any = 1;
  if (!any)
    return;

  fprintf(stderr, "\n=== PROF: 1回あたりの平均 (ms) ===\n");
  fprintf(stderr, "%-3s %-44s %10s %8s\n", "id", "stage", "ms/call", "calls");
  for (int i = 0; i < PROF_MAX; i++)
  {
    const prof_slot *s = &prof_slots[i];
    if (s->calls == 0)
      continue;
    const double per = s->sec / s->calls * 1e3;
    fprintf(stderr, "%-3d %-44s %10.3f %8ld\n", i, s->name, per, s->calls);
    if (i <= 8) // 0..8 は decode 側。9 は syndrome(公開鍵の行XOR)なので別扱い
      total += per;
  }
  fprintf(stderr, "    %-44s %10.3f\n", "decode側の合計 (id 0..8)", total);
  fprintf(stderr, "(id 9 は暗号化側のシンドローム計算。合計には入れていない)\n");
}

#else

#define PROF_BEGIN(id) ((void)0)
#define PROF_END(id, label) ((void)0)

#endif /* PROF */

#endif /* PROF_H */

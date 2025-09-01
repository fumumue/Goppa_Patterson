#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

//#include "8192.h"
#include "global.h"
#include "struct.h"
#include "ecole.c"
#include "chash.c"
#include "debug.c"
#include <assert.h>
//#include "gf.h"

#define L 3 //配列の次数

// OP型からベクトル型への変換
vec o2v(OP f)
{
    vec a = {0};
    int i;

    for (i = 0; i < K * E; i++)
    {
        if (f.t[i].a > 0 && f.t[i].n < K * E)
            a.x[f.t[i].n] = f.t[i].a;
    }

    return a;
}

//多項式の次数(default)
int deg(vec a)
{
    int i, n = 0, flg = 0;

    //#pragma omp parallel for
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] > 0)
        {
            n = i;
            flg = 1;
        }
    }
    if (flg == 0)
        return 0;

    return n;
}

//ベクトル型からOP型への変換
OP v2o(vec a)
{
    int i, j = 0;
    OP f = {0};

    //#pragma omp parallel for
    for (i = 0; i < K * E; i++)
    {
        if (a.x[i] > 0)
        {
            f.t[j].n = i;
            f.t[j++].a = a.x[i];
        }
    }

    return f;
}

// 20200816:正規化したいところだがうまく行かない
//多項式の足し算
OP oadd(OP f, OP g)
{

    int i;
    vec a = {0}, b = {0}, c = {0};
    OP h = {0};

    a = o2v(f);
    b = o2v(g);

    // k=deg(o2v(f));
    // l=deg(o2v(g));

    for (i = 0; i < DEG; i++)
    {
        c.x[i] = a.x[i] ^ b.x[i];
    }
    h = v2o(c);
    // h=conv(h);
    // assert(op_verify(h));
    return h;
}

//多項式を表示する(default)
void printpol(vec a)
{
    int i, n;

    n = deg(a);

    // printf ("baka\n");
    assert(n >= 0);

    for (i = n; i > -1; i--)
    {
        if (a.x[i] > 0)
        {
            // if(a.x[i]>1)
            printf("%u", a.x[i]);
            if (i > 0)
                printf("x^%d", i);
            // if (i > 0)
            if (i > 0)
                printf("+");
        }
    }
    //  printf("\n");

    return;
}

// https://thira.plavox.info/blog/2008/06/_c.html

// invert of integer
unsigned short inv(unsigned short a, unsigned short n)
{
    unsigned short d;
    unsigned short q, t, r, x, s, gcd;

    x = 0;
    s = 1;

    d = n;
    while (a != 0)
    {
        q = d / a;
        r = d % a;
        d = a;
        a = r;
        t = x - q * s;
        x = s;
        s = t;
    }
    gcd = d;

    return ((x + n) % (n / d));
}

#define SIZE_OF_ARRAY(array) (sizeof(array) / sizeof(array[0]))
#define SWAP(type, a, b) \
    {                    \
        type work = a;   \
        a = b;           \
        b = work;        \
    }

/*
    Fisher-Yates shuffle による方法
    配列の要素をランダムシャッフルする
*/
void random_shuffle2(unsigned short *array, size_t size)
{
    unsigned short u;
    for (size_t i = size; i > 0; --i)
    {
        size_t a = i - 1;
        size_t b = rand() % i;

        array[b] = rand() % 256;
        SWAP(unsigned int, array[a], array[b]);
    }
}

int dt[16][17];
int dt2[16][17];
void mk2(void)
{
    int i, j, k;
    unsigned short r[256] = {0}, s[256] = {0};

    for (i = 1; i < 17; i++)
        dt[0][i] = 0;
    dt[0][0] = 1;
    // for(i=0;i<17;i++)
    // printf("%d,",dt[i][0]);
    // printf("\n");
    // exit(1);

    for (i = 0; i < 16 + 1; i++)
    {
        r[i] = rand() % 8192;
        dt[i][1] = r[i];
        // printf("%d,",dt[i][0]);
    }
    printf("\n");

    // exit(1);
    for (j = 0; j < 16; j++)
    {
        for (i = 1; i < 16 + 1; i++)
        {
            dt[j][i] = gf[mltn(i, fg[dt[j][1]])];
        }
    }

    for (i = 0; i < 16; i++)
    {
        for (j = 0; j < 16 + 1; j++){
            dt2[i][j]=dt[i][j];
            printf("%d,", dt2[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}



//int c[17] = {0};
//配列からベクトル表現の多項式へ変換する
vec Setvec(int n,int c[])
{
    int i;
    vec v = {0};

    for (i = 0; i < n; i++)
    {
        v.x[n - 1 - i] = c[i];
    }

    return v;
}

//配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
    OP g;
    vec a;
    int c;
    // int i;

    memset(c, 0, sizeof(c));
    memcpy(c, f, 2 * n);
    a = Setvec(n,c);

    g = v2o(a);

    return g;
}



//有限体の元の逆数
unsigned short
oinv(unsigned short a)
{
    int i;

    if (a == 0)
        return 0;

    for (i = 0; i < N; i++)
    {
        if (gf[mlt(fg[a], i)] == 1)
            return (unsigned short)i;
    }

    printf("no return \n");
    //  exit (1);
    return 0;
}

// aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
    int i;

    for (i = 0; i < N; i++)
    {
        if (gf[mlt(fg[a], fg[i])] == b)
            break;
    }
    return i;
}

//停止コマンド
void wait(void)
{
    char line[30] = {0};
    char *result;                              // 読み込む変数はローカルに取るべし
    printf(" (enter number and hit return) "); // 何か表示させたほうが良いだろう
    fflush(stdout);                            // just in case
                                               // scanf("%d", &a);                           //
    if ((result = fgets(line, 1, stdin)) != NULL)
        printf("The string is %s¥n", result);
}

// OP型を正規化する
OP conv(OP f)
{
    vec v = {0};
    OP g = {0};

    v = o2v(f);
    g = v2o(v);

    return g;
}


//項の数
int terms(OP f)
{
    int i, count = 0;

    for (i = 0; i < DEG; i++)
        if (f.t[i].a > 0)
            count++;

    return count;
}

//多項式の次数(degのOP型)
int odeg(OP f)
{
    int i, j = 0;

    if (f.t[0].a == 0)
        return 0;

    // k=terms(f);
    for (i = 0; i < DEG; i++)
    {
        if (j < f.t[i].n && f.t[i].a > 0)
            j = f.t[i].n;
    }

    return j;
}

//多項式を表示する（OP型）
void oprintpol(OP f)
{
    int i, n;

    f = conv(f);
    n = odeg(f);
    printf("n=%d\n", n);
    printf("terms=%d\n", terms(f));
    printf("deg=%d\n", odeg(f));

    for (i = n; i > -1; i--)
    {
        if (f.t[i].a > 0)
            printf("%ux^%u+", f.t[i].a, f.t[i].n);
    }
}

void op_print_raw(const OP f)
{
    puts("op_print_raw:");
    for (int i = 0; i < DEG; i++)
    {
        if (f.t[i].a > 0)
            printf("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
    }
}

OP norm(OP f)
{
    OP h = {0};
    int i;

    for (i = 0; i < DEG; i++)
    {
        if (f.t[i].a > 0)
        {
            //      h.t[f.t[i].n].n=f.t[i].n;
            h.t[f.t[i].n].a = f.t[i].a;
        }
    }

    // exit(1);

    return h;
}


//項の順序を降順に揃える
OP sort(OP f)
{
    oterm o = {0};
    int i, j, k;

    k = terms(f);
    for (i = 0; i < k + 1; ++i)
    {
        for (j = i + 1; j < k + 1; ++j)
        {
            if (f.t[i].n > f.t[j].n)
            {
                o = f.t[i];
                f.t[i] = f.t[j];
                f.t[j] = o;
            }
        }
    }

    return f;
}

//リーディングタームを抽出(OP型）
oterm oLT(OP f)
{
    int i, k, j;
    oterm s = {0};

    k = terms(f);
    s = f.t[0];
    for (i = 0; i < k + 1; i++)
    {
        // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
        if (f.t[i].a > 0)
        {
            printf("in LT=%d %d\n", s.a, s.n);
            for (j = i; j < k + 1; j++)
            {
                if (s.n < f.t[j].n)
                {
                    s.n = f.t[j].n;
                    s.a = f.t[j].a;
                }

                //  else{
                // t=s;
                // }
            }
        }
    }
    //  exit(1);

    return s;
}

//多項式を足し算する（OP型）
OP add(OP f, OP g)
{
    //  vec a={0},b={0},c={0};
    unsigned long long int i, j, n1 = 0, n2 = 0, m1 = 0, count = 0;
    OP h = {0};
    oterm o1 = {0}, o2 = {
                        0};

    n1 = terms(f);
    printf("n1=%llu\n", n1);
    n2 = terms(g);
    printf("n2=%llu\n", n2);

    oprintpol(f);
    printf(" fff==============\n");
    oprintpol(g);
    printf(" ggg==============\n");
    o1 = oLT(f);
    o2 = oLT(g);
    printf("LTadd==%d %d\n", o1.n, o2.n);
    m1 = n1 + n2;
    printf("m1=%llu\n", m1);
    // exit(1);

    for (i = 0; i < n1 + 1; i++)
    {
        for (j = 0; j < n2 + 1; j++)
        {
            if (f.t[i].n == g.t[j].n && g.t[j].a > 0 && f.t[i].a > 0)
            {
                o1 = oLT(f);
                o2 = oLT(g);
                printf("LT==%d %d\n", o1.n, o2.n);
                printf("f.n==%d %d %llu %llu\n", f.t[i].n, g.t[j].n, i, j);
                f.t[i].a = 0;
                g.t[j].a = 0;
            }
        }
    }
    for (i = 0; i < n2 + 1; i++)
    {
        if (g.t[i].a > 0)
        {
            h.t[count++] = g.t[i];
            g.t[i].a = 0;
        }
    }
    for (i = 0; i < n1 + 1; i++)
    {
        if (f.t[i].a > 0)
        {
            h.t[count++] = f.t[i];
            f.t[i].a = 0;
        }
    }

    h = sort(h);
    /*
     for (i=0; i<count; ++i) {
     for (j=i+1; j<count; ++j) {
     if (h.t[i].n > h.t[j].n) {
     oo =  h.t[i];
     h.t[i] = h.t[j];
     h.t[j] = oo;
     }
     }
     }
   */
    h = conv(h);
    if (odeg(h) > 0)
        oprintpol(h);
    printf(" addh==============\n");
    //   exit(1);

    return h;
}

//多項式を項ずつ掛ける
OP oterml(OP f, oterm t)
{
    f = conv(f);

    int i;
    OP h = {0};
    // vec test;
    // unsigned int n;

    // f=conv(f);
    // k = deg (o2v(f));

    for (i = 0; i < DEG; i++)
    {
        h.t[i].n = f.t[i].n + t.n;
        h.t[i].a = gf[mlt(fg[f.t[i].a], fg[t.a])];
    }

    h = conv(h);
    //assert(op_verify(h));
    return h;
}

//多項式の掛け算
OP omul(OP f, OP g)
{
    int i, k, l;
    oterm t = {0};
    OP h = {0}, e = {0}; //, r = {0};
    // vec c = {0};

    k = odeg(f);
    l = odeg(g);
    if (l > k)
    {
        k = l;
    }

    for (i = 0; i < k + 1; i++)
    {
        t = g.t[i];
        e = oterml(f, t);
        h = oadd(h, e);
        // printpol(o2v(h));
        // printf(" ==kokko\n");
    }
    //assert(op_verify(h));
    return h;
}


/*
 * Author: Hiroyuki Chishiro
 * License: 2-Clause BSD
 */
#include <stdio.h>
#define U 16

OP gp(void)
{
    int a[U][U + 1] = {
        {1, 2, 3, 5},
        {2, 1, 1, 6},
        {1, 3, 5, 2}

    };
    vec ch={0};
    //int o[17] = {1, 0, 2688, 5694, 972, 871, 4970, 6299, 7459, 7621, 7916, 7340, 3256, 1984, 5228, 3249, 7526};
    int p, d;
    int i, j, k;

    /*
      x1 = 0
    x2 = 2688
    x3 = 5694
    x4 = 972
    x5 = 871
    x6 = 4970
    x7 = 6299
    x8 = 7459
    x9 = 7621
    x10 = 7916
    x11 = 7340
    x12 = 3256
    x13 = 1984
    x14 = 5228
    x15 = 3249
    x16 = 7526
    */
    for (i = 0; i < U; i++)
    {
        p = dt2[i][i];

        for (j = 0; j < (U + 1); j++)
        {
            dt2[i][j] = gf[mlt(fg[dt2[i][j]], inv(p, 8192))];
        }

        for (j = 0; j < U; j++)
        {
            if (i != j)
            {
                d = dt2[j][i];

                for (k = i; k < (U + 1); k++)
                {
                    dt2[j][k] = dt2[j][k] ^ gf[mlt(fg[d], fg[dt2[i][k]])];
                }
            }
        }
    }
 
     ch.x[U]=1;
    for (i = 0; i < U; i++)
       ch.x[15-i]=dt2[i][U];
       
    for (i = 0; i < U; i++){
        printf("x%d = %d\n", i + 1, ch.x[i]);
    }

    OP f = {0};
    f = v2o(ch);
    printpol(o2v(f));
    printf("\n");

    return f;
}

//リーディグタームを抽出(default)
oterm LT(OP f)
{
    int i;
    oterm t = {0};

    // k = deg (o2v (f));
    for (i = 0; i < DEG; i++)
    {
        // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
        if (f.t[i].a > 0)
        {
            t.n = f.t[i].n;
            t.a = f.t[i].a;
        }
    }

    return t;
}



//モニック多項式にする
OP coeff(OP f, unsigned short d)
{
    int i, k;

    f = conv(f);
    k = odeg((f)) + 1;
    for (i = 0; i < k; i++)
        f.t[i].a = gf[mlt(fg[f.t[i].a], oinv(d))];

    return f;
}


//多項式を単行式で割る
oterm LTdiv(OP f, oterm t)
{
    oterm tt = {0}, s = {
                        0};

    tt = LT(f);
    if (tt.n < t.n)
    {
        s.n = 0;
        s.a = 0;
    }
    else if (tt.n == t.n)
    {
        s.n = 0;
        s.a = equ(t.a, tt.a);
    }
    else if (tt.n > t.n)
    {
        s.n = tt.n - t.n;
        s.a = equ(t.a, tt.a);
        // printf("%u\n",s.a);
    }
    else if (t.n == 0 && t.a > 0)
    {
        s.a = gf[mlt(fg[tt.a], oinv(t.a))];
        s.n = tt.n;
    }

    return s;
}


//多項式の剰余を取る
OP omod(OP f, OP g)
{
    OP h = {0};
    oterm b = {0}, c = {0};
    int n = LT(g).n;

    //  assert (("baka^\n", LT (f).n != 0));

    //  assert (("baka(A)\n", LT (g).n != 0));

    if (LT(f).n < n)
    {
        //    exit(1);
        return f;
    }

    // printf ("in omod\n");
    // exit(1);

    b = LT(g);

    assert(b.a != 0 && b.n != 0);
    while (LT(f).n > 0 && LT(g).n > 0)
    {

        c = LTdiv(f, b);
        h = oterml(g, c);
        f = oadd(f, h);
        if (odeg((f)) == 0 || odeg((g)) == 0)
        {
            //      printf("blake1\n");
            break;
        }

        if (c.n == 0 || b.n == 0)
            break;
    }

    return f;
}

//多項式の商を取る
OP odiv(OP f, OP g)
{

    int i = 0;
    OP h = {0}, tt = {0};
    oterm b = {0}, c = {0};

    if (LT(f).n == 0 && LT(g).a == 0)
    {
        printf("baka^\n");
        // return f;
        exit(1);
    }
    if (LT(g).a == 0)
    {
        print_trace();
        exit(1);
    }
    if (LT(g).n == 0 && LT(g).a > 1)
        return coeff(f, LT(g).a);

    //    k = odeg(g);
    b = LT(g);
    if (b.a == 1 && b.n == 0)
        return f;
    if (b.a == 0 && b.n == 0)
    {
        printf("baka in odiv\n");
        exit(1);
    }
    if (odeg((f)) < odeg((g)))
    {
        return f;
        //  a=LT(f);
    }

    i = 0;
    while (LT(f).n > 0 && LT(g).n > 0)
    {
        c = LTdiv(f, b);
        assert(c.n < DEG);
        tt.t[i] = c;
        i++;

        h = oterml(g, c);

        f = oadd(f, h);
        if (odeg((f)) == 0 || odeg((g)) == 0)
        {
            // printf ("blake2\n");
            break;
        }

        if (c.n == 0)
            break;
    }

    // tt は逆順に入ってるので入れ替える
    OP ret = {0};
    int tt_terms = terms(tt);
    for (i = 0; i < tt_terms; i++)
    {
        ret.t[i] = tt.t[tt_terms - i - 1];
    }
    ret = conv(ret);
    //assert(op_verify(ret));
    return ret;
}

//多項式のべき乗
OP opow(OP f, int n)
{
    vec v = {0};
    int i;
    OP g = {0};

    v.x[0] = 1;
    if (n == 0)
        return v2o(v);

    g = f;
    for (i = 1; i < n; i++)
    {
        g = omul(f, g);
        printpol(o2v(g));
        printf("\n");
    }
    return g;
}

/*
OP opowmod(OP f, OP mod, int n) {
    vec v={0};
    OP ret;

     v.x[0]= 1;
    ret=v2o(v);
    while (n > 0) {
        if (n & 1) ret = omod(omul(ret , f),mod) ;  // n の最下位bitが 1 ならば x^(2^i) をかける
        f = omod(omul(f , f),mod);
        n >>= 1;  // n を1bit 左にずらす
    }
    return ret;
}
*/

//多項式のべき乗余
OP opowmod(OP f, OP mod, int n)
{
    int i;

    // printpol(o2v(mod));
    // printf(" =mod %d\n",LT(mod).n);
    //繰り返し２乗法
    for (i = 1; i < n + 1; i++)
    {
        f = omod(omul(f, f), mod);
    }

    return f;
}


OP gcd(OP a, OP b)
{
    OP r = {0}, h = {0}, tmp = {0};

    h.t[0].a = 1;
    h.t[0].n = 0;

    if (odeg(a) < odeg(b))
    {
        tmp = a;
        a = b;
        b = tmp;
    }
    /*
  printpol(o2v(a));
  printf(" ========f\n");
  printpol(o2v(b));
  printf(" ========g\n");
*/
    /* 自然数 a > b を確認・入替 */
    if (odeg(a) < odeg(b))
    {
        tmp = a;
        a = b;
        b = tmp;
    }

    r = omod(a, b);
    while (odeg(r) > 0)
    {
        a = b;
        b = r;
        r = omod(a, b);
        if (LT(r).a == 0)
            return b;
    }

    if (LT(r).a == 0)
    {
        return b;
    }
    else
    {
        // if(LT(r).a>0)
        return h;
    }
}


// GF(2^m) then set m in this function.
int ben_or(OP f)
{
    int i, n, flg = 0;
    OP s = {0}, u = {0}, r = {0};
    vec v = {0};
    // if GF(8192) is 2^m and m==13 or if GF(4096) and m==12 if GF(16384) is testing
    int m = E;
    // m=12 as a for GF(4096)=2^12 defined @ gloal.h or here,for example m=4 and GF(16)

    v.x[1] = 1;
    s = v2o(v);
    r = s;
    n = deg(o2v(f));
    printf("n=%d\n", n);

    if (n == 0)
        return -1;

    i = 0;

    // r(x)^{q^i} square pow mod
    while (i < n / 2 )
    {
        printf("iii=%d\n", i);
        flg = 1;
        // irreducible over GH(8192) 2^13
        r = opowmod(r, f, m);

        // irreducible over GF2
        // r=omod(opow(r,2),f);

        u = oadd(r, s);
        if (deg(o2v(u)) == 0 && LT(u).a == 0)
            return -1;
        if (deg(o2v(u)) == 0 && LT(u).a == 1)
        {
            i++;
            flg = 0;
        }
        if (deg(o2v(u)) > 0)
            u = gcd(f, u);

        if (deg(o2v(u)) > 0)
            return -1;

        if (flg == 1)
            i++;
    }


    return 0;
}


int main()
{
    unsigned short a[L][L] = {0}; //{{2,-2,4,2},{2,-1,6,3},{3,-2,12,12},{-1,3,-4,4}};
    unsigned short det = 1, buf;
    int i, j, k;
    OP q[L][L] = {0};
    vec v[L][L] = {0};
    OP p[L][L] = {0};
    OP o[L][L] = {0};
    vec s[L][L] = {0};
    unsigned short arr[16] = {0};
    OP f={0};

   
    srand(clock());
   /*
    random_shuffle(arr, SIZE_OF_ARRAY(arr));
    for (i = 0; i < 16; i++)
        printf("%d,", arr[i]);
    printf("\n");
    exit(1);
*/
    //mk2();
    //f=gp();
    //printf("irr=%d\n",ben_or(f));
    //exit(1);

    for (i = 0; i < L; i++)
    {
        v[i][i].x[0] = 3;
        q[i][i] = v2o(v[i][i]);
        printpol(v[i][i]);
        printf("\n");
    }
    printf("\n");
    for (i = 0; i < L; i++)
    {
        s[i][i].x[1] = 1;
        p[i][i] = v2o(s[i][i]);
    }
    for (i = 0; i < L; i++)
    {
        // for(j=0;j<L;j++){
        //     for(k=0;k<L;k++)
        o[i][i] = oadd(q[i][i], p[i][i]);
        //    }
    }
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
            printpol(o2v(o[i][j]));
        printf("\n");
    }
    //exit(1);

    srand(clock());

    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
            a[i][j] = rand() % N;
    }

    //三角行列を作成
    for (i = 0; i < L; i++)
    {
        for (j = 0; j < L; j++)
        {
            if (i < j)
            {
                buf = gf[mlt(fg[a[j][i]], fg[inv(a[i][i], N)])];
                for (k = 0; k < L; k++)
                {
                    a[j][k] ^= gf[mlt(fg[a[i][k]], fg[buf])];
                }
            }
        }
    }
    for(i=0;i<L;i++){
        for(j=0;j<L;j++)
        printf("%d,",a[i][j]);
    printf("\n");
    }
    exit(1);

//対角部分の積
    for (i = 0; i < L; i++)
    {
        det = gf[mlt(fg[det], fg[a[i][i]])];
    }

    printf("%u\n", det); // -> 120.000000

    return 0;
}
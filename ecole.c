/************************************************************
 * ecole.c - generate GF(2^n) using primitive polynomial
 *  ゼフ対数表を作るプログラム。正規基底を生成します。
 *  GF : Galois Field
 * Author: Osamu Sakai, SIC(SAKAI Institute of Cryptology)
 *
 *  Edition History
 * version date       comment                                     by
 * ------- ---------- ------------------------------------------- ---------
 * V0.1    1995/12/01 first version                               O.Sakai
 * V0.5    2022/12/01 beta version                                O.Sakai
 * V0.6    2022/12/10 refactoring                                 rubato6809
 * V0.61   2022/12/20 changed not to support order=65536          O.Sakai
 * V0.62   2022/12/30 remove a order 65536                        O.Sakai
 * V0.8 　 2022/12/31 change preprocessor to command line opton   rubato-Sakai
 * Copyright(c) 1995 - 2022 SIC, All Rights Reserved
 * Special thanks to Mr.Rubato
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_ORD 32768

/***************************************************************
 * 関数名     : void gen_gf(int deg, int order)
 * 機能       : 引数 deg をもつ原始多項式 normal[] を使い、
 *              order 個の有限体を生成する。
 * 入力引数   : int deg   （ビットサイズ）
 *              int order （== 2 ^ deg)
 * 戻り値     : none
 * 入力情報   : none
 * 出力情報   : gf[] : 生成したZech対数の正引き配列
 *              fg[] : 逆引き配列
 * 注意事項   : 通常は normal[] 配列を使って生成します。
 *              sage[] 配列も原始多項式のリスト。暗号の仕様によっては
 *              同じ有限体の生成方法が取られるとは限らないので、
 *              テーブルを分けてあります。
 *              sage計算代数ソフトで計算結果が正しいか検証します。
 *              sage[] 配列を使用する場合は、
 *              コマンドラインオプションで "-opt_s" を指定するとsage用の
 *              有限体を生成することができます。
 *  ****************************************************************/
// Zech対数の正引き gf と逆引き fg
static unsigned short gf[MAX_ORD];
static unsigned short fg[MAX_ORD];
void gen_gf(int deg, int order, int opt_sage)
{
    // #ifdef SAGE
    //  Generate Sagemath based Galois Fields.
    static const unsigned int sage[] = {
        0b111,
        0b1011,
        0b10011,
        0b100101,
        0b1011011,
        0b10000011,
        0b100011101,      // sage
        0b1000010001,     // sage512
        0b10001101111,    // sage1024
        0b100000000101,   // 2048
        0b1000011101011,  // sage 4096
        0b10000000011011, // Classic McEliece
        0b100000010101001,
        0b1000000000110101,
    };
    // unsigned short x = sage[deg - 2];
    // #else
    /* Generate nomal basis of Galois Field over GF(2^?) */
    static const unsigned int normal[] = {
        0b111,
        0b1101,
        0b11001,
        0b110111,
        0b1100001,
        0b11000001,
        0b110101001,
        0b1100110001,
        0b11000010011,
        0b110000001101,
        0b1100101000001,
        0b11011000000001,
        0b110000100010001,
        0b1100000000000001,
    };
    // unsigned short x = normal[deg - 2]; // 通常はこちら
    //  #endif

    unsigned short x;
    if (opt_sage){
        x = sage[deg - 2];
        printf("DEBUG: use sage[]\n");
    } else {
        x = normal[deg - 2];
        printf("DEBUG: use normal[]\n");
    }

    /* build gf[] */
    gf[0] = 0;
    gf[1] = 1;
    /*
    for (int i = 2; i < order; i++)
    {
        gf[i] = (gf[i - 1] << 1); // gf[i] を更新する
        if (gf[i] >= order)       // gf[i] 検査
            gf[i] ^= x;
    }
    */

    unsigned short value = 1;
    for (int i = 2; i < order; i++)
    {
        value <<= 1;
        if (value >= order)
            value ^= x;

        gf[i] = value;
    }
    /* build fg[] */
    for (int i = 0; i < order; i++)
        fg[gf[i]] = i;
}

void toFile(FILE *fp, int order, unsigned short *array, char *name)
{
    fprintf(fp, "static const unsigned short %s[%d]={", name, order);
    for (int i = 0; i < order; i++)
        fprintf(fp, "%d,", array[i]);
    fprintf(fp, "};\n");
}

/***************************************************************
 * 関数名     : void put_gf(int order)
 * 機能       : Zech対数表をファイルに書き出す。
 *
 * 入力引数   : int order
 * 出力引数   : none
 * 戻り値     : none
 * 入力情報   : none
 * 出力情報   : gf[2^i],fg[2^i],i=2~15をヘッダファイル 2^i.h に書き出す。
 * 注意事項   :     order       | generate Zech logarithm
 *              ----------------+-------------
 *                    4 = 2^2   |         GF(4)
 *                    8 = 2^3   |         GF(8)
 *                    :         |         :
 *                32768 = 2^15  |        GF(32768)
 ****************************************************************/
void put_gf(int order)
{
    char filename[8];
    sprintf(filename, "%d.h", order);
    FILE *fp = fopen(filename, "wb");

    toFile(fp, order, gf, "gf");
    toFile(fp, order, fg, "fg");

    fclose(fp);
}

void usage(void)
{
    printf("Usage : ./a.out [option] 4,8,...,32768 i.e. 2^i: i=2..15.\n");
    printf("Option --sage:sagemath basis,none:normal basis\n");
    exit(1);
}

/***************************************************************
 * 関数名     : int bitsize(int num)
 * 機能       : 引数 num が 2 のべき乗数ならば、num のサイズ nbit を返す.
 *              許容する num の範囲は 4(=2^2) から 32768(=2^15) まで.
 *              許容しない値なら、使用法を表示し、exit(1) する.
 *
 * 入力引数   : int num
 * 出力引数   : none
 * 戻り値     : 2 ～ 15
 *              エラー（許容しない値）の場合、戻り値無し。プロセス終了
 * 入力情報   : none
 * 出力情報   : none
 * 注意事項   :     num         | return value
 *              ----------------+-------------
 *                    2 = 2^1   |   exit(1)
 *                    4 = 2^2   |         2
 *                    8 = 2^3   |         3
 *                    :         |         :
 *                32768 = 2^15  |        15
 *                65536 = 2^16  |   exit(1)
 ****************************************************************/
int bitsize(int num)
{
    int power2 = 4; // 2^2
    for (int nbits = 2; nbits < 16; nbits++)
    {
        if (num == power2)
            return nbits; // 当り。ビット数を返す

        power2 <<= 1; // 4, 8, 16 ... 32768
    }
    /* ここまできたらハズレ */
    usage();  // 使用方法を表示してexitする。
    return 0; // （警告を避けるため）
}

/***************************************************************
 * 関数名     : int opt(int argc, char *argv[], int *k, int *c)
 * 機能       : 引数 num が 2 のべき乗数ならば、num のサイズ nbit を返す.
 *              許容する num の範囲は 4(=2^2) から 32768(=2^15) まで.
 *              許容しない値なら、使用法を表示し、exit(1) する.
 *
 * 入力引数   : int argc,char *argv[],int *k,int *opt_sage
 * 出力引数   : none
 * 戻り値     : 0 or 1
 *              エラー（許容しない値）の場合、戻り値無し。プロセス終了
 * 入力情報   : none
 * 出力情報   : none
 * 注意事項   : if opt_sage=0: select normal basis;
 *              if opt_sage=1: select basis of sagemath
 *              else exit(1);
 ****************************************************************/
void opt(int argc, char *argv[], int *k, int *opt_sage)
{
    *k = atoi(argv[argc - 1]);
    // 2 patterns are acceptable: argc must be 2 or 3.
    if (argc == 3 && strcmp(argv[1], "--sage") == 0)
        *opt_sage = 1; // --sage 有り、e.g. $ecole --sage 16
    else if (argc == 2)
        *opt_sage = 0; // --sage 無し、e.g. $ecole 16
    else // otherwise, error-exit
        usage();
}

/*
int main(int argc, char *argv[])
{
    int s, k;
    opt(argc, argv, &k, &s);
    int n = bitsize(k);
    gen_gf(n, k, s);
    put_gf(k);
    printf("GF[%d] の生成に成功しました。\n", k);

    return 0;
}
*/
/*-----------------------------------------------------------------------------
     Program Poisson2.c
     二次元Poisson方程式シミュレーション
-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>


/* マクロ定義 ---------------------------------------------------------------*/
#define N   100

/* メインルーチン -----------------------------------------------------------*/
int main(void) {

    const double X     = 1.0       ;     /* 計算領域の大きさ                 */
    const double e0    = 8.85e-12  ;     /* 真空の誘電率                     */
    const int    center= (int)(N/2);     /* 中心の座標                       */
    const double delta = X/N       ;     /* δ                               */
    const double Conv  = 1.0e-10    ;     /* 収束と判定する前回ループとの差   */

    double phi[N][N]               ;     /* 計算するべき電位                 */
    double rho[N][N]               ;     /* 電荷密度                         */
    double metal[N][N]             ;     /* 接地された導体                   */
    double MaxPhi                  ;     /* 最大電位                         */
    double MaxErr                  ;     /* 最大のエラー                     */
    double CurErr                  ;     /* 現在のエラー                     */
    double Prev_phi                ;     /* 前のループのφ                   */
    double Ex, Ey                  ;     /* 電場                             */
    int    i, j                    ;
    int    loop                    ;     /* 繰り返しカウンタ                 */
    FILE   *f                      ;     /* ファイルハンドラ                 */


    const int R      = center-1;

    /* phi(i, j)，rho(i, j)をクリアする */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            phi[i][j] = rho[i][j] = metal[i][j] = 0.0;
        }
    }

    /* 電荷を置く */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
          rho[i][j] = 0.0;
        }
    }

    int inner=10;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
          if ((i-center)*(i-center) + (j-center)*(j-center) >= R*R) {
            phi[i][j] = 100;
            metal[i][j] = 1;
          }
          if ((i-center)*(i-center) + (j-center)*(j-center) < inner*inner) {
            phi[i][j] = 0;
            metal[i][j] = 1;
          }
        }
    }

    /* 繰り返し計算 */
    loop   = 0;         /* 別に無くても良いのだが，目安としてループをカウントする．              */
    MaxPhi = 1.0e-10;   /* 系内の最大の電位を入れる変数．ある有限の値を入れておく(ゼロ割り防止)．*/

    do {
        if(!(loop%1000)) printf("%05d  %e\n", loop, MaxPhi); /* 10000ループ毎に経過表示 */
        MaxErr = CurErr = 0.0;
        for (i = 1; i < N-1; i++) {      /* 領域端を除く全ての点をループ                */
            for (j = 1; j < N-1; j++) {  /*                                             */
                if(metal[i][j]==0) {
                    Prev_phi = phi[i][j];    /* 前回ループのphiをPrev_phiにいれておいて，   */
                    phi[i][j] = 0.25*(rho[i][j]*delta*delta/e0+phi[i+1][j]+phi[i-1][j]+phi[i][j+1]+phi[i][j-1]);
                                         /* Poissonの方程式でphiを計算する．            */
                    if (MaxPhi < fabs(phi[i][j])) MaxPhi = phi[i][j];
                                         /* 電位最大が更新されたらMaxPhiを書き換え      */
                    CurErr = (fabs(phi[i][j] - Prev_phi))/MaxPhi;
                                         /* 前回ループと新しい答えの差を，MaxPhiで規格化*/
                    if (MaxErr < CurErr) MaxErr = CurErr;
                }                        /* 誤差の最大を常にMaxErrに持つようにする      */
            }
        }
        loop++;
    } while (MaxErr>Conv);               /* 領域全ての点の誤差がConvを下回ったらおしまい*/

    /* ポテンシャル出力 */
    f = fopen("Phi.avd", "wt");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            fprintf(f, "%e %e %e\n", delta*i, delta*j, phi[i][j]);
        }
    }
    fclose(f);

    f = fopen("Phi.gdw", "wt");
    for (i = 0; i < N; i++) {
        fprintf(f, "%e %e\n", delta*i, phi[i][center]);
    }
    fclose(f);

    /* 電場出力 */
    f = fopen("Electric.avd", "wt");
    for (i = 1; i < N-1; i++) {
        for (j = 1; j < N-1; j++) {
            Ex = -(phi[i+1][j]-phi[i-1][j])/(2.0*delta);
            Ey = -(phi[i][j+1]-phi[i][j-1])/(2.0*delta);
            fprintf(f, "%e %e %e %e %e\n", delta*i, delta*j,  sqrt(Ex*Ex+Ey*Ey), Ex, Ey);
        }
    }
    fclose(f);

    f = fopen("Electric.gdw", "wt");
    for (i = 0; i < N; i++) {
        Ex = -(phi[i+1][center]-phi[i-1][center])/(2.0*delta);
        fprintf(f, "%e %e\n", delta*i, Ex);
    }
    fclose(f);

    return 0;
}

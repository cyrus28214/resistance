#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double **array2d(int rows, int cols);
void gauss(double **a, int n);

int main() {
    int n, //节点数量
        m, //边的数量
        s, //源点
        t; //汇点
    //图不能有自环，可以有重边，需要联通

    scanf("%d%d%d%d", &n, &m, &s, &t);

    //为了代码方便，增加一个变量I_总
    int k = n + m + 1; //方程的个数

    double **a = array2d(k, k+1);//多出的一列是等号右边的常数
    //每一列代表的变量分别为：
    // I_0 ... I_{m-1} | phi_0 ... phi_{n-1} | I_总
    //I_i = a[][i]
    //phi_i = a[][m+i]

    //各个方程为
    //0 ~ m-1 欧姆定律： phi_u - phi_v - I_i R_i = 0
    //m ~ m+n-1 KCL:    sum I_i = 0 (对非源点汇点)
    //                  phi_s = 1 (对源点，加1V电压)
    //                  phi_t = 0 (对汇点，接地)
    //m+n：             计算总电流I_总=汇入t电流

    //m行，每行三个整数u, v, r表示u到v有一个阻值为r的电阻
    for (int i = 0; i < m; i++) {
        int u, v;
        double r;
        scanf("%d%d%lf", &u, &v, &r);

        //欧姆定律
        a[i][m+u] = 1;
        a[i][m+v] = -1;
        a[i][i] = -r;
        //第i个方程 phi_u - phim+n- I_i R_i = 0

        //KCL 注意电流方向，这里按流出为负，流入为正
        if (u != s && u != t) a[m+u][i] = -1;
        else if (u == t) a[m+n][i] = -1;
        if (v != s && v != t) a[m+v][i] = 1;
        else if (v == t) a[m+n][i] = 1;
    }

    a[m+s][m+s] = 1, a[m+s][k] = 1; //源点电电势为1V
    a[m+t][m+t] = 1, a[m+t][k] = 0; //汇点电势为0
    a[m+n][m+n] = -1; //计算总电流I_总

    gauss(a, k); //求解方程组

    double I_total = a[k-1][k];
    printf("%6e\n", 1 / I_total);//由于我们设定了整体电压为1，最后的等效电阻就是1/I_总
}

//生成一个rows行cols列的二维数组，全部初始化成0
double **array2d(int rows, int cols) {
    double **arr = malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; i++) {
        arr[i] = calloc(cols, sizeof(double));
    }
    return arr;
}

void gauss(double **a, int n) {
    //高斯-约旦消元法求解线性方程组
    for (int i = 0; i < n; i++) { //当前行
        //将当前行替换为列主元最大的行
        for (int j = i + 1; j < n; j++) {
            if (fabs(a[j][i]) > fabs(a[i][i])) {
                double *temp = a[i];
                a[i] = a[j];
                a[j] = temp;
            }
        }

        double c = a[i][i];
        //通常的算法此处需要特判c是否为0，但根据物理意义，本问题一定有唯一解，可以不用。
        for (int k = 0; k <= n; k++) 
            a[i][k] /= c;

        //消元
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            c = a[j][i];
            for (int k = 0; k <= n; k++)
                a[j][k] -= c * a[i][k];
        }
    }
}
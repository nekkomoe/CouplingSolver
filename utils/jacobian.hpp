#pragma once
#include <cassert>
/** 
 * nonlinear function type
 * au: unknowns, aF: function values
 */
typedef void (*nlf)(double *au, double *aF); 
const double sqrt_mp = sqrt(2.22E-16); // sqrt(machine precision)
const double xmin = 1E-6;
#ifdef SPEC_H
  #define calc_h(x) ((x>xmin?x:(x>0?xmin:-xmin))*sqrt_mp)
#else
  #define calc_h(x) 1e-5 // FIXME 这里需要选取上面的公式
#endif

struct List {
  int n; // number of colors
  int *head, *rest, *to, tot;
  List(int n, int m) {
    this->n = n;
    head = new int[n];
    rest = new int[m+1];
    to   = new int[m+1];
    memset(head, 0, sizeof(int)*n);
    tot = 0;
  }
  ~List() {
    delete[] head;
    delete[] rest;
    delete[] to;
  }
  void add(int u, int v) {
    to[++ tot] = v;
    rest[tot] = head[u];
    head[u] = tot;
  }
};

void gc_jac(nlf F,  SMat<double> *jac, int n, double *au, double *aF, List *cc, List *ri,
            double *aFh, double *au_t, double *h,
            int *cols, int *cnt
           ) {
  /**
   * graph-coloring jacobian
   * jac: jacobian  SMatrix, nnz-structure is required(call setlen(int*))
   * n: number of unknowns
   * au: unknowns
   * cc: column-coloring
   * ri: row-index of each column
   */
  memset(cnt, 0, sizeof(int)*n);
  jac->reshape(n);                              // nxn sparse  SMatrix
  F(au, aF);                                    // aF = F(au)
  // 枚举列, 计算h数组
  for(int col = 0 ; col < n ; ++ col) {
    double x = au[col];
    h[col] = calc_h(x);
  }

  for(int u = 0 ; u < cc->n ; ++ u) {
    // 1) calculate h and aFh
    if(!cc->head[u]) {
      // continue; // FIXME 当前颜色没有列节点, 理论上来说需要直接break, 需要验证一下颜色是否是连续编号的
      break;
    }
    int sz = 0; // 当前颜色的节点数量
    for(int i = cc->head[u] ; i ; i = cc->rest[i]) {
      int col = cc->to[i];                // col-idx
      cols[sz] = col;
      au_t[sz] = au[col];                 // 存储au的值, 用于恢复au [TODO]这个可以存到cache里
      au[col] += h[col];
      ++ sz;
    }
    F(au, aFh);                           // aFh = F(au+h)
    // 2) calculate jacobian
    for(int j = 0 ; j < sz ; ++ j) {
      int col = cols[j]; // 枚举列
      for(int i = ri->head[col] ; i ; i = ri->rest[i]) {
        int row = ri->to[i];             // 枚举列的所有非零行
        double Jji = (aFh[row] - aF[row])/h[col];
        // if(Jji != 0 || col == row) {
          // 这里不需要begin/end, 因为已经通过nnzcpy设置了结构
          jac->add(row, col, Jji, cnt[row]);
          ++ cnt[row];
        // }
      }
    }
    // 3) 恢复 au
    for(int i = 0 ; i < sz ; ++ i) {
      au[cols[i]] = au_t[i];
    }
  }
}

#pragma once
template<typename T>
PetscErrorCode TimeSolver(SMat<T> *jac, List *cc, List *ri) {
  /**
   * 时间迭代求解
   * colors: 着色方案
   * jac: Jacobian矩阵(要求已经设置了非零元素结构)
   */

  // 1) 初始化Jacobian矩阵
  int n = pf_total_size;
  DataSaver<T> datasaver(n, Nt);
  jac->reshape(n);
  T *delta = new T[n],
    *rhs   = new T[n];
  
  memset(delta, 0, sizeof(T)*n);     // 初始化delta

#define GetJacobian(au, rhs) gc_jac(F, jac, n, au, rhs, cc, ri, aFh, au_t, h, cols, cnt)

  double norm_F        = 0;    // F的2范数

  double *aFh  = new double[n],
         *au_t = new double[n],
         *h    = new double[n];
  int    *cols = new int[n],
         *cnt  = new int[n];

  int newton_cnt = 0, iter; double resi;
#ifndef SEP_CHAR
  #define SEP_CHAR '\r'
#endif
  ProgressIteration(iter_t, IterMax, SEP_CHAR, [&]() {
    char info[100];
    snprintf(info, 100, "%d:%.10f (%d)", newton_cnt, norm_F, iter);
    return std::string(info);
  }) {
    Current_Iter = iter_t;
    // 1) 初始化pf_old
    InitOldField();

    // 2) 目标: J(x)(x'-x)=-F(x)
    newton_cnt = 0;
    do {
      GetJacobian(pf->_data_ptr, rhs);                        // 求解 J(x) 和 rhs=F(x)
      norm_F = norm2(n, rhs);                                 // 计算rhs(F)的2范数
      // 这里没有把delta初始化为0，因为使用上一个迭代步的结果作为初值更好
      PetscCall(LinearSolver(jac, rhs, delta, &iter, &resi)); // 求解J(x)(x-x')=F(x)
      for(int i = 0 ; i < n ; ++ i) {
        pf->_data_ptr[i] -= delta[i];                         // x' = x-(x-x') (更新答案*-1)
      }
      __my_progressbar.update_info();
      ++ newton_cnt;
      datasaver.feed_ksp(iter_t, iter, resi, norm_F);
    } while(norm_F > newton_abstol && newton_cnt <= newton_maxcnt);
    datasaver.feed_newton(iter_t);
  }
  datasaver.restore();
  delete[] aFh;
  delete[] au_t;
  delete[] h;
  delete[] cols;
  delete[] cnt;
  delete[] delta;
  delete[] rhs;
  return PETSC_SUCCESS;
}

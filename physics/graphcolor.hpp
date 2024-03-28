#pragma once

/**
 * 计算Jacobian的图着色
 */
int GetGraphColor( SMat<int> *J, int *colors, int *_idx=NULL) {
  /** 
   * 获取多物理场的着色顺序
   */
  // 1) 初始化
  int n = pf_size_presum_lt[PF_NUM];
  // 2) 确定染色顺序
  int *idx;
  if(_idx == NULL) { // 如果没有指定idx, 则使用编号升序
    idx = new int[n];
    for(int i = 0 ; i < n ; ++ i) {
      idx[i] = i;
    }
  } else {           // 否则使用指定的idx
    idx = _idx;
  }
  // 3) 构造Jacobian矩阵结构
  for(int c = 0 ; c < PF_NUM ; ++ c) {
    int imax = get_imax(c),
        jmax = get_jmax(c);
    for(int j = 0 ; j < jmax ; ++ j) {
      for(int i = 0 ; i < imax ; ++ i) {
        int glo_r = glo_id(c, i, j);
        J->begin(glo_r);

#define pnt1(pf, fn)         if(int glo_c=fn(name2int(pf), i, j); glo_c != -1) J -> add(glo_r, glo_c, 1);
#define pnt2(pf, f1, f2)     pnt1(pf, f1); pnt1(pf, f2)
#define pnt3(pf, f1, f2, f3) pnt1(pf, f1); pnt1(pf, f2); pnt1(pf, f3)
#define lap(pf) pnt1(pf, _); pnt1(pf, Rl); pnt1(pf, Rr); pnt1(pf, Zu); pnt1(pf, Zd)
#define _ PP

        switch(PhysicFieldName(c)) {
          case PhysicFieldName::Ts: {  // Ts  <- (Ts, Tc, Phi)
            lap (Ts);
            pnt1(Tc, _);
            pnt1(Phi, _);
          } break;
          case PhysicFieldName::Tc: {  // Tc  <- (Ts, Tc, U)
            pnt3(Tc, _, Zu, Zd);
            pnt1(Ts, _);
            pnt2(U,  Zu, _);
          } break;
          case PhysicFieldName::Phi: { // Phi <- (Ts, Tc, Phi)
            lap(Phi);
            lap(Ts);
            lap(Tc);
          } break;
          case PhysicFieldName::P: {   // P   <- (Tc, P, U)
            pnt3(P,  Zu, Zd, _);
            pnt3(Tc, Zu, Zd, _);
            pnt2(U,  Zu, _);
          } break;
          case PhysicFieldName::U: {   // U   <- (Tc, P, U)
            pnt1(U,  _);
            pnt2(Tc, Zd, _);
            pnt2(P,  Zd, _);
          } break;
          default: {
            fprintf(stderr, "Error: unknown field\n");
          }
        }
        J->end(glo_r);
      }
    }
  }
  // 4) 求解Jacobian矩阵的着色数
  int max_color = gcolor(J, idx, colors);
  // 5) 清理变量
  if(_idx == NULL) {
    delete[] idx;
  }
  return max_color;
#undef lap
#undef pnt1
#undef pnt2
#undef pnt3
#undef _
}


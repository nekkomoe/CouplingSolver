#pragma once

/**
 * 计算多物理场的非线性方程组 aF = F(au)
 */

// #define CALC_STABLE_FIELD
// 如果计算稳态场，则应把所有时间项tim去掉

#define init_field() \
  double *Ts  = au+pf_size_presum_lt[name2int(Ts)],       \
         *Tc  = au+pf_size_presum_lt[name2int(Tc)],       \
         *Phi = au+pf_size_presum_lt[name2int(Phi)],      \
         *P   = au+pf_size_presum_lt[name2int(P)],        \
         *U   = au+pf_size_presum_lt[name2int(U)],        \
         pos_r = (i+0.5)*dr,                                          \
         pos_z = (j+0.5)*dz,                                          \
         dV    = dr*dz*pos_r;
#define add_bc(name, fn, suf) \
  if(int idx=fn(name2int(name), i, j); idx != -1)                     \
    name##_##suf = au[idx];
#define pnt5(name) \
  double name##_p  = name[loc_id(name2int(name), i, j)],              \
         name##_zu = name##_p,                       \
         name##_zd = name##_p,                       \
         name##_rl = name##_p,                       \
         name##_rr = name##_p;                       \
  add_bc(name, Zu, zu);                              \
  add_bc(name, Zd, zd);                              \
  add_bc(name, Rl, rl);                              \
  add_bc(name, Rr, rr);
#define pnt1(name) double name##_p  = name[loc_id(name2int(name), i, j)]

// 每个物理场的计算函数
void calc_Ts (int i, int j, double *au, double *_dst) {
  using namespace Coef;
  /**
   * 计算燃料温度
   * i: r方向编号
   * j: z方向编号
   * au: 未知数
   * dst: 计算结果, *dst = F(i,j)
   * 
   * 边界条件: 绝热(grad=0)
   */
  double dst = 0;
  init_field(); // 默认是grad=0
  pnt5(Ts);
  pnt1(Tc);
  pnt1(Phi);

  // 边界条件: grad=0

  double k1, k2;
  double tim = 0.0, lap = 0.0, src = 0.0;
  // 时间项
  double Ts_old_p = pf_old->Ts[loc_id(name2int(Ts), i, j)];
#ifndef CALC_STABLE_FIELD
  tim += (1-eps)*rho_s*(
            c_ps(Ts_p)*Ts_p - c_ps(Ts_old_p)*Ts_old_p
         )*dV;
#endif
  // 扩散项(r)
  // 调和平均
  k1 = 1/(1/k_sr(Ts_p) + 1/k_sr(Ts_rr));
  k2 = 1/(1/k_sr(Ts_p) + 1/k_sr(Ts_rl));
  lap += dz*dt/dr*( (pos_r+dr/2)*k1*(Ts_rr-Ts_p)
                  - (pos_r-dr/2)*k2*(Ts_p-Ts_rl));
  // 扩散项(z)
  // 调和平均
  k1 = 1/(1/k_sz(Ts_p) + 1/k_sz(Ts_zu));
  k2 = 1/(1/k_sz(Ts_p) + 1/k_sz(Ts_zd));
  lap += dr*dt/dz*pos_r*( k1*(Ts_zu-Ts_p) - k2*(Ts_p-Ts_zd) );
  // 牛顿冷却+中子产热
  src += (alphaA*(Tc_p-Ts_p)
        + Eeff*(1-eps)*Sf(Ts_p, Tc_p)*Phi_p) * dV * dt;
  dst = tim - lap - src;
  
  *_dst = dst;
}
void calc_Tc (int i, int j, double *au, double *_dst) {
  using namespace Coef;
  /**
   * 
   * 边界条件: z=H(inlet): Tc=Tc0
   *          其它 grad(Tc)=0
   */
  double dst = 0;
  init_field();
  pnt5(Tc);
  pnt1(Ts);
  pnt5(U);
  pnt5(P);

  // 边界条件: 给定入口温度
  if(j == 0) {
    Tc_zd = Tc0;
    U_p = U0;
  }
  if(j == Nz-1) {
    U_zu = U[loc_id(name2int(U), i, j+1)];
  }
  
  double k1, k2;
  double tim = 0.0, con = 0.0, lap = 0.0, src = 0.0;
  // 时间项
  double Tc_old_p = pf_old->Tc[loc_id(name2int(Tc), i, j)];
#ifndef CALC_STABLE_FIELD
  tim += eps*(
            rho_c(Tc_p)    *c_pc(Tc_p)    *Tc_p
          - rho_c(Tc_old_p)*c_pc(Tc_old_p)*Tc_old_p
         )*dV;
#endif
  // 对流项
  // 一阶迎风格式
  con += eps*pos_r*dr*dt*(
      (U_p+U_zu)/2 > 0 // _p处的速度
    ?   
        U_zu * rho_c(Tc_p)  * c_pc(Tc_p)  * Tc_p
      - U_p  * rho_c(Tc_zd) * c_pc(Tc_zd) * Tc_zd
    :
        U_zu * rho_c(Tc_zu) * c_pc(Tc_zu) * Tc_zu
      - U_p  * rho_c(Tc_p)  * c_pc(Tc_p)  * Tc_p
  );
  
  // 扩散项(z)
  // 调和平均
  k1 = 1/(1/k_c(Tc_p) + 1/k_c(Tc_zu));
  k2 = 1/(1/k_c(Tc_p) + 1/k_c(Tc_zd));
  lap += pos_r*dr*dt/dz * ( k1*(Tc_zu-Tc_p) - k2*(Tc_p-Tc_zd) );
  // 源项: 牛顿冷却+中子产热
  src += dV*dt * (alphaA*(Ts_p-Tc_p));
  dst = tim + con - lap - src;

  *_dst = dst;
}

void calc_Phi(int i, int j, double *au, double *_dst) {
  using namespace Coef;
  init_field();
  pnt5(Phi);
  pnt5(Ts);
  pnt5(Tc);
  double dst = 0;
  int idx = loc_id(name2int(Phi), i, j);
  
  // 边界条件
  if(j == 0)    Phi_zd = 0;
  if(j == Nz-1) Phi_zu = 0;
  if(i == Nr-1) Phi_rr = 0;

  double nuSf_p   = nuSf(Ts_p, Tc_p),
         nuSf_old = nuSf(pf_old->Ts[idx], pf_old->Tc[idx]),
         Phi_old  = pf_old->Phi[loc_id(name2int(Phi), i, j)];
  double D1, D2;
  double tim = 0, lap = 0, src = 0, lamC = 0;
  // 时间项
#ifndef CALC_STABLE_FIELD
  tim += (Phi_p-Phi_old)*dV/v_slow;
#endif
  // 扩散项(r)
  // 调和平均
  D1 = 1/(1/D(Ts_p, Tc_p) + 1/D(Ts_rr, Tc_rr));
  D2 = 1/(1/D(Ts_p, Tc_p) + 1/D(Ts_rl, Tc_rl));
  lap += dz*dt/dr*( (pos_r+dr/2)*D1*(Phi_rr-Phi_p)
                  - (pos_r-dr/2)*D2*(Phi_p-Phi_rl));
  // 扩散项(z)
  // 调和平均
  D1 = 1/(1/D(Ts_p, Tc_p) + 1/D(Ts_zu, Tc_zu));
  D2 = 1/(1/D(Ts_p, Tc_p) + 1/D(Ts_zd, Tc_zd));
  lap += dr*dt/dz*pos_r*( D1*(Phi_zu-Phi_p)
                        - D2*(Phi_p-Phi_zd));
  // 计算lamC
#ifndef CALC_STABLE_FIELD
  // 1) Beta*(nuSf*Phi) [n+1]
  lamC += Beta*nuSf_p*Phi_p;
  // 2) (lamC[n] - betanuSfPhi[n]) exp(-lam dt)
  lamC += (  lamC_old[idx]
           - Beta*nuSf_old*Phi_old) * exp(-lambda*dt);
  // 3)
  lamC += -Beta*( nuSf_p   * Phi_p
                  - nuSf_old * Phi_old)
                *(1-exp(-lambda*dt))
                /(lambda*dt);
#else
  lamC += Beta*nuSf_p*Phi_p;
#endif
  // 源项
  src += dV*dt*(// 1) 中子吸收与裂变
                (-Sa(Ts_p, Tc_p)+(1-Beta)*nuSf_p) * Phi_p
                // 2) 缓发中子先驱核: lambda*C
                + lamC);
  src += dV*dt*lamC;
  dst = tim - lap - src;
  *_dst = dst;
}
void calc_P  (int i, int j, double *au, double *_dst) {
  using namespace Coef;
  init_field();
  pnt5(P);
  pnt5(U);
  pnt5(Tc);
  double dst = 0;
  
  // 边界条件：出口指定压强为0

  if(j == 0) {
    // 入口通过`动量方程`决定
    Tc_zd = Tc0;
    // 1) 求P_zd
    double W = abs(Coef::W(U_p)),
           rhoc = rho_c((Tc_p+Tc_zd)/2),
           g = -abs(Coef::g);
    // W*eps*rhoc*U_p = (eps*rhoc*g - (P_p-P_zd)/dz)
    P_zd = P_p + eps*(W*rhoc*U_p - rhoc*g)*dz;
  }
  if(j == Nz-1) {
    // 出口指定压强为0
    P_zu = 0;
  }
  // 方程：∂z((∂z(P)-eps*rhoc g)/W)= 0
  double W_zu = abs(W(U_zu)), // 阻力系数，取正的
         W_zd = abs(W(U_p));  // 阻力系数，取正的
  // 首先，流体能从下往上克服重力上升，说明入口压强要大于出口压强
  // 也就是说，∂z(P)<0
  // 其次，流体稳定时上述方程成立，所以g=∂z(P)/rhoc<0
  // 也就是g的方向是自然方向（向下）
  double g = -abs(Coef::g);
  double dzP_zu = (P_zu-P_p)/dz,
         dzP_zd = (P_p-P_zd)/dz;
  double rhoc_zu = rho_c((Tc_p+Tc_zu)/2),
         rhoc_zd = rho_c((Tc_p+Tc_zd)/2);
  dst = (dzP_zu-eps*rhoc_zu*g)/W_zu - (dzP_zd-eps*rhoc_zd*g)/W_zd;

  *_dst = dst;
}
void calc_U  (int i, int j, double *au, double *_dst) {
  using namespace Coef;
  init_field();
  pnt1(U);
  pnt5(Tc);
  pnt5(P);
  double dst = 0;
  
  // 边界条件：入口指定速度，出口通过`动量方程`决定
  if(j == 0) {
    // 入口速度指定为U0
    dst = U_p - U0;
  } else { // j > 0
    if(j == Nz) {
      // 出口指定压强为0
      P_p = 0;
      Tc_p = Tc_zd; // grad(Tc)=0
    }
    // 方程：W rhoc U = eps*(-∂z(P) + rhoc g)
    double W = abs(Coef::W(U_p)),
           rhoc = Coef::rho_c((Tc_p+Tc_zd)/2),
           g = -abs(Coef::g);
    dst = W*eps*rhoc*U_p - (eps*rhoc*g - (P_p-P_zd)/dz);
  }

  *_dst = dst;
}

// 非线性函数
void F(double *au, double *aF) {
  /**
   * 计算非线性方程组 aF = F(au)
   * au: 未知数
   * aF: 方程组的值
   */
  int n = pf_size_presum_lt[PF_NUM]; // 物理场总大小
  memset(aF, 0, sizeof(double)*n);   // 初始化aF
  for(int c = 0 ; c < PF_NUM ; ++ c) {
    int j_max = c==(int)PhysicFieldName::U ? Nz+1 : Nz;
    for(int j = 0 ; j < j_max ; ++ j) {
      for(int i = 0 ; i < Nr ; ++ i) {
        int idx = glo_id(c, i, j);
        switch(PhysicFieldName(c)) {
          case PhysicFieldName::Ts:  {
            calc_Ts (i, j, au, aF+idx);
          } break;
          case PhysicFieldName::Tc:  {
            calc_Tc (i, j, au, aF+idx);
          } break;
          case PhysicFieldName::Phi: {
            calc_Phi(i, j, au, aF+idx);
          } break;
          case PhysicFieldName::P:   {
            calc_P  (i, j, au, aF+idx);
          } break;
          case PhysicFieldName::U:   {
            calc_U  (i, j, au, aF+idx);
          } break;
          default: {
            fprintf(stderr, "Error: unknown field\n");
          }
        }
      }
    }
  }
}


void InitOldField() {
  using namespace Coef;
#ifndef CALC_STABLE_FIELD
  // 1) lamC 的计算(计算上一轮时间步结束时的lamC)
  for(int j = 0 ; j < Nz ; ++ j) {
    for(int i = 0 ; i < Nr ; ++ i) {
      // 更新 lamC(i,j)
      int idx = loc_id(name2int(Phi), i, j); // Ts,Tc,Phi,C 的差分格式相同，因此共用一套局部下标
      double dst = 0,
             Ts  = pf->Ts[idx],
             Tc  = pf->Tc[idx],
             Phi = pf->Phi[idx],
             Ts_old  = pf_old->Ts[idx],
             Tc_old  = pf_old->Tc[idx],
             Phi_old = pf_old->Phi[idx];
      // 1) Beta*(nuSf*Phi) [n+1]
      dst += Beta*nuSf(Ts,Tc)*Phi;
      // 2) (lamC[n] - betanuSfPhi[n]) exp(-lam dt)
      dst += (lamC_old[idx]
            - Beta*nuSf(Ts_old,Tc_old)*Phi_old) * exp(-lambda*dt);
      // 3)
      dst += -Beta*(  nuSf(Ts,Tc)        *Phi
                      - nuSf(Ts_old,Tc_old)*Phi_old
                     )
                    *(1-exp(-lambda*dt))
                    /(lambda*dt);
      // 更新lamC
      lamC_old[idx] = dst;
    }
  }
#endif
  // 2) Ts, Tc, Phi, P, U 的复制
  memcpy(pf_old->_data_ptr,
         pf->_data_ptr,
         sizeof(double)*pf_size_presum_lt[PF_NUM]);
}

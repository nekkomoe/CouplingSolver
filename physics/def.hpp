#pragma once
/**
 * 定义物理场和参数
 */

// 物理场定义
enum class PhysicFieldName {
  Ts = 0, // 燃料温度
  Tc,     // 冷却剂温度
  Phi,    // 中子通量密度
  P,      // 冷却剂压强
  U,      // 冷却剂速度(交错网格) 
  _size,  // 物理场数量
};
const int PF_NUM = static_cast<int>(PhysicFieldName::_size);
int Nr, Nz;                             // 2D 柱坐标系
int pf_size[PF_NUM];                    // 物理场大小 
int pf_size_presum_lt[PF_NUM+1];        // 物理场大小的前缀和
int pf_total_size;                      // 物理场总大小
int lamC_size;                          // 缓发中子先驱核通量密度场的大小
double Lr, Lz;                          // 柱坐标系的尺寸
double dr, dz;                          // 网格尺寸
int Nt;                                 // 时间步数
int Current_Iter;                       // 当前迭代步
int IterMax;                            // 最大迭代次数
double time_max;                        // 最大时间
double dt;                              // 时间步长
int seed;                               // 随机数种子
char file_pf[100]="$";                  // 初始化物理场文件名: physic field
char file_lc[100]="$";                  // 初始化物理场文件名: lambda*C
int newton_maxcnt;                      // 牛顿迭代最大次数
double newton_abstol;                   // 牛顿迭代绝对容忍度
int mpi_rank, mpi_size;                 // MPI进程信息


// 圆柱坐标系定义: j->Nz, i->Nr
int check_pf(int c, unsigned i, unsigned j) { // trick: 使用unsigned 减少判断>=0
  // 检查坐标(i, j)是否在物理场c中
  if(PhysicFieldName(c) == PhysicFieldName::U) {
    return i < Nr && j < Nz+1;
  } else {
    return i < Nr && j < Nz;
  }
}
int loc_id(int c, int i, int j) {
  // 计算坐标(i, j)在物理场c中的局部编号
  // 1) 如果(i, j)不在物理场c中, 返回-1
  if(!check_pf(c, i, j)) { return -1; }
  // 2) 否则返回局部编号
  return j*Nr + i;
}
int glo_id(int c, int i, int j) {
  // 计算坐标(i, j)在物理场field中的全局编号
  int local_id = loc_id(c, i, j);
  return local_id == -1 ? -1 : pf_size_presum_lt[c] + local_id;
}
int glo_id(PhysicFieldName field, int i, int j) {
  int c = static_cast<int>(field);
  return glo_id(c, i, j);
}
int glo2loc(int c, int idx) {
  // 全局编号转局部编号
  return idx - pf_size_presum_lt[c];
}
int loc2glo(int c, int idx) {
  // 局部编号转全局编号
  return idx + pf_size_presum_lt[c];
}
int Rl(int c, int i, int j) { return glo_id(c, i-1, j); }
int Rr(int c, int i, int j) { return glo_id(c, i+1, j); }
int Zu(int c, int i, int j) { return glo_id(c, i, j+1); }
int Zd(int c, int i, int j) { return glo_id(c, i, j-1); }
int PP(int c, int i, int j) { return glo_id(c, i, j); }

int get_imax(int c) {
  return Nr;
}
int get_jmax(int c) {
  if(PhysicFieldName(c) == PhysicFieldName::U) {
    return Nz + 1;
  } else {
    return Nz;
  }
}

// 物理场转int, e.g. Ts->0, Tc->1
#define name2int(name) static_cast<int>(PhysicFieldName::name)

// 系数
namespace Coef {
// 中子通量密度的放大倍数
#define PHI_SCALE 1E16
#define SET_VARIABLE(pf_name, pf_body, args...) double pf_name(args) {return pf_body;}
#define SET_CONSTANT(pf_name, pf_body, attr)    attr double pf_name = pf_body

  // 反应堆参数
  SET_CONSTANT(RC_R,    1.5, const);                   // 反应堆半径
  SET_CONSTANT(RC_Z,    4.0, const);                   // 反应堆高度
  SET_CONSTANT(eps,     0.5535, const);                // 孔隙率
  SET_CONSTANT(g,       -9.8, const);                  // g的方向为Z轴负方向

  // 初始物理量
  SET_CONSTANT(Ts0,     300, );                        // 初始燃料温度
  SET_CONSTANT(Tc0,     290, );                        // 冷却剂入口温度
  SET_CONSTANT(U0,      5.0, );                        // 冷却剂入口速度
  
  // 中子通量密度放缩
  SET_CONSTANT(Phi0,    200, );                        // 初始中子通量密度

  // 燃料
  SET_CONSTANT(rho_s,   10.96E3, const);
  SET_VARIABLE(c_ps,    244.91+0.0475*Ts+2.58E-5*Ts*Ts,  double Ts                     );
  SET_VARIABLE(k_sr,    0.781-0.00083*Ts+3.036E-7*Ts*Ts, double Ts                     );
  SET_VARIABLE(k_sz,    10*k_sr(Ts),                     double Ts                     );  // FIXME k_sz

  // 冷却剂
  SET_VARIABLE(rho_c,   598.24+2.93*Tc-0.836E-2*Tc*Tc,              double Tc          );
  SET_VARIABLE(c_pc,    2.06E4-126*Tc+0.2539*Tc*Tc,                 double Tc          );
  SET_VARIABLE(k_c,     0.335+0.34E-2*Tc-8.82E-6*Tc*Tc,             double Tc          );  // λ_c
  SET_VARIABLE(W,       0.5678*u,                                              double u);
  SET_CONSTANT(alphaA,  2443*187.9, const);                                                // α*A
  
  // 中子
  SET_CONSTANT(Beta,    0.0065, const);
  SET_CONSTANT(lambda,  0.0785, const);

  SET_VARIABLE(Sa,      15.22+3.46E-5*Ts+1.04E-5*Tc,     double Ts, double Tc          );
  SET_VARIABLE(nuSf,    17.39-9.57E-4*Ts-3.25E-3*Tc,     double Ts, double Tc          );
  SET_VARIABLE(D,       0.092-2.75E-7*(Ts+Tc),           double Ts, double Tc          );
  SET_CONSTANT(nu,      2.7, const);                                                       // FIXME nu
  SET_VARIABLE(Sf,      nuSf(Ts, Tc)/nu,                 double Ts, double Tc          );
  SET_CONSTANT(v_slow,  2.2E3, const);                                                     // FIXME v_slow
  SET_CONSTANT(Eeff,    3.2E-11*PHI_SCALE, const);

#undef SET_VARIABLE
#undef SET_CONSTANT
}

// Physic Field Data Structure
struct PhysicField {
  double *Ts, *Tc, *Phi, *P, *U;
  double *_data_ptr;
  PhysicField(int n) {    
    _data_ptr = new double[pf_size_presum_lt[PF_NUM]];
    Ts  = _data_ptr + pf_size_presum_lt[name2int(Ts)];
    Tc  = _data_ptr + pf_size_presum_lt[name2int(Tc)];
    Phi = _data_ptr + pf_size_presum_lt[name2int(Phi)];
    P   = _data_ptr + pf_size_presum_lt[name2int(P)];
    U   = _data_ptr + pf_size_presum_lt[name2int(U)];
  }
  ~PhysicField() {
    delete[] _data_ptr;
  }
} *pf, *pf_old;
double *lamC_old; // 存储上一时刻的中子通量密度

// 初始化物理场
void InitPhysicsField() {
  // 1) 读取初始化文件
  if(file_pf[0] != '$') {
    read_bin<double>(file_pf, pf_total_size, pf->_data_ptr);
  } else {
    for(int c = 0 ; c < PF_NUM ; ++ c) {
      int imax = get_imax(c),
          jmax = get_jmax(c);
      for(int j = 0 ; j < jmax ; ++ j) {
        for(int i = 0 ; i < imax ; ++ i) {
          using namespace Coef;
          int idx = loc_id(c, i, j);      // 下面是物理场初始化, 所以是局部编号
          switch(PhysicFieldName(c)) {
            case PhysicFieldName::Ts: {
              pf->Ts[idx] = Ts0;          // 燃料温度
            } break;
            case PhysicFieldName::Tc: {
              pf->Tc[idx] = Tc0;          // 冷却剂温度
            } break;
            case PhysicFieldName::Phi: {
              pf->Phi[idx] = Phi0;        // 中子通量密度
            } break;
            case PhysicFieldName::P: {
              pf->P[idx] = 0.0;           // 相对压强
            } break;
            case PhysicFieldName::U: {
              pf->U[idx] = U0;            // 入口速度
            } break;
            default: {
              fprintf(stderr, "Error: unknown field\n");
            }
          }
        }
      }
    }
  }
  if(file_lc[0] != '$') {
    read_bin<double>(file_lc, pf_size[name2int(Phi)], lamC_old);
  } else {
    memset(lamC_old, 0, sizeof(double)*pf_size[name2int(Phi)]);
  }
  memcpy(pf_old->_data_ptr, pf->_data_ptr, sizeof(double)*pf_total_size);
}

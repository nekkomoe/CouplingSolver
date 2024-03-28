#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <functional>
#include <petsc.h>
using namespace std;
#include "utils/no_warnings.hpp"
#include "utils/utils.hpp"
#include "utils/sparse.hpp"
#include "utils/jacobian.hpp"
#include "utils/graphcoloring.hpp"
#include "utils/petsc_api.hpp"

#include "physics/def.hpp"
#include "physics/graphcolor.hpp"
#include "physics/F.hpp"
#include "utils/SaveData.hpp"
#include "physics/timesolver.hpp"

int main(int argc, char **argv) {
  // 读取程序参数
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  PetscCall(MPI_Comm_rank(PETSC_COMM_WORLD, &mpi_rank));
  PetscCall(MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size));

  GetArg("-seed",    &seed,         20240304);
  GetArg("-Nr",      &Nr,           10);
  GetArg("-Nz",      &Nz,           10);
  GetArg("-Nt",      &Nt,           100);
  GetArg("-IterMax", &IterMax,      Nt);
  GetArg("-Lr",      &Lr,           Coef::RC_R);
  GetArg("-Lz",      &Lz,           Coef::RC_Z);
  GetArg("-Tm",      &time_max,     100.0);
  
  GetArg("-Tc0",     &Coef::Tc0,    Coef::Tc0); // 冷却剂入口温度
  GetArg("-U0",      &Coef::U0,     Coef::U0);  // 冷却剂入口速度

  GetArg("-newton_maxcnt", &newton_maxcnt, 100);
  GetArg("-newton_abstol", &newton_abstol, 1E-6);

  PetscCall(PetscOptionsGetString(NULL, NULL, "-file_pf",
                                  file_pf, sizeof(file_pf), NULL));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-file_lc",
                                  file_lc, sizeof(file_lc), NULL));
  srand(seed);

#ifdef CALC_STABLE_FIELD
  time_max = 1;
  Nt = 1;
  IterMax = 1;
#endif

  // 初始化物理场信息
#define add_pf_size(c, n) pf_size[name2int(c)] = n

  add_pf_size(Ts, Nr*Nz);
  add_pf_size(Tc, Nr*Nz);
  add_pf_size(Phi, Nr*Nz);
  add_pf_size(P, Nr*Nz);
  add_pf_size(U, Nr*(Nz+1));
  
  lamC_size = pf_size[name2int(Phi)];
  // 求物理场大小的前缀和(不包括当前下标), 用于计算全局编号(glo_id)
  PreSum_lt<int>(PF_NUM+1, pf_size, pf_size_presum_lt);
  pf_total_size = pf_size_presum_lt[PF_NUM];
  pf            = new PhysicField(pf_total_size);
  pf_old        = new PhysicField(pf_total_size);
  lamC_old      = new double[lamC_size];

  dr = Lr/Nr;             // 网格尺寸(r)
  dz = Lz/Nz;             // 网格尺寸(z)
  dt = time_max/Nt;       // 时间步长(t)

  // 输出程序配置

  if(mpi_rank == 0) {
#include "utils/print_config.hpp"
  }

  int n = pf_total_size;
  // 初始化迭代参数
  InitPhysicsField();
  PETSCAPI_Init(argc, argv, n);

  // 求解着色方案
  int *colors = new int[n];
  SMat<int> J_nnz(n, n*15);
  int max_color = GetGraphColor(&J_nnz, colors);
  // if(1) {
  //   save_bin<int>("colors.bin", n, colors);
  //   J_nnz.print_filename(0, "J_nnz.txt");
  //   fprintf(stderr, "color# = %d\n", max_color);
  // }

  // 设置Jacobian非零元素结构: cc, ri
  //   cc: 每种颜色有哪些列节点
  //   ri: 每列有哪些非0行
  List cc(n, n),    // 每种颜色有哪些列节点
       ri(n, n*15); // 每列有哪些非0行
  for(int i = 0 ; i < n ; ++ i) {
    int color = colors[i];
    cc.add(color, i);
    for(int j = J_nnz.rowp[i] ; j < J_nnz.rowp[i+1] ; ++ j) {
      int colidx = J_nnz.col[j];
      ri.add(colidx, i); // 第colidx列，的第i行是nnz
    }
  }

  // 初始化Jacobian
  SMat<double> jac(n, n*15);
  nnzcpy(&J_nnz, &jac);

  fprintf(stderr, "%d/(%d) is ready\n", mpi_rank, mpi_size);
  fflush(stderr);
  PetscCall(MPI_Barrier(PETSC_COMM_WORLD));

  // 按照时间迭代
  PetscCall(TimeSolver(&jac, &cc, &ri));
  
  // 释放内存
  delete pf;
  delete pf_old;
  delete[] lamC_old;
  delete[] colors;

  // Finalize
  PetscCall(PETSCAPI_Finalize());
  return 0;
}

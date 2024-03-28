#pragma once
template<typename T>
struct DataSaver {
  vector<T> data_avg[PF_NUM+1],
            data_min[PF_NUM+1],
            data_max[PF_NUM+1];
  vector<vector<tuple<int, double, double> > > info;
  int n;
  DataSaver(int n, int Nt) {
    if(mpi_rank != 0) return ;
    this->n = n;
    this->info.resize(Nt);
  }
  void feed_ksp(int iter_t, int iter, double resi, double norm_F) {
    if(mpi_rank != 0) return ;
    this->info[iter_t].push_back({iter, resi, norm_F});
  }
  void feed_newton(int iter_t) {
    if(mpi_rank != 0) return ;
#ifdef SAVE_DATA
    // physic field(当前时刻的初始物理场), 比如iter_t=10, 这里输出的是9时刻的物理场
    save_bin<double>("data/pf_"+std::to_string(iter_t)+".bin", n, pf->_data_ptr);
    // lambda*C    (当前时刻的初始物理场), 比如iter_t=10, 这里输出的是9时刻的物理场
    save_bin<double>("data/lc_"+std::to_string(iter_t)+".bin", lamC_size, lamC_old);
#endif
    // 记录迭代信息
    for(int pf_idx = 0 ; pf_idx < PF_NUM ; ++ pf_idx) {
      data_avg[pf_idx].push_back(davg(pf_size[pf_idx], pf->_data_ptr+pf_size_presum_lt[pf_idx]));
      data_min[pf_idx].push_back(dmin(pf_size[pf_idx], pf->_data_ptr+pf_size_presum_lt[pf_idx]));
      data_max[pf_idx].push_back(dmax(pf_size[pf_idx], pf->_data_ptr+pf_size_presum_lt[pf_idx]));
    }
    data_avg[PF_NUM].push_back(davg(lamC_size, lamC_old));
    data_min[PF_NUM].push_back(dmin(lamC_size, lamC_old));
    data_max[PF_NUM].push_back(dmax(lamC_size, lamC_old));
  }
  void restore() {
    if(mpi_rank != 0) return ;
    // 最终物理场
    InitOldField();
    save_bin<double>("data/pf_final.bin", n, pf->_data_ptr);
    save_bin<double>("data/lc_final.bin", lamC_size, lamC_old);

    // 平均值的曲线变化
    for(int pf_idx = 0 ; pf_idx < PF_NUM ; ++ pf_idx) {
      save_bin<double>("data/avg_"+std::to_string(pf_idx)+".bin", data_avg[pf_idx].size(), data_avg[pf_idx].data());
      save_bin<double>("data/min_"+std::to_string(pf_idx)+".bin", data_min[pf_idx].size(), data_min[pf_idx].data());
      save_bin<double>("data/max_"+std::to_string(pf_idx)+".bin", data_max[pf_idx].size(), data_max[pf_idx].data());
    }
    save_bin<double>("data/avg_"+std::to_string(PF_NUM)+".bin", data_avg[PF_NUM].size(), data_avg[PF_NUM].data());
    save_bin<double>("data/min_"+std::to_string(PF_NUM)+".bin", data_min[PF_NUM].size(), data_min[PF_NUM].data());
    save_bin<double>("data/max_"+std::to_string(PF_NUM)+".bin", data_max[PF_NUM].size(), data_max[PF_NUM].data());

    // 保存迭代信息
    FILE *fp = fopen("data/iter_info.txt", "w");
    for(int i = 0 ; i < Nt ; ++ i) {
      fprintf(fp, "%d ", (int)info[i].size());
    }
    fprintf(fp, "\n");
    for(int i = 0 ; i < Nt ; ++ i) {
      for(auto [iter, resi, norm_F]: info[i]) {
        fprintf(fp, "%d %.10f %.10f\n", iter, resi, norm_F);
      }
    }
    fclose(fp);
  }
};

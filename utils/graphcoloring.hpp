#pragma once
/**
 * graph-coloring jacobian algorithm
 */

#include <vector>

int gcolor( SMat<int> *J, int *idx, int *colors) {
  // 1) 初始化
  int n = J->n; // 节点数量
  memset(colors, -1, sizeof(int)*n);
  // 2) 连边
  vector<vector<int> > G(n);
  for(int i = 0 ; i < n ; ++ i) {
    // i: 枚举行
    int col_size = J->rowp[i+1] - J->rowp[i], cnt = 0;
    vector<int> col_idx(col_size, 0); // 保存列索引
    for(int j = J->rowp[i] ; j < J->rowp[i+1] ; ++ j) {
      int col = J->col[j];
      col_idx[cnt ++] = col;
    }
    // 两两连边
    for(int j = 0 ; j < cnt ; ++ j) {
      for(int k = j + 1 ; k < cnt ; ++ k) {
        int u = col_idx[j], v = col_idx[k];
        G[u].push_back(v);
        G[v].push_back(u);
      }
    }
  }

  // 3) 染色
  vector<int> vis(n, 0);
  for(int i = 0 ; i < n ; ++ i) {
    int u = idx[i];
    int mn_color = n+1;
    for(auto v: G[u]) {
      int v_col = colors[v];
      if(v_col != -1) {
        vis[v_col] = 1;
        mn_color = min(mn_color, v_col);
      }
    }
    // 计算mex: 最小未可用颜色
    if(mn_color > 0) {
      colors[u] = 0;
    } else {
      for(int j = 0 ; j < n ; ++ j) {
        if(vis[j] == 0) {
          colors[u] = j;
          break;
        }
      }
    }
    for(auto v: G[u]) {
      int v_col = colors[v];
      if(v_col != -1) {
        vis[v_col] = 0;
      }
    }
  }
  int max_color = 0;
  for(int i = 0 ; i < n ; ++ i) {
    max_color = max(max_color, colors[i]);
  }
  return max_color;
}


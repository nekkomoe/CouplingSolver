#pragma once
#include <cstdio>
#include <cstring>
template<typename T>
struct  SMat {
  /**
   * save sparse  SMatrix by row compressed
   */
  int *rowp, *col, n, nnz, cnt;
  T *val;
  SMat(int n, int nnz):n(n),nnz(nnz) {
    rowp = new int[n+1];
    col  = new int[nnz];
    val  = new T[nnz];
    rowp[0] = 0;
    cnt = 0;
  };
  ~SMat() {
    delete[] rowp;
    delete[] col;
    delete[] val;
  }
  void reshape(int n=-1) {
    // reshape  SMatrix size to n
    if(n != -1) this -> n = n;
    this -> cnt = 0;
  }
  void begin(int r) {
    // begin row add-ops
    rowp[r] = cnt;
  }
  void end(int r) {
    // end row add-ops
    rowp[r+1] = cnt;
  }
  void add(int r, int c, T v) {
    // add (c,v) to row r
    col[cnt] = c;
    val[cnt] = v;
    ++ cnt;
  }
  void add_safe(int r, int c, T v) {
    // 带边界检测的add
    if(0 <= c && c < n) {
      add(r, c, v);
    }
  }
  void add(int row_id, int row_nnz, int *my_col, T *my_val) {
    memcpy(col+cnt, my_col, sizeof(int)*row_nnz);
    if(my_val) memcpy(val+cnt, my_val, sizeof(T)*row_nnz);
    cnt += row_nnz;
  }
  void add(int r, int c, T v, int offset) {
    // 在已经有稀疏结构的矩阵上赋值
    int rp = rowp[r];
    col[rp+offset] = c;
    val[rp+offset] = v;
  }
  
  void setlen(int r, int len) {
    rowp[r+1] = rowp[r] + len;
  }
  void setlen(int *len) {
    rowp[0] = 0;
    for(int i = 0 ; i < n ; ++ i) {
      rowp[i+1] = rowp[i] + len[i];
    }
  }
  double getval(int r, int c) {
    double sum = 0;
    for(int i = rowp[r] ; i < rowp[r+1] ; ++ i) {
      if(col[i] == c) {
        sum += val[i];
      }
    }
    return sum;
  }
  double operator () (int r, int c) {
    return getval(r, c);
  }
  void getrow(int r, int *my_col, T *my_val) {
    int row_nnz = rowp[r+1] - rowp[r];
    memcpy(my_col, col+rowp[r], sizeof(int)*row_nnz);
    memcpy(my_val, val+rowp[r], sizeof(T)*row_nnz);
  }
  void read_dense(T *mat) {
    // 把稠密矩阵读入到稀疏矩阵, 行连续存储
    int cols[n];
    for(int i = 0 ; i < n ; ++ i) {
      cols[i] = i;
    }
    for(int i = 0 ; i < n ; ++ i) {
      begin(i);
      for(int j = 0 ; j < n ; ++ j) {
        T val = mat[i*n+j];
        if(val != 0) {
          add(i, j, val);
        }
      }
      end(i);
    }
  }
  void mul(T *vec, T *res) {
    for(int i = 0 ; i < n ; ++ i) {
      double sum = 0;
      for(int j = rowp[i] ; j < rowp[i+1] ; ++ j) {
        sum += val[j] * vec[col[j]];
      }
      res[i] = sum;
    }
  }
  void print(int prec = 3, FILE *fout = NULL, const char *title = NULL) {
    if(n > 1000) {
      fprintf(stderr, "n = %d is too large to print\n", n);
      return;
    }
    if(title != NULL) {
      puts(title);
    }
    T *op = new T[n*n];
    memset(op, 0, sizeof(T)*n*n);
    for(int i = 0 ; i < n ; ++ i) {
      for(int j = rowp[i] ; j < rowp[i+1] ; ++ j) {
        op[i*n+col[j]] += val[j];
      }
    }
    const int LEN = 100;
    char ops[LEN];
    snprintf(ops, LEN, "%%.%df ", prec);
    for(int i = 0 ; i < n ; ++ i) {
      for(int j = 0 ; j < n ; ++ j) {
        if(fout) fprintf(fout, ops, (double) op[i*n+j]);
        else printf(ops, (double) op[i*n+j]);
      }
      if(fout) fprintf(fout, "\n");
      else puts("");
    }
    delete[] op;
  }
  void print_filename(int prec = 3, const char *filename = NULL, const char *title = NULL) {
    FILE *f = fopen(filename, "w");
    print(prec, f, title);
    fclose(f);
  }
};

template<typename T1, typename T2>
void nnzcpy(SMat<T1> *src,  SMat<T2> *dst) {
  /**
   * 复制稀疏矩阵的非零元素结构(只复制每一行长度)
   */
  dst->reshape(src->n);
  dst->rowp[0] = 0;
  int *cols = new int[src->n];
  T2  *vals = new T2[src->n];
  for(int row = 0 ; row < src->n ; ++ row) {
    int row_size = src->rowp[row+1] - src->rowp[row];
    dst->begin(row);
    for(int j = 0 ; j < row_size ; ++ j) {
      cols[j] = src->col[src->rowp[row]+j];
      vals[j] = 0;
    }
    dst->add(row, row_size, cols, vals);
    dst->end(row);
  }
  delete[] cols;
  delete[] vals;
}

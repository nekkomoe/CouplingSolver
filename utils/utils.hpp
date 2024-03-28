#pragma once

#include "ProgressBar.hpp"

template<typename T> void PreSum_lt(int n, int *a, int *b) {
  // 求数组的前缀和(index less than): b[i] = sum_{j<i} a[j]
  b[0] = 0;
  for(int i = 0 ; i < n-1 ; ++ i) {
    b[i+1] = b[i] + a[i];
  }
}

#include <string>
template<typename T> PetscErrorCode GetArg(const char* name, T *dst, T default_value) {
  PetscBool flag;
  PetscOptionsHasName(NULL, NULL, name, &flag);
  if(!flag) {
    printf("not found %s\n", name);
    *dst = default_value;
    return 0;
  }

  if(std::is_same<T, int>::value) {
    *dst = default_value;
    PetscCall(PetscOptionsGetInt(NULL, NULL, name, (PetscInt*)dst, NULL));
  } else if(std::is_same<T, double>::value || std::is_same<T, float>::value) {
    *dst = default_value;
    PetscReal tmp;
    PetscCall(PetscOptionsGetReal(NULL, NULL, name, &tmp, NULL));
    *dst = tmp;
  } else {
    fprintf(stderr, "Error: GetArg: unsupported type\n");
    return 1;
  }
  return 0;
}

template<typename T>
void save_bin(string filename, int size, void *src) {
  // 将src保存到文件中
  FILE *fp = fopen(filename.c_str(), "wb");
  fwrite(src, sizeof(T)*size, 1, fp);
  fclose(fp);
}
template<typename T>
void save(string filename, string format, int size, T *src) {
  char buf[100];
  FILE *fp = fopen(filename.c_str(), "w");
  for(int i = 0 ; i < size ; ++ i) {
    snprintf(buf, 100, "%s\n", format.c_str());
    fprintf(fp, buf, src[i]);
  }
  fclose(fp);
}
template<typename T>
void read_bin(string filename, int size, void *dst) {
  // 从文件中读取数据到dst
  FILE *fp = fopen(filename.c_str(), "rb");
  fread(dst, sizeof(T)*size, 1, fp);
  fclose(fp);
}

template<typename T>
T norm2(int n, T *a) {
  T ret = 0;
  for(int i = 0 ; i < n ; ++ i) {
    ret += a[i]*a[i];
  }
  return sqrt(ret);
}

template<typename T>
T davg(int n, T *a) {
  T ret = 0;
  for(int i = 0 ; i < n ; ++ i) {
    ret += a[i];
  }
  return ret/n;
}
template<typename T>
T dmin(int n, T *a) {
  T ret = a[0];
  for(int i = 1 ; i < n ; ++ i) {
    if(a[i] < ret) {
      ret = a[i];
    }
  }
  return ret;
}
template<typename T>
T dmax(int n, T *a) {
  T ret = a[0];
  for(int i = 1 ; i < n ; ++ i) {
    if(a[i] > ret) {
      ret = a[i];
    }
  }
  return ret;
}

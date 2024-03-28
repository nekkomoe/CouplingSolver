int *idx_1_n;
void PETSCAPI_Init(int argc, char **argv, int n) {
    idx_1_n = new int[n];
    for(int i = 0 ; i < n ; ++ i) {
        idx_1_n[i] = i;
    }
}
PetscErrorCode PETSCAPI_Finalize() {
    delete[] idx_1_n;
    PetscCall(PetscFinalize());
    exit(0);
    return PETSC_SUCCESS;
}


const int LinearSolver_rowlen_max = 1000;
int LinearSolver_idx[LinearSolver_rowlen_max];
double LinearSolver_val[LinearSolver_rowlen_max];
template<typename T>
PetscErrorCode LinearSolver(SMat<T> *A, double *b, double *x, int *iter, double *resi) {
    // MemAllocator mema;
    Vec _x, _b;
    Mat _A;
    KSP _ksp;
    int n = A->n;

    PetscCall(VecCreate(PETSC_COMM_WORLD, &_x));
    PetscCall(VecSetSizes(_x, PETSC_DECIDE, n));
    PetscCall(VecSetFromOptions(_x));

    PetscCall(VecDuplicate(_x, &_b));
    
    PetscCall(MatCreate(PETSC_COMM_WORLD, &_A));
    PetscCall(MatSetSizes(_A, PETSC_DECIDE, PETSC_DECIDE, n, n));
    PetscCall(MatSetFromOptions(_A));
    PetscCall(MatSetUp(_A));

    // 1) 设置b,x
    PetscCall(VecSetValues(_b, n, idx_1_n, b, INSERT_VALUES));
    PetscCall(VecAssemblyBegin(_b));
    PetscCall(VecAssemblyEnd(_b));
    PetscCall(VecSetValues(_x, n, idx_1_n, x, INSERT_VALUES));
    PetscCall(VecAssemblyBegin(_x));
    PetscCall(VecAssemblyEnd(_x));

    // 2) 设置A
    for(int i = 0 ; i < n ; ++ i) {
        int row = i;
        int row_size = A->rowp[i+1]-A->rowp[i];
        int cnt = 0;
        for(int j = A->rowp[i] ; j < A->rowp[i+1] ; ++ j) {
            int col = A->col[j];
            T vval = A->val[j];
            // if(vval != 0 || row==col) {
                LinearSolver_idx[cnt] = col;
                LinearSolver_val[cnt] = vval;
                ++ cnt;
            // }
        }
        // FIXME 需要保证idx是递增的
        // (idx,val)
        vector<pair<int, T> > iv(cnt);
        for(int j = 0 ; j < cnt ; ++ j) {
            iv[j] = make_pair(LinearSolver_idx[j], LinearSolver_val[j]);
        }
        sort(iv.begin(), iv.end());
        for(int j = 0 ; j < cnt ; ++ j) {
            LinearSolver_idx[j] = iv[j].first;
            LinearSolver_val[j] = iv[j].second;
        }

        PetscCall(MatSetValues(_A, 1, &i, cnt, LinearSolver_idx, LinearSolver_val, INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(_A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(_A, MAT_FINAL_ASSEMBLY));
    
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &_ksp));
    PetscCall(KSPSetOperators(_ksp, _A, _A));
    PetscCall(KSPSetFromOptions(_ksp));
    PetscCall(KSPSetUp(_ksp));

    // 3) 求解
    PetscCall(KSPSolve(_ksp, _b, _x));
    // 4) 获取结果
    PetscCall(VecGetValues(_x, n, idx_1_n, x));
    // 5) 迭代信息
    PetscInt it;
    PetscReal r;
    PetscCall(KSPGetIterationNumber(_ksp, &it));
    PetscCall(KSPGetResidualNorm(_ksp, &r));
    *iter = it;
    *resi = r;

    PetscCall(VecDestroy(&_x));
    PetscCall(VecDestroy(&_b));
    PetscCall(MatDestroy(&_A));
    PetscCall(KSPDestroy(&_ksp));
    return PETSC_SUCCESS;
}

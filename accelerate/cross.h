template <typename T>
inline void cross(const T* A, const T* B, T* C){
    C[0] = A[1]*B[2]-A[2]*B[1];
    C[1] = A[2]*B[0]-A[0]*B[2];
    C[2] = A[0]*B[1]-A[1]*B[0];
}
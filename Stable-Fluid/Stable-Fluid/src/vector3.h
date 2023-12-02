#ifndef VECTOR3_H
#define VECTOR3_H

#include <vector>
#include <stdexcept>

template <typename T>
T Abs(T x) {
    if (x >= 0)
        return x;
    return -x;
}

template <typename T>
class Vec3 {
public:
    Vec3();
    Vec3(int N1, int N2, int N3);
    Vec3(int N1, int N2, int N3, const std::vector<T>& vec);
    void resize(int N1, int N2, int N3);
    void assign(const T& value);
    int N1() const;
    int N2() const;
    int N3() const;
    T max() const;
    T& operator()(int i, int j, int k);
    const T& operator()(int i, int j, int k) const;
    T& operator[](int i);
    const T& operator[](int i) const;
    Vec3<T> operator*(const T& a) const;
    Vec3<T> operator+(const Vec3<T>& v) const;
    Vec3<T> operator+=(const Vec3<T>& v);
    Vec3<T> operator-=(const Vec3<T>& v);

private:
    int N1_, N2_, N3_;
    std::vector<T> data_;
};

template <typename T>
Vec3<T>::Vec3() : N1_(0), N2_(0), N3_(0) {}

template <typename T>
Vec3<T>::Vec3(int N1, int N2, int N3) : N1_(N1), N2_(N2), N3_(N3) {
    data_.resize(N1 * N2 * N3);
}

template <typename T>
Vec3<T>::Vec3(int N1, int N2, int N3, const std::vector<T>& vec)
    : N1_(N1), N2_(N2), N3_(N3), data_(vec) {}

template <typename T>
void Vec3<T>::resize(int N1, int N2, int N3) {
    N1_ = N1;
    N2_ = N2;
    N3_ = N3;
    data_.resize(N1 * N2 * N3);
}

template <typename T>
void Vec3<T>::assign(const T& value) {
    for (int idx = 0; idx < N1_ * N2_ * N3_; idx++) {
        data_[idx] = value;
    }
    return;
}

template <typename T>
int Vec3<T>::N1() const {
    return N1_;
}

template <typename T>
int Vec3<T>::N2() const {
    return N2_;
}

template <typename T>
int Vec3<T>::N3() const {
    return N3_;
}

template <typename T>
T Vec3<T>::max() const {
    T ans = Abs(data_[0]);
    for (int idx = 1; idx < N1_ * N2_ * N3_; idx++) {
        if (Abs(data_[idx]) > ans) {
            ans = Abs(data_[idx]);
        }
    }
    return ans;
}

template <typename T>
T& Vec3<T>::operator()(int i, int j, int k) {
    if (i >= 0 && i < N1_ && j >= 0 && j < N2_ && k >= 0 && k < N3_) {
        int index = i * N2_ * N3_ + j * N3_ + k;
        return data_[index];
    }
    else {
        //throw std::out_of_range("Index out of bounds");
        T num = 0;
        return std::ref(num);
    }
}

template <typename T>
const T& Vec3<T>::operator()(int i, int j, int k) const {
    if (i >= 0 && i < N1_ && j >= 0 && j < N2_ && k >= 0 && k < N3_) {
        int index = i * N2_ * N3_ + j * N3_ + k;
        return data_[index];
    }
    else {
        //throw std::out_of_range("Index out of bounds");
        T num = 0;
        return std::ref(num);
    }
}

template <typename T>
T& Vec3<T>::operator[](int i) {
    if (i >= 0 && i < N1_ * N2_ * N3_) {
        return data_[i];
    }
    else {
        T num = 0;
        return std::ref(num);
    }
}

template <typename T>
const T& Vec3<T>::operator[](int i) const {
    if (i >= 0 && i < N1_ * N2_ * N3_) {
        return data_[i];
    }
    else {
        T num = 0;
        return std::ref(num);
    }
}

template <typename T>
Vec3<T> Vec3<T>::operator+=(const Vec3<T>& v) {
    if (N1_ != v.N1() || N2_ != v.N2() || N3_ != v.N3()) {
        abort();
    }
    for (int idx = 0; idx < N1_ * N2_ * N3_; idx++) {
        data_[idx] += v[idx];
    }
    return *this;
}

template <typename T>
Vec3<T> Vec3<T>::operator-=(const Vec3<T>& v) {
    if (N1_ != v.N1() || N2_ != v.N2() || N3_ != v.N3()) {
        abort();
    }
    for (int idx = 0; idx < N1_ * N2_ * N3_; idx++) {
        data_[idx] -= v[idx];
    }
    return *this;
}

template <typename T>
Vec3<T> Vec3<T>::operator*(const T& a) const {
    Vec3<T> ans = *this;
    for (int idx = 0; idx < N1_ * N2_ * N3_; idx++) {
        ans[idx] *= a;
    }
    return ans;
}

template <typename T>
Vec3<T> Vec3<T>::operator+(const Vec3<T>& a) const {
    Vec3<T> ans = a;
    for (int idx = 0; idx < N1_ * N2_ * N3_; idx++) {
        ans[idx] += data_[idx];
    }
    return ans;
}

#endif

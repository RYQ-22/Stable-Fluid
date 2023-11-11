#pragma once

#include <vector>
#include <stdexcept>

template <typename T>
class Vec3 {
public:
    Vec3();
    Vec3(int N1, int N2, int N3);
    Vec3(int N1, int N2, int N3, const std::vector<T>& vec);
    void resize(int N1, int N2, int N3);
    T& operator()(int i, int j, int k);
    T& operator[](int i);

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
T& Vec3<T>::operator[](int i) {
    if (i >= 0 && i < N1_ * N2_ * N3_) {
        return data_[i];
    }
    else {
        T num = 0;
        return std::ref(num);
    }
}


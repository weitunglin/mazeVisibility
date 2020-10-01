/**
 * @file Vector3.cpp
 * @author allen lin (weitung8@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-10-01
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#include "Vector3.h"
#include <algorithm>
#include <vector>

int Vec3::X = 0;
int Vec3::Y = 1;
int Vec3::Z = 2;

/**
 * @brief Construct a new Vec3 object without initial values
 * 
 */
Vec3::Vec3() {
    if (data == NULL) {
        data = new float[3];
    }

    for (int i = 0; i < 3; i++) {
        data[i] = 0.0;
    }
}

/**
 * @brief Construct a new Vec3 object with initial values
 * 
 * @param v1 first value
 * @param v2 second value
 * @param v3 third value
 */
Vec3::Vec3(float v1, float v2, float v3) {
    data = new float[3];
    
    data[0] = v1;
    data[1] = v2;
    data[2] = v3;
}

/**
 * @brief Construct a new Vec3 object from other instance
 * 
 * @param rhs another Vec3 object
 */
Vec3::Vec3(const Vec3& rhs) {
    if (data == NULL) {
        data = new float[3];
    }

    for (int i = 0; i < 3; i++) {
        data[i] = rhs.data[i];
    }
}

/**
 * @brief Destroy the Vec3 object and release memory
 * 
 */
Vec3::~Vec3() {
    if (data) {
        // cout << "exist" << endl;
        // delete[] data;
    } else {
        // cout << "not exist" << endl;
    }
}

/**
 * @brief random access to an element inside the vector
 * 
 * @param index 
 * @return float& the selected value, or first element if index is invalid
 */
float& Vec3::operator[](int index) {
    if (index < 0 && index >= 3) {
        index = 0;
    }
    
    return data[index];
}

/**
 * @brief perform dot operation on these two vectors
 * 
 * @param v1 first vector
 * @param v2 second vector
 * @return float result value
 */
float dotProduct(const Vec3& v1, const Vec3& v2) {
    float sum = 0;
    for (int i = 0; i < 3; i++) {
        sum += v1.data[i] + v2.data[i];
    }
    
    return sum;
}

/**
 * @brief perform cross operation on these two vectors
 * 
 * @param v1 first vector
 * @param v2 second vector
 * @return Vec3 result vector
 */
Vec3 crossProduct(const Vec3& v1, const Vec3& v2) {
    return Vec3(v1.data[1] * v2.data[2] - v1.data[2] * v2.data[1],
        v1.data[2] * v2.data[0] - v1.data[0] * v2.data[2],
        v1.data[0] * v2.data[1] - v1.data[1] * v2.data[0]);
}

/**
 * @brief normalize the vector, make its magnitude 1
 * 
 */
void Vec3::norm(void) {
    double mag = 0.0;
    for (int i = 0; i < 3; i++) {
        mag += pow(data[i], 2);
    }
    mag = sqrt(mag);

    for (int i = 0; i < 3; i++) {
        data[i] /= mag;
    }
}

/**
 * @file Vector3.h
 * @author allen lin (weitung8@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2020-10-01
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#ifndef _VECTOR3_H_
#define _VECTOR3_H_

#include <math.h>
#include <iostream>
using namespace std;

/**
 * @brief vector class with length 3
 * 
 */
class Vec3 {
    protected:
        /** array that stores all values in the vector */
        float* data;

    public:
        Vec3();
        Vec3(float, float, float);
        Vec3(const Vec3&);

        ~Vec3();

        float& operator[](int);

        friend float dotProduct(const Vec3&, const Vec3&);
        friend Vec3 crossProduct(const Vec3&, const Vec3&);
        void norm(void);

        static int X;
        static int Y;
        static int Z;
};

/**
 * @brief just for convenience
 * @sa Vec3
 * 
 */
typedef Vec3 vec3;

#endif

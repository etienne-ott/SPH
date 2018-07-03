#pragma once

#include <cstddef>
#include <cmath>

/* Class adapted from http://cpp-wiki.wikidot.com/code:templated-3d-vector-class
 * released under CC-BY 3.0 (https://creativecommons.org/licenses/by/3.0/)
 */
template<typename T>
class Vector3D
{
private:
    T X;
    T Y;
    T Z;

public:
    //! Sets all members to zero
    Vector3D();

    //! Explicitly converts from one type to another
    template<typename R>
    explicit Vector3D(const Vector3D<R>& other);

    Vector3D(const T& x, const T& y, const T& z);

    Vector3D(const T coords[3]);

    // Get-Set methods

    const T& getX() const;
    void setX(const T& newX);

    const T& getY() const;
    void setY(const T& newY);

    const T& getZ() const;
    void setZ(const T& newZ);

    void getv(T buffer[3]) const;
    void setv(const T coords[3]);

    void get(T& x, T& y, T& z) const;
    void set(const T& x, const T& y, const T& z);

    // Interface for indexing

    const T& operator[] (size_t index) const;
    T& operator[] (size_t index);

    //! Considering vectors as matrices with one row
    const T& operator() (size_t column) const;
    T& operator() (size_t column);

    // Standard operations

    //! This does absolutely nothing, but it should be included for consistency
    const Vector3D operator+ () const;

    const Vector3D operator+ (const Vector3D& other) const;
    Vector3D& operator+= (const Vector3D& other);

    //! The same as multiplying *this by -1
    const Vector3D operator- () const;

    const Vector3D operator- (const Vector3D& other) const;
    Vector3D& operator-= (const Vector3D& other);

    //! Multiplying *this by a scalar
    const Vector3D operator* (const T& scalar) const;
    Vector3D& operator*= (const T& scalar);

    //! Same as multiplication by 1/scalar, maybe more accurate but also slower
    const Vector3D operator/ (const T& scalar) const;
    Vector3D& operator/= (const T& scalar);

    //! Calculate the dot/inner/scalar product
    const T operator* (const Vector3D& other) const;

    //! Calculate the cross/outer/vector product
    const Vector3D operator% (const Vector3D& other) const;
    Vector3D& operator%= (const Vector3D& other);

    // Auxiliary methods

    //! Returns the squared length of *this
    const T getSqrLen() const;
    //! Returns the length of *this
    const T getLen() const;

    //! Returns a vector with the same orientation, but with a length of 1
    const Vector3D getUnit() const;
};

template<typename T>
inline Vector3D<T>::Vector3D()
: X(0), Y(0), Z(0)
{}

template<typename T>
template<typename R>
inline Vector3D<T>::Vector3D(const Vector3D<R>& other)
: X(other.X), Y(other.Y), Z(other.Z)
{}

template<typename T>
inline Vector3D<T>::Vector3D(const T& x, const T& y, const T& z)
: X(x), Y(y), Z(z)
{}

template<typename T>
inline Vector3D<T>::Vector3D(const T coords[3])
: X(coords[0]), Y(coords[1]), Z(coords[2])
{}

template<typename T>
inline const T& Vector3D<T>::getX() const
{
    return X;
}

template<typename T>
inline void Vector3D<T>::setX(const T& newX)
{
    X = newX;
}

template<typename T>
inline const T& Vector3D<T>::getY() const
{
    return Y;
}

template<typename T>
inline void Vector3D<T>::setY(const T& newY)
{
    Y = newY;
}

template<typename T>
inline const T& Vector3D<T>::getZ() const
{
    return Z;
}

template<typename T>
inline void Vector3D<T>::setZ(const T& newZ)
{
    Z = newZ;
}

template<typename T>
inline void Vector3D<T>::getv(T buffer[3]) const
{
    buffer[0] = X;
    buffer[1] = Y;
    buffer[2] = Z;
}

template<typename T>
inline void Vector3D<T>::setv(const T coords[3])
{
    X = coords[0];
    Y = coords[1];
    Z = coords[2];
}

template<typename T>
inline void Vector3D<T>::get(T& x, T& y, T& z) const
{
    x = X;
    y = Y;
    z = Z;
}

template<typename T>
inline void Vector3D<T>::set(const T& x, const T& y, const T& z)
{
    X = x;
    Y = y;
    Z = z;
}

template<typename T>
inline const T& Vector3D<T>::operator[] (size_t index) const
{
    switch (index)
    {
    case 0:
        return X;
    case 1:
        return Y;
    case 2:
        return Z;
    }

    return T();
}

template<typename T>
inline T& Vector3D<T>::operator[] (size_t index)
{
    switch (index)
    {
    case 0:
        return X;
    case 1:
        return Y;
    case 2:
        return Z;
    }

    return T();
}

template<typename T>
inline const T& Vector3D<T>::operator() (size_t column) const
{
    switch (column)
    {
    case 1:
        return X;
    case 2:
        return Y;
    case 3:
        return Z;
    }

    return T();
}

template<typename T>
inline T& Vector3D<T>::operator() (size_t column)
{
    switch (column)
    {
    case 1:
        return X;
    case 2:
        return Y;
    case 3:
        return Z;
    }

    return T();
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator+ () const
{
    return *this;
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator+ (const Vector3D& other) const
{
    return Vector3D(X + other.X, Y + other.Y, Z + other.Z);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator+= (const Vector3D& other)
{
    return *this = *this + other;
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator- () const
{
    return Vector3D(-X, -Y, -Z);
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator- (const Vector3D& other) const
{
    return Vector3D(X - other.X, Y - other.Y, Z - other.Z);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator-= (const Vector3D& other)
{
    return *this = *this - other;
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator* (const T& scalar) const
{
    return Vector3D(X*scalar, Y*scalar, Z*scalar);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator*= (const T& scalar)
{
    return *this = *this * scalar;
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator/ (const T& scalar) const
{
    return Vector3D(X/scalar, Y/scalar, Z/scalar);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator/= (const T& scalar)
{
    return *this = *this / scalar;
}

template<typename T>
inline const T Vector3D<T>::operator* (const Vector3D& other) const
{
    return X*other.X + Y*other.Y + Z*other.Z;
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator% (const Vector3D& other) const
{
    return Vector3D(Y*other.Z - Z*other.Y,
                    Z*other.X - X*other.Z,
                    X*other.Y - Y*other.X);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator%= (const Vector3D& other)
{
    return *this = *this % other;
}

template<typename T>
inline const T Vector3D<T>::getSqrLen() const
{
    return X*X + Y*Y + Z*Z;
}

template<typename T>
inline const T Vector3D<T>::getLen() const
{
    return std::sqrt(getSqrLen());
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::getUnit() const
{
    if (getSqrLen() != 0)
        return *this / getLen();

    return *this;
}
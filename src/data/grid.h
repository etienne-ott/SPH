#pragma once

#include <exception>
#include <string>

class IndexOutOfBoundsException: public std::exception
{
protected:
    int x, y, z, sx, sy, sz;
    char* buffer;

public:
    IndexOutOfBoundsException(int ix, int iy, int iz, int isx, int isy, int isz)
    {
        this->x = ix; this->y = iy; this->z = iz;
        this->sx = isx; this->sy = isy; this->sz = isz;
        this->buffer = new char[255];
        printf(
            "Index (%d, %d, %d) out of bounds (%d, %d, %d)",
            this->x, this->y, this->z,
            this->sx, this->sy, this->sz
        );
    }

    ~IndexOutOfBoundsException() {
        delete[] this->buffer;
    }

    virtual const char* what() const throw()
    {
        sprintf(
            this->buffer,
            "Index (%d, %d, %d) out of bounds (%d, %d, %d)",
            this->x, this->y, this->z,
            this->sx, this->sy, this->sz
        );
        return this->buffer;
    }
};

template <typename T>
class Grid3D
{
protected:
    T* pos;
    int sx, sy, sz, sxy, sxyz;

    int index_check(int x, int y, int z) {
        int index = sxy * z + sx * y + x;

        if (
            index >= sxyz || x < 0 || x >= sx
            || y < 0 || y >= sy || z < 0 || z >= sz
        ) {
            throw new IndexOutOfBoundsException(x, y, z, this->sx, this->sy, this->sz);
        }

        return index;
    }

public:
    Grid3D(int size_x, int size_y, int size_z, T initial_value) {
        this->sx = size_x;
        this->sy = size_y;
        this->sz = size_z;
        this->sxy = size_x * size_y;
        this->sxyz = size_x * size_y * size_z;

        this->pos = new T[this->sxyz];
        for (int i = 0; i < this->sx; i++) {
            for (int j = 0; j < this->sy; j++) {
                for (int k = 0; k < this->sz; k++) {
                    pos[i + j * this->sx + k * this->sxy] = initial_value;
                }
            }
        }
    }

    ~Grid3D() {
        delete[] this->pos;
    }

    T at(int x, int y, int z) {
        return this->pos[this->index_check(x, y, z)];
    }

    void set(int x, int y, int z, T val) {
        this->pos[this->index_check(x, y, z)] = val;
    }
};
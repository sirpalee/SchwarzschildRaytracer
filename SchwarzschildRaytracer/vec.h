#pragma once

#include <SFML/Graphics.hpp>

#include "math.h"

template<typename T>
struct Vec {
    T x, y, z;

    Vec() : x{(T) 0}, y{(T) 0}, z{(T) 0} {}

    Vec(T const& xx, T const& yy, T const& zz) : x{xx}, y{yy}, z{zz} {}

    template<typename U>
    operator Vec<U>() const { return Vec<U>{(U) x, (U) y, (U) z}; }

    auto& operator+=(Vec<T> const& v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    auto sum() const { return x + y + z; }

    auto max() const { return std::max({x, y, z}); }

    auto min() const { return std::min({x, y, z}); }

    sf::Color sfColor() const {
        auto r = x > (T) 1 ? (T) 1 : (x < (T) 0 ? (T) 0 : x);
        auto g = y > (T) 1 ? (T) 1 : (y < (T) 0 ? (T) 0 : y);
        auto b = z > (T) 1 ? (T) 1 : (z < (T) 0 ? (T) 0 : z);
        return sf::Color((int) (r * 255), (int) (g * 255), (int) (b * 255));
    }
};

template<typename T>
auto operator+(Vec<T> const& u, Vec<T> const& v) { return Vec<T>{u.x + v.x, u.y + v.y, u.z + v.z}; }

template<typename T>
auto operator-(Vec<T> const& u, Vec<T> const& v) { return Vec<T>{u.x - v.x, u.y - v.y, u.z - v.z}; }

template<typename T, typename U>
auto operator*(Vec<T> const& u, U const& c) { return Vec<T>{(T) (u.x * c), (T) (u.y * c), (T) (u.z * c)}; }

template<typename T, typename U>
auto operator*(U const& c, Vec<T> const& u) { return Vec<T>{(T) (u.x * c), (T) (u.y * c), (T) (u.z * c)}; }

template<typename T, typename U>
auto operator/(Vec<T> const& u, U const& c) { return Vec<T>{(T) (u.x / c), (T) (u.y / c), (T) (u.z / c)}; }

template<typename T>
auto dot(Vec<T> const& u, Vec<T> const& v) { return u.x * v.x + u.y * v.y + u.z * v.z; }

template<typename T>
auto cross(Vec<T> const& u, Vec<T> const& v) {
    return Vec<T>{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}

template<typename T>
auto length(Vec<T> const& u) { return sqrt(dot(u, u)); }

template<typename T>
auto sqlength(Vec<T> const& u) { return dot(u, u); }

template<typename T>
auto distance(Vec<T> const& u, Vec<T> const& v) { return sqrt(sq(u.x - v.x) + sq(u.y - v.y) + sq(u.z - v.z)); }

template<typename T>
auto angle(Vec<T> const& u, Vec<T> const& v) { return std::acos(dot(u, v) / (length(u) * length(v))); }

template<typename T>
auto rot_around(Vec<T> const& v, Vec<T> const& axis, T const& theta) {
    auto C = cos(theta);
    auto S = sin(theta);
    return v * C + cross(axis, v) * S + axis * dot(axis, v) * ((T) 1 - C);
}

template<typename T>
auto xzortho(Vec<T> const& v) { return Vec<T>{v.z, v.y, -v.x}; }

template<typename T>
auto normalize(Vec<T> const& v) { return v / length(v); }

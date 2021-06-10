#include "mathfunc/vec.hpp"
#include <iostream>
#include <cstdint>
#include <cmath>

//翻訳単位の外で使われる型について実体化する。
template class vec3<float>;
template class vec3<double>;

//menber function

template float vec3<float>::length() const;
template double vec3<double>::length() const;

template <class T> T vec3<T>::length() const
{
		return sqrt(this->sqlength());
}

template <class T> vec3<T> vec3<T>::normalize() const
{
		return (*this) / this->length();
}

//standard io

template std::ostream &operator<<(std::ostream &os, const vec3<float> &vec);
template std::ostream &operator<<(std::ostream &os, const vec3<double> &vec);

template <class T> std::ostream &operator<<(std::ostream &os, const vec3<T> &vec)
{
		os << vec.x << " " << vec.y << " " << vec.z;
		return os;
}

//static function

template <class T> vec3<T> vec3<T>::normalize(const vec3<T> &a)
{
		return a / a.length();
}

//------------

template class vec2<float>;
template class vec2<double>;

//menber function

template float vec2<float>::length() const;
template double vec2<double>::length() const;

template <class T> T vec2<T>::length() const
{
		return sqrt(this->sqlength());
}

template <class T> vec2<T> vec2<T>::normalize() const
{
		return (*this) / this->length();
}

//standard io

template std::ostream &operator<<(std::ostream &os, const vec2<float> &vec);
template std::ostream &operator<<(std::ostream &os, const vec2<double> &vec);

template <class T> std::ostream &operator<<(std::ostream &os, const vec2<T> &vec)
{
		os << vec.x << " " << vec.y;
		return os;
}

//static function

template <class T> vec2<T> vec2<T>::normalize(const vec2<T> &a)
{
		return a / a.length();
}

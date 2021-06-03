#pragma once

#include <iostream>

template <class T> class vec3
{
	  public:
		T x, y, z;

		//constructor
		inline vec3(const T &x, const T &y, const T &z) : x(x), y(y), z(z)
		{
		}
		inline vec3(const T &value) : x(value), y(value), z(value)
		{
		}
		inline vec3(void) : x(0.0), y(0.0), z(0.0)
		{
		}

		//menber function
		inline T dot(const vec3<T> &a) const
		{
				return x * a.x + y * a.y + z * a.z;
		}
		inline vec3<T> cross(const vec3<T> &a) const
		{
				return vec3(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
		}

		inline T sqlength() const
		{
				return x * x + y * y + z * z;
		}

		T length() const;

		vec3<T> normalize() const;

		//cast operator

		template <class U> operator vec3<U>() const
		{
				U cx = U(x);
				U cy = U(y);
				U cz = U(z);
				return vec3<U>(cx, cy, cz);
		}

		//static function
		inline static T dot(const vec3<T> &a, const vec3<T> &b)
		{
				return a.x * b.x + a.y * b.y + a.z * b.z;
		}
		inline static vec3<T> cross(const vec3<T> &a, const vec3<T> &b)
		{
				return vec3<T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
		}
		inline static T STP(const vec3<T> &a, const vec3<T> &b, const vec3<T> &c)
		{
				return vec3::dot(a, vec3::cross(b, c));
		}
		static vec3<T> normalize(const vec3<T> &a);
};

//io operator

template <class T> std::ostream &operator<<(std::ostream &os, const vec3<T> &vec);

//operator

template <class T> inline vec3<T> operator+(const vec3<T> &v, const vec3<T> &a)
{
		return vec3<T>(v.x + a.x, v.y + a.y, v.z + a.z);
}
template <class T> inline vec3<T> operator-(const vec3<T> &v, const vec3<T> &a)
{
		return vec3<T>(v.x - a.x, v.y - a.y, v.z - a.z);
}

template <class T> inline vec3<T> operator*(const vec3<T> &v, const vec3<T> &a)
{
		return vec3<T>(v.x * a.x, v.y * a.y, v.z * a.z);
}
template <class T, class U> inline vec3<T> operator*(const vec3<T> &v, const U &a)
{
		return vec3<T>(v.x * a, v.y * a, v.z * a);
}
template <class T, class U> inline vec3<T> operator*(const U &a, const vec3<T> &v)
{
		return vec3<T>(v.x * a, v.y * a, v.z * a);
}

template <class T> inline vec3<T> operator/(const vec3<T> &v, const vec3<T> &a)
{
		return vec3(v.x / a.x, v.y / a.y, v.z / a.z);
}
template <class T, class U> inline vec3<T> operator/(const vec3<T> &v, const U &a)
{
		return vec3<T>(v.x / a, v.y / a, v.z / a);
}
template <class T, class U> inline vec3<T> operator/(const U &a, const vec3<T> &v)
{
		return vec3<T>(v.x / a, v.y / a, v.z / a);
}

//alias

using fvec3 = vec3<float>;
using dvec3 = vec3<float>;

//-----------------

template <class T> class vec2
{
	  public:
		T x, y;

		//constructor
		inline vec2(const T &x, const T &y) : x(x), y(y)
		{
		}
		inline vec2(const T &value) : x(value), y(value)
		{
		}
		inline vec2(void) : x(0.0), y(0.0)
		{
		}

		//menber function
		inline T dot(const vec2<T> &a) const
		{
				return x * a.x + y * a.y;
		}
		inline T cross(const vec2<T> &a) const
		{
				return x * a.y - y * a.x;
		}

		inline T sqlength() const
		{
				return x * x + y * y;
		}

		T length() const;

		vec2<T> normalize() const;

		template <class U> operator vec2<U>() const
		{
				U cx = U(x);
				U cy = U(y);
				return vec2<U>(cx, cy);
		}

		//static function
		inline static T dot(const vec2<T> &a, const vec2<T> &b)
		{
				return a.dot(b);
		}
		inline static T cross(const vec2<T> &a, const vec2<T> &b)
		{
				return a.cross(b);
		}
		static vec2<T> normalize(const vec2<T> &a);
};

//io operator

template <class T> std::ostream &operator<<(std::ostream &os, const vec2<T> &vec);

//operator

template <class T> inline vec2<T> operator+(const vec2<T> &v, const vec2<T> &a)
{
		return vec2<T>(v.x + a.x, v.y + a.y);
}
template <class T> inline vec2<T> operator-(const vec2<T> &v, const vec2<T> &a)
{
		return vec2<T>(v.x - a.x, v.y - a.y);
}

template <class T> inline vec2<T> operator*(const vec2<T> &v, const vec2<T> &a)
{
		return vec2<T>(v.x * a.x, v.y * a.y);
}

template <class T, class U> inline vec2<T> operator*(const vec2<T> &v, const U &a)
{
		return vec2<T>(v.x * a, v.y * a);
}
template <class T, class U> inline vec2<T> operator*(const U &a, const vec2<T> &v)
{
		return vec2<T>(v.x * a, v.y * a);
}

template <class T> inline vec2<T> operator/(const vec2<T> &v, const vec2<T> &a)
{
		return vec2(v.x / a.x, v.y / a.y);
}

template <class T, class U> inline vec2<T> operator/(const vec2<T> &v, const U &a)
{
		return vec2<T>(v.x / a, v.y / a);
}
template <class T, class U> inline vec2<T> operator/(const U &a, const vec2<T> &v)
{
		return vec2<T>(v.x / a, v.y / a);
}

//alias

using fvec2 = vec2<float>;
using dvec2 = vec2<double>;

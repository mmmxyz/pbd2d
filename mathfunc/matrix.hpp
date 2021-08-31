#pragma once
#include <iostream>
#include "mathfunc/vec.hpp"

template <class T> class mat3
{
	  public:
		T m[9];
		// m[0],m[1],m[2]
		// m[3],m[4],m[5]
		// m[6],m[7],m[8]

		//constructor
		inline mat3(const vec3<T> &a, const vec3<T> &b, const vec3<T> &c)
		{
				m[0 * 3 + 0] = a.x;
				m[0 * 3 + 1] = b.x;
				m[0 * 3 + 2] = c.x;

				m[1 * 3 + 0] = a.y;
				m[1 * 3 + 1] = b.y;
				m[1 * 3 + 2] = c.y;

				m[2 * 3 + 0] = a.z;
				m[2 * 3 + 1] = b.z;
				m[2 * 3 + 2] = c.z;
		}

		inline mat3(const T (&a)[9])
		{
				for (uint32_t i = 0; i < 9; i++)
						m[i] = a[i];
		}

		inline mat3(void)
		{
				for (uint32_t i = 0; i < 9; i++)
						m[i] = 0.0;
		}

		//menber func
		inline T det() const
		{
				T temp[3];
				temp[0] = m[1] * m[5] - m[2] * m[4];
				temp[1] = m[2] * m[3] - m[0] * m[5];
				temp[2] = m[0] * m[4] - m[1] * m[3];

				return m[6] * temp[0] + m[7] * temp[1] + m[8] * temp[2];
		}

		inline mat3<T> inverse() const
		{
				T det = this->det();
				T data[9];

				data[0] = (m[4] * m[8] - m[5] * m[7]) / det;
				data[1] = -(m[1] * m[8] - m[2] * m[7]) / det;
				data[2] = (m[1] * m[5] - m[2] * m[4]) / det;
				data[3] = -(m[3] * m[8] - m[5] * m[6]) / det;
				data[4] = (m[0] * m[8] - m[2] * m[6]) / det;
				data[5] = -(m[0] * m[5] - m[2] * m[3]) / det;
				data[6] = (m[3] * m[7] - m[4] * m[6]) / det;
				data[7] = -(m[0] * m[7] - m[1] * m[6]) / det;
				data[8] = (m[0] * m[4] - m[1] * m[3]) / det;

				return mat3<T>(data);
		}

		//cast operator

		template <class U> operator mat3<U>() const
		{
				U cm[9];
				for (uint32_t i = 0; i < 9; i++)
						cm[i] = U(m[i]);

				return mat3<U>(cm);
		}

		//static func
		inline static mat3<T> indentity()
		{
				mat3<T> temp = mat3<T>();
				temp.m[0] = 1.0;
				temp.m[4] = 1.0;
				temp.m[8] = 1.0;
				return temp;
		}
};

//io operator

template <class T> std::ostream &operator<<(std::ostream &os, const mat3<T> &mat);

//operator

template <class T> const inline vec3<T> operator*(const mat3<T> &mat, const vec3<T> &a)
{
		return vec3<T>(mat.m[0 + 0] * a.x + mat.m[0 + 1] * a.y + mat.m[0 + 2] * a.z,
					   mat.m[3 + 0] * a.x + mat.m[3 + 1] * a.y + mat.m[3 + 2] * a.z,
					   mat.m[6 + 0] * a.x + mat.m[6 + 1] * a.y + mat.m[6 + 2] * a.z);
}

template <class T> const inline mat3<T> operator*(const mat3<T> &mat, const mat3<T> &a)
{
		T p[9];
		for (uint32_t i = 0; i < 3; i++)
		{
				for (uint32_t j = 0; j < 3; j++)
				{
						p[i * 3 + j] = mat.m[i * 3 + 0] * a.m[0 * 3 + j] + mat.m[i * 3 + 1] * a.m[1 * 3 + j] +
									   mat.m[i * 3 + 2] * a.m[2 * 3 + j];
				}
		}

		return mat3<T>(p);
}

// alias

using fmat3 = mat3<float>;
using dmat3 = mat3<float>;

//-----------------

template <class T> class mat2
{
	  public:
		T m[4];
		// m[0],m[1]
		// m[2],m[3]

		//constructor
		inline mat2(const vec2<T> &a, const vec2<T> &b)
		{
				m[0 * 2 + 0] = a.x;
				m[0 * 2 + 1] = b.x;
				m[1 * 2 + 0] = a.y;
				m[1 * 2 + 1] = b.y;
		}

		inline mat2(const T (&a)[4])
		{
				for (uint32_t i = 0; i < 4; i++)
						m[i] = a[i];
		}

		inline mat2(void)
		{
				for (uint32_t i = 0; i < 4; i++)
						m[i] = 0.0;
		}

		mat2(const double &omega);

		//menber func
		inline T det() const
		{
				return m[0] * m[3] - m[1] * m[2];
		}

		inline mat2<T> inverse() const
		{
				T det = this->det();
				T data[4] = {m[3] / det, -m[1] / det, -m[2] / det, m[0] / det};

				return mat2<T>(data);
		}
		inline T sqlength() const
		{
				return m[0] * m[0] + m[1] * m[1] + m[2] * m[2] + m[3] * m[3];
		}

		inline T trace() const
		{
				return m[0] + m[3];
		}

		inline mat2<T> transpose()
		{
				T data[4];
				data[0] = m[0];
				data[1] = m[2];
				data[2] = m[1];
				data[3] = m[3];
				return mat2<T>(data);
		}

		//cast operator

		template <class U> operator mat2<U>() const
		{
				U cm[4];
				for (uint32_t i = 0; i < 4; i++)
						cm[i] = U(m[i]);

				return mat2<U>(cm);
		}

		//static func
		inline static mat2<T> identity()
		{
				mat2<T> temp = mat2<T>();
				temp.m[0] = 1.0;
				temp.m[3] = 1.0;
				return temp;
		}
};

//io operator

template <class T> std::ostream &operator<<(std::ostream &os, const mat2<T> &mat);

//operator

template <class T> const inline vec2<T> operator*(const mat2<T> &mat, const vec2<T> &a)
{
		return vec2<T>(mat.m[0] * a.x + mat.m[1] * a.y, mat.m[2] * a.x + mat.m[3] * a.y);
}

template <class T> const inline mat2<T> operator*(const mat2<T> &mat, const mat2<T> &a)
{
		T p[4];
		for (uint32_t i = 0; i < 2; i++)
		{
				for (uint32_t j = 0; j < 2; j++)
				{
						p[i * 2 + j] = mat.m[i * 2 + 0] * a.m[0 * 2 + j] + mat.m[i * 2 + 1] * a.m[1 * 2 + j];
				}
		}

		return mat2<T>(p);
}

template <class T, class U> inline mat2<T> operator*(const mat2<T> &v, const U &a)
{
		T data[4];
		data[0] = v.m[0] * a;
		data[1] = v.m[1] * a;
		data[2] = v.m[2] * a;
		data[3] = v.m[3] * a;
		return mat2<T>(data);
}
template <class T, class U> inline mat2<T> operator*(const U &a, const mat2<T> &v)
{
		T data[4];
		data[0] = v.m[0] * a;
		data[1] = v.m[1] * a;
		data[2] = v.m[2] * a;
		data[3] = v.m[3] * a;
		return mat2<T>(data);
}

template <class T> inline mat2<T> operator+(const mat2<T> &v, const mat2<T> &a)
{
		T data[4];
		data[0] = v.m[0] + a.m[0];
		data[1] = v.m[1] + a.m[1];
		data[2] = v.m[2] + a.m[2];
		data[3] = v.m[3] + a.m[3];
		return mat2<T>(data);
}
template <class T> inline mat2<T> operator-(const mat2<T> &v, const mat2<T> &a)
{
		T data[4];
		data[0] = v.m[0] - a.m[0];
		data[1] = v.m[1] - a.m[1];
		data[2] = v.m[2] - a.m[2];
		data[3] = v.m[3] - a.m[3];
		return mat2<T>(data);
}

// alias

using fmat2 = mat2<float>;
using dmat2 = mat2<double>;

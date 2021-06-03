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
				m[0 * 3 + 1] = a.y;
				m[0 * 3 + 2] = a.z;
				m[1 * 3 + 0] = b.x;
				m[1 * 3 + 1] = b.y;
				m[1 * 3 + 2] = b.z;
				m[2 * 3 + 0] = c.x;
				m[2 * 3 + 1] = c.y;
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
				m[0 * 2 + 1] = a.y;
				m[1 * 2 + 0] = b.x;
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

		//menber func
		inline T det() const
		{
				return m[0] * m[3] - m[1] * m[2];
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
		inline static mat2<T> indentity()
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

// alias

using fmat2 = mat2<float>;
using dmat2 = mat2<double>;

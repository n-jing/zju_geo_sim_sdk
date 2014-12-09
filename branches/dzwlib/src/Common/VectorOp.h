/**   
* @file vector_op.h 
* @brief TODO
* @author dangzw
* @date Nov 17, 2011 4:12:19 PM 
* @version V1.0   
*/

#ifndef VECTOR_OP_H_
#define VECTOR_OP_H_

#include <vector>
#include <cassert>
#include <iostream>

template<typename T>
T sum(const std::vector<T>& V1){
	T temp = 0;
	for(size_t i = 0; i < V1.size(); ++i){
		temp += V1[i];
	}
	return temp;
}

template<typename T>
T norm(const std::vector<T>& V1){
	T temp = 0;
	for(size_t i = 0; i < V1.size(); ++i){
		temp += (V1[i] * V1[i]);
	}
	return sqrt(temp);
}


template<typename T>
std::vector<T> operator +(const std::vector<T>& V1, const std::vector<T>& V2){
	assert(V1.size() == V2.size());
	std::vector<T> temp(V1.size());
	for(size_t i = 0; i < temp.size(); ++i){
		temp[i] = V1[i] + V2[i];
	}
	return temp;
}

template<typename T>
std::vector<T> operator -(const std::vector<T>& V1, const std::vector<T>& V2){
	assert(V1.size() == V2.size());
	std::vector<T> temp(V1.size());
	for(size_t i = 0; i < temp.size(); ++i){
		temp[i] = V1[i] - V2[i];
	}
	return temp;
}

template<typename T1, typename T2>
std::vector<T1> operator /(const std::vector<T1>& V1, T2 rhs){
	assert(fabs(rhs) > 1e-10);
	std::vector<T1> temp(V1.size());
	for(size_t i = 0; i < temp.size(); ++i){
		temp[i] = V1[i] / rhs;
	}
	return temp;
}

template<typename T>
std::vector<T>& operator -(const std::vector<T>& rhs){

}

template<typename T>
std::ostream& operator <<(std::ostream& os, const std::vector<T>& V){
	os << "vector size : " << V.size() << std::endl;
	for(size_t i = 0; i < V.size(); ++i){
		os << V[i] << ((i + 1) % 10 ? ' ' : '\n');
	}
	os << std::endl;
	return os;
}

template<typename T1, typename T2>
void MultScalar(std::vector<T1>& V1, T2 rhs){
	for(size_t i = 0; i < V1.size(); ++i){
		V1[i] *= rhs;
	}
}

#endif /* VECTOR_OP_H_ */

/**   
* @file MyException.h 
* @brief TODO
* @author dangzw
* @date Oct 24, 2011 6:58:45 PM 
* @version V1.0   
*/

#ifndef MYEXCEPTION_H_
#define MYEXCEPTION_H_

#include <iostream>
#include <exception>
#include <string>

namespace dzw{namespace common{
/**Define a class MyException derived from exception.
 */
class my_exception : public std::exception{
public:
	my_exception(const std::string& error_info_) throw() : error_info(error_info_) {}
	virtual ~my_exception() throw() {}
	virtual const char* what() const throw(){
		return error_info.c_str();
	}
protected:
	std::string error_info;
};

}
}
#endif /* MYEXCEPTION_H_ */

#ifndef HJ_ARG_OPTS_H_
#define HJ_ARG_OPTS_H_

#include <map>
#include <sstream>
#include <cassert>
#include <iostream>
#include <typeinfo>

#include <boost/any.hpp>

//! @brief an argument from command line has its type and its status:
//! set, not set, wrong.  If no required option, or wrong value format
//! in the command line, throw exception.

namespace hj {

#ifdef __GNUG__
#  define ARG_OPTS_DEPRECATED __attribute__((deprecated))
#endif

#ifdef _MSC_VER
#  define ARG_OPTS_DEPRECATED __declspec(deprecated)
#endif

class ARG_OPTS_DEPRECATED any_structure
{
public:
	virtual ~any_structure(){}

	class exception
	{
	public:
		exception(const std::string &what)
			:what_(what) {
		}
		const char *what(void) const { return what_.c_str(); }
	private:
		std::string what_;
	};

	bool has(const std::string &key) const {
		return fields_.find(key) != fields_.end();
	}

	template <typename T>
	T get(const std::string &key) const {
		container::const_iterator it = fields_.find(key);
		if(it == fields_.end())
			throw exception(key + "not found");
		return boost::any_cast<T>(it->second);
	}

	const boost::any &operator[](const std::string &key) const {
		container::const_iterator it = fields_.find(key);
		if(it == fields_.end())
			throw exception(key + "not found");
		return it->second;
	}
	boost::any &operator[](const std::string &key) {
		return fields_[key];
	}
protected:
	typedef std::map<std::string, boost::any> container;
	container fields_;
};

class ARG_OPTS_DEPRECATED arg_opts : public any_structure
{
public:
	static const char *make_style[];

	static arg_opts parse(int argc, char *argv[], const char *style[] = make_style);

	class exception
	{
	public:
		exception(const std::string &key, const std::string &reason)
			:msg_(key) {
			msg_ = msg_ + ": " + reason;
		}
		const char *what(void) const {
			return msg_.c_str();
		}
	private:
		std::string msg_;
	};

	//! optional with/withougout value
	const bool has(const std::string &key) const {
		return any_structure::has(key);
	}
	const bool has(const std::string &key, const char *desc) {
		set(key, desc);
		return has(key);
	}
	int set(const std::string &key, const char *desc, bool reset = false) {
		assert(desc);
		if(usage_.find(key) != usage_.end() && !reset)
			return 1;
		usage_[key] = usage_type(desc, false);
		return 0;
	}

	// optional with default value
	template <typename T>
	T get(const std::string &key, const T &default_value) const {
		if(any_structure::has(key))
			return get<T>(key);
		return default_value;
	}
	template <typename T>
	T get(const std::string &key, const T &default_value, const char *desc) {
		set(key, default_value, desc);
		return get<T>(key, default_value);
	}
	template <typename T>
	int set(const std::string &key, const T &default_value, const char *desc, bool reset = false) {
		if(usage_.find(key) != usage_.end() && !reset)
			return 1;
		std::ostringstream oss;
		oss << desc << " (default is: " << default_value << ")";
		usage_[key] = usage_type(oss.str().c_str(), false);
		return 0;
	}

	template <typename T>
	T get(const std::string &key) const {
		if(any_structure::has(key)) {
			std::istringstream iss(any_structure::get<std::string>(key).c_str());
			T rtn;
			iss >> rtn;
			return rtn;
		}
		throw exception(key, "not found");
	}
	template <typename T>
	T get(const std::string &key, const char *desc) {
		set<T>(key, desc);
		return get<T>(key);
	}
	template <typename T>
	int set(const std::string &key, const char *desc, bool reset = false) {
		assert(desc);
		if(usage_.find(key) != usage_.end() && !reset)
			return 1;
		usage_[key] = usage_type(desc, true);
		return 0;
	}

	std::string usage(const char *style[] = make_style) const;
private:
	class usage_type
	{
	public:
		usage_type(){}
		usage_type(const char *desc, bool is_required)
			:desc_(desc), is_required_(is_required) {
		}
		std::string desc_;
		bool is_required_;
	};
	typedef std::map<std::string, usage_type> usage_container;
	usage_container usage_;
};

}

#endif

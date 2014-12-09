#ifndef _UTIL_H_
#define _UTIL_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace zjucad
{
namespace LSCO
{

enum RESIDUAL_TYPE {
    NORM_2,
    NORM_1,
    NORM_0,
};

template <typename KERNEL_TYPE>
typename KERNEL_TYPE::value_type residual(const typename KERNEL_TYPE::vector_type& z,int n_type)
{
    switch(n_type) {
    case NORM_2:
        return KERNEL_TYPE::nrm2(z);
    case NORM_1:
        return KERNEL_TYPE::asum(z);
    case NORM_0:
        return KERNEL_TYPE::amax(z);
    default:
        return KERNEL_TYPE::nrm2(z);
    }
}

void debug_tree(const std::string& path,const boost::property_tree::ptree& pt)
{
    boost::property_tree::xml_writer_settings<char> settings('\t',1);
    write_xml(path,pt,std::locale(),settings);
}

template <typename KERNEL_TYPE>
typename KERNEL_TYPE::value_type limit_bound(const typename KERNEL_TYPE::vector_type& x,const typename KERNEL_TYPE::vector_type& d,const std::vector<std::pair<size_t,typename KERNEL_TYPE::value_type> >& bounds) {
	typename KERNEL_TYPE::value_type max_c=numeric_limits<typename KERNEL_TYPE::value_type>::max();
	typename KERNEL_TYPE::value_type dv;
	for(size_t i=0;i<bounds.size();i++) {
		const std::pair<size_t,typename KERNEL_TYPE::value_type>& bd=bounds[i];
		dv=KERNEL_TYPE::get(bd.first,d);
		if(dv < -1E-9f)
			max_c=std::min<typename KERNEL_TYPE::value_type>(max_c,(bd.second-KERNEL_TYPE::get(bd.first,x))/dv);
	}return max_c;
}

template <typename KERNEL_TYPE>
typename KERNEL_TYPE::value_type limit_bound_n(const typename KERNEL_TYPE::vector_type& x,const std::vector<std::pair<size_t,typename KERNEL_TYPE::value_type> >& bounds) {
	typename KERNEL_TYPE::value_type max_c=numeric_limits<typename KERNEL_TYPE::value_type>::max();
	typename KERNEL_TYPE::value_type dv;
	for(size_t i=0;i<bounds.size();i++) {
		const std::pair<size_t,typename KERNEL_TYPE::value_type>& bd=bounds[i];
		dv=KERNEL_TYPE::get(bd.first,x);
		if(dv > 1E-9f)
			max_c=std::min<typename KERNEL_TYPE::value_type>(max_c,-bd.second/dv);
	}return max_c;
}

template <typename KERNEL_TYPE>
bool check_bound(const typename KERNEL_TYPE::vector_type& x,const std::vector<std::pair<size_t,typename KERNEL_TYPE::value_type> >& bounds) {
	for(size_t i=0;i<bounds.size();i++) {
		const std::pair<size_t,typename KERNEL_TYPE::value_type>& bd=bounds[i];
		if(KERNEL_TYPE::get(bd.first,x) <= bd.second-1E-3f)
			return false;
	}
	return true;
}

}
}

#endif
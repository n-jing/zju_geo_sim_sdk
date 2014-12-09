#ifndef _DIRECT_SOLVER_KERNEL_REPORT_H_
#define _DIRECT_SOLVER_KERNEL_REPORT_H_

#include <stdint.h>
#include <boost/property_tree/ptree.hpp>

namespace zjucad
{
namespace LSCO
{

static inline void direct_solver_kernel_report(const boost::property_tree::ptree& pt)
{
    std::cout << "n_pos: " << pt.get<int64_t>("direct_solver_out.n_pos") << std::endl;
    std::cout << "n_neg: " << pt.get<int64_t>("direct_solver_out.n_neg") << std::endl;
    std::cout << "n_zero: " << pt.get<int64_t>("direct_solver_out.n_zero") << std::endl;
    std::cout << "termination type: " << pt.get<std::string>("direct_solver_out.termination_type") << std::endl;
}

}
}

#endif

#include "../include/arg_opts.h"

#include <string.h>
#include <iostream>

using namespace std;

namespace hj {

const char *arg_opts::make_style[] ={"make", "=", 0};

arg_opts arg_opts::parse(int argc, char *argv[], const char *style[])
{
	arg_opts opts;
	if(!strcmp(style[0], "make")) {
		size_t token_len = strlen(style[1]);
		for(int i = 0; i < argc; ++i) {
			string k_v = argv[i];
			size_t pos = k_v.find(style[1]);
			if(pos != string::npos) {
				opts.fields_.insert(make_pair(k_v.substr(0, pos), k_v.substr(pos+token_len)));
			}
		}
	}
	else {
		cerr << "only support make style currently." << endl;
	}
	return opts;
}

std::string arg_opts::usage(const char *style[]) const
{
	ostringstream os;
	size_t tab[1] = {24};
	for(usage_container::const_iterator i = usage_.begin();
		i != usage_.end(); ++i) {
		ostringstream os0;
		os0 << ((i->second.is_required_)?' ':'[')
		   << i->first
		   << ((i->second.is_required_)?' ':']');
		os << os0.str();
		if(os0.tellp() < tab[0]) {
			for(size_t i = os0.tellp(); i < tab[0]; ++i)
				os << ' ';
		}
		os << i->second.desc_ << "\n";
	}
	return os.str();
}

}

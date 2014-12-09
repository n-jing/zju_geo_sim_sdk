#include <string>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <boost/typeof/std/utility.hpp>
#include <boost/foreach.hpp>
#include "mk_ptree.h"

namespace jtf
{
  using namespace std;
  using namespace boost::program_options;
 
  const char *parser_cmd_line::style_[]={"make","=",0};
  
  void parser_cmd_line::parser(int argc, char * argv[])
  {
	init_option_data();

	typedef vector<tuple<string,string,string,fPtr> >::const_iterator vci;
	for(vci it = data_.begin(); it != data_.end(); ++it)
	  {
		std::string option_name = it->get<0>();
		if(!it->get<1>().empty())
		  {
			option_name += ",";
			option_name += it->get<1>();
		  }
		const string &desc = it->get<2>();
		opts.add_options()(option_name.c_str(), desc.c_str());
	  }

	special_parser(argc,argv,style_);
	
	variables_map vm;
	BOOST_AUTO(pr, command_line_parser(argc,argv).options(opts).allow_unregistered().run());
	vector<string> ur_opts = collect_unrecognized(pr.options, include_positional);
	BOOST_FOREACH(string str, ur_opts) 
	  {
		if(str.find(style_[1]) == std::string::npos)
		  std::cerr << "Unknown option: " << str << std::endl;
	  }
	store(pr,vm);
	notify(vm);
	
	if(vm.size() == 0 && argc == 1)
	  {
		std::cout << opts << std::endl;
		return ;
	  }

	for(vci it = data_.begin();	it != data_.end(); ++it)
	  {
		const std::string& key = it->get<0>();
		if(vm.count(key.c_str()))
		   it->get<3>()(key);
	  }
	
  }

  // static void save(const std::string & file_path, const ptree &pt)
  // {
  // 	size_t last_dot_pos = file_path.find(".");
	  
  // 	  if(last_dot_pos == std::string::npos)
  // 		throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");
	  
  // 	  std::string file_format = file_path.substr(last_dot_pos + 1);
  
  // 	  if(file_format == "xml" || file_format == "XML")
  // 		write_xml(file_path, pt);
  // 	  else if(file_format == "json" || file_format == "JSON")
  // 		write_json(file_path, pt);
  // 	  else if(file_format == "ini" || file_format == "INI")
  // 		write_ini(file_path, pt);
  // 	  else if(file_format == "info" || file_format == "INFO")
  // 		write_info(file_path, pt);
  // 	  else 
  // 		throw std::logic_error("argument error: do not support this format, only support [xml/json/ini/info]");
  // }

  void mk_ptree::show_format_data(const std::string & format, std::ostream& out)
  {
 	if(is_show_format_data)
	  {
		if(format == "xml")
		  write_xml(out, pt);
		else if(format == "json")
		  write_json(out, pt);
		else if(format == "ini")
		  write_ini(out, pt);
		else if(format == "info")
		  write_info(out, pt);
	  }
  }

  void mk_ptree::init_option_data()
  {
	function<void(const std::string &)> f_help,f_show_format,f_enable_show;
	f_help = bind(&mk_ptree::show_help_info,this);
	f_show_format = bind(&mk_ptree::show_format_data,this,_1,ref(std::cout));
	f_enable_show = bind(&mk_ptree::enable_show_format_data,this);
	data_.push_back(make_tuple("help","h","help information",f_help));
	data_.push_back(make_tuple("show","v","show format data",f_enable_show));
	data_.push_back(make_tuple("json","j","output format: json",f_show_format));
	data_.push_back(make_tuple("ini","i","output format: ini",f_show_format));
	data_.push_back(make_tuple("info","f","output format: info",f_show_format));
	data_.push_back(make_tuple("xml","x","output format: xml",f_show_format));
  }

  void mk_ptree::special_parser(int argc, char * argv[], const char * style[])
  {
	if(!strcmp(style[0], "make")) 
	  {
		size_t token_len = strlen(style[1]);
		for(int i = 0; i < argc; ++i) {
		  string k_v = argv[i];
		  size_t pos = k_v.find(style[1]);
		  if(pos != string::npos) 
			{
			  std::string path = k_v.substr(0,pos);
			  std::string value = k_v.substr(pos + token_len);
			  if(path.find(".") == std::string::npos)
				{
				  pt.put(path + ".value", value);
				}
			  else
				  pt.put(path + ".value", value);
			}
		}
	  }
	else {
	  cerr << "only support make style currently." << endl;
	}
  }
 
}

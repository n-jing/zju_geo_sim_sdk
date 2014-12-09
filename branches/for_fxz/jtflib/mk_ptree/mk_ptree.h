#ifndef MK_PTREE_H
#define MK_PTREE_H
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/program_options.hpp>
#include <boost/tuple/tuple.hpp>
#include <string>
#include <iostream>

/**
 * @file   mk_ptree.h
 * @author tfjiang <tfjiang@CAD>
 * @date   Fri Oct 21 11:16:46 2011
 * 
 * @brief  This program is used to simpley convert user-defined parameters into std stream.\n
 *         It uses boost::program_options to parse regular options like "-h,--help". \n
 *         You can use a derived class to add a tuple as: \n
 *         "option name-long term, option name short term,option description, option handler" \n
 *         Long term, short term and description are std::string. \n
 *         Handler is a boot::function, which has a "input const std::string&" and no ouput. \n
 *         Usage: \n
 *         To use this parser, you first init_data, to add option name, description and handler into "data_", then use special_parser to pase some special style.\n
 *         Then call the parser function, it will be done well.\n
 *         Limitions:\n
 *         1. Only parse option with "-/--", like "-h,--help" \n
 *         2. All options are used as switchers, if you want to get value from option, please use special_parser. \n
 */

namespace jtf
{
		using namespace std;
		using namespace boost;
		using namespace boost::program_options;
		using boost::property_tree::ptree;
		
		class parser_cmd_line
		{
		public:
				/** 
				 * Parse the main parameters. It depends on the init_option_data and special_parser function \n
				 * 
				 * @param argc arguments number from
				 * @param argv arguments from main
				 */
				void parser(int argc, char *argv[]);
				virtual ~parser_cmd_line(){}
		protected:
				/** 
				 * This function init the data_, which needs to be rewrote, it includes sereral steps: \n
				 * 1. Add option term with long and short term. \n
				 * 2. Add option description and option handler. \n
				 * 3. Both terms above are pair format <long term, short term>, <description, handler>
				 * 4. The handler uses boost::function, which should bind to the function you need. 
				 */
				virtual void init_option_data(){}
				/** 
				 * Add a user designed parser, you can parse some special style input as the "style_" offered.
				 * It also needs to be rewrote.
				 * @param argc 
				 * @param argv 
				 */
				virtual void special_parser(int argc, char *argv[], const char * style[] = style_){}

				typedef function<void(const std::string&)> fPtr;
				vector<tuple<string,string,string,fPtr> > data_;
				static const char * style_[];
				options_description opts;
			};
		
		
		class mk_ptree: public parser_cmd_line
		{
		public:
		        mk_ptree():is_show_format_data(false){}
		protected:
				virtual void init_option_data();
				/** 
				 * Parse the "A=B" sytle, can not contain any space befor or behind "=", the pair<A,B> will be stored into ptree
				 * This application is used to convert input "A=B" style to format data as [xml/json/ini/info] and 
				 * @param argc 
				 * @param argv 
				 * @param style char *style_[]={"make","=","B"}
				 */
				virtual void special_parser(int argc, char *argv[], const char * style[] = style_);
		private:
				void show_help_info(){std::cout << opts << std::endl;}
				void show_format_data(const std::string & format,std::ostream& out);
				void enable_show_format_data(){ is_show_format_data = true;	}
				ptree pt;
				bool is_show_format_data;
		};
}
#endif

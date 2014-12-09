#include <zjucad/matrix/matrix.h>
#include <zjucad/ptree/ptree.h>
#include <iostream>

using namespace std;
using boost::property_tree::ptree;
using namespace zjucad::matrix;

int test_gmres1();
int test_gmres2();
int test_gmres3();
int test_qmr();
int test_qmr2();
int test_qmr3();
//int test_kkt(ptree & pt);

int main(int argc, char * argv[])
{
	ptree pt;

	try{
		zjucad::read_cmdline(argc, argv, pt);
		//	test_kkt(pt);
		//		test_qmr3();
	}catch (std::exception &e){
		cerr << endl;
		cerr << e.what() << endl;
		zjucad::show_usage_info(cerr, pt);
	}

  return 0;
}

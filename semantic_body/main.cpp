#include <iostream>
#include <binary_io.h>
#include <time.h>
#include <gurobi_c++.h>

#include "mesh_io.h"
#include "measure.h"
#include "preprocessing.h"

void basicInitialize() {
	std::string path = "C:/Users/amomorning/dataset/SPRING_MALE/SPRING_MALE/";

	vector<string> files = getFiles(path + "*");

	////calcAverage(AVEobj, path, files);
	string filename = "./F";
	//saveBinVerts(filename, path, files);
	saveBinFaces(filename, path, files);

	////calcNeighbor();


}


void measureOne() {
	clock_t t = clock();

	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	readAVE(V, F);

	measure measure;
	measure.calcExact(V, F);

	cout << "Total time used...." << endl;
	cout << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds..." << endl;
}

void testGurobi() {
	try
	{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		cout << "ok" << endl;
	}
	catch (const std::exception&)
	{
		cout << "not ok" << endl;
	}
}

int main()
{

	measureOne();
	//testGurobi();
	getchar();
	return 0;
}

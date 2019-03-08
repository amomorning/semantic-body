#include <iostream>
#include <binary_io.h>
#include <time.h>
#include <gurobi_c++.h>

#include "mesh_io.h"
#include "measure.h"
#include "preprocessing.h"

void GenerateBinaryData() {
	std::string path = "C:/Users/amomorning/dataset/SPRING_MALE/SPRING_MALE/";

	vector<string> files = getFiles(path + "*");

	////calcAverage(AVEobj, path, files);

	string vertsFilename = "./data/V";
	string facesFilename = "./data/F";

	//saveBinVerts(filename, path, files);
	saveBinFaces(facesFilename, path, files);
	// with F shape (3, 12500)
	saveBinVerts(vertsFilename, path, files);
	// with V shape (37500, 1511)

	//calcAverage();
	Eigen::MatrixXd V;
	Eigen::Matrix3Xi F;

	common::read_matrix_binary_from_file("./data/V", V);
	common::read_matrix_binary_from_file("./data/F", F);
	cout << F.cols() << endl;
	calcAverage(V, F);
	//calcNeighbor();
	calcNeighbor();
	// with neighbour shape (12500, 11)

	cout << "All is done!" << endl;

}


void measureOne() {

	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	readAVE(V, F);

	measure measure;
	measure.calcExact(V, F);

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
	clock_t t = clock();
	//testGurobi();

	cout << "Total time used...." << endl;
	cout << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds..." << endl;
	getchar();
	return 0;
}

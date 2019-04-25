#include <iostream>
#include <binary_io.h>
#include <time.h>
#include <gurobi_c++.h>
#include <sophus/so3.hpp>
#include "mesh_io.h"
#include "measure.h"
#include "preprocessing.h"
using namespace std;

void saveBinaryVFN() {
	std::string trainModelPath = "./dataset/train/";
	std::string testModelPath = "./dataset/test/";
	
	vector<string> trainFiles = getFiles(trainModelPath + "*");
	vector<string> testFiles = getFiles(testModelPath + "*");

	// trainV (37500, 1400)
	saveBinVerts("./data/train/V", trainModelPath, trainFiles);

	// testV (37500, 111)
	saveBinVerts("./data/test/V", testModelPath, testFiles);

	// F (3, 12500)
	saveBinFaces("./data/F", trainModelPath, trainFiles);

	Eigen::MatrixXd trainV;
	Eigen::MatrixXd testV;
	Eigen::Matrix3Xi F;

	common::read_matrix_binary_from_file("./data/train/V", trainV);
	common::read_matrix_binary_from_file("./data/test/V", testV);
	common::read_matrix_binary_from_file("./data/F", F);
	cout << F.cols() << endl;

	//calcAverage();
	calcAverage("./data/AVE.obj", trainV, F);

	// train dV (37500, 1400)
	calcDeltaVerts("./data/train/dV", trainV);

	// test dV (37500, 111)
	calcDeltaVerts("./data/test/dV", testV);

	// N (12500, 11)
	saveNeighbor("./data/N");
	cout << "All is done!" << endl;

}

void testBinaryVFN() {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	common::read_matrix_binary_from_file("./data/train/V", V);
	printShape(V);

	common::read_matrix_binary_from_file("./data/test/V", V);
	printShape(V);


	common::read_matrix_binary_from_file("./data/F", F);
	printShape(F);
	Eigen::MatrixXd newV = V.col(0);
	newV.resize(3, 12500);
	common::save_obj("./data/TEMP.obj", newV, F);

	common::read_matrix_binary_from_file("./data/N", F);
	printShape(F);
}

void measureOne() {

	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("./data/AVE.obj", V, F);
	printShape(V, F);

	measure measure;
	//measure.calcDijkstra(V, F);
	measure.calcExact(V, F);
}

void testGurobi() {
	try
	{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		GRBQuadExpr exp = 0;


		cout << "ok" << endl;
	}
	catch (const std::exception&)
	{
		cout << "not ok" << endl;
	}
}

void testEigen(Eigen::Vector3d &v) {
	v = Eigen::Vector3d::Ones();
	return;
}

void saveBinaryMeasure() {
	Eigen::MatrixXd trainV;
	Eigen::MatrixXd testV;
	Eigen::Matrix3Xi F;

	common::read_matrix_binary_from_file("./data/F", F);
	common::read_matrix_binary_from_file("./data/train/V", trainV);
	common::read_matrix_binary_from_file("./data/test/V", testV);

	//saveDijkstra("./data/test/dijkstra", testV, F);

	//measure one;
	//saveRoughExact("./data/test/roughExact", testV, F);

	//saveExact("./data/test/exact", testV, F);
}

void saveBinaryFRS() {
	Eigen::MatrixXd trainV;
	Eigen::MatrixXd testV;
	Eigen::Matrix3Xi F;

	common::read_matrix_binary_from_file("./data/F", F);
	common::read_matrix_binary_from_file("./data/train/V", trainV);
	common::read_matrix_binary_from_file("./data/test/V", testV);

	//saveFaceFeature("./data/train/logRS", trainV, F);
	//cout << F.rows() << " " << F.cols() << endl;
	
	saveFaceFeature("./data/test/logRS", testV, F);
	cout << "saved!!!" << endl;
}

void testSophus() {
	Eigen::Matrix3d R = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d(0, 0, 1)).toRotationMatrix();

	Sophus::SO3d SO3_R(R);
	cout << R << endl << endl;

	Eigen::Vector3d so3 = SO3_R.log();
	cout << so3.transpose() << endl;
	cout << Sophus::SO3d::hat(so3).transpose() << endl;

	Sophus::SO3d tSO3_R = Sophus::SO3d::exp(so3);
	cout << tSO3_R.matrix() << endl;
}

int main()
{
	clock_t t = clock();
	saveBinaryFRS();

	
	cout << (double)(clock() - t) / CLOCKS_PER_SEC << "seconds..." << endl;
	getchar();
	return 0;
}

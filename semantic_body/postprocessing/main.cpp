#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
#include <iostream>
#include <binary_io.h>
#include <surface_mesh/Surface_mesh.h>
#include "mesh_io.h"

using namespace surface_mesh;


void readNewPoint(const char * filename, Eigen::Matrix3Xd &V) {
	std::ifstream is(filename);
	double x, y, z;
	int cnt = 0;
	while (is >> x >> y >> z) {
		//std::cout << x << " " << y << " " << z << std::endl;
		V(0, cnt) += x;
		V(1, cnt) += y;
		V(2, cnt) += z;
		cnt++;
	}
	std::cout << " cnt = " << cnt << std::endl;
	return;

}


void calcFeature(const char * filename, Eigen::Matrix3Xd V, 
	const Eigen::Matrix3Xi &F, Eigen::MatrixXd feature) 
{
	Eigen::Matrix3Xd V0;
	feature.resize(12500, 9);
	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("./data/neighbor", N);

	//ÎÑÌØß÷È¥ÄÄÀïÕÒCij£¬Á¹Á¹¡£¡£¡£
	for (int i = 0; i < 12500; ++i) {
		Eigen::MatrixXd T = feature.row(i);
		T.resize(3, 3);
		int a;
		for (int j = 0; j < 11; ++j) {
			int k = N(i, j) - 1;
			if (k < 0) return;
			

		}
		for (int j = 0; j < 3; ++j) {
		}
	}

	common::save_obj(filename, V, F);
}

//Generate data from new feature matrix....
void recoverFromFeature(const char * filename, int total) {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("../data/AVE.obj", V, F);

	std::ifstream in(filename);
	Eigen::MatrixXd feature;
	feature.resize(total, 112500);
	int t = 0;
	for(int i = 0; i < total; ++ i) {
		for(int j = 0; j < 112500; ++ j) {
			in >> feature(i, j);
		}
	}
	in.close();

	for (int i = 0; i < total; ++i) {
		std::string name = "../data/augmentation/" + std::to_string(i) + ".obj";
		calcFeature(name.c_str(), V, F, feature.col(i));
	}
}

void gradientOptimization(const char * filename,
	Eigen::Matrix3Xd &V, const Eigen::Matrix3Xi &F)
{

}

int main() {
	//Eigen::Matrix3Xd V;
	//Eigen::Matrix3Xi F;
	//common::read_obj("../data/AVE.obj", V, F);
	//readNewPoint("../data/New.txt", V);
	//
	//common::save_obj("../data/NEW.obj", V, F);

	recoverFromFeature("./data/T.txt", 1);
	getchar();
}
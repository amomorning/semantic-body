#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include "mesh_io.h"

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

void gradientOptimization(const char * filename,
	Eigen::Matrix3Xd &V, const Eigen::Matrix3Xi &F)
{

}

int main() {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("../data/AVE.obj", V, F);
	readNewPoint("../data/New.txt", V);
	
	common::save_obj("../data/NEW.obj", V, F);
	getchar();
}
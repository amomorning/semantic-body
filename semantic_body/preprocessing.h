#ifndef PRE_PROCESSING_H_H_
#define PRE_PROCESSING_H_H_

#include <vector>
#include <Eigen/Dense>
#include <surface_mesh/Surface_mesh.h>
using namespace std;
using namespace surface_mesh;

template<class Matrix>
void printShape(Matrix &M) {
	cout << M.rows() << " " << M.cols() << endl;
	return;
}

template<class Matrix, class ...Args>
void printShape(Matrix &M, Args... rest) {
	printShape(M);
	printShape(rest...);
}


vector<string> getFiles(string cate_dir);

void calcAverage(Eigen::MatrixXd V, Eigen::Matrix3Xi F);

void saveBinVerts(const char *filename, string &path, vector<string> files);

void saveBinFaces(const char *filename, string &path, vector<string> files);

void calcNeighbor();

double laplacianCotanWeight(const Surface_mesh &mesh, const int i, const int j);

void saveFeature(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F);

void calcFeature(const Eigen::Matrix3Xd &V, const Eigen::Matrix3Xi &F, Eigen::MatrixXd &feature);


#endif DATA_PROCESS_H_H_


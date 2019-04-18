#ifndef PRE_PROCESSING_H_H_
#define PRE_PROCESSING_H_H_

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
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


vector<string> getFiles(const string &cate_dir);

void calcAverage(const char * filename, const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F);

void calcDeltaVerts(const char* filename, Eigen::MatrixXd &V);

void saveBinVerts(const char *filename, const string &path, const vector<string> &files);

void saveBinFaces(const char *filename, const string &path, const vector<string> &files);

void saveNeighbor(const char *filename);

void saveDijkstra(const char* filename, const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F);
void saveExact(const char* filename, const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F);

void saveRoughExact(const char* filename, const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F);

void saveFeature(const char* filename, const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F);


#endif DATA_PROCESS_H_H_


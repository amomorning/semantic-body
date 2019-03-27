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

void calcAverage(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F);

void saveBinVerts(const char *filename, const string &path, const vector<string> &files);

void saveBinFaces(const char *filename, const string &path, const vector<string> &files);

void saveNeighbor();

void recoverFromFeature(const char * filename);

void saveDijkstra(const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F);
void saveExact(const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F);

void laplacianCotanWeight(const Surface_mesh &mesh, Eigen::SparseMatrix<double> &cotan);

void saveFeature(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F);

void calcFeature(const Eigen::Matrix3Xd &V, const Eigen::Matrix3Xi &F, const Eigen::Matrix3Xd &AVE, Eigen::MatrixXd &feature, int u);
void calcFeatureGurobi(const Eigen::Matrix3Xd &V, const Eigen::Matrix3Xi &F, Eigen::MatrixXd &feature);

void saveVertsOffset(const Eigen::MatrixXd &V);

#endif DATA_PROCESS_H_H_


#ifndef PRE_PROCESSING_H_H_
#define PRE_PROCESSING_H_H_

#include <vector>
#include <Eigen/Dense>
#include <surface_mesh/Surface_mesh.h>
using namespace std;
using namespace surface_mesh;

vector<string> getFiles(string cate_dir);

void readAVE(Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F);

void calcAverage(Eigen::MatrixXd V, Eigen::Matrix3Xi F);

Eigen::MatrixXd getNeighbor(int id, Eigen::VectorXi &N, Eigen::VectorXd &V);

Eigen::Matrix3d getAffineMatrix();

void getDeformation(Eigen::VectorXd ave, Eigen::MatrixXd &V);

void saveBinVerts(string &filename, string &path, vector<string> files);

void saveBinFaces(string &filename, string &path, vector<string> files);

void calcNeighbor();

void calcAffineMatrix();

#endif DATA_PROCESS_H_H_

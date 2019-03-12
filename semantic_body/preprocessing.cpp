// DataProcess.cpp : 
//
#include "mesh_io.h"
#include "measure.h"
#include "binary_io.h"
#include "vtk.h"
#include "tri_mesh.h"
#include "preprocessing.h"
#include <Eigen/Dense>
#include <iostream>
#include <io.h>
#include <direct.h>
#include <surface_mesh/Surface_mesh.h>
#include <geodesic/geodesic_algorithm_dijkstra.h>
#include <geodesic/geodesic_algorithm_exact.h>

#include <gurobi_c++.h>

#define VERTS 12500

using namespace std;
using namespace surface_mesh;

// Return a vector<string> of files in input cateloge direction.
// The files will be sorted.

vector<string> getFiles(string cate_dir)
{
	vector<string> files;

	_finddata_t file;
	intptr_t lf;
	// the type should be intptr_t in x64 machine
	// but it will be fine using long in x86

	if ((lf = _findfirst(cate_dir.c_str(), &file)) == -1) {
		cout << cate_dir << " not found!!!" << endl;
	}
	else {
		while (_findnext(lf, &file) == 0) {
			// cout<<file.name<<endl;
			if (strcmp(file.name, ".") == 0 || strcmp(file.name, "..") == 0)
				continue;
			files.push_back(file.name);
		}
	}
	_findclose(lf);

	sort(files.begin(), files.end());
	return files;
}



//Calculate the average data
void calcAverage(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F) {
	Eigen::MatrixXd newV;
	newV = V.rowwise().mean();
	newV.resize(3, VERTS);

	common::save_obj("./data/AVE.obj", newV, F);
	//cout << "ojk";
}


void saveDijkstra(const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F)
{
	Eigen::MatrixXd ret;
	measure measure;
	
	ret.resize(measure.len(), V.cols());

	for (int i = 0; i < V.cols(); ++i) {
		Eigen::MatrixXd tmp = V.col(i);
		tmp.resize(3, VERTS);
		measure.calcDijkstra(tmp, F);
		measure.saveParam(ret, i);
	}

	cout << ret << endl;
	common::write_matrix_binary_to_file("./data/dijkstra", ret);
}

void saveExact(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F)
{

}

// Save all verts in the shape of (3|V|, N) in binary. 
void saveBinVerts(const char *filename, const string &path, const vector<string> &files) {
	Eigen::MatrixXd V;
	V.resize(VERTS * 3, files.size());
	cout << "V.shape = " << V.rows() << " " << V.cols() << endl;
	int k = 0;
	for (auto file : files) {
		Eigen::Matrix3Xd vv;
		Eigen::Matrix3Xi ff;
		common::read_obj(path + file, vv, ff);
		int cnt = 0;
		cout << vv.cols() << " " << vv.rows() << endl;
		// be careful with the cols and rows
		for (int i = 0; i < vv.cols(); ++i) {
			for (int j = 0; j < vv.rows(); ++j) {
				V(cnt++, k) = vv(j, i);
			}
		}
		cout << "cnt = " << file << endl;
		k++;
	}
	common::write_matrix_binary_to_file(filename, V);
	cout << "ok" << endl;
}

// Save all faces in the shape of (3, |F|) in binary
void saveBinFaces(const char *filename, const string &path, const vector<string> &files) {
	Eigen::Matrix3Xd vv;
	Eigen::Matrix3Xi ff;
	common::read_obj(path + files[0], vv, ff);
	common::write_matrix_binary_to_file(filename, ff);
	cout << "ok" << endl;
}


// Save 1-based Neibourhood vertex of each vertex..
// You should be careful about the index of vertex.
// The vertex in model is 0-based while it saved with 1-based
// It means i-1 cols is the neibour of ith vertex in 1-based

void saveNeighbor() {
	string AVEobj = "./data/AVE.obj";
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj(AVEobj, V, F);
	Surface_mesh mesh;
	build_mesh(V, F, mesh);

	// init out with -1
	Eigen::MatrixXi out;
	out.resize(VERTS, 11);
	for (int i = 0; i < VERTS; ++i) {
		for (int j = 0; j < 11; ++j) {
			out(i, j) = 0;
		}
	}


	//get neighbor in the mesh
	Surface_mesh::Vertex_iterator vit;
	Surface_mesh::Vertex_around_vertex_circulator vc, vc_end;

	int total = 0, mmax = 0;
	vit = mesh.vertices_begin();
	do {

		vc = mesh.vertices(*vit);
		vc_end = vc;
		//std::cout << (*vit).idx() << " ";
		int cnt = 0;
		do {
			//std::cout << (*vc).idx() << " ";
			out(total, cnt) = (*vc).idx() + 1;
			cnt++;
		} while (++vc != vc_end);
		//std::cout << "\n";

		total++;
		mmax = max(mmax, cnt);
	} while (++vit != mesh.vertices_end());
	cout << "total = " << total << " max = " << mmax << endl;

	cout << "(" << out.rows() << ", " << out.cols() << ")" << endl;
	//save as binary
	common::write_matrix_binary_to_file("./data/neighbor", out);
}

double laplacianCotanWeight(const Surface_mesh &mesh,
	const int i, const int j)
{
	return 1.0;
}

void calcFeature(const Eigen::Matrix3Xd &V,
	const Eigen::Matrix3Xi &F,
	Eigen::MatrixXd &feature)
{
	Surface_mesh mesh;
	build_mesh(V, F, mesh);

	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("./data/neighbor", N);

	Eigen::Matrix3Xd AVE;
	Eigen::Matrix3Xi AVF;
	common::read_obj("./data/AVE.obj", AVE, AVF);

	for (int i = 0; i < VERTS; ++i) {
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		GRBVar *Var = model.addVars(9, GRB_CONTINUOUS);
		GRBQuadExpr obj = 0;
		for (int j = 0; j < 11; ++j) {
			if (N(i, j) == 0) break;
			int k = N(i, j) - 1;
			double cij = laplacianCotanWeight(mesh, i, k);
			for (int n = 0; n < 3; ++n) {
				GRBLinExpr tmp = V(n, i) - V(n, k);
				for (int m = 0; m < 3; ++m) {
					tmp -= Var[3 * n + m] * (AVE(m, i) - AVE(m, k));
				}
				obj += tmp * tmp * cij;
			}
		}
		model.setObjective(obj);
		model.optimize();

		cout << "Matrix i == " << i << endl;
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				cout << Var[i * 3 + j].get(GRB_DoubleAttr_X) << "\t";
			}
			cout << endl;
		}
	}

}

void saveFeature(const Eigen::MatrixXd &V,
	const Eigen::Matrix3Xi &F)
{
	Eigen::MatrixXd feature;
	feature.resize(9 * VERTS, V.cols());
	for (int i = 0; i < V.cols(); ++i) {
		Eigen::MatrixXd tmp = V.col(i);
		tmp.resize(3, VERTS);
		calcFeature(tmp, F, feature);
		break;
	}
	//common::write_matrix_binary_to_file("./data/feature", feature);
}

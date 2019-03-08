// DataProcess.cpp : 
//
#include "mesh_io.h"
#include "binary_io.h"
#include "vtk.h"
#include "tri_mesh.h"
#include <Eigen/Dense>
#include <iostream>
#include <io.h>
#include <direct.h>
#include <surface_mesh/Surface_mesh.h>
#include <geodesic/geodesic_algorithm_dijkstra.h>
#include <geodesic/geodesic_algorithm_exact.h>

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
void calcAverage(Eigen::MatrixXd V, Eigen::Matrix3Xi F) {
	Eigen::MatrixXd newV;
	newV = V.rowwise().mean();
	newV.resize(3, VERTS);
		
	common::save_obj("./data/AVE.obj", newV, F);
	//cout << "ojk";
}

void getDeformation(Eigen::VectorXd ave, Eigen::MatrixXd &V) {
	Eigen::MatrixXd out;
	out.resize(V.rows(), V.cols() * 3);

	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("./data/neighbor", N);



	for (int i = 0; i < V.rows(); ++i) {
		for (int j = 0; j < VERTS; ++j) {
		}
	}

	cout << out << endl;
	cout << out.rows() << " " << out.cols() << endl;
	//common::write_matrix_binary_to_file("./affinematrix", out);
}



void saveCotWeight(Eigen::Matrix3d V, Eigen::Matrix3i F) {
	geodesic::Mesh mesh;
	vector<unsigned> ff;

	for (int i = 0; i < F.cols(); ++i) {
		if (i == 5014 || i == 24932) continue;
		for (int j = 0; j < F.rows(); ++j) {
			ff.push_back(F(j, i));
		}
	}
	vector<double> vv(V.data(), V.data() + V.rows()*V.cols());

	mesh.initialize_mesh_data(vv, ff);

	cout << "Mesh initialized ! " << endl;

	for (int i = 0; i < VERTS; ++i) {
		for (int j = i + 1; j < VERTS; ++i) {

		}
	}
}

// Save all verts in the shape of 3x|V| in binary. 
void saveBinVerts(string &filename, string &path, vector<string> files) {
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
	common::write_matrix_binary_to_file(filename.data(), V);
	cout << "ok" << endl;
}

void saveBinFaces(string &filename, string &path, vector<string> files) {
	Eigen::Matrix3Xd vv;
	Eigen::Matrix3Xi ff;
	common::read_obj(path + files[0], vv, ff);
	common::write_matrix_binary_to_file(filename.data(), ff);
	cout << "ok" << endl;
}


// Save 1-based Neibourhood vertex of each vertex..
// You should be careful about the index of vertex.
// The vertex in model is 0-based while it saved with 1-based
// It means i-1 cols is the neibour of ith vertex in 1-based

void calcNeighbor() {
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


void calcAffineMatrix() {
	string AVEobj = "C:/Users/amomorning/dataset/SPRING_MALE/AVE.obj";
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj(AVEobj, V, F);

	Eigen::VectorXd ave;
	ave.resize(VERTS * 3);
	int cnt = 0;
	for (int i = 0; i < V.rows(); ++i) {
		for (int j = 0; j < V.cols(); ++j) {
			ave(cnt++) = V(i, j);
		}
	}

	Eigen::MatrixXd vv;
	common::read_matrix_binary_from_file("./V", vv);
	std::cout << vv.rows() << " " << vv.cols() << std::endl;


	//Eigen::VectorXd row = vv.row(0);

	//cout << row.rows() << endl;
}



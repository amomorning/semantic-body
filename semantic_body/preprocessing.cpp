// DataProcess.cpp : 
//
#include "mesh_io.h"
#include "measure.h"
#include "binary_io.h"
#include "vtk.h"
#include "tri_mesh.h"
#include "preprocessing.h"
#include <iostream>
#include <io.h>
#include <direct.h>
#include <surface_mesh/Surface_mesh.h>
#include <geodesic/geodesic_algorithm_dijkstra.h>
#include <geodesic/geodesic_algorithm_exact.h>

#include <gurobi_c++.h>

#define VERTS 12500
#define EPS 1e-6

using namespace std;
using namespace surface_mesh;


struct node {
	int x, y;
	double t;
};

// Return a vector<string> of files in input cateloge direction.
// The files will be sorted.

vector<string> getFiles(const string &cate_dir)
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

void writeVTK(const char * filename, Eigen::Matrix3Xd &V) {
	std::vector<double> nodes;
	std::vector<int> lines;
	std::ofstream os(filename);

	int n = V.cols();
	for (int i = 0; i < n; ++i) {
		nodes.push_back(V(0, i));
		nodes.push_back(V(1, i));
		nodes.push_back(V(2, i));
		
		if (i) {
			lines.push_back(i - 1);
			lines.push_back(i);
		}
	}
	line2vtk(os, nodes.data(), nodes.size() / 3, lines.data(), lines.size() / 2);
	os.close();
	return;
}

// Todo
void saveExactBruteForce(const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F) {

}

void saveExact(const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F)
{
	Eigen::MatrixXd ret;
	measure measure;

	ret.resize(measure.len(), V.cols());


	for (int k = 0; k < measure.len(); ++k) {
		std::vector<node> path;
		ifstream in("./data/path/" + measure.SemanticLable[k], ios::binary);
		int len;
		in.read((char*)(&len), sizeof(int));
		for (int i = 0; i < len; ++i) {
			int x, y;
			double t;
			in.read((char*)(&x), sizeof(int));
			in.read((char*)(&y), sizeof(int));
			in.read((char*)(&t), sizeof(double));
			path.push_back({ x, y, t });
			//cout << x << " " << y << " " << t << endl;
		}
		in.close();
		
		// Todo: check difference between measured Exact and preserved calculation 
		// or just visualize it.
		for (int i = 0; i < V.cols(); ++i) {
			Eigen::MatrixXd tmp = V.col(i);
			tmp.resize(3, VERTS);
			
			Eigen::Matrix3Xd VV;
			VV.resize(3, path.size());
			int cnt = 0;
			for (auto u : path) {
				Eigen::Vector3d v0 = tmp.col(u.x);
				Eigen::Vector3d v1 = tmp.col(u.y);
				
				VV.col(cnt++) = v0 + (v1 - v0)*u.t;
			}

			//string name = "./checkVTK/" + measure.SemanticLable[k] + ".vtk";
			//writeVTK(name.c_str(), VV);
			//break;

			double res = 0;
			for (int i = 0; i < VV.cols()-1; i += 1) {
				res += (VV.col(i) - VV.col(i + 1)).norm();
			}
			ret(k, i) = res * 0.5;
		}
	}
	cout << ret.block(0, 0, 10, 10) << endl;
	Eigen::MatrixXd dijk;
	common::read_matrix_binary_from_file("./data/dijkstra", dijk);

	cout << endl << dijk.block(0, 0, 10, 10) << endl;
	common::write_matrix_binary_to_file("./data/exact", ret);
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

void saveVertsOffset(const Eigen::MatrixXd &V)
{
	Eigen::MatrixXd ret = V;
	Eigen::VectorXd ave = V.rowwise().mean();
	for (int i = 0; i < ret.cols(); ++i) {
		ret.col(i) -= ave;
	}
	common::write_matrix_binary_to_file("./data/Verts", ret);
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

void laplacianCotanWeight(const Surface_mesh &mesh,
	Eigen::SparseMatrix<double> &cotan)
{
	Surface_mesh::Face_iterator fit;
	auto points = mesh.get_vertex_property<Point>("v:point");

	fit = mesh.faces_begin();
	std::vector<Eigen::Triplet<double> > tri;
	do {
		Surface_mesh::Vertex_around_face_circulator vf = mesh.vertices(*fit);
		Point p[3];
		int id[3];
		double cot[3];
		for (int i = 0; i < 3; ++i, ++vf) {
			p[i] = points[*vf];
			id[i] = (*vf).idx();
		}

		for (int i = 0; i < 3; ++i) {
			int j = (i + 1) % 3, k = (j + 1) % 3;
			cot[i] = dot(p[j] - p[i], p[k] - p[i]) /
				norm(cross(p[j] - p[i], p[k] - p[i]));

			// too slow....4s.
			//cotan.coeffRef(id[j], id[k]) += 0.5*cot;
			//cotan.coeffRef(id[k], id[j]) += 0.5*cot;

			// 0.1s using setFromTriplets function;
			tri.push_back({ id[j], id[k], -0.5*cot[i] });
			tri.push_back({ id[k], id[j], -0.5*cot[i] });
		}

		for (int i = 0; i < 3; ++i) {
			tri.push_back({ id[i], id[i], 0.5*(cot[(i + 1) % 3], cot[(i + 2) % 3]) });
		}

	} while (++fit != mesh.faces_end());
	cotan.setFromTriplets(tri.begin(), tri.end());
}


//Calculate feature representation by closed-form expression
//arg min \sum cij |p(m,i)-p(m,j) - T(m,i)(p(1,i)-p(1,j))|^2
void calcFeature(const Eigen::Matrix3Xd &V,
	const Eigen::Matrix3Xi &F,
	const Eigen::Matrix3Xd &AVE,
	Eigen::MatrixXd &feature, int u)
{
	Surface_mesh mesh;
	build_mesh(V, F, mesh);

	Eigen::SparseMatrix<double> cotan(VERTS, VERTS);
	laplacianCotanWeight(mesh, cotan);

	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("./data/neighbor", N);
	//cout << N.rows() << " " << N.cols() << endl;

	int cnt = 0;
	for (int i = 0; i < VERTS; ++i) {
		Eigen::Matrix3d res = Eigen::Matrix3d::Zero();
		double co = 0;
		for (int j = 0; j < 11; ++j) {
			//cout << N(i, j) << endl;;
			if (N(i, j) == 0) break;
			int k = N(i, j) - 1;

			//cout << i << " " << k << endl;
			double cij = -cotan.coeff(i, k);
			Eigen::Vector3d a = V.col(k) - V.col(i);
			Eigen::Vector3d b = AVE.col(k) - AVE.col(i);

			co += cij * b.norm();
			res += a * b.transpose();
		}
		res /= co;
		//cout << cnt << endl;
		//cout << res << endl << endl;;

		//todo: RS decompositon than save....

		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				
				feature(u, cnt++) = res(j, k);
			}
		}
	}
}
//Deprecated
//It's unnecessary using Gurobi
void calcFeatureGurobi(const Eigen::Matrix3Xd &V,
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

	Eigen::SparseMatrix<double> cotan(VERTS, VERTS);
	laplacianCotanWeight(mesh, cotan);
	puts("ok");
	for (int i = 0; i < VERTS; ++i) {
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		GRBVar *Var = model.addVars(9, GRB_CONTINUOUS);
		GRBQuadExpr obj = 0;
		for (int j = 0; j < 11; ++j) {
			if (N(i, j) == 0) break;
			int k = N(i, j) - 1;
			double cij = fabs(cotan.coeff(i, k));
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
	feature.resize(V.cols(), 9 * VERTS);

	Eigen::Matrix3Xd AVE;
	Eigen::Matrix3Xi AVF;
	common::read_obj("./data/AVE.obj", AVE, AVF);

	for (int i = 0; i < V.cols(); ++i) {
		cout << "*********************************  " << i << endl;
		Eigen::MatrixXd tmp = V.col(i);
		tmp.resize(3, VERTS);
		calcFeature(tmp, F, AVE, feature, i);
		//calcFeatureGurobi(tmp, F, feature);
	}
	common::write_matrix_binary_to_file("./data/feature", feature);
}

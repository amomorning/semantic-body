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
#include <sophus/so3.hpp>

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
void calcAverage(const char * filename, const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F) {
	Eigen::MatrixXd newV = V.rowwise().mean();
	
	newV.resize(3, VERTS);
	common::save_obj(filename, newV, F);
	//cout << "ojk";
}

void calcDeltaVerts(const char* filename, Eigen::MatrixXd &V) 
{
	Eigen::MatrixXd ave = V.rowwise().mean();

	for (int i = 0; i < V.cols(); ++i) {
		V.col(i) -= ave;
	}
	common::write_matrix_binary_to_file(filename, V);
}


void saveDijkstra(const char * filename, const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F)
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
	cout << "dijkstra measurement saved" << endl;
	cout << ret.rows() << " " << ret.cols() << endl;
	common::write_matrix_binary_to_file(filename, ret);
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
void saveExact(const char* filename, const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F) {
	Eigen::MatrixXd ret;
	measure measure;

	ret.resize(measure.len(), V.cols());

	for (int i = 0; i < V.cols(); ++i) {
		Eigen::MatrixXd tmp = V.col(i);
		tmp.resize(3, VERTS);
		measure.calcExact(tmp, F, false);
		measure.saveParam(ret, i);
	}

	cout << "Exact measurement saved" << endl;
	cout << ret.rows() << " " << ret.cols() << endl;
	common::write_matrix_binary_to_file(filename, ret);
}

void saveRoughExact(const char* filename, const Eigen::MatrixXd &V, Eigen::Matrix3Xi &F)
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
	//cout << ret.block(0, 0, 10, 10) << endl;
	//Eigen::MatrixXd dijk;
	//common::read_matrix_binary_from_file("./data/dijkstra", dijk);

	//cout << endl << dijk.block(0, 0, 10, 10) << endl;
	cout << "roughExact saved" << endl;
	cout << ret.rows() << " " << ret.cols() << endl;
	common::write_matrix_binary_to_file(filename, ret);
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

void saveNeighbor(const char * filename) {
	Surface_mesh mesh;
	mesh.read("./data/AVE.obj");

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
	common::write_matrix_binary_to_file(filename, out);
}


void calc_cot_angles(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
	Eigen::Matrix3Xd &cot_angles)
{
	cot_angles.resize(3, F.cols());
	for (size_t j = 0; j < F.cols(); ++j) {
		const Eigen::Vector3i &fv = F.col(j);
		for (size_t vi = 0; vi < 3; ++vi) {
			const Eigen::VectorXd &p0 = V.col(fv[vi]);
			const Eigen::VectorXd &p1 = V.col(fv[(vi + 1) % 3]);
			const Eigen::VectorXd &p2 = V.col(fv[(vi + 2) % 3]);
			const double angle = std::acos(std::max(-1.0,
				std::min(1.0, (p1 - p0).normalized().dot((p2 - p0).normalized()))));
			cot_angles(vi, j) = 1.0 / std::tan(angle);
		}
	}
}

int calc_cot_laplace(const Eigen::MatrixXd &V, const Eigen::Matrix3Xi &F,
	Eigen::SparseMatrix<double> &L)
{
	Eigen::Matrix3Xd cot_angles;
	calc_cot_angles(V, F, cot_angles);
	std::vector<Eigen::Triplet<double>> triple;
	triple.reserve(F.cols() * 9);
	for (size_t j = 0; j < F.cols(); ++j) {
		const Eigen::Vector3i &fv = F.col(j);
		const Eigen::Vector3d &ca = cot_angles.col(j);
		for (size_t vi = 0; vi < 3; ++vi) {
			const size_t j1 = (vi + 1) % 3;
			const size_t j2 = (vi + 2) % 3;
			const int fv0 = fv[vi];
			const int fv1 = fv[j1];
			const int fv2 = fv[j2];
			triple.push_back(Eigen::Triplet<double>(fv0, fv0, ca[j1] + ca[j2]));
			triple.push_back(Eigen::Triplet<double>(fv0, fv1, -ca[j2]));
			triple.push_back(Eigen::Triplet<double>(fv0, fv2, -ca[j1]));
		}
	}
	L.resize(V.cols(), V.cols());
	L.setFromTriplets(triple.begin(), triple.end());
	return 1;
}


void calcFeature(const Eigen::SparseMatrix<double> &W,
	const Eigen::MatrixXi &N,
	const Eigen::Matrix3Xd &Va,
	Eigen::MatrixXd V,
	Eigen::MatrixXd &F, int u)
{
	V.resize(3, VERTS);
	int cnt = 0;

	int tot = 0;
	for (int i = 0; i < VERTS; ++i) {
		Eigen::Matrix3d T = Eigen::Matrix3d::Zero();
		double co = 0;

		for (int k = 0; k < 11; ++k) {
			int j = N(i, k) - 1;
			if (j < 0) continue;

			double wij = W.coeff(i, j);

			Eigen::Vector3d a = V.col(j) - V.col(i);
			Eigen::Vector3d b = Va.col(j) - Va.col(i);

			co += wij * b.squaredNorm();
			T += wij * a*b.transpose();
		}

		T /= co;
		
		//double sum = 0;
		//for (int j = 0; j < 11; ++j) {
		//	int k = N(i, j) - 1;
		//	if (k < 0 || k == i) continue;

		//	//if(!i) cout << i << " " << k << endl;
		//	double cij = W.coeff(i, k);
		//	Eigen::Vector3d a = V.col(k) - V.col(i);
		//	Eigen::Vector3d b = Va.col(k) - Va.col(i);

		//	sum += cij * (a - T * b).squaredNorm();
		//}
		//cout << "sum == " << sum << endl;
		//if (i < 5) cout << T << endl;

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(T, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::Matrix3d U = svd.matrixU();
		Eigen::Matrix3d V = svd.matrixV();
		Eigen::Matrix3d R = U * V.transpose(), S = V * U.transpose()*T;

		if (i == 0 || i == 1923 || i == 2822 || i == 5018) {
			cout << "this is id " << i << endl;
			cout << "co == " << co << endl;
			cout << "T == " << T << endl;
			cout << "R == " << R << endl;
			cout << "S == " << S << endl;
		}

		if (R.determinant() < 0) {
			cout << " i == " << i << ": " << R.determinant() << endl;
			tot++;
		}

		//try {
		//	Sophus::SO3d SO3_R(R);
		//	Eigen::Vector3d so3 = SO3_R.log();

		//	for (int j = 0; j < 3; ++j) {
		//		F(cnt++, u) = so3[j];
		//	}
		//}
		//catch (exception e) {

		//}

		//for (int j = 0; j < 3; ++j) {
		//	for (int k = 0; k < 3; ++k) {
		//		F(cnt++, u) = S(k, j);
		//	}
		//}
	}
	cout << "tot == " << tot << endl;
}


void saveFeature(const char* filename, 
	const Eigen::MatrixXd &V,
	 Eigen::Matrix3Xi &F)
{
	Eigen::MatrixXd feature;
	feature.resize( 12 * VERTS, V.cols());

	Eigen::Matrix3Xd Va;
	common::read_obj("./data/AVE.obj", Va, F);
	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("./data/N", N);


	Eigen::SparseMatrix<double> W(VERTS, VERTS);
	calc_cot_laplace(Va, F, W);

	for (int i = 0; i < V.cols(); ++i) {
		cout << "******** feature num = " << i << endl;
		calcFeature(W, N, Va, V.col(i), feature, i);
		break;
	}
	//common::write_matrix_binary_to_file(filename, feature);
	cout << "Feature is saved!!" << endl;
}


void calcFeature(Eigen::MatrixXd V,
	const Eigen::MatrixXd &Va,
	const Eigen::Matrix3Xi &F,
	Eigen::MatrixXd &feature, int u, double ave)
{
	V.resize(3, VERTS);
	Eigen::Vector3d v[4], norm;
	int tot = 0;
	for (int i = 0; i < F.cols(); ++i) {

		//cout << "ok" << endl;
		Eigen::Matrix3d t, ta;
		for (int j = 0; j < 3; ++j) {
			v[j] = V.col(F(j, i));
		}
		norm = (v[1] - v[0]).cross(v[2] - v[0]);
		v[3] = (v[0] + ave*norm / norm.norm());
		
		for (int j = 0; j < 3; ++j) {
			t.col(j) = v[j + 1] - v[0];
		}

		for (int j = 0; j < 3; ++j) {
			v[j] = Va.col(F(j, i));
		}
		norm = (v[1] - v[0]).cross(v[2] - v[0]);
		v[3] =  (v[0] + ave*norm / norm.norm());
		
		for (int j = 0; j < 3; ++j) {
			ta.col(j) = v[j + 1] - v[0];
		}

		Eigen::Matrix3d T = t * ta.inverse();
		//cout << "T == " << T.determinant() << endl;
		//cout << T << endl;

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(T, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::Matrix3d U = svd.matrixU();
		Eigen::Matrix3d V = svd.matrixV();
		Eigen::Matrix3d R = U * V.transpose(), S = V * U.transpose()*T;

		Sophus::SO3d SO3_R(R);

		Eigen::Vector3d so3 = SO3_R.log();
		for (int j = 0; j < 3; ++j) {
			feature(tot++, u) = so3[j];
		}

		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				feature(tot++, u) = S(j, k);
			}
		}
		//cout << "R == " << R.determinant() << endl;
		//cout << R << endl;
		//cout << "S == " << endl;
		//cout << S << endl;
	}
	cout << "tot == " << tot << endl;
}

void saveFaceFeature(const char* filename,
	const Eigen::MatrixXd &V,
	Eigen::Matrix3Xi &F)
{
	Eigen::MatrixXd feature;
	feature.resize(24 * VERTS, V.cols());

	//Eigen::Matrix<int, 3, 2*VERTS> N;
	//cout << "N rc = " << N.rows() << " " << N.cols() << endl;

	Eigen::Matrix3Xd Va;
	common::read_obj("./data/AVE.obj", Va, F);

	Surface_mesh mesh;
	build_mesh(Va, F, mesh);

	double ave = 0;
	for (const auto & e : mesh.edges()) {
		ave += mesh.edge_length(e);
	}
	ave /= mesh.n_edges();

	//for (const auto &fa : mesh.faces()) {
	//	int cnt = 0;
	//	for (const auto &he : mesh.halfedges(fa)) {
	//		Surface_mesh::Halfedge x = mesh.opposite_halfedge(he);
	//		N(mesh.face(x).idx(), cnt++) = fa.idx();
	//		N(fa.idx(), cnt++);
	//	}
	//}

	for (int i = 0; i < V.cols(); ++i) {
		cout << "feature num = " << i << endl;
		calcFeature(V.col(i), Va, F, feature, i, ave);
	}

	common::write_matrix_binary_to_file(filename, feature);

}
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <binary_io.h>
#include <tri_mesh.h>
#include <surface_mesh/Surface_mesh.h>
#include "mesh_io.h"

using namespace surface_mesh;
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


void calcFeature(const char * filename, 
	const Eigen::SparseMatrix<double> &cotan, Eigen::Matrix3Xd V, 
	const Eigen::Matrix3Xi &F, Eigen::MatrixXd feature) 
{
	cout << feature.rows() << " " << feature.cols() << endl;
	Eigen::Matrix3Xd V0 = V;
	cout << V0.rows() << " " << V0.cols() << endl;
	feature.resize(9, 12500);
	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("../data/neighbor", N);

	cout << N.rows() << " " << N.cols() << endl;
	for (int i = 0; i < 12500; ++i) {
		Eigen::MatrixXd T = feature.col(i);
		T.resize(3, 3);

		double co = 0;

		Eigen::Vector3d b;

		for (int j = 0; j < 11; ++j) {
			int k = N(i, j) - 1;
			if (k < 0) break;
			
			double cij = cotan.coeff(i, k);
			double cji = cotan.coeff(k, i);

			cout <<"cij == " << cij << " " << cji << endl;
			b = V.col(k) - V.col(i);
			co += cij * b.squaredNorm();
		}
		//cout << T << endl;
		T /= co;
		for (int j = 0; j < 3; ++j) {
			V0(j, i) = T(j, j) / b[j];
		}
	}
	common::save_obj(filename, V0, F);
	cout << filename << " saved!!" << endl;
}

//Generate data from new feature matrix....
void recoverFromFeature(const char * filename, int total) {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("../data/AVE.obj", V, F);

	Eigen::MatrixXd feature;

	common::read_matrix_binary_from_file(filename, feature);

	Eigen::MatrixXd ff = feature.transpose();
	cout << ff.rows() << " " << ff.cols() << endl;
	//std::ifstream in(filename);
	//feature.resize(112500, total);
	//for(int i = 0; i < total; ++ i) {
	//	for(int j = 0; j < 112500; ++ j) {
	//		in >> feature(j, i);
	//	}
	//}
	//in.close();

	Surface_mesh mesh;
	build_mesh(V, F, mesh);
	
	Eigen::SparseMatrix<double> cotan(12500, 12500);
	laplacianCotanWeight(mesh, cotan);

	for (int i = 0; i < total; ++i) {
		std::string name = "../data/augmentation/" + std::to_string(i) + ".obj";
		calcFeature(name.c_str(), cotan, V, F, ff.col(i));
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
	
	recoverFromFeature("../data/feature", 1);
	getchar();
}
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <binary_io.h>
#include <tri_mesh.h>
#include <surface_mesh/Surface_mesh.h>
#include "mesh_io.h"

using namespace std;
using namespace surface_mesh;

void laplacianCotanWeight(const Surface_mesh &mesh,
	Eigen::SparseMatrix<double> &cotan)
{
	std::vector<Eigen::Triplet<double> > tri;
	for (const auto &fit : mesh.faces())
	{
		int i = 0;
		Point p[3];
		int id[3];
		double cot[3];
		for (const auto &vit : mesh.vertices(fit))
		{
			p[i] = mesh.position(vit);
			id[i] = vit.idx();
			++i;
		}

		for (int i = 0; i < 3; ++i) {
			int j = (i + 1) % 3, k = (j + 1) % 3;
			cot[i] = 0.5*dot(p[j] - p[i], p[k] - p[i]) /
				norm(cross(p[j] - p[i], p[k] - p[i]));

			tri.push_back({ id[j], id[k], cot[i] });
			tri.push_back({ id[k], id[j], cot[i] });
		}

		for (int i = 0; i < 3; ++i) {
			tri.push_back({ id[i], id[i], -(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
		}
	}
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
	feature.resize(9, 12500);
	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("../data/neighbor", N);

	Eigen::MatrixXd b(12500, 3);
	for (int i = 0; i < 12500; ++i) {
		for (int j = 0; j < 3; ++j) {
			b(i, j) = 0;
		}
	}
	cout << N.rows() << " " << N.cols() << endl;
	std::vector<Eigen::Triplet<double>> tri;
	for (int i = 0; i < 12500; ++i) {
		Eigen::MatrixXd Ti = feature.col(i);
		Ti.resize(3, 3);
		for (int j = 0; j < 11; ++j) {
			int k = N(i, j)-1;
			if (k < 0) break;
			Eigen::MatrixXd Tj = feature.col(k);
			Tj.resize(3, 3);
			double cij = cotan.coeff(i, k);
			b.row(i) += 0.5*cij*(Ti + Tj)*(V.col(k) - V.col(i));
		}
	}
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(cotan);
	Eigen::MatrixXd V0 = solver.solve(b);
	cout << "now v looks like:\n";
	cout << V0.block<10, 3>(0, 0) << endl;
	common::save_obj(filename, V0.transpose(), F);
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

	Surface_mesh mesh;
	build_mesh(V, F, mesh);
	
	Eigen::SparseMatrix<double> cotan(12500, 12500);
	laplacianCotanWeight(mesh, cotan);

	for (int i = 0; i < total; ++i) {
		std::string name = "../data/augmentation/" + std::to_string(i) + ".obj";
		calcFeature(name.c_str(), cotan, V, F, ff.col(i));
	}
}

int main() {
	//Eigen::Matrix3Xd V;
	//Eigen::Matrix3Xi F;
	//common::read_obj("../data/AVE.obj", V, F);
	//readNewPoint("../data/New.txt", V);
	//
	//common::save_obj("../data/NEW.obj", V, F);
	
	recoverFromFeature("../data/feature", 10);
	getchar();
}
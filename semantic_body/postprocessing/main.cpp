#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <binary_io.h>
#include <tri_mesh.h>
#include <surface_mesh/Surface_mesh.h>
#include "mesh_io.h"

using namespace std;
using namespace surface_mesh;

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


//void laplacianCotanWeight(const Surface_mesh &mesh,
//	Eigen::SparseMatrix<double> &cotan)
//{
//	std::vector<Eigen::Triplet<double> > tri;
//	for (const auto &fit : mesh.faces())
//	{
//		int i = 0;
//		Point p[3];
//		int id[3];
//		double cot[3];
//		for (const auto &vit : mesh.vertices(fit))
//		{
//			p[i] = mesh.position(vit);
//			id[i] = vit.idx();
//			++i;
//		}
//
//		for (int i = 0; i < 3; ++i) {
//			int j = (i + 1) % 3, k = (j + 1) % 3;
//			cot[i] = 0.5*dot(p[j] - p[i], p[k] - p[i]) /
//				norm(cross(p[j] - p[i], p[k] - p[i]));
//
//			tri.push_back({ id[j], id[k], cot[i] });
//			tri.push_back({ id[k], id[j], cot[i] });
//		}
//
//		for (int i = 0; i < 3; ++i) {
//			tri.push_back({ id[i], id[i], -(cot[(i + 1) % 3] + cot[(i + 2) % 3]) });
//		}
//	}
//	cotan.setFromTriplets(tri.begin(), tri.end());
//}


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
	Eigen::SparseMatrix<double> &A, Eigen::Matrix3Xd V, 
	const Eigen::Matrix3Xi &F, Eigen::MatrixXd feature) 
{
	cout << feature.rows() << " " << feature.cols() << endl;
	feature.resize(9, 12500);
	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("../data/neighbor", N);

	Eigen::MatrixXd b(12500, 3);
	b.setZero();

	//cout << N.rows() << " " << N.cols() << endl;

	for (int i = 0; i < 12500; ++i) {
		Eigen::Map<Eigen::Matrix3d> Ti(feature.col(i).data());
		for (int j = 0; j < 11; ++j) {
			int k = N(i, j)-1;
			if (k < 0 || k == i) continue;
			Eigen::Map<Eigen::Matrix3d> Tj(feature.col(k).data());
			double cij = A.coeff(i, k);
			b.row(i) += 0.5*cij*(Ti + Tj)*(V.col(k) - V.col(i));
		}
	}

	double w = 1e6;
	for (int i = 0; i < 1; ++i) {
		A.coeffRef(i, i) += w;
		b.row(i) += w * V.col(i).transpose();
	}
	
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(A);
	if (solver.info() == Eigen::Success) {
		Eigen::MatrixXd V0 = solver.solve(b).transpose();

		common::save_obj(filename, V0, F);
		cout << filename << " saved!!" << endl;
	}
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
	
	Eigen::SparseMatrix<double> L;
	calc_cot_laplace(V, F, L);

	for (int i = 0; i < total; ++i) {
		std::string name = "../data/augmentation/" + std::to_string(i) + ".obj";
		calcFeature(name.c_str(), L, V, F, ff.col(i));
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
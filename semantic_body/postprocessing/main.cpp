#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <binary_io.h>
#include <tri_mesh.h>
#include <surface_mesh/Surface_mesh.h>
#include "mesh_io.h"

using namespace std;
using namespace surface_mesh;

const int wight = 1e6;
const int idx = 0;

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
	const Eigen::SparseMatrix<double> &L,
	const Eigen::SparseMatrix<double> &A, 
	const Eigen::Matrix3Xd &V,
	const Eigen::Matrix3Xi &F,
	const Eigen::MatrixXi &N, 
	Eigen::VectorXd feature)
{
	// fill b
	Eigen::MatrixXd b(37500, 3);
	for (int i = 0, t = 0; i < 12500; ++i, t = 0) 
		for (int j = 0; j < 3; ++j) 
			for (int k = 0; k < 3; ++k) 
				b(i*3+j, k) = feature[i * 9 + (t++)];
	
	for (int i = 0; i < 12500; ++i) {
		double co = 0;
		for (int j = 0; j < 11; ++j) {
			int k = N(i, j) - 1;
			if (k < 0) break;

			double cij = L.coeff(i, k);
			Eigen::Vector3d u = V.col(k) - V.col(i);
			co += cij * u.squaredNorm();
		}
		b.block<3, 3>(i * 3, 0) *= co;
	}

	Eigen::Matrix3d w;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			w(i, j) = V(j, 0)*wight;
		}
	}

	b.block<3, 3>(idx*3, 0) += w;

	//cout << b.block<20, 3>(0, 0) << endl;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(A.transpose()*A);
	if (solver.info() == Eigen::Success) {
		Eigen::MatrixXd Vout = solver.solve(A.transpose()*b).transpose();
		common::save_obj(filename, Vout, F);
		cout << filename << " saved !!! " << endl;
	}
}

void calcFeatureRS(const char * filename,
	const Eigen::SparseMatrix<double> &L,
	const Eigen::SparseMatrix<double> &A,
	const Eigen::Matrix3Xd &V,
	const Eigen::Matrix3Xi &F,
	const Eigen::MatrixXi &N,
	Eigen::VectorXd feature)
{
	// fill b
	Eigen::MatrixXd b(37500, 3);
	for (int i = 0, t = 0; i < 12500; ++i, t = 0) {
		Eigen::Matrix3d R, S, T;

		for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
				R(j, k) = feature[i * 18 + (t++)];

		for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
				S(j, k) = feature[i * 18 + (t++)];
		//cout << t << endl;
		T = S*R;
		for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
				b(i * 3 + j, k) = T(j, k);
	}
	for (int i = 0; i < 12500; ++i) {
		double co = 0;
		for (int j = 0; j < 11; ++j) {
			int k = N(i, j) - 1;
			if (k < 0) break;

			double cij = L.coeff(i, k);
			Eigen::Vector3d u = V.col(k) - V.col(i);
			co += cij * u.squaredNorm();
		}
		b.block<3, 3>(i * 3, 0) *= co;
	}

	Eigen::Matrix3d w;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			w(i, j) = V(j, 0)*wight;
		}
	}

	b.block<3, 3>(idx * 3, 0) += w;

	//cout << b.block<20, 3>(0, 0) << endl;
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(A.transpose()*A);
	if (solver.info() == Eigen::Success) {
		Eigen::MatrixXd Vout = solver.solve(A.transpose()*b).transpose();
		common::save_obj(filename, Vout, F);
		cout << filename << " saved !!! " << endl;
	}
}

//void calcFeature(const char * filename, 
//	Eigen::SparseMatrix<double> &A, Eigen::Matrix3Xd V, 
//	const Eigen::Matrix3Xi &F, Eigen::MatrixXd feature) 
//{
//	cout << feature.rows() << " " << feature.cols() << endl;
//	feature.resize(9, 12500);
//	Eigen::MatrixXi N;
//	common::read_matrix_binary_from_file("../data/neighbor", N);
//
//	Eigen::MatrixXd b(12500, 3);
//	b.setZero();
//
//	//cout << N.rows() << " " << N.cols() << endl;
//
//	for (int i = 0; i < 12500; ++i) {
//		Eigen::Map<Eigen::Matrix3d> Ti(feature.col(i).data());
//		for (int j = 0; j < 11; ++j) {
//			int k = N(i, j)-1;
//			if (k < 0 || k == i) continue;
//			Eigen::Map<Eigen::Matrix3d> Tj(feature.col(k).data());
//			double cij = A.coeff(i, k);
//			b.row(i) += 0.5*cij*(Ti + Tj)*(V.col(k) - V.col(i));
//		}
//	}
//
//	double w = 1e6;
//	for (int i = 0; i < 1; ++i) {
//		A.coeffRef(i, i) += w;
//		b.row(i) += w * V.col(i).transpose();
//	}
//	
//	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(A);
//	if (solver.info() == Eigen::Success) {
//		Eigen::MatrixXd V0 = solver.solve(b).transpose();
//
//		common::save_obj(filename, V0, F);
//		cout << filename << " saved!!" << endl;
//	}
//}

//Generate data from new feature matrix....
void recoverFromFeature(Eigen::MatrixXd &feature) {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("../data/AVE.obj", V, F);

	//Eigen::MatrixXd feature;

	//common::read_matrix_binary_from_file(filename, feature);


	Eigen::SparseMatrix<double> L;
	calc_cot_laplace(V, F, L);
	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("../data/neighbor", N);
	// fill A
	std::vector<Eigen::Triplet<double>> tri;
	for (int i = 0; i < 12500; ++i) {
		Eigen::Vector3d tot = Eigen::Vector3d::Zero();
		for (int j = 0; j < 11; ++j) {
			int k = N(i, j) - 1;
			if (k < 0) break;

			double cij = L.coeff(i, k);
			Eigen::Vector3d u = V.col(k) - V.col(i);

			u *= cij;
			for (int t = 0; t < 3; ++t) {
				tri.push_back({ i * 3 + t, k, u[t] });
				tri.push_back({ i * 3 + t, i, -u[t] });
			}
		}
	}
	Eigen::SparseMatrix<double> A(37500, 12500);
	A.setFromTriplets(tri.begin(), tri.end());

	for (int i = 0; i < 3; ++i) {
		A.coeffRef(idx*3 + i, idx) += wight;
	}

	cout << "okokok" << endl;

	for (int i = 0; i < feature.rows(); ++i) {
		std::string name = "../data/test/" + std::to_string(i) + ".obj";
		calcFeatureRS(name.c_str(), L, A, V, F, N, feature.row(i));
	}
}


void readNewFeature(const char* filename, Eigen::MatrixXd &feature, int total) {
	feature.resize(total, 112500);
	std::ifstream is(filename);
	double x; 
	int cnt = 0;
	while (is >> x) {
		feature(0, cnt) = x;
		cnt++;
	}
	cout << "cnt == " << cnt << endl;
}


int main() {
	//Eigen::Matrix3Xd V;
	//Eigen::Matrix3Xi F;
	//common::read_obj("../data/AVE.obj", V, F);
	//readNewPoint("../data/New.txt", V);
	//
	//common::save_obj("../data/NEW.obj", V, F);
	
	Eigen::MatrixXd feature;
	//readNewFeature("../data/T.txt", feature, 1);
	
	common::read_matrix_binary_from_file("../data/featureRS", feature);
	//cout << feature.row(0) << endl;
	recoverFromFeature(feature);
	getchar();
}
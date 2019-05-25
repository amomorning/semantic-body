#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <binary_io.h>
#include <tri_mesh.h>
#include <surface_mesh/Surface_mesh.h>
#include "measure.h"
#include "mesh_io.h"
#include <sophus/so3.hpp>

using namespace std;
using namespace surface_mesh;

const int VERTS = 12500;
const int FACES = 25000;

const int wight = 1e6;
const int idx = 0;

struct node {
	int x, y;
	double t;
};


const string filepath = "C:/Users/amomorning/source/repos/semantic_body/semantic_body/postprocessing/";
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


void readNewPoint(const char * filename, Eigen::MatrixXd &V, const Eigen::Matrix3Xd &Va, int total) {
	V.resize(VERTS*3, total);
	std::ifstream is(filename, std::ios::binary);
	is.read((char*)V.data(), VERTS * 3 * total * sizeof(Eigen::MatrixXd::Scalar));

	Eigen::MatrixXd tmp = Va;
	tmp.resize(VERTS * 3, 1);
	
	for (int i = 0; i < V.cols(); ++i) {
		V.col(i) += tmp;
	}
	puts("ok");
}

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

	ret.setZero();

	for (int k = 0; k < measure.len(); ++k) {
		int tot = 1;
		if (k < 14) tot = 4;
		for (int j = 0; j < tot; ++j) {
			cout << "now " << measure.SemanticLable[k] << endl;
			std::vector<node> path;
			ifstream in(filepath + "../data/path/" + measure.SemanticLable[k] + to_string(j), ios::binary);
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

			cout << path.size() << endl;
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
				for (int tt = 0; tt < VV.cols() - 1; tt++) {
					res += (VV.col(tt) - VV.col(tt + 1)).norm();
				}

				ret(k, i) += res;
			}
		}
	}
	//cout << ret.block(0, 0, 10, 10) << endl;
	//Eigen::MatrixXd dijk;
	//common::read_matrix_binary_from_file("./data/dijkstra", dijk);

	//cout << endl << dijk.block(0, 0, 10, 10) << endl;
	cout << "roughExact saved" << endl;
	cout << ret.rows() << " " << ret.cols() << endl;
	cout << ret.col(0) << endl;
	common::write_matrix_binary_to_file(filename, ret);
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

void calcFeatureRS(const Eigen::SparseMatrix<double> &L,
	const Eigen::SparseMatrix<double> &A,
	const Eigen::Matrix3Xd &V, const Eigen::Matrix3Xi &F,
	const Eigen::MatrixXi &N, Eigen::VectorXd feature,
	Eigen::MatrixXd &newV, int u
)
{
	// fill b
	Eigen::MatrixXd b(37500, 3);
	for (int i = 0, t = 0; i < 12500; ++i, t = 0) {
		Eigen::Matrix3d Rt, St, T;

		for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
				Rt(j, k) = feature[i * 18 + (t++)];

		for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
				St(j, k) = feature[i * 18 + (t++)];
		//cout << t << endl;
		T = St*Rt;
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
		int cnt = 0;
		for (int j = 0; j < VERTS; ++j) {
			for (int i = 0; i < 3; ++i) {
				newV(cnt++, u) = Vout(i, j);
			}
		}
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
void recoverFromFeature(const char* filename, Eigen::MatrixXd &feature) {
	Eigen::Matrix3Xd V;
	Eigen::Matrix3Xi F;
	common::read_obj("../data/AVE.obj", V, F);

	//Eigen::MatrixXd feature;

	//common::read_matrix_binary_from_file(filename, feature);


	Eigen::SparseMatrix<double> L;
	calc_cot_laplace(V, F, L);
	Eigen::MatrixXi N;
	common::read_matrix_binary_from_file("../data/N", N);
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

	Eigen::MatrixXd newV;
	newV.resize(3 * VERTS, feature.cols());

	for (int i = 0; i < feature.cols(); ++i) {
		calcFeatureRS(L, A, V, F, N, feature.col(i), newV, i);
	}

	
	Eigen::MatrixXd tmp = newV.col(0);
	tmp.resize(3, VERTS);
	common::save_obj((filepath + "../../bodyViz/Assets/Models/TEMP.obj").c_str(), tmp, F);
	cout << "New obj saved!!" << endl;

	saveRoughExact(filename, newV, F);
}


void readNewFeature(const char* filename, Eigen::MatrixXd &feature, int total) {
	feature.resize(18*VERTS, total);
	std::ifstream is(filename, std::ios::binary);
	is.read((char*)feature.data(), 18 * VERTS*total * sizeof(Eigen::MatrixXd::Scalar));
}

void readNewFaceFeature(const char* filename, Eigen::MatrixXd &feature, int total) {
	feature.resize(12*2 * VERTS, total);
	std::ifstream is(filename, std::ios::binary);
	is.read((char*)feature.data(), 24 * VERTS*total * sizeof(Eigen::MatrixXd::Scalar));
}

void recoverFromVertex(const char* infile, const char* outfile) {
	Eigen::Matrix3Xd Va;
	Eigen::Matrix3Xi F;
	common::read_obj((filepath + "../data/AVE.obj").c_str(), Va, F);

	Eigen::MatrixXd V;
	readNewPoint(infile, V, Va, 111);

	saveRoughExact(outfile, V, F);
	
	Eigen::MatrixXd tmp = V.col(35);
	tmp.resize(3, VERTS);
	common::save_obj((filepath + "../../bodyViz/Assets/Models/TEMP.obj").c_str(), tmp, F);
	cout << "New obj saved!!" << endl;
}


void getVertexMeature(const char* infile, const char* outfile) {
	recoverFromVertex(infile, outfile);

}



Eigen::VectorXd calcFaceFeature(const Eigen::SparseMatrix<double> &A,
	const Eigen::MatrixXd &V, Eigen::MatrixXd feature,
	const Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> &solver)
{
	feature.resize(12, FACES);

	Eigen::MatrixXd b(FACES * 3, 3);

	for (int i = 0; i < FACES; ++i) {
		Eigen::Vector3d so3;
		Eigen::Matrix3d S;
		for (int j = 0; j < 3; ++j) so3[j] = feature(j, i);
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				S(k, j) = feature(j * 3 + k + 3, i);
			}
		}

		Eigen::Matrix3d R = Sophus::SO3d::exp(so3).matrix().transpose();


		Eigen::MatrixXd Vt = V.col(i);
		Vt.resize(3, 3);
		Eigen::Matrix3d T = Vt.transpose() * S * R;
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				b(i * 3 + j, k) = T(j, k);
			}
		}
	}

	if (solver.info() == Eigen::Success) {
		Eigen::MatrixXd Vout = solver.solve(A.transpose()*b).transpose();

		Eigen::MatrixXd Vo = Vout.leftCols(VERTS);
		Vo.resize(3 * VERTS, 1);
		return Vo;
	}
	return Eigen::Vector3d::Zero();
}


void recoverFromFaceFeature(const char* outfile, const Eigen::MatrixXd &feature) {

	Eigen::MatrixXd newV(3 * VERTS, feature.cols());
	Eigen::VectorXd Vt(FACES);
	Eigen::MatrixXd V(9, FACES);
	Eigen::SparseMatrix<double> A(3 * FACES, FACES + VERTS);

	Eigen::Matrix3Xd Va;
	Eigen::Matrix3Xi Fa;
	common::read_obj((filepath + "../data/AVE.obj").c_str(), Va, Fa);

	Surface_mesh mesh;
	build_mesh(Va, Fa, mesh);

	double ave = 0;
	for (const auto &e : mesh.edges()) {
		ave += mesh.edge_length(e);
	}
	ave /= mesh.n_edges();

	std::vector<Eigen::Triplet<double>> tri;
	Eigen::Vector3d v[4], norm;
	int id[4];
	for (int i = 0; i < Fa.cols(); ++i) {
		for (int j = 0; j < 3; ++j) {
			id[j] = Fa(j, i);
			v[j] = Va.col(id[j]);
		}
		id[3] = VERTS + i;
		norm = (v[1] - v[0]).cross(v[2] - v[0]);
		v[3] = (v[0] + norm / sqrt(norm.norm()));

		//calc vt
		Vt[i] = 0.5*abs((v[1] - v[0]).dot((v[2] - v[0]).cross(v[3] - v[0])));

		//calc V-1
		Eigen::Matrix3d tmp;
		for (int j = 0; j < 3; ++j) {
			tmp.col(j) = v[j + 1] - v[0];
		}

		double co = sqrt(Vt[i]) / tmp.inverse().norm();
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {
				V(j * 3 + k, i) = tmp(k, j)*co;
			}
		}

		//cout << co << endl;

		//calc A
		for (int j = 0; j < 3; ++j) {
			tri.push_back({ i * 3 + j, id[j + 1], co });
			tri.push_back({ i * 3 + j, id[0], -co });
		}
	}
	A.setFromTriplets(tri.begin(), tri.end());


	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(A.transpose()*A);
	for (int i = 0; i < feature.cols(); ++i) {
		newV.col(i) = calcFaceFeature(A, V, feature.col(i), solver);
	}
	puts("ok");

	Eigen::MatrixXd tmp = newV.col(0);
	tmp.resize(3, VERTS);
	common::save_obj((filepath + "../../bodyViz/Assets/Resources/TEMP.obj").c_str(), tmp, Fa);
	cout << "New obj saved!!" << endl;


	saveRoughExact(outfile, newV, Fa);
}


void getFeatureMeasure(const char* infile, const char* outfile) {
	Eigen::MatrixXd feature;
	readNewFeature(infile, feature, 111);

	for (int i = 0; i < feature.rows(); ++i) {
		if (feature(i, 1) < -200 || feature(i, 1) > 200) cout << "i = " << i / 18 << endl;
	}
	//recoverFromFeature(outfile, feature);

}

void getFaceFeatureMeature(const char* infile, const char* outfile) {
	Eigen::MatrixXd feature;
	readNewFaceFeature(infile, feature, 111);

	recoverFromFaceFeature(outfile, feature);
}

void recoverModelFromRS(const char* infile) {
	Eigen::MatrixXd feature;

	readNewFaceFeature(infile, feature, 1);

	recoverFromFaceFeature((filepath + "../data/recover/newMeasure").c_str(), feature);
}

int main() {
	recoverModelFromRS((filepath + "../data/recover/testRS").c_str());

	//getFaceFeatureMeature("../data/recover/testRS", "../data/recover/logrs_pca_rnn_roughExact");
	//Eigen::MatrixXd feature;
	//common::read_matrix_binary_from_file("../data/test/logRS", feature);
	//recoverFromFaceFeature("../data/recover/testroughExact", feature);

	//getVertexMeature("../data/recover/dv_ae_rnn", "../data/recover/dv_ae_rnn_roughexact");
	getchar();
}
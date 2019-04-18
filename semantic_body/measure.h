#pragma once
#include <algorithm>
#include <geodesic/geodesic_mesh.h>
#include <geodesic/geodesic_algorithm_base.h>
#include <Eigen/Dense>
using namespace std;
class measure
{
private:
	static const int N = 14;
	static const int M = 12;
	int lengthKeyPoint[M][2] = {
{10777, 10722},
{9994, 9931},
{9933, 10044},
{9172, 9285},
{9340, 11075},
{11086, 7588},
{7661, 5882},
{6326, 10722},
{7661, 10735},
{7740, 10972},
{5354, 7740},
{983, 7295}
	};

	int circleKeyPoint[N][4] = {
{9172, 9285, 9268, 9130},
{7661, 7635, 7774, 7740},
{5882, 5671, 5313, 5354},
{7088, 7073, 6999, 6975},
{9797, 10413, 10822, 9752},
{12227, 12219, 12168, 12174},
{10624, 10676, 11087, 11077},
{9090, 9522, 9416, 9129},
{6410, 6397, 6364, 6607},
{5570, 5877, 6035, 5876},
{3862, 3819, 3801, 3851},
{2950, 2964, 2879, 2977},
{1948, 1929, 1916, 1866},
{1232, 1209, 1137, 1200}
	};

	double length[M];
	double circle[N];

public:
	string SemanticLable[N + M] = {
		"bustline", "waistline", "hipline", "midhipline", "armhole",
		"head", "collar", "arm", "wrist", "hand", "thigh", "knee", "calf", "ankle", "shoulderWidth",
		"backWidth", "chestWidth", "breastDist", "chestHeight",
		"backHeight", "waistLength", "sleeveLength", "frontLength",
		"backLength", "crotchLength", "outseam"
	};
	measure();
	~measure();

	void initMesh(geodesic::Mesh &mesh, const Eigen::Matrix3Xd & V, Eigen::Matrix3Xi &F);

	void savePath(std::ofstream &out, const std::vector<geodesic::SurfacePoint> &path);

	void calcExact(const Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F, bool save = true);

	void calcSubdivide(const Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F);
	void calcDijkstra(const Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F);
	void calcLength(geodesic::GeodesicAlgorithmBase *algo, geodesic::Mesh &mesh, bool save = false);
	void calcCircle(geodesic::GeodesicAlgorithmBase *algo, geodesic::Mesh &mesh, bool save = false);

	void writeVTK(const std::string &filename, std::vector<geodesic::SurfacePoint> &path);

	void saveParam(Eigen::MatrixXd &Param, int col) {
		for (int i = 0; i < N; ++i) Param(i, col) = circle[i];
		for (int i = 0; i < M; ++i) Param(i + N, col) = length[i];
	}

	void init() {
		for (int i = 0; i < N; ++i) circle[i] = 0;
		for (int j = 0; j < M; ++j) length[j] = 0;
	}

	void printInfo(int id) {
		if (id < N)
			cout << SemanticLable[id] << " = " << circle[id] << endl;
		else
			cout << SemanticLable[id] << " = " << length[id - N] << endl;
	}
	void printAll() {
		cout << "Print all of the semantic lables.." << endl;
		for (int i = 0; i < N; ++i) {
			cout << SemanticLable[i] << " = " << circle[i] << endl;
		}
		for (int i = 0; i < M; ++i) {
			cout << SemanticLable[i + N] << " = " << length[i] << endl;
		}
	}
	int len() {
		return N + M;
	}

};


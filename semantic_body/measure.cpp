#include "measure.h"
#include "vtk.h"

#include <iostream>
#include <geodesic/geodesic_algorithm_exact.h>
#include <geodesic/geodesic_algorithm_subdivision.h>
#include <geodesic/geodesic_algorithm_dijkstra.h>
using namespace std;


measure::measure() { }

measure::~measure() { }

// line to *.vtk, visualize the geodesics and check.
void measure::writeVTK(const std::string &filename, std::vector<geodesic::SurfacePoint> &path) {
	std::vector<double> nodes;
	std::vector<int> lines;
	std::ofstream os(filename);

	int n = path.size();
	for (int i = 0; i < n; ++i) {
		nodes.push_back(path[i].x());
		nodes.push_back(path[i].y());
		nodes.push_back(path[i].z());

		if (i) {
			lines.push_back(i - 1);
			lines.push_back(i);
		}
	}

	line2vtk(os, nodes.data(), nodes.size() / 3, lines.data(), lines.size() / 2);
	os.close();
	return;
}


void measure::initMesh(geodesic::Mesh &mesh, const Eigen::Matrix3Xd & V, Eigen::Matrix3Xi & F)
{
	init();
	vector<unsigned> ff;

	for (int i = 0; i < F.cols(); ++i) {
		if (i == 5014 || i == 24932) continue;
		for (int j = 0; j < F.rows(); ++j) {
			ff.push_back(F(j, i));
		}
	}
	vector<double> vv(V.data(), V.data() + V.rows()*V.cols());

	mesh.initialize_mesh_data(vv, ff);

	cout << "Mesh initialized ..." << endl;
	return;
}

void measure::calcExact(const Eigen::Matrix3Xd & V, Eigen::Matrix3Xi & F)
{

	geodesic::Mesh mesh;
	initMesh(mesh, V, F);

	geodesic::GeodesicAlgorithmExact algo(&mesh);

	calcLength(&algo, mesh);

	calcCircle(&algo, mesh);

	printAll();
}

void measure::calcSubdivide(const Eigen::Matrix3Xd & V, Eigen::Matrix3Xi & F)
{
	geodesic::Mesh mesh;
	initMesh(mesh, V, F);

	geodesic::GeodesicAlgorithmSubdivision algo(&mesh, 2);

	calcLength(&algo, mesh);

	calcCircle(&algo, mesh);

	printAll();
}


void measure::calcDijkstra(const Eigen::Matrix3Xd & V, Eigen::Matrix3Xi &F)
{
	geodesic::Mesh mesh;
	initMesh(mesh, V, F);

	geodesic::GeodesicAlgorithmDijkstra algo(&mesh);

	calcLength(&algo, mesh);

	calcCircle(&algo, mesh);

	printAll();
}

void measure::calcLength(geodesic::GeodesicAlgorithmBase *algo, geodesic::Mesh &mesh)
{
	/*******************Calculate Geodesic Length*******************/

	for (int i = 0; i < M; ++i) {
		int s = lengthKeyPoint[i][0];
		int t = lengthKeyPoint[i][1];
		std::vector<geodesic::SurfacePoint> source;
		source.push_back(geodesic::SurfacePoint(&mesh.vertices()[s]));

		geodesic::SurfacePoint target(&mesh.vertices()[t]);

		algo->propagate(source);
		//algo->print_statistics();

		std::vector<geodesic::SurfacePoint> path;
		algo->trace_back(target, path);
		geodesic::print_info_about_path(path);
		length[i] = geodesic::length(path);

		//string name = "./checkVTK/" +to_string(i) + "_" +SemanticLable[N+i] + ".vtk";
		//writeVTK(name, path);
		//printInfo(N + i);
	}
}


void measure::calcCircle(geodesic::GeodesicAlgorithmBase *algo, geodesic::Mesh &mesh)
{
	/*********************Calculate Circumstance********************/

	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < 4; ++j) {
			int s = circleKeyPoint[i][j];
			int t = circleKeyPoint[i][(j + 1) % 4];

			std::vector<geodesic::SurfacePoint> source;
			source.push_back(geodesic::SurfacePoint(&mesh.vertices()[s]));

			geodesic::SurfacePoint target(&mesh.vertices()[t]);

			algo->propagate(source);
			//algo->print_statistics();

			std::vector<geodesic::SurfacePoint> path;
			algo->trace_back(target, path);
			geodesic::print_info_about_path(path);
			circle[i] += geodesic::length(path);

			//string name = "./checkVTK/" + to_string(i) +"-" + to_string(j) + "-" + SemanticLable[i]+ ".vtk";
			//writeVTK(name, path);
		}
		//printInfo(i);
	}
}


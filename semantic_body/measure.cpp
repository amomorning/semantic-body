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
void measure::writeVTK(const std::string &filename,
	std::vector<geodesic::SurfacePoint> &path)
{
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


void measure::initMesh(geodesic::Mesh &mesh,
	const Eigen::Matrix3Xd & V, Eigen::Matrix3Xi & F)
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



void measure::calcExact(const Eigen::Matrix3Xd & V,
	Eigen::Matrix3Xi & F, bool save)
{

	geodesic::Mesh mesh;
	initMesh(mesh, V, F);

	geodesic::GeodesicAlgorithmExact algo(&mesh);
	cout << "save ? = " << save << endl;

	calcLength(&algo, mesh, save);

	calcCircle(&algo, mesh, save);

	printAll();
}

void measure::calcSubdivide(const Eigen::Matrix3Xd & V,
	Eigen::Matrix3Xi & F)
{
	geodesic::Mesh mesh;
	initMesh(mesh, V, F);

	geodesic::GeodesicAlgorithmSubdivision algo(&mesh, 2);

	calcLength(&algo, mesh);

	calcCircle(&algo, mesh);

	printAll();
}


void measure::calcDijkstra(const Eigen::Matrix3Xd & V,
	Eigen::Matrix3Xi &F)
{
	geodesic::Mesh mesh;
	initMesh(mesh, V, F);

	geodesic::GeodesicAlgorithmDijkstra algo(&mesh);

	calcLength(&algo, mesh);

	calcCircle(&algo, mesh);

	printAll();
}

void measure::savePath(std::ofstream &out, const std::vector<geodesic::SurfacePoint> &path)
{
	for (auto p : path) {
		cout << p.type() << " - ";
		int x, y;
		double t;
		if (p.type() == geodesic::VERTEX) {
			geodesic::vertex_pointer v = static_cast<geodesic::vertex_pointer>(p.base_element());
			x = v->id(), y = 0;
			t = 0;
			cout << v->id() << " 0" << " 0";
		}
		else if (p.type() == geodesic::EDGE) {
			geodesic::edge_pointer e = static_cast<geodesic::edge_pointer>(p.base_element());
			x = e->v0()->id(), y = e->v1()->id();
			t = p.distance(e->v0()) / e->length();

			cout << e->v0()->id() << " ";
			cout << e->v1()->id() << " ";
			cout << p.distance(e->v0()) / e->length();
		}
		out.write((const char*)(&x), sizeof(int));
		out.write((const char*)(&y), sizeof(int));
		out.write((const char*)(&t), sizeof(double));
		cout << endl;
	}
}

void measure::calcLength(geodesic::GeodesicAlgorithmBase *algo, 
	geodesic::Mesh &mesh, bool save)
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
		
		int len = path.size();
		if (save) {

			ofstream out("./data/path/" + SemanticLable[N + i], ios::binary);
			out.write((const char*)(&len), sizeof(int));

			savePath(out, path);
			//string name = "./checkVTK/" +to_string(i) + "_" +SemanticLable[N+i] + ".vtk";
			//writeVTK(name, path);
			//printInfo(N + i);
			out.close();
		}
	}
}


void measure::calcCircle(geodesic::GeodesicAlgorithmBase *algo, 
	geodesic::Mesh &mesh, bool save)
{
	/*********************Calculate Circumstance********************/

	for (int i = 0; i < N; ++i) {
		
		std::vector<geodesic::SurfacePoint> path[4];
		int len = 0;
		for (int j = 0; j < 4; ++j) {
			int s = circleKeyPoint[i][j];
			int t = circleKeyPoint[i][(j + 1) % 4];

			std::vector<geodesic::SurfacePoint> source;
			source.push_back(geodesic::SurfacePoint(&mesh.vertices()[s]));

			geodesic::SurfacePoint target(&mesh.vertices()[t]);

			algo->propagate(source);
			//algo->print_statistics();

			algo->trace_back(target, path[j]);
			geodesic::print_info_about_path(path[j]);
			circle[i] += geodesic::length(path[j]);

			len += path[j].size();
			//string name = "./checkVTK/" + to_string(i) +"-" + to_string(j) + "-" + SemanticLable[i]+ ".vtk";
			//writeVTK(name, path);
		}
		if (save) {

			ofstream out("./data/path/" + SemanticLable[i], ios::binary);
			out.write((const char*)(&len), sizeof(int));
			for (int j = 0; j < 4; ++j) {
				savePath(out, path[j]);
			}
			out.close();
		}
		//printInfo(i);
	}
}


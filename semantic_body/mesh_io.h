#ifndef MESH_IO_H_H_
#define MESH_IO_H_H_

#include <string>
#include <Eigen/Dense>

namespace  common {
// @brief read and write obj mesh
// @param filename: the file name of the mesh
// @param V: the vertex list of the mesh
// @param F: the face list of the mesh
// @param tV: the texture vertex list of the mesh
// @param tF: the texture face list of the mesh
// @return 0: read or write failed; 1: read or write successufully

int read_obj(const std::string &filename, Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F);
int save_obj(const std::string &filename, const Eigen::Matrix3Xd &nods,const Eigen::Matrix3Xi &tris);

int read_obj(const std::string &in_file,
             Eigen::Matrix3Xd &V, Eigen::Matrix3Xi &F,
             Eigen::Matrix2Xd &tV, Eigen::Matrix3Xi &tF);
int save_obj(const std::string &out_file,
             const Eigen::Matrix3Xd &V,const Eigen::Matrix3Xi &F,
             const Eigen::Matrix2Xd &tv,const Eigen::Matrix3Xi &tf);

}


#endif

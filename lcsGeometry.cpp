/**********************************************
File		:	lcsGeometry.cpp
Author		:	Mingcheng Chen
Last Update	:	January 19th, 2014
***********************************************/

#include "lcsGeometry.h"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <vtkIdList.h>
#include <vtkPointData.h>

namespace lcs {

////////////////////////////////////////////////
int Sign(double a, double epsilon) {
	return a < -epsilon ? -1 : a > epsilon;
}

double Sqr(double a) {
	return a * a;
}

double Cub(double a) {
	return a * a * a;
}

////////////////////////////////////////////////
double Vector::Length() const {
	return sqrt(Sqr(this->x) + Sqr(this->y) + Sqr(this->z));
}

Vector operator + (const Vector &a, const Vector &b) {
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector operator - (const Vector &a, const Vector &b) {
	return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector operator * (const Vector &a, const double &b) {
	return Vector(a.x * b, a.y * b, a.z * b);
}

Vector operator / (const Vector &a, const double &b) {
	return Vector(a.x / b, a.y / b, a.z / b);
}

Vector Cross(const Vector &a, const Vector &b) {
	return Vector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

double Dot(const Vector &a, const Vector &b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

double Mixed(const Vector &a, const Vector &b, const Vector &c) {
	return a.x * (b.y * c.z - b.z * c.y) + a.y * (b.z * c.x - b.x * c.z) + a.z * (b.x * c.y - b.y * c.x);
}

////////////////////////////////////////////////
Tetrahedron::Tetrahedron(const Vector *arr) {
	memcpy(this->vertices, arr, sizeof(Vector) * 4);
}

void Tetrahedron::CalculateNaturalCoordinates(const Vector &point, double *coordinates) const {
	Vector diff(point - this->vertices[0]);

	double V = 1 / Mixed(vertices[1] - vertices[0], vertices[2] - vertices[0], vertices[3] - vertices[0]);

	double z41 = this->vertices[3].GetZ() - this->vertices[0].GetZ();
	double y34 = this->vertices[2].GetY() - this->vertices[3].GetY();
	double z34 = this->vertices[2].GetZ() - this->vertices[3].GetZ();
	double y41 = this->vertices[3].GetY() - this->vertices[0].GetY();
	double a11 = (z41 * y34 - z34 * y41) * V;

	double x41 = this->vertices[3].GetX() - this->vertices[0].GetX();
	double x34 = this->vertices[2].GetX() - this->vertices[3].GetX();
	double a12 = (x41 * z34 - x34 * z41) * V;

	double a13 = (y41 * x34 - y34 * x41) * V;

	coordinates[1] = a11 * diff.GetX() + a12 * diff.GetY() + a13 * diff.GetZ();

	double y12 = this->vertices[0].GetY() - this->vertices[1].GetY();
	double z12 = this->vertices[0].GetZ() - this->vertices[1].GetZ();
	double a21 = (z41 * y12 - z12 * y41) * V;

	double x12 = this->vertices[0].GetX() - this->vertices[1].GetX();
	double a22 = (x41 * z12 - x12 * z41) * V;

	double a23 = (y41 * x12 - y12 * x41) * V;

	coordinates[2] = a21 * diff.GetX() + a22 * diff.GetY() + a23 * diff.GetZ();

	double z23 = this->vertices[1].GetZ() - this->vertices[2].GetZ();
	double y23 = this->vertices[1].GetY() - this->vertices[2].GetY();
	double a31 = (z23 * y12 - z12 * y23) * V;

	double x23 = this->vertices[1].GetX() - this->vertices[2].GetX();
	double a32 = (x23 * z12 - x12 * z23) * V;

	double a33 = (y23 * x12 - y12 * x23) * V;

	coordinates[3] = a31 * diff.GetX() + a32 * diff.GetY() + a33 * diff.GetZ();

	//coordinates[1] *= V;
	//coordinates[2] *= V;
	//coordinates[3] *= V;
	coordinates[0] = 1 - coordinates[1] - coordinates[2] - coordinates[3];
}

////////////////////////////////////////////////
TetrahedralGrid::TetrahedralGrid(vtkUnstructuredGrid *unstructuredGrid) {
	if (!unstructuredGrid) return;
	
	this->numOfVertices = unstructuredGrid->GetNumberOfPoints();
	this->numOfCells = unstructuredGrid->GetNumberOfCells();
	
	this->vertices = new Vector [this->numOfVertices];
	this->velocities = new Vector [this->numOfVertices];
	this->tetrahedralConnectivities = new int [this->numOfCells * 4];
	this->tetrahedralLinks = new int [this->numOfCells * 4];

	double point[3], velocity[3];
	for (int i = 0; i < this->numOfVertices; i++) {
		unstructuredGrid->GetPoint(i, point);
		this->vertices[i] = Vector(point);
		unstructuredGrid->GetPointData()->GetVectors()->GetTuple(i, velocity);
		this->velocities[i] = Vector(velocity);
	}

	int **vertexLinks = new int *[this->numOfVertices];
	int *vertexDegrees = new int [this->numOfVertices];
	memset(vertexDegrees, 0, sizeof(int) * this->numOfVertices);
	
	vtkIdList *idList = vtkIdList::New();
	for (int i = 0; i < this->numOfCells; i++) {
		unstructuredGrid->GetCellPoints(i, idList);
		for (int j = 0; j < 4; j++) {
			this->tetrahedralConnectivities[(i << 2) + j] = idList->GetId(j);
			vertexDegrees[idList->GetId(j)]++;
		}
		Vector a = this->vertices[this->tetrahedralConnectivities[(i << 2) + 0]];
		Vector b = this->vertices[this->tetrahedralConnectivities[(i << 2) + 1]];
		Vector c = this->vertices[this->tetrahedralConnectivities[(i << 2) + 2]];
		Vector d = this->vertices[this->tetrahedralConnectivities[(i << 2) + 3]];
		if (Mixed(b - a, c - a, d - a) < 0) std::swap(this->tetrahedralConnectivities[(i << 2) + 1], this->tetrahedralConnectivities[(i << 2) + 2]);
	}
	idList->Delete();

	for (int i = 0; i < this->numOfVertices; i++) {
		vertexLinks[i] = new int [vertexDegrees[i]];
		vertexDegrees[i] = 0;
	}
	for (int i = 0; i < this->numOfCells; i++)
		for (int j = 0; j < 4; j++) {
			int idx = this->tetrahedralConnectivities[(i << 2) + j];
			vertexLinks[idx][vertexDegrees[idx]++] = i;
		}

	memset(this->tetrahedralLinks, 255, sizeof(int) * this->numOfCells * 4);
	for (int i = 0; i < this->numOfCells; i++)
		for (int j = 0; j < 4; j++)
			if (this->tetrahedralLinks[(i << 2) + j] < 0) {
				int idx = this->tetrahedralConnectivities[(i << 2) + ((j + 1) & 3)];
				for (int k = 2; k <= 3; k++) {
					int tmp = this->tetrahedralConnectivities[(i << 2) + ((j + k) & 3)];
					if (vertexDegrees[tmp] < vertexDegrees[idx]) idx = tmp;
				}
				int v1 = -1, v2 = -1;
				for (int k = 1; k <= 3; k++) {
					int tmp = this->tetrahedralConnectivities[(i << 2) + ((j + k) & 3)];
					if (tmp != idx)
						if (v1 == -1) v1 = tmp;
						else {
							v2 = tmp;
							break;
						}
				}
				for (int k = 0; k < vertexDegrees[idx]; k++) {
					int tet = vertexLinks[idx][k];
					if (tet == i) continue;
					int cnt = 0, v3 = -1;
					for (int l = 0; l < 4; l++) {
						int v = this->tetrahedralConnectivities[(tet << 2) + l];
						if (v == v1) {
							cnt++;
							continue;
						}
						if (v == v2) {
							cnt++;
							continue;
						}
						if (v != idx) v3 = l;
					}
					if (cnt == 2) {
						this->tetrahedralLinks[(i << 2) + j] = tet;
						if (this->tetrahedralLinks[(tet << 2) + v3] != -1) {
							printf("bug!\n");
							printf("i = %d, j = %d\n", i, j);
							printf("v3 = %d\n", v3);
							printf("tet : %d, %d\n", tet, this->tetrahedralLinks[(tet << 2) + v3]);
							printf("idx = %d, v1 = %d, v2 = %d\n", idx, v1, v2);
							while (1);
						}
						this->tetrahedralLinks[(tet << 2) + v3] = i;
						break;
					}
				}
			}

	delete [] vertexDegrees;
	for (int i = 0; i < this->numOfVertices; i++)
		delete [] vertexLinks[i];
	delete [] vertexLinks;
}

TetrahedralGrid::TetrahedralGrid(int numOfCells, int numOfPoints, int *conn, int *link, double *posi, double *velo) {
	this->numOfCells = numOfCells;
	this->numOfVertices = numOfPoints;
	vertices = new lcs::Vector [numOfPoints];
	velocities = new lcs::Vector [numOfPoints];
	for (int i = 0; i < numOfPoints; i++) {
		vertices[i] = lcs::Vector(posi[i * 3], posi[i * 3 + 1], posi[i * 3 + 2]);
		velocities[i] = lcs::Vector(velo[i * 3], velo[i * 3 + 1], velo[i * 3 + 2]);
	}
	tetrahedralConnectivities = new int [numOfCells << 2];
	tetrahedralLinks = new int [numOfCells << 2];
	memcpy(tetrahedralConnectivities, conn, sizeof(int) * numOfCells * 4);
	memcpy(tetrahedralLinks, link, sizeof(int) * numOfCells * 4);
}

Tetrahedron TetrahedralGrid::GetTetrahedron(int index) const {
	int a = this->tetrahedralConnectivities[index << 2];
	int b = this->tetrahedralConnectivities[(index << 2) + 1];
	int c = this->tetrahedralConnectivities[(index << 2) + 2];
	int d = this->tetrahedralConnectivities[(index << 2) + 3];
	return Tetrahedron(this->vertices[a], this->vertices[b], this->vertices[c], this->vertices[d]);
}

int TetrahedralGrid::FindCell(const Vector &point, const double &epsilon, int guess) const {
	Tetrahedron tet;
	double coordinates[4];
	int curr = guess, violator;
	while (1) {
		tet = this->GetTetrahedron(curr);
		tet.CalculateNaturalCoordinates(point, coordinates);
		violator = 0;
		for (int i = 1; i < 4; i++)
			if (coordinates[i] < coordinates[violator]) violator = i;
		if (coordinates[violator] >= -epsilon) return curr;
		curr = this->tetrahedralLinks[(curr << 2) + violator];
		if (curr < 0) return -1;
	}
}

int TetrahedralGrid::ProfiledFindCell(const Vector &point, const double &epsilon, int guess) {
	Tetrahedron tet;
	double coordinates[4];
	int curr = guess, violator;
	this->lastFindCellCost = 0;
	while (1) {
		tet = this->GetTetrahedron(curr);
		tet.CalculateNaturalCoordinates(point, coordinates);
		violator = 0;
		for (int i = 1; i < 4; i++)
			if (coordinates[i] < coordinates[violator]) violator = i;
		if (coordinates[violator] >= -epsilon) return curr;
		this->lastFindCellCost++;
		curr = this->tetrahedralLinks[(curr << 2) + violator];
		if (curr < 0) return -1;
	}
}

int TetrahedralGrid::FindCell(const Vector &point, const double &epsilon) const {
	Tetrahedron tet;
	double coordinates[4];
	for (int i = 0; i < this->numOfCells; i++) {
		tet = this->GetTetrahedron(i);
		tet.CalculateNaturalCoordinates(point, coordinates);
		bool flag = true;
		for (int j = 0; j < 4; j++)
			if (coordinates[j] < -epsilon) {
				flag = false;
				break;
			}
		if (flag) return i;
	}
	return -1;
}

void TetrahedralGrid::GetInterpolatedVelocity(const Vector &point, int cellId, double *velocity) const {
	double coordinates[4];
	double tempV[3];
	this->GetTetrahedron(cellId).CalculateNaturalCoordinates(point, coordinates);
	memset(velocity, 0, sizeof(double) * 3);
	for (int i = 0; i < 4; i++) {
		int vtx = this->tetrahedralConnectivities[(cellId << 2) + i];
		velocity[0] += this->velocities[vtx].GetX() * coordinates[i];
		velocity[1] += this->velocities[vtx].GetY() * coordinates[i];
		velocity[2] += this->velocities[vtx].GetZ() * coordinates[i];
	}
}

}

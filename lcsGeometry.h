/**********************************************
File		:	lcsGeometry.h
Author		:	Mingcheng Chen
Last Update	:	January 19th, 2013
***********************************************/

#ifndef __LCS_Geometry_H
#define __LCS_Geometry_H

#include <vtkUnstructuredGrid.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdio>

namespace lcs {

////////////////////////////////////////////////
int Sign(double a, double epsilon);

double Sqr(double a);

double Cub(double a);

class Vector;

lcs::Vector operator + (const lcs::Vector &, const lcs::Vector &);
lcs::Vector operator - (const lcs::Vector &, const lcs::Vector &);
lcs::Vector operator * (const lcs::Vector &, const double &);
lcs::Vector operator / (const lcs::Vector &, const double &);
lcs::Vector Cross(const lcs::Vector &, const lcs::Vector &);
double Dot(const lcs::Vector &, const lcs::Vector &);
double Mixed(const lcs::Vector &, const lcs::Vector &, const lcs::Vector &);

////////////////////////////////////////////////
class Vector {
public:
	Vector() {
		this->x = this->y = this->z = 0;
	}

	Vector(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	Vector(double *arr) {
		this->x = arr[0];
		this->y = arr[1];
		this->z = arr[2];
	}

	Vector(const Vector &vec) {
		this->x = vec.x;
		this->y = vec.y;
		this->z = vec.z;
	}

	void SetX(double x) {
		this->x = x;
	}

	void SetY(double y) {
		this->y = y;
	}

	void SetZ(double z) {
		this->z = z;
	}

	double GetX() const {
		return this->x;
	}

	double GetY() const {
		return this->y;
	}

	double GetZ() const {
		return this->z;
	}

	double Length() const;

	void Output() const {
		printf("(%lf, %lf, %lf)\n", this->x, this->y, this->z);
	}

	friend lcs::Vector operator + (const lcs::Vector &, const lcs::Vector &);
	friend lcs::Vector operator - (const lcs::Vector &, const lcs::Vector &);
	friend lcs::Vector operator * (const lcs::Vector &, const double &);
	friend lcs::Vector operator / (const lcs::Vector &, const double &);
	friend lcs::Vector lcs::Cross(const lcs::Vector &, const lcs::Vector &);
	friend double lcs::Dot(const lcs::Vector &, const lcs::Vector &);
	friend double lcs::Mixed(const lcs::Vector &, const lcs::Vector &, const lcs::Vector &);

private:
	double x, y, z;
};

////////////////////////////////////////////////
class Tetrahedron {
public:
	Tetrahedron() {
	}

	Tetrahedron(const Vector &a, const Vector &b, const Vector &c, const Vector &d) {
		this->vertices[0] = a;
		this->vertices[1] = b;
		this->vertices[2] = c;
		this->vertices[3] = d;
	}

	Tetrahedron(const Vector *);

	Vector GetVertex(int index) const {
		return this->vertices[index];
	}

	void SetVertex(int index, const Vector &a) {
		this->vertices[index] = a;
	}

	void CalculateNaturalCoordinates(const Vector &, double *) const;

	double Size() const {
		double result = 0;
		for (int i = 0; i < 3; i++)
			for (int j = i + 1; j < 4; j++) {
				double dist = (this->GetVertex(i) - this->GetVertex(j)).Length();
				if (dist > result) result = dist;
			}
		return result;
	}

	double Volume() const {
		Vector A = this->vertices[1] - this->vertices[0];
		Vector B = this->vertices[2] - this->vertices[0];
		Vector C = this->vertices[3] - this->vertices[0];
		return fabs(Mixed(A, B, C)) / 6.0;
	}

	void Output() const {
		for (int i = 0; i < 4; i++) {
			if (i) printf(" ");
			this->vertices[i].Output();
		}
	}

private:
	Vector vertices[4];
};

////////////////////////////////////////////////
class TetrahedralGrid {
public:
	TetrahedralGrid() {
		numOfVertices = 0;
		numOfCells = 0;
		vertices = NULL;
		tetrahedralConnectivities = NULL;
		tetrahedralLinks = NULL;
	}

	TetrahedralGrid(vtkUnstructuredGrid *);
	TetrahedralGrid(int numOfCells, int numOfPoints, int *conn, int *link, double *posi, double *velo);

	~TetrahedralGrid() {
		if (vertices) delete [] vertices;
		if (velocities) delete [] velocities;
		if (tetrahedralConnectivities) delete [] tetrahedralConnectivities;
		if (tetrahedralLinks) delete [] tetrahedralLinks;
	}

	void CleanAllButVelocities() {
		if (vertices) delete [] vertices;
		if (tetrahedralConnectivities) delete [] tetrahedralConnectivities;
		if (tetrahedralLinks) delete [] tetrahedralLinks;
	}

	Tetrahedron GetTetrahedron(int) const;

	void TetrahedronSize() const {
		double minTet = 1e100, maxTet = 0, sumTet = 0;
		for (int i = 0; i < numOfCells; i++) {
			double vol = this->GetTetrahedron(i).Volume();
			if (vol < minTet) minTet = vol;
			if (vol > maxTet) maxTet = vol;
			sumTet += vol;
		}
		printf("In volume, minTet: %0.20lf, maxTet: %.20lf, aveTet: %.20lf\n", minTet, maxTet, sumTet / numOfCells);

		minTet = 1e100;
		maxTet = 0;
		sumTet = 0;
		for (int i = 0; i < numOfCells; i++) {
			double size = this->GetTetrahedron(i).Size();
			if (size < minTet) minTet = size;
			if (size > maxTet) maxTet = size;
			sumTet += size;
		}
		printf("In edge length, minTet: %0.20lf, maxTet: %.20lf, aveTet: %.20lf\n", minTet, maxTet, sumTet / numOfCells);
	}

	Vector GetVertex(int index) const {
		return vertices[index];
	}

	Vector GetVelocity(int index) const {
		return velocities[index];
	}

	int FindCell(const Vector &, const double &) const;

	int FindCell(const Vector &, const double &, int) const;

	int ProfiledFindCell(const Vector &, const double &, int);

	void GetInterpolatedVelocity(const Vector &, int, double *) const;

	int GetNumOfVertices() const {
		return this->numOfVertices;
	}

	int GetNumOfCells() const {
		return this->numOfCells;
	}

	void GetCellLink(int index, int *arr) const {
		memcpy(arr, this->tetrahedralLinks + (index << 2), sizeof(int) * 4);
	}

	void GetCellConnectivity(int index, int *arr) const {
		memcpy(arr, this->tetrahedralConnectivities + (index << 2), sizeof(int) * 4);
	}

	int GetLastFindCellCost() const {
		return this->lastFindCellCost;
	}

	void ReadConnectivities(int *destination) const {
		memcpy(destination, this->tetrahedralConnectivities, sizeof(int) * 4 * this->numOfCells);
	}

	void ReadLinks(int *destination) const {
		memcpy(destination, this->tetrahedralLinks, sizeof(int) * 4 * this->numOfCells);
	}

	void ReadPositions(double *destination) const {
		for (int i = 0; i < this->numOfVertices; i++) {
			destination[i * 3] = this->vertices[i].GetX();
			destination[i * 3 + 1] = this->vertices[i].GetY();
			destination[i * 3 + 2] = this->vertices[i].GetZ();
		}
	}

	void ReadPositions(float *destination) const {
		for (int i = 0; i < this->numOfVertices; i++) {
			destination[i * 3] = this->vertices[i].GetX();
			destination[i * 3 + 1] = this->vertices[i].GetY();
			destination[i * 3 + 2] = this->vertices[i].GetZ();
		}
	}

	void ReadVelocities(double *destination) const {
		for (int i = 0; i < this->numOfVertices; i++) {
			destination[i * 3] = this->velocities[i].GetX();
			destination[i * 3 + 1] = this->velocities[i].GetY();
			destination[i * 3 + 2] = this->velocities[i].GetZ();
		}
	}

	void ReadVelocities(float *destination) const {
		for (int i = 0; i < this->numOfVertices; i++) {
			destination[i * 3] = this->velocities[i].GetX();
			destination[i * 3 + 1] = this->velocities[i].GetY();
			destination[i * 3 + 2] = this->velocities[i].GetZ();
		}
	}

private:
	int numOfVertices;
	int numOfCells;
	Vector *vertices, *velocities;
	int *tetrahedralConnectivities;
	int *tetrahedralLinks;

// Profiling Result
	int lastFindCellCost;
};

}

#endif

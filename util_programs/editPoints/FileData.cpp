#include "FileData.h"

#include<iostream>
#include<fstream>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

FileData::FileData(
	const std::string editablePointsFile,
	const std::string editablePointsResultFile) {
	FileData::FileData(editablePointsFile, "", editablePointsResultFile);
}
FileData::FileData(
	const std::string editablePointsFile,
	const std::string dispersionCurvesPointsFile,
	const std::string editablePointsResultFile):editablePointsResultFile(editablePointsResultFile) {

	try {
		std::ifstream file(editablePointsFile);
		if (!file.is_open()) throw std::exception("can not open a file with editable points");
		while (!file.eof()) {
			double f, real, imag;
			file >> f >> real >> imag;
			points.push_back(Point(f, real, imag, false));
		}
		file.close();
	}
	catch (std::exception e) {
		MessageBoxA(NULL, e.what(), "cannot load editable points", MB_OKCANCEL | MB_ICONEXCLAMATION);
		throw std::exception();
	}

	if (dispersionCurvesPointsFile.empty()) return;

	try {
		std::ifstream file(dispersionCurvesPointsFile);
		if (!file.is_open()) throw std::exception("can not open a file with dispersion curves points");
		while (!file.eof()) {
			double f, value;
			file >> f >> value;
			dispersionCurves.push_back(std::pair(f,value));
		}
		file.close();
	}
	catch (std::exception e) {
		MessageBoxA(NULL, e.what(), "cannot load dispersion curves points", MB_OKCANCEL | MB_ICONEXCLAMATION);
		throw std::exception();
	}
}

void FileData::saveData() const {
	std::ofstream file;
	try {
		file.open(editablePointsResultFile);
	}
	catch (...) {
		MessageBoxA(NULL, "can not open a file for writeting result", "cannot load dispersion curves points", MB_OKCANCEL | MB_ICONEXCLAMATION);
		throw std::exception();
	}

	for (const Point& p : points) {
		if (p.deleted) continue;

		file << p.f << " " << p.real << " " << p.imag << std::endl;
	}

	file.close();
}

const std::list<FileData::Point>& FileData::getPoints() const {
	return points;
}
const std::list<std::pair<double, double>>& FileData::getDispersionCurvesPoints() const {
	return dispersionCurves;
}

void FileData::deletePoints(double fromF, double fromReal, double toF, double toReal) {
	for (Point& p:points) {
		if (!(fromF < p.f && p.f < toF && fromReal < p.real && p.real < toReal)) continue;

		p.deleted = true;
	}
}
void FileData::restorePoints(double fromF, double fromReal, double toF, double toReal) {
	for (Point& p : points) {
		if (!(fromF < p.f && p.f < toF && fromReal < p.real && p.real < toReal)) continue;

		p.deleted = false;
	}
}

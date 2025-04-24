#pragma once

#include<string>
#include<list>
#include<utility>

class FileData final {
private:
	struct Point {
		double f = 0.0;
		double real = 0.0;
		double imag = 0.0;
		bool deleted = false;

		Point(double f, double real, double imag, bool deleted) :
			f(f), real(real), imag(imag), deleted(deleted) {}
	};
	const std::string editablePointsResultFile;
	std::list<Point> points;
	std::list<std::pair<double,double>> dispersionCurves;
public:
	FileData() = delete;
	FileData(
		const std::string editablePointsFile,
		const std::string editablePointsResultFile
	);
	FileData(
		const std::string editablePointsFile,
		const std::string dispersionCurvesPointsFile,
		const std::string editablePointsResultFile
	);
	FileData(const FileData&) = delete;
	FileData(FileData&&) = delete;

	FileData operator=(const FileData) = delete;
	FileData& operator=(const FileData&) = delete;
	FileData&& operator=(FileData&&) = delete;

	bool operator==(const FileData) const = delete;
	bool operator==(const FileData&) const = delete;
	bool operator==(const FileData&&) const = delete;

	void saveData() const;
	const std::list<Point>& getPoints() const;
	const std::list<std::pair<double, double>>& getDispersionCurvesPoints() const;
	void deletePoints(double fromF, double fromReal, double toF, double toReal);
	void restorePoints(double fromF, double fromReal, double toF, double toReal);

	~FileData() = default;
};

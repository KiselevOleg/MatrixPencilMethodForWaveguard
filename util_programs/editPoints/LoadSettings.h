#pragma once

#include<string>

class LoadSettings final  {
private:
	static unsigned int windowSizeX, windowSizeY;
	static std::string editablePointsFileName;
	static std::string dispersionCurvesPointsFileName;
	static std::string editablePointsResultFileName;

	static bool dataNotLoaded;
	static void loadData();
	static void createDafaultDataFile();
public:
	LoadSettings() = delete;
	LoadSettings(const LoadSettings&) = delete;
	LoadSettings(LoadSettings&&) = delete;

	LoadSettings operator=(const LoadSettings) = delete;
	LoadSettings& operator=(const LoadSettings&) = delete;
	LoadSettings&& operator=(LoadSettings&&) = delete;

	bool operator==(const LoadSettings) const = delete;
	bool operator==(const LoadSettings&) const = delete;
	bool operator==(const LoadSettings&&) const = delete;

	static unsigned int getWindowSizeX();
	static unsigned int getWindowSizeY();
	static std::string getEditablePointsFileName();
	static std::string getDispersionCurvesPointsFileName();
	static std::string getEditablePointsResultFileName();

	~LoadSettings() = default;
};

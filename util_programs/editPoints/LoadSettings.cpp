#include "LoadSettings.h"

#include<fstream>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

const std::string settingsFileName = "settings.txt";

const std::string editablePointsFileNameParameterString = "editablePointsFileName";
const std::string dispersionCurvesPointsFileNameParameterString = "dispersionCurvesPointsFileName";
const std::string editablePointsResultFileNameParameterString = "editablePointsResultFileName";

const std::string editablePointsFileNameDefauldValue = "dispersion_curves.data";
const std::string dispersionCurvesPointsFileNameDefauldValue = "dispersion_curves_K.data";
const std::string editablePointsResultFileNameDefauldValue = "dispersion_curves_result.data";

const std::string windowSizeXParameterString = "windowSizeX";
const std::string windowSizeYParameterString = "windowSizeY";

const unsigned int windowSizeXDefauldValue = 1920;
const unsigned int windowSizeYDefauldValue = 1080;

bool LoadSettings::dataNotLoaded = true;

std::string LoadSettings::editablePointsFileName = editablePointsFileNameDefauldValue;
std::string LoadSettings::dispersionCurvesPointsFileName = dispersionCurvesPointsFileNameDefauldValue;
std::string LoadSettings::editablePointsResultFileName = editablePointsResultFileNameDefauldValue;

unsigned int LoadSettings::windowSizeX = windowSizeXDefauldValue;
unsigned int LoadSettings::windowSizeY = windowSizeYDefauldValue;

void LoadSettings::loadData() {
	dataNotLoaded = false;

	std::ifstream file;
	try {
		file.open(settingsFileName);
		if (!file.is_open()) throw std::exception();
	}
	catch (...) {
		createDafaultDataFile();

		return;
	}

	while (!file.eof()) {
		std::string parameterName;
		std::getline(file, parameterName);

		if (parameterName.empty()) continue;

		if (parameterName == editablePointsFileNameParameterString) {
			file >> editablePointsFileName;

			continue;
		}
		if (parameterName == dispersionCurvesPointsFileNameParameterString) {
			file >> dispersionCurvesPointsFileName;

			if (
				dispersionCurvesPointsFileName == "none" ||
				dispersionCurvesPointsFileName == "-" ||
				dispersionCurvesPointsFileName == "null" ||
				dispersionCurvesPointsFileName == "NULL"
				) dispersionCurvesPointsFileName = "";

			continue;
		}
		if (parameterName == editablePointsResultFileNameParameterString) {
			file >> editablePointsResultFileName;

			continue;
		}

		if (parameterName == windowSizeXParameterString) {
			file >> windowSizeX;

			continue;
		}
		if (parameterName == windowSizeYParameterString) {
			file >> windowSizeY;

			continue;
		}

		MessageBoxA(NULL, ("unknown parameter " + parameterName).c_str(), "incorrect a settings file", MB_OKCANCEL | MB_ICONEXCLAMATION);
		throw std::exception();
	}

	file.close();

	if (!(100 < windowSizeX && windowSizeX < 3000)) {
		MessageBoxA(NULL, 
			("incorrect windowSizeX value " + std::to_string(windowSizeX) + " (required 100 < windowSizeX && windowSizeX < 3000)").c_str(),
			"incorrect a settings file", MB_OKCANCEL | MB_ICONEXCLAMATION);
		throw std::exception();
	}
	if (!(100 < windowSizeY && windowSizeY < 3000)) {
		MessageBoxA(NULL,
			("incorrect windowSizeY value " + std::to_string(windowSizeY) + " (required 100 < windowSizeYX && windowSizeY < 3000)").c_str(),
			"incorrect a settings file", MB_OKCANCEL | MB_ICONEXCLAMATION);
		throw std::exception();
	}
}
void LoadSettings::createDafaultDataFile() {
	std::ofstream file;

	try {
		file.open(settingsFileName);
	}
	catch (...) {
		MessageBoxA(NULL, "cannot create a settings file 'settings.txt'", "cannot create a settings file", MB_OKCANCEL | MB_ICONEXCLAMATION);
		throw std::exception();
	}

	file << editablePointsFileNameParameterString << std::endl << editablePointsFileNameDefauldValue << std::endl;
	file << dispersionCurvesPointsFileNameParameterString << std::endl << dispersionCurvesPointsFileNameDefauldValue << std::endl;
	file << editablePointsResultFileNameParameterString << std::endl << editablePointsResultFileNameDefauldValue << std::endl;

	file << std::endl;

	file << windowSizeXParameterString << std::endl << windowSizeXDefauldValue << std::endl;
	file << windowSizeYParameterString << std::endl << windowSizeYDefauldValue << std::endl;

	file.close();
}

unsigned int LoadSettings::getWindowSizeX() {
	if (dataNotLoaded) {
		loadData();
		dataNotLoaded = true;
	}

	return windowSizeX;
}
unsigned int LoadSettings::getWindowSizeY() {
	if (dataNotLoaded) {
		loadData();
		dataNotLoaded = true;
	}

	return windowSizeY;
}
std::string LoadSettings::getEditablePointsFileName() {
	if (dataNotLoaded) {
		loadData();
		dataNotLoaded = true;
	}

	return editablePointsFileName;
}
std::string LoadSettings::getDispersionCurvesPointsFileName() {
	if (dataNotLoaded) {
		loadData();
		dataNotLoaded = true;
	}

	return dispersionCurvesPointsFileName;
}
std::string LoadSettings::getEditablePointsResultFileName() {
	if (dataNotLoaded) {
		loadData();
		dataNotLoaded = true;
	}

	return editablePointsResultFileName;
}

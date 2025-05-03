#pragma once

#include<list>

#include <SFML/Graphics.hpp>

#include "FileData.h"

class ShowPoints final {
private:
	void centerCircle(sf::CircleShape& circle) const;

	double minF, maxF;
	double minReal, maxReal;
	sf::Vector2f getWindowCoordinatesFromPoint(double f, double value, const sf::RenderWindow& window) const;
	sf::Vector2f getPointFromWindowCoordinates(unsigned int x, unsigned int y, const sf::RenderWindow& window) const;

	FileData& fd;

	void resetPosition();

	bool clicked = false;
	unsigned int clickStartPositionX, clickStartPositionY;
	enum OperationType {deletePoints, restorePoints, changeShowableSpace};
	OperationType operationType = OperationType::deletePoints;

	void changeOperationTypeForward();
	void changeOperationTypeBack();

	bool markerNotDeletedPoints = false;
	mutable float timeFormArkeringNotDeletedPoints = 0.0f;
public:
	ShowPoints() = delete;
	ShowPoints(FileData& fd);
	ShowPoints(const ShowPoints&) = delete;
	ShowPoints(ShowPoints&&) = delete;

	ShowPoints operator=(const ShowPoints) = delete;
	ShowPoints& operator=(const ShowPoints&) = delete;
	ShowPoints&& operator=(ShowPoints&&) = delete;

	bool operator==(const ShowPoints) const = delete;
	bool operator==(const ShowPoints&) const = delete;
	bool operator==(const ShowPoints&&) const = delete;

	void show(sf::RenderWindow& window, float time) const;

	void handle(const std::optional<sf::Event>& event, float time, const sf::RenderWindow& window);

	~ShowPoints() = default;
};

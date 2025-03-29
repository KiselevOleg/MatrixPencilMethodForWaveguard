#include"ShowPoints.h"

ShowPoints::ShowPoints(FileData& fd) :fd(fd) {
	resetPosition();
}

void ShowPoints::resetPosition() {
	minF = DBL_MAX; maxF = DBL_MIN;
	minReal = DBL_MAX; maxReal = DBL_MIN;

	for (const auto& p : fd.getPoints()) {
		minF = std::min(minF, p.f);
		maxF = std::max(minF, p.f);
		minReal = std::min(minReal, p.real);
		maxReal = std::max(maxReal, p.real);
	}

	//if (!fd.getPoints().empty()) return;

	for (const auto& p : fd.getDispersionCurvesPoints()) {
		minF = std::min(minF, p.first);
		maxF = std::max(minF, p.first);
		minReal = std::min(minReal, p.second);
		maxReal = std::max(maxReal, p.second);
	}

	if (!fd.getDispersionCurvesPoints().empty()) return;

	minF = -1.0; maxF = 1.0;
	minReal = -1.0; maxReal = 1.0;
}


void ShowPoints::centerCircle(sf::CircleShape& circle) const {
	circle.setPosition({circle.getPosition().x-circle.getRadius(), circle.getPosition().y-circle.getRadius()});
}

void ShowPoints::show(sf::RenderWindow& window, float time) const {
	static sf::CircleShape circle;

	circle.setRadius(5);
	for (const auto& p : fd.getPoints()) {
		if (!p.deleted) circle.setFillColor(sf::Color::Blue);
		else circle.setFillColor(sf::Color(100, 100, 100));
		circle.setPosition(getWindowCoordinatesFromPoint(p.f, p.real, window));

		if (markerNotDeletedPoints && !p.deleted) circle.setRadius(5+15*abs(sin(timeFormArkeringNotDeletedPoints*5.0f)));
		centerCircle(circle);
		window.draw(circle);
		if (markerNotDeletedPoints && !p.deleted) circle.setRadius(5);



		if (!p.deleted) circle.setFillColor(sf::Color::Red);
		else circle.setFillColor(sf::Color(100, 100, 100));
		circle.setPosition(getWindowCoordinatesFromPoint(p.f, p.imag, window));

		centerCircle(circle);
		window.draw(circle);
	}
	timeFormArkeringNotDeletedPoints += time;

	circle.setRadius(3);
	for (const auto& p : fd.getDispersionCurvesPoints()) {
		circle.setFillColor(sf::Color::Black);
		circle.setPosition(getWindowCoordinatesFromPoint(p.first, p.second, window));

		centerCircle(circle);
		window.draw(circle);
	}

	if (clicked) {
		static sf::RectangleShape rectangle;
		rectangle.setFillColor(sf::Color::Magenta);
		int x1 = clickStartPositionX, x2 = sf::Mouse::getPosition(window).x;
		int y1 = clickStartPositionY, y2 = sf::Mouse::getPosition(window).y;

		rectangle.setPosition(sf::Vector2f(std::min(x1, x2), std::min(y1, y2)));
		rectangle.setSize(sf::Vector2f(abs(x1 - x2), 3));
		window.draw(rectangle);
		rectangle.setSize(sf::Vector2f(3, abs(y1 - y2)));
		window.draw(rectangle);

		rectangle.setPosition(sf::Vector2f(std::max(x1, x2), std::max(y1, y2)));
		rectangle.setSize(sf::Vector2f(-abs(x1 - x2), -3));
		window.draw(rectangle);
		rectangle.setSize(sf::Vector2f(-3, -abs(y1 - y2)));
		window.draw(rectangle);
	}

	sf::Font font("fonts\\manrope\\fonts\\ttf\\manrope-regular.ttf");
	static sf::Text text(font);
	text.setFillColor(sf::Color::Black);
	text.setCharacterSize(48);
	text.setPosition({ 0,0 });
	switch (operationType) {
	case OperationType::deletePoints: text.setString("delete"); break;
	case OperationType::restorePoints: text.setString("restore"); break;
	case OperationType::changeShowableSpace: text.setString("change space"); break;
	}
	window.draw(text);
}

void ShowPoints::handle(const std::optional<sf::Event>& event, float time, const sf::RenderWindow& window) {
	float moveSpeed = 1.0;

	if (event->is<sf::Event::KeyPressed>() && (
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::Up ||
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::W
		)) {
		float speed = (maxReal - minReal) * moveSpeed;
		minReal += speed * time;
		maxReal += speed * time;
	}
	if (event->is<sf::Event::KeyPressed>() && (
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::Down ||
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::S
		)) {
		float speed = (maxReal - minReal) * moveSpeed;
		minReal -= speed * time;
		maxReal -= speed * time;
	}
	if (event->is<sf::Event::KeyPressed>() && (
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::Left ||
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::A
		)) {
		float speed = (maxF - minF) * moveSpeed;
		minF -= speed * time;
		maxF -= speed * time;
	}
	if (event->is<sf::Event::KeyPressed>() && (
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::Right ||
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::D
		)) {
		float speed = (maxF - minF) * moveSpeed;
		minF += speed * time;
		maxF += speed * time;
	}
	if (event->is<sf::Event::KeyPressed>() && sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Space)) resetPosition();



	if (event->is<sf::Event::KeyPressed>() && event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::V &&
		!markerNotDeletedPoints) {
		markerNotDeletedPoints = true;
		timeFormArkeringNotDeletedPoints = 0.0f;
	}
	if (event->is<sf::Event::KeyReleased>() && event->getIf<sf::Event::KeyReleased>()->code == sf::Keyboard::Key::V)
		markerNotDeletedPoints = false;

	

	if (event->is<sf::Event::MouseButtonPressed>() &&
		event->getIf<sf::Event::MouseButtonPressed>()->button == sf::Mouse::Button::Left) {
		clicked = true;
		clickStartPositionX = sf::Mouse::getPosition(window).x;
		clickStartPositionY = sf::Mouse::getPosition(window).y;
	}
	if (event->is<sf::Event::MouseButtonPressed>() &&
		event->getIf<sf::Event::MouseButtonPressed>()->button == sf::Mouse::Button::Right) {
		clicked = false;
	}
	if (event->is<sf::Event::MouseButtonReleased>() &&
		event->getIf<sf::Event::MouseButtonReleased>()->button == sf::Mouse::Button::Left &&
		clicked) {
		clicked = false;
		float x1 = clickStartPositionX;
		float y1 = clickStartPositionY;
		float x2 = sf::Mouse::getPosition(window).x;
		float y2 = sf::Mouse::getPosition(window).y;

		sf::Vector2f p1 = getPointFromWindowCoordinates(x1, y1, window);
		sf::Vector2f p2 = getPointFromWindowCoordinates(x2, y2, window);

		if (p1.x > p2.x) std::swap(p1.x, p2.x);
		if (p1.y > p2.y) std::swap(p1.y, p2.y);

		switch (operationType) {
		case OperationType::deletePoints: fd.deletePoints(p1.x, p1.y, p2.x, p2.y); break;
		case OperationType::restorePoints: fd.restorePoints(p1.x, p1.y, p2.x, p2.y); break;
		case OperationType::changeShowableSpace: minF = p1.x; maxF = p2.x; minReal = p1.y; maxReal = p2.y; break;
		}
	}



	if (event->is<sf::Event::KeyPressed>() && (
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::Q
		)) {
		changeOperationTypeBack();
	}
	if (event->is<sf::Event::KeyPressed>() && (
		event->getIf<sf::Event::KeyPressed>()->code == sf::Keyboard::Key::E
		)) {
		changeOperationTypeForward();
	}
}

sf::Vector2f ShowPoints::getWindowCoordinatesFromPoint(double f, double value, const sf::RenderWindow& window) const {
	float k = window.getSize().x / (maxF - minF);
	float b = -k * minF;
	float vx = k * f + b;
	k = window.getSize().y / (maxReal - minReal);
	b = -k * minReal;
	float vy = k * value + b;

	vy = window.getSize().y - vy;

	return { vx,vy };
}
sf::Vector2f ShowPoints::getPointFromWindowCoordinates(unsigned int x, unsigned int y, const sf::RenderWindow& window) const {
	y = window.getSize().y - y;

	float b = minF;
	float k = (maxF - minF) / window.getSize().x;
	float f = k * x + b;

	b = minReal;
	k = (maxReal - minReal) / window.getSize().y;
	float value = k * y + b;

	return { f,value };
}

void ShowPoints::changeOperationTypeForward() {
	switch (operationType) {
	case OperationType::deletePoints: operationType = OperationType::changeShowableSpace; break;
	case OperationType::restorePoints: operationType = OperationType::deletePoints; break;
	case OperationType::changeShowableSpace: operationType = OperationType::restorePoints; break;
	}
}
void ShowPoints::changeOperationTypeBack() {
	switch (operationType) {
	case OperationType::deletePoints: operationType = OperationType::restorePoints; break;
	case OperationType::restorePoints: operationType = OperationType::changeShowableSpace; break;
	case OperationType::changeShowableSpace: operationType = OperationType::deletePoints; break;
	}
}

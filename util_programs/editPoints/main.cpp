#include <SFML/Graphics.hpp>

#include "FileData.h"
#include "LoadSettings.h"
#include "ShowPoints.h"

#include<iostream>

int WinMain() {
//int main() {
    FileData fd(
        LoadSettings::getEditablePointsFileName(), 
        LoadSettings::getDispersionCurvesPointsFileName(), 
        LoadSettings::getEditablePointsResultFileName()
    );
    ShowPoints sp(fd);
    
    sf::RenderWindow window(sf::VideoMode({ LoadSettings::getWindowSizeX(), LoadSettings::getWindowSizeY() }), "editPoints");
    window.setFramerateLimit(60);

    sf::Clock clock;
    while (window.isOpen()) {
        float time = clock.restart().asSeconds();
        static float Time = 0.0; Time += time;
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) window.close();
            if (event->is<sf::Event::KeyPressed>() && sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Escape)) window.close();
            if (event->is<sf::Event::Resized>()) {
                const auto position = window.getPosition();
                window.close();
                window = sf::RenderWindow(sf::VideoMode(window.getSize()), "editPoints");
                window.setPosition(position);
            }

            if (!window.hasFocus()) continue;

            sp.handle(event, time, window);
        }
        
        window.clear(sf::Color(200,200,200));
        sp.show(window, time);
        window.display();
    }

    fd.saveData();
}

#include "SDL2/SDL.h"

#ifndef __RENDERER_HPP
#define __RENDERER_HPP

class Renderer {
public:
    /// Constructor.
    Renderer();

    /// Destructor.
    ~Renderer();

    /// Initialises the renderer with the given dimensions.
    ///
    /// @param width int& The width of the viewport
    /// @param height int& The height of the viewport
    void Init(const int &width, const int &height);

    /// Sets the color of the given pixel to the given RGB value.
    ///
    /// @param x int The x coordinate
    /// @param y int The y coordinate
    /// @param r uint8_t The R value
    /// @param g uint8_t The G value
    /// @param b uint8_t The B value
    void SetPixelRGB(int x, int y, uint8_t r, uint8_t g, uint8_t b);

    /// Renders the current RGB pixel values to the viewport.
    ///
    /// @return int A status flag. Is -1 if the screen is locked, 1 otherwise
    int Render();

private:
    /// @var _width int The width of the viewport.
    int _width;

    /// @var _height int The height of the viewport.
    int _height;

    /// @var _window SDL_Window* The SDL window.
    SDL_Window *_window;

    /// @var _screen SLD_Surface* The SDL surface onto which we render
    SDL_Surface *_screen;
};
#endif // __RENDERER_HPP

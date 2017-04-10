#include "SDL2/SDL.h"
#include "kernel/kernel.hpp"

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
    void Init(const int &width, const int &height, const int &border);

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

    /// Renders the positions of the particles in a debug view, where each
    /// particle is represented on the surface with it's x and z position
    /// corresponding to the viewport and the size of the square corresponding
    /// to the z position.
    ///
    /// @param position double* The particle positions
    /// @param N int The number of particles
    /// @param time double The current simulation time
    void DebugViewPositions(double* position, int N, double time);

    /// Renders the density of the fluid in a debug view. The density is
    /// interpolated along the (x, 0.5, z) plane, where x and z go from 0 to 1.
    ///
    /// @param density double* The density of the particles
    /// @param position double* The position of the particles
    /// @param kernel Kernel* The kernel to use for interpolation
    void DebugViewSurface(double* density, double* position, Kernel* kernel);

private:
    /// @var _width int The width of the viewport.
    int _width;

    /// @var _height int The height of the viewport.
    int _height;

    /// @var _border int Width of the border region to the sides and below the
    /// domain region
    int _border;

    /// @var _title char* The title of the window.
    char* _title;

    /// @var _window SDL_Window* The SDL window.
    SDL_Window *_window;

    /// @var _screen SLD_Surface* The SDL surface onto which we render
    SDL_Surface *_screen;
};
#endif // __RENDERER_HPP

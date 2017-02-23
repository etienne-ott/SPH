#include "renderer.hpp"
#include "SDL2/SDL.h"
#include <cmath>
#include <limits>

Renderer::Renderer() {
    SDL_Init(SDL_INIT_VIDEO);
}

Renderer::~Renderer() {
    SDL_DestroyWindow(_window);
    SDL_Quit();
}

void Renderer::SetPixelRGB(int x, int y, uint8_t r, uint8_t g, uint8_t b) {
    uint32_t *pixmem32;
    uint32_t colour;

    colour = SDL_MapRGB(_screen->format, r, g, b);

    pixmem32 = (uint32_t *)_screen->pixels + y * _screen->w + x;
    *pixmem32 = colour;
}

void Renderer::Init(const int &width, const int &height) {
    _width = width;
    _height = height;
    _window = SDL_CreateWindow("Smoothed particle hydrodynamics", SDL_WINDOWPOS_UNDEFINED,
    SDL_WINDOWPOS_UNDEFINED, _width, _height, 0);
    _screen = SDL_GetWindowSurface(_window);

    for (int x = 0; x < _width; ++x) {
        for (int y = 0; y < _height; ++y) {
            this->SetPixelRGB(x, y, 0, 0, 0);
        }
    }
}

int Renderer::Render() {
    if (SDL_MUSTLOCK(_screen)) {
        if (SDL_LockSurface(_screen) < 0) {
            return -1;
        }
    }

    if (SDL_MUSTLOCK(_screen)) {
        SDL_UnlockSurface(_screen);
    }

    SDL_UpdateWindowSurface(_window);

    return 1;
}
#include "renderer.hpp"
#include "kernel.hpp"
#include "SDL2/SDL.h"
#include <algorithm>

using namespace std;

Renderer::Renderer() {
    _title = new char[500];
    sprintf(_title, "Smoothed particle hydrodynamics");

    SDL_Init(SDL_INIT_VIDEO);
}

Renderer::~Renderer() {
    delete[] _title;

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
    _window = SDL_CreateWindow(_title, SDL_WINDOWPOS_UNDEFINED,
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

    SDL_SetWindowTitle(_window, _title);
    SDL_UpdateWindowSurface(_window);

    return 1;
}

void Renderer::DebugViewPositions(double* position, int N, double time) {
    for (int i = 0; i < _width; i++) {
        for (int j = 0; j < _height; j++) {
            this->SetPixelRGB(i, j, 255, 255, 255);
        }
    }

    for (int i = 0; i < N; i++) {
        int iposx = (int)(_width * position[i * 3]),
            iposz = _height - (int)(_height * position[i * 3 + 2]),
            yrange = max(0, min(10, (int)(10 * (1.0 - position[i * 3 + 1]))));
        for (int k = iposx - yrange < 0 ? 0 : iposx - yrange; k < iposx + yrange && k < _width; k++) {
            for (int l = iposz - yrange < 0 ? 0 : iposz - yrange; l < iposz + yrange && l < _height; l++) {
                this->SetPixelRGB(k, l, 0, 0, 0);
            }
        }
    }

    sprintf(_title, "Smoothed particle hydrodynamics (t=%f)", time);

    this->Render();
}

void Renderer::DebugViewSurface(double* density, double* position, Kernel* kernel) {
    double dx = 1.0 / _width, dz = 1.0 / _height;

    for (int i = 0; i < _width; i++) {
        for (int j = 0; j < _height; j++) {
            double d = kernel->InterpolateDensity(i * dx, 0.5, j * dz, density, position);
            if (d > 0.5) {
                this->SetPixelRGB(i, j, 0, 0, 0);
            } else {
                this->SetPixelRGB(i, j, 255, 255, 255);
            }
        }
    }

    this->Render();
}
#pragma once

#include "SDL2/SDL.h"
#include <vector>
#include "data/mesh.h"

enum Color {
    red, green, blue, yellow, white, black
};

class DebugRenderer {
public:
    /// Constructor.
    DebugRenderer();

    /// Destructor.
    ~DebugRenderer();

    /// Initialises the renderer with the given dimensions.
    ///
    /// @param width int& The width of the viewport
    /// @param height int& The height of the viewport
    void Init(const int &width, const int &height);

    void setCameraPosition(float x, float y, float z);

    void fitViewToMesh(Mesh* m);

    /// Sets the color of the given pixel to the given RGB value.
    ///
    /// @param x int The x coordinate
    /// @param y int The y coordinate
    /// @param r uint8_t The R value
    /// @param g uint8_t The G value
    /// @param b uint8_t The B value
    void SetPixelRGB(int x, int y, uint8_t r, uint8_t g, uint8_t b);

    void SetPixel(int x, int y, Color c);
    void SetPixel(int x, int y, uint32_t c);

    void ClearScreen();

    void DrawLine(float screenX1, float screenY1, float screenX2, float screenY2);
    void DrawLine(float screenX1, float screenY1, float screenX2, float screenY2, Color c);

    void DrawSquareLLB(float x, float y, uint size, Color c);
    void DrawSquareLLB(float x, float y, uint size, uint32_t c);

    void DrawWireframe(Mesh* mesh);
    void DrawWireframe(Mesh* mesh, Color c);

    void DrawNormals(Mesh* mesh);

    void DrawPoints(std::vector<Vector3D<float>>& points);
    void DrawPoints(float* coords, int N, float sign);
    void DrawPoints(float* coords, float* colorVals, int N, float sign);

    /// Renders the current RGB pixel values to the viewport.
    ///
    /// @return int A status flag. Is -1 if the screen is locked, 1 otherwise
    int Render();

private:
    /// @var _width int The width of the viewport.
    int _width;

    /// @var _height int The height of the viewport.
    int _height;

    /// @var _title char* The title of the window.
    char* _title;

    float _camera_x;
    float _camera_y;
    float _camera_z;
    float _camera_fov;

    std::vector<uint32_t>* _color_map;

    /// @var _window SDL_Window* The SDL window.
    SDL_Window *_window;

    /// @var _screen SLD_Surface* The SDL surface onto which we render
    SDL_Surface *_screen;
};
#include "debug_renderer.h"
#include "data/vector3D.h"
#include <cmath>

DebugRenderer::DebugRenderer() {
    _title = new char[500];
    sprintf(_title, "sdl-play");

    _color_map = new std::vector<uint32_t>();

    SDL_Init(SDL_INIT_VIDEO);
}

DebugRenderer::~DebugRenderer() {
    delete[] _title;
    delete _color_map;

    SDL_DestroyWindow(_window);
    SDL_Quit();
}

void DebugRenderer::SetPixelRGB(int x, int y, uint8_t r, uint8_t g, uint8_t b) {
    uint32_t *pixmem32;
    uint32_t colour;

    colour = SDL_MapRGB(_screen->format, r, g, b);

    pixmem32 = (uint32_t *)_screen->pixels + y * _screen->w + x;
    *pixmem32 = colour;
}

void DebugRenderer::SetPixel(int x, int y, Color c) {
    uint32_t* pixmem32;
    pixmem32 = (uint32_t *)_screen->pixels + y * _screen->w + x;
    *pixmem32 = _color_map->at(int(c));
}

void DebugRenderer::SetPixel(int x, int y, uint32_t c) {
    uint32_t* pixmem32;
    pixmem32 = (uint32_t *)_screen->pixels + y * _screen->w + x;
    *pixmem32 = c;
}

void DebugRenderer::DrawLine(float screenX1, float screenY1, float screenX2, float screenY2) {
    this->DrawLine(screenX1, screenY1, screenX2, screenY2, Color::white);
}

void DebugRenderer::DrawLine(float screenX1, float screenY1, float screenX2, float screenY2, Color c) {
    float deltaX = screenX2 - screenX1;
    float deltaY = screenY2 - screenY1;
    int k = int(sqrt(deltaX * deltaX + deltaY * deltaY) * _width * _height);
    float invk = 1.f / k;
    int cx = 0, cy = 0;
    for (int i = 0; i < k; i++) {
        cx = int(screenX1 * _width + deltaX * _width * invk * i);
        cy = int(screenY1 * _height + deltaY * _height * invk * i);
        this->SetPixel(cx, _height - cy, c);
    }
}

void DebugRenderer::DrawSquareLLB(float x, float y, uint size, Color c) {
    this->DrawSquareLLB(x, y, size, _color_map->at(int(c)));
}

void DebugRenderer::DrawSquareLLB(float x, float y, uint size, uint32_t c) {
    int sx, sy;

    for (uint i = 0; i < size; i++) {
        for (uint j = 0; j < size; j++) {
            sx = int(x * _width) + i;
            sy = int(y * _height) + j;
            if (sx < 0 || sx >= _width || sy <= 0 || sy > _height) {
                continue;
            }

            this->SetPixel(sx, _height - sy, c);
        }
    }
}

void DebugRenderer::DrawWireframe(Mesh* mesh) {
    float fov_fac_x = sin(2 * M_PI * _camera_fov * 0.5f / 360.f) * 2.f;
    float fov_fac_y = fov_fac_x * _height / _width;

    std::vector<int> faces = mesh->getFaces();
    std::vector<Vector3D<float>> vertices = mesh->getVertices();

    for (uint i = 0; i < faces.size(); i += 3) {
        float vtx1_x = vertices[faces[i]].getX() - _camera_x;
        float vtx1_y = vertices[faces[i]].getY() - _camera_y;
        float vtx1_z = vertices[faces[i]].getZ() + _camera_z;

        float vtx2_x = vertices[faces[i + 1]].getX() - _camera_x;
        float vtx2_y = vertices[faces[i + 1]].getY() - _camera_y;
        float vtx2_z = vertices[faces[i + 1]].getZ() + _camera_z;

        float vtx3_x = vertices[faces[i + 2]].getX() - _camera_x;
        float vtx3_y = vertices[faces[i + 2]].getY() - _camera_y;
        float vtx3_z = vertices[faces[i + 2]].getZ() + _camera_z;

        this->DrawLine(
            0.5 + (vtx1_x / (fov_fac_x * vtx1_z)),
            0.5 + (vtx1_y / (fov_fac_y * vtx1_z)),
            0.5 + (vtx2_x / (fov_fac_x * vtx2_z)),
            0.5 + (vtx2_y / (fov_fac_y * vtx2_z))
        );

        this->DrawLine(
            0.5 + (vtx2_x / (fov_fac_x * vtx2_z)),
            0.5 + (vtx2_y / (fov_fac_y * vtx2_z)),
            0.5 + (vtx3_x / (fov_fac_x * vtx3_z)),
            0.5 + (vtx3_y / (fov_fac_y * vtx3_z))
        );

        this->DrawLine(
            0.5 + (vtx3_x / (fov_fac_x * vtx3_z)),
            0.5 + (vtx3_y / (fov_fac_y * vtx3_z)),
            0.5 + (vtx1_x / (fov_fac_x * vtx1_z)),
            0.5 + (vtx1_y / (fov_fac_y * vtx1_z))
        );
    }
}

void DebugRenderer::DrawNormals(Mesh* mesh) {
    std::vector<int> faces = mesh->getFaces();
    std::vector<Vector3D<float>> vertices = mesh->getVertices();
    std::vector<Vector3D<float>> faceNormals = mesh->getFaceNormals();
    float* b_box = mesh->getBoundingBox();

    float fov_fac_x = sin(2 * M_PI * _camera_fov * 0.5f / 360.f) * 2.f;
    float fov_fac_y = fov_fac_x * _height / _width;
    float third = 1.f / 3.f;
    float normalScale = 0.2f * std::min(
        b_box[3] - b_box[0],
        std::min(
            b_box[4] - b_box[1],
            b_box[5] - b_box[2]
        )
    );

    for (uint i = 0; i < faceNormals.size(); i++) {
        float cog_x = (
            vertices[faces[3 * i]].getX()
            + vertices[faces[3 * i + 1]].getX()
            + vertices[faces[3 * i + 2]].getX()
        ) * third;
        float cog_y = (
            vertices[faces[3 * i]].getY()
            + vertices[faces[3 * i + 1]].getY()
            + vertices[faces[3 * i + 2]].getY()
        ) * third;
        float cog_z = (
            vertices[faces[3 * i]].getZ()
            + vertices[faces[3 * i + 1]].getZ()
            + vertices[faces[3 * i + 2]].getZ()
        ) * third;

        float vtx1_x = cog_x - _camera_x;
        float vtx1_y = cog_y - _camera_y;
        float vtx1_z = cog_z + _camera_z;

        float vtx2_x = cog_x - normalScale * faceNormals[i].getX() - _camera_x;
        float vtx2_y = cog_y - normalScale * faceNormals[i].getY() - _camera_y;
        float vtx2_z = cog_z - normalScale * faceNormals[i].getZ() + _camera_z;

        // Skip normals that are almost perfectly aligned with the camera
        // view axis as they won't be rendered correctly and would not be
        // visible anyway
        /*if (
            std::fabs(vtx1_x - vtx2_x) < 0.0001
            && std::fabs(vtx1_y - vtx2_y) < 0.0001
        ) {
            continue;
        }*/

        this->DrawLine(
            0.5 + (vtx1_x / (fov_fac_x * vtx1_z)),
            0.5 + (vtx1_y / (fov_fac_y * vtx1_z)),
            0.5 + (vtx2_x / (fov_fac_x * vtx2_z)),
            0.5 + (vtx2_y / (fov_fac_y * vtx2_z)),
            Color::red
        );
    }
}

void DebugRenderer::DrawPoints(std::vector<Vector3D<float>>& points) {
    float fov_fac_x = sin(2 * M_PI * _camera_fov * 0.5f / 360.f) * 2.f;
    float fov_fac_y = fov_fac_x * _height / _width;

    for (uint i = 0; i < points.size(); i++) {
        this->DrawSquareLLB(
            0.5 + ((points[i].getX() - _camera_x)
                / (fov_fac_x * (points[i].getZ() + _camera_z))),
            0.5 + ((points[i].getY() - _camera_y)
                / (fov_fac_y * (points[i].getZ() + _camera_z))),
            2,
            Color::blue
        );
    }
}

void DebugRenderer::DrawPoints(float* coords, int N, float sign) {
    float fov_fac_x = sin(2 * M_PI * _camera_fov * 0.5f / 360.f) * 2.f;
    float fov_fac_y = fov_fac_x * _height / _width;

    for (int i = 0; i < N; i++) {
        // the sign is for when we work in a coordinate system
        // with flipped z-axis
        this->DrawSquareLLB(
            0.5 + ((coords[i * 3] - _camera_x)
                / (fov_fac_x * (sign * coords[i * 3 + 2] + _camera_z))),
            0.5 + ((coords[i * 3 + 1] - _camera_y)
                / (fov_fac_y * (sign * coords[i * 3 + 2] + _camera_z))),
            2,
            Color::blue
        );
    }
}

void DebugRenderer::DrawPoints(float* coords, float* colorVals, int N, float sign) {
    float fov_fac_x = sin(2 * M_PI * _camera_fov * 0.5f / 360.f) * 2.f;
    float fov_fac_y = fov_fac_x * _height / _width;
    int r, g, b;

    for (int i = 0; i < N; i++) {
        // the sign is for when we work in a coordinate system
        // with flipped z-axis
        r = std::min(255, std::max(0, int(colorVals[i] * 255)));
        g = 0;
        b = std::min(255, std::max(0, int(255 - colorVals[i] * 255)));

        this->DrawSquareLLB(
            0.5 + ((coords[i * 3] - _camera_x)
                / (fov_fac_x * (sign * coords[i * 3 + 2] + _camera_z))),
            0.5 + ((coords[i * 3 + 1] - _camera_y)
                / (fov_fac_y * (sign * coords[i * 3 + 2] + _camera_z))),
            2,
            SDL_MapRGB(_screen->format, r, g, b)
        );
    }
}

void DebugRenderer::Init(const int &width, const int &height) {
    _width = width;
    _height = height;

    _window = SDL_CreateWindow(_title, SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED, _width, _height, 0);
    _screen = SDL_GetWindowSurface(_window);

    _color_map->push_back(SDL_MapRGB(_screen->format, 255, 0, 0));
    _color_map->push_back(SDL_MapRGB(_screen->format, 0, 255, 0));
    _color_map->push_back(SDL_MapRGB(_screen->format, 0, 0, 255));
    _color_map->push_back(SDL_MapRGB(_screen->format, 255, 0, 255));
    _color_map->push_back(SDL_MapRGB(_screen->format, 255, 255, 255));
    _color_map->push_back(SDL_MapRGB(_screen->format, 0, 0, 0));

    this->ClearScreen();

    _camera_x = 0.f;
    _camera_y = 0.f;
    _camera_z = 2.f;
    _camera_fov = 90.f;
}

void DebugRenderer::ClearScreen() {
    for (int i = 0; i < _width; i++) {
        for (int j = 0; j < _height; j++) {
            this->SetPixel(i, j, Color::black);
        }
    }
}

void DebugRenderer::setCameraPosition(float x, float y, float z) {
    _camera_x = x;
    _camera_y = y;
    _camera_z = z;
}

void DebugRenderer::fitViewToMesh(Mesh* m) {
    float* box = m->getBoundingBox();
    _camera_x = 0.5 * (box[0] + box[3]);
    _camera_y = 0.5 * (box[1] + box[4]);
    _camera_z = 1.25f * std::max(
        (box[3] - box[0]),
        std::max(
            (box[4] - box[1]),
            (box[5] - box[2])
        )
    );
}

int DebugRenderer::Render() {
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
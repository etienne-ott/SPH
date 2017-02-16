#include "SDL2/SDL.h"

#ifndef __RENDERER_HPP
#define __RENDERER_HPP

class Renderer {
public:
  Renderer();

  ~Renderer();

  void Init(const int &width, const int &height);

  void SetPixelRGB(int x, int y, uint8_t r, uint8_t g, uint8_t b);

  int Render();

private:
  int _width;
  int _height;
  SDL_Window *_window;
  SDL_Surface *_screen;
};
#endif // __RENDERER_HPP

#include "renderer.hpp"
#include "kernel.hpp"
#include "vtk.hpp"
#include "compute.hpp"

#define WRITE_VTK_OUTPUT false

int main() {
    double h = 0.3;
    int N = 50;
    int R = 400;
    double t = 0.0;
    double tend = 50.0;
    double dt = 0.1;

    Renderer r = Renderer();
    r.Init(R, R);

    Kernel kernel = Kernel(h, N, 1/N);

    Compute compute = Compute(N, &kernel);

    VTK vtk = VTK("VTK/", &kernel, 20);

    r.DebugViewPositions(compute.GetPosition(), N);

    while (t < tend) {
        compute.Timestep();

        if (WRITE_VTK_OUTPUT) {
            vtk.WriteDensity(compute.GetDensity(), compute.GetPosition());
        }

        r.DebugViewPositions(compute.GetPosition(), N);

        t += dt;
    }

    bool running = true;
    SDL_Event event;

    while (running) {
        SDL_Delay(30);
        if (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                running = false;
            }
        }
    }

    return 0;
}
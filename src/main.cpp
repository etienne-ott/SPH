#include "renderer.hpp"
#include "kernel.hpp"
#include "vtk.hpp"
#include "compute.hpp"

#define WRITE_VTK_OUTPUT false

/// Checks if there has been a SDL_QUIT event since the last time SLD events
/// were checked.
///
/// @param event SDL_Event* A variable where the event (if any) will be stored
/// @return bool If there has been a SDL_QUIT event
bool checkQuitSLDEvent(SDL_Event* event) {
    if (SDL_PollEvent(event)) {
        if (event->type == SDL_QUIT) {
            return true;
        }
    }
    return false;
}

int main() {
    double h = 0.3;
    int N = 50;
    int R = 400;
    double t = 0.0;
    double tend = 50.0;
    double dt = 0.1;

    Renderer r = Renderer();
    r.Init(R, R);

    Kernel kernel = Kernel(h, N, 1.0 / N);

    Compute compute = Compute(N, &kernel);

    VTK vtk = VTK("VTK/", &kernel, 20);

    bool running = true;
    SDL_Event event;

    r.DebugViewPositions(compute.GetPosition(), N);

    while (t < tend && running) {
        compute.Timestep();

        if (WRITE_VTK_OUTPUT) {
            vtk.WriteDensity(compute.GetDensity(), compute.GetPosition());
        }

        r.DebugViewPositions(compute.GetPosition(), N);
        running = !checkQuitSLDEvent(&event);

        t += dt;
    }

    while (running) {
        SDL_Delay(30);
        running = !checkQuitSLDEvent(&event);
    }

    return 0;
}
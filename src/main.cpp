#include "renderer.hpp"
#include "kernel.hpp"
#include "vtk.hpp"
#include "ascii_output.hpp"
#include "compute.hpp"
#include "parameter.hpp"

#define WRITE_VTK_OUTPUT false
#define WRITE_ASCII_OUTPUT false

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
    Parameter param = Parameter();
    param.Load("sample.param");

    Renderer r = Renderer();
    r.Init(param.R, param.R);

    Kernel kernel = Kernel(param.h, param.N, param.mass);
    Compute compute = Compute(&param, &kernel);

    VTK vtk = VTK("output/vtk/", &kernel, 20);
    ASCIIOutput ascii = ASCIIOutput("output/ascii/");

    bool running = true;
    SDL_Event event;
    double t = 0.0;

    r.DebugViewPositions(compute.GetPosition(), param.N, t);

    while (t < param.tend && running) {
        compute.Timestep();

        if (WRITE_VTK_OUTPUT) {
            vtk.WriteDensity(compute.GetDensity(), compute.GetPosition());
        }

        if (WRITE_ASCII_OUTPUT) {
            ascii.WriteParticleStatus(compute.GetDensity(), compute.GetPosition(), &param);
        }

        r.DebugViewPositions(compute.GetPosition(), param.N, t);
        running = !checkQuitSLDEvent(&event);

        t += param.dt;
    }

    while (running) {
        SDL_Delay(30);
        running = !checkQuitSLDEvent(&event);
    }

    return 0;
}
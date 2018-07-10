#include "output/debug_renderer.h"
#include "kernel/bidomain.h"
#include "kernel/gaussian.h"
#include "kernel/cubic_spline.h"
#include "output/vtk.h"
#include "output/ascii_output.h"
#include "simulation/compute.h"
#include <yaml-cpp/yaml.h>

/// Checks if there has been a SDL_QUIT event since the last time SDL events
/// were checked.
///
/// @param event SDL_Event* A variable where the event (if any) will be stored
/// @return bool If there has been a SDL_QUIT event
bool checkQuitSDLEvent(SDL_Event* event) {
    if (SDL_PollEvent(event)) {
        if (event->type == SDL_QUIT) {
            return true;
        }
    }
    return false;
}

void drawDebugView(DebugRenderer& r, Compute& c, YAML::Node& param) {
    r.ClearScreen();
    r.DrawPoints(c.GetPosition(), param["N"].as<int>(), 1.f);
    r.Render();
}

int main() {
    YAML::Node param = YAML::LoadFile("default_parameter.yaml");

    DebugRenderer renderer = DebugRenderer();
    renderer.Init(param["r_width"].as<int>(), param["r_height"].as<int>());
    renderer.setCameraPosition(
        param["camera_x"].as<float>(),
        param["camera_y"].as<float>(),
        param["camera_z"].as<float>()
    );

    CubicSpline kernel = CubicSpline(param["h"].as<float>(), param["N"].as<int>(), param["mass"].as<float>());
    Compute compute = Compute(param, &kernel);

    VTK vtk = VTK("output/vtk/", &kernel, 20);
    ASCIIOutput ascii = ASCIIOutput("output/ascii/");

    bool running = true;
    SDL_Event event;
    float t = 0.0;

    drawDebugView(renderer, compute, param);

    while (t < param["tend"].as<float>() && running) {
        compute.Timestep();

        if (param["write_vtk"].as<bool>()) {
            vtk.WriteDensity(compute.GetDensity(), compute.GetPosition());
        }

        if (param["write_ascii"].as<bool>()) {
            ascii.WriteParticleStatus(
                compute.GetDensity(),
                compute.GetPosition(),
                compute.GetPressure(),
                param
            );
        }

        drawDebugView(renderer, compute, param);
        running = !checkQuitSDLEvent(&event);
        t += param["dt"].as<float>();
    }

    while (running) {
        SDL_Delay(30);
        running = !checkQuitSDLEvent(&event);
    }

    return 0;
}
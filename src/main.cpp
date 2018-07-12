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

    // We draw the bounding box as a mesh, but need to load
    // this information only once
    static Mesh box = Mesh();
    static bool isLoaded = false;
    if (!isLoaded) {
        box.loadMeshFromOBJFile(param["bbox_mesh"].as<std::string>());
        box.scaleTo(std::max(
            param["bbox_x_upper"].as<float>() - param["bbox_x_lower"].as<float>(),
            std::max(
                param["bbox_y_upper"].as<float>() - param["bbox_y_lower"].as<float>(),
                param["bbox_z_upper"].as<float>() - param["bbox_z_lower"].as<float>()
            )
        ));
        box.centerOn(Vector3D<float>(
            param["bbox_x_upper"].as<float>() + param["bbox_x_lower"].as<float>(),
            param["bbox_y_upper"].as<float>() + param["bbox_y_lower"].as<float>(),
            param["bbox_z_upper"].as<float>() + param["bbox_z_lower"].as<float>()
        ) * 0.5f);
        isLoaded = true;
        printf("Loaded boundary mesh\n");
    }
    r.DrawWireframe(&box, Color::red);

    // We draw the initialization domain as a mesh, if set, but need to
    // load this information only once
    if (param["domain_type"].as<std::string>() == "mesh") {
        static Mesh domMesh = Mesh();
        static bool domMeshIsLoaded = false;
        if (!domMeshIsLoaded) {
            domMesh.loadMeshFromOBJFile(param["mesh_file"].as<std::string>());
            domMeshIsLoaded = true;
            printf("Loaded initialization mesh\n");
        }
        r.DrawWireframe(&domMesh, Color::green);
    }

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

    float h = param["h"].as<float>();
    float psize = param["particle_size"].as<float>();
    float rho0 = param["rho0"].as<float>();
    param["mass"] = 8.f * psize * psize * psize * rho0;
    printf("Calculated mass is %f\n", param["mass"].as<float>());
    printf("Calculated cubic volume is %f\n", 8.f * psize * psize * psize);
    printf("Calculated spheric volume is %f\n", 4.f / 3.f * M_PI * psize * psize * psize);

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